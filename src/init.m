%% MODEL SETUP & INITIAL CONDITIONS

% load colorbar
load ocean
TINY = 1e-16;

% initialise grid coordinates
h = D./Nz;                       % grid spacing [m]
W = D*Nx/Nz;                     % domain width [m]
x = linspace(h/2,W-h/2,Nx);      % x-coordinate vector
z = linspace(h/2,D-h/2,Nz);      % z-coordinate vector
[X,Z] = meshgrid(x,z);           % coordinate arrays

% set ghosted index arrays
if strcmp(BC_VP{2},'periodic')
    icx = [Nx,1:Nx,1];
    icz = [1,1:Nz,Nz];
else
    icx = [1,1:Nx,Nx];
    icz = [1,1:Nz,Nz];
end

if bnchm; mms; else  % construct manufactured solution if running benchmark [AP: Where is mms defined?]

% initialise smooth random noise
rng(15); 
smth = smth*Nz^2*1e-4;
rp   = randn(Nz,Nx);
for i = 1:round(smth*2)
    rp = rp + diffus(rp,ones(size(rp))/8,1,[1,2],BC_VP);
    rp = rp - mean(mean(rp));
end
rp = rp./max(abs(rp(:)));

% set basic initial conditions
switch finit  % initial porosity
    case 'linear'
        f = f0 + (f1-f0) .* Z/D;
    case 'layer'
        f = f0 + (f1-f0) .* (1+erf((Z/D-zlay)/wlay))/2;
end
switch Tinit  % initial temperature
    case 'linear'
        T = T0 + (T1-T0) .* Z/D;
    case 'layer'
        T = T0 + (T1-T0) .* (1+erf((Z/D-zlay)/wlay))/2;
end
switch Cinit  % initial salinity
    case 'linear'
        C = C0 + (C1-C0) .* Z/D;
    case 'layer'
        C = C0 + (C1-C0) .* (1+erf((Z/D-zlay)/wlay))/2;
end

% update initial condition within structures
if ~any(isnan(unit(:)))
    for i = 1:size(unit,3)
        if ~isnan(fstruct(i)); f = f + unit(:,:,i).*fstruct(i); end
        if ~isnan(Tstruct(i)); T = T + unit(:,:,i).*Tstruct(i); end
        if ~isnan(Cstruct(i)); C = C + unit(:,:,i).*Cstruct(i); end
    end
end

% add temperature loaded from array
if exist('fArray','var')
    f = f + fArray;
end
if exist('TArray','var')
    T = T + TArray;
end
if exist('CArray','var')
    C = C + CArray;
end

% prepare treatment for water and air
if exist('wat','var')
    rp(wat==1) = 0;
    f (wat==1) = 0.5;
    T (wat==1) = T_wat;
    C (wat==1) = C_wat;
    wat_surf   = diff(wat,1)>0;
    wat_base   = diff(wat,1)<0;
    [jbed, ibed] = find(diff(wat(icz,icx))<0);
    sub2ind     = sub2ind(size(wat),jbed,ibed);
else
    wat = zeros(Nz,Nx);
end
if exist('air','var') 
    rp(air==1) = 0;
    f (air==1) = 1.0;
    T (air==1) = T_air;
    C (air==1) = 0;
else
    air = zeros(Nz,Nx);
end

% smoothing applied to initial fields to avoid sharp contrasts
for i=1:smth
    f = f + diffus(f,ones(size(f))/8,1,[1,2],{'',''});
    T = T + diffus(T,ones(size(T))/8,1,[1,2],{'',''});
    C = C + diffus(C,ones(size(C))/8,1,[1,2],{'',''});
    wat = wat + diffus(wat,ones(size(wat))/8,1,[1,2],{'',''});
    air = air + diffus(air,ones(size(air))/8,1,[1,2],{'',''});
end

% add smooth random perturbations
f = f + df.*rp;
T = T + dT.*rp;
C = C + dC.*rp;

% enforce bounds on porosity
f = max(1e-3,min(1-1e-3,f));

% adjust boundary layer to top boundary conditions
T = T + (Ttop-T).*exp(-max(0,Z)/h);
C = C + (Ctop-C).*exp(-max(0,Z)/h);

% adjustment for treatment of air and water
f = (1-wat-air).*f + wat.*0.5   + air.*1.0  ;
T = (1-wat-air).*T + wat.*T_wat + air.*T_air;
C = (1-wat-air).*C + wat.*C_wat + air.*0.0  ;

res_T = zeros(Nz,Nx);  % residual for temperature equation
res_C = zeros(Nz,Nx);  % residual for salinity equation
res_V = zeros(Nz,Nx);  % residual for vapour equation

% initialise vapour phase
Plith = 2600.*grav.*Z;
Vq = vapour(T,C,Plith);
V  = Vq;
phsr_V = 0.*V;

end

% store initial fields
fin = f;
Tin = T;
Cin = C;
Vin = V;

% initialise evolving values 
if wat_evolve
    T_wat = mean(T(wat>=1),'all');
    C_wat = mean(C(wat>=1),'all');
    V_wat = mean(V(wat>=1),'all');
end

% get permeability [m2]
k = k0 * f(icz,icx).^n;  % Kozeny-Carman relationship

% get Darcy coefficient [m2/Pas]
K  = k/mu;
Kz = (K(1:end-1,:)+K(2:end,:))./2;
Kx = (K(:,1:end-1)+K(:,2:end))./2;

% get iterative step size for p-solution
dtau = (h/2)^2./K(2:end-1,2:end-1);

% update density difference
rho  = rhol0.*(1 - aT.*T(icz,icx) ...
                 - aC.*C(icz,icx) ...
                 - aV.*V(icz,icx) );% .* (1-air(icz,icx));
Drho = rho - mean(rho,2);
% Drho(air(icz,icx)+wat(icz,icx)>=0.5) = 0;  % set air and water to zero
Drhoz = (Drho(1:end-1,:)+Drho(2:end,:))./2;
if bnchm; Drhoz = Drho_mms(:,:,step+1); end

% get dimensionless numbers
Pr0  = 1./geomean(K(:))./rhol0./kT;
ScC0 = 1./geomean(K(:))./rhol0./kC;
ScV0 = 1./geomean(K(:))./rhol0./kV;
LeC0 = kT./kC;
LeV0 = kT./kV;
RbC0 = abs(aT.*(Ttop-Tbot)) ./ abs(aC.*(Ctop-Cbot));
RbV0 = abs(aT.*(Ttop-Tbot)) ./ (aV.*(max(V(:))-min(V(:))));
RaT0 = rhol0 .* max(0,-aT.*(Ttop-Tbot))   .* grav .* geomean(K(:)) .* D ./ kT;
RaC0 = rhol0 .* max(0,-aC.*(Ctop-Cbot))   .* grav .* geomean(K(:)) .* D ./ kC;
RaV0 = rhol0 .* aV.*(max(V(:))-min(V(:))) .* grav .* geomean(K(:)) .* D ./ kV;
Ra   = (RaT0+RaC0+RaV0).*ones(size(T));

% prepare solution & residual arrays for VP solver
w = zeros(Nz+1,Nx+2);   % vertical Darcy speed
u = zeros(Nz+2,Nx+1);   % horizontal Darcy speed
p = zeros(Nz+2,Nx+2);   % pore fluid pressure
res_p = zeros(Nz,Nx)./dtau;  % residual for pressure equation

% initialise timing parameters
dTdt = 0.*T;
dCdt = 0.*C;
dVdt = 0.*V;
step = 0;
time = 0;


%% PLOT INITIAL CONDITIONS

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',14};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};

axh = 6.00; axw = 7.50; %   Height and width of axis
ahs = 1.50; avs = 1.00; %   Horzontal and vertial distance between axis
axb = 1.50; axt = 0.90; %   Bottom and top;Size of page relative to axis
axl = 1.75; axr = 0.90; %   Right and left; spacing of axis to page

% prepare and plot figure for mechanical solution fields
fh1 = figure(1); clf; colormap(ocean);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 2*axw + 1*ahs + axr;
set(fh1,UN{:},'Position',[3 3 fw fh]);
set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh1,'Color','w','InvertHardcopy','off');
set(fh1,'Resize','off');
ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(2) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(3) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(4) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

set(fh1, 'CurrentAxes', ax(1))
imagesc(x,z,f); axis equal tight; box on; cb = colorbar; hold on
if ~any(isnan(unit(:)))
    for i = 1:size(unit,3)
        contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
    end
end
clim([min(f(:)),max(f(:))])
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Initial Porosity [vol]',TX{:},FS{:})

set(fh1, 'CurrentAxes', ax(2))
imagesc(x,z,C); axis equal tight;  box on; cb = colorbar; hold on
if ~any(isnan(unit(:)))
    for i = 1:size(unit,3)
        contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
    end
end
clim([min(C(:)),max(C(:))])
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:});title('Initial Salinity [wt]',TX{:},FS{:})

set(fh1, 'CurrentAxes', ax(3))
imagesc(x,z,T); axis equal tight;  box on; cb = colorbar; hold on
if ~any(isnan(unit(:)))
    for i = 1:size(unit,3)
        contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
    end
end
clim([min(T(:)),max(T(:))])
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:});title('Initial Temperature [C]',TX{:},FS{:})

set(fh1, 'CurrentAxes', ax(4))
imagesc(x,z,V); axis equal tight; box on; cb = colorbar; hold on
if ~any(isnan(unit(:)))
    for i = 1:size(unit,3)
        contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
    end
end
clim([min(V(:)-eps),max(V(:)+eps)])
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Initial Vapour [wt]',TX{:},FS{:})
drawnow

% print figure to file
if svout
    print(fh1,[outdir,'/',runID,'/',runID,'_init_',int2str(step/nout)],'-dpng','-r200')
        save([outdir,'/',runID,'/',runID,'_',int2str(step/nout)],'x','z','u','w','p','f','T','C','V','dTdt','dCdt','dVdt','K','Drho','time','Ra','unit');
end