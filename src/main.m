%% CHIXCULUB: IMPACT HYDROTHERMAL CIRCULATION MODEL

fprintf('\n\n')
fprintf('*************************************************************\n');
fprintf('*****  RUN CHIXCULUB MODEL | %s  **********\n',datetime('now'));
fprintf('*************************************************************\n');
fprintf('\n   run ID: %s \n\n',runID);

% create output directory
if ~isfolder([outdir,'/',runID])
    mkdir([outdir,'/',runID]);
end

% save input parameters & options
parfile = [outdir,'/',runID,'/',runID,'_par'];
save(parfile);


%% MODEL SETUP & INITIAL CONDITIONS

% initialise grid coordinates
h = D./N;     % grid spacing
x = linspace(-h/2,D+h/2,N+2);
z = linspace(-h/2,D+h/2,N+2);
[X,Z] = meshgrid(x,z);

% initialise smooth random noise
rng(5);
rn = rand(N+2,N+2) - 0.5;
for i=1:smth
   rn(2:end-1,2:end-1) = rn(2:end-1,2:end-1) ...
                       + diff(rn(:,2:end-1),2,1)./8 ...
                       + diff(rn(2:end-1,:),2,2)./8;
    rn([1 end],:) = rn([2 end-1],:);
    rn(:,[1 end]) = rn(:,[end-1 2]);
end
rn = rn./max(abs(rn(:)));

% set initial condition
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

% add linear structures (faults, aquifers, etc.)
% get indicator functions
indstruct = zeros([size(f),length(zstruct)]);
for i = 1:length(zstruct)
    indstruct(:,:,i) = abs(Z-D/2)<=hstruct(i)/2 & abs(X-D/2)<=wstruct(i)/2;
    indstruct(:,:,i) = imrotate(indstruct(:,:,i),astruct(i),'crop');
    indstruct(:,:,i) = circshift(indstruct(:,:,i),-round(D/2-zstruct(i))/D*N,1);
    indstruct(:,:,i) = circshift(indstruct(:,:,i),-round(D/2-xstruct(i))/D*N,2);
    indstruct(z>zstruct(i)+(hstruct(i)+abs(cosd(astruct(i))*wstruct(i)))/2 | z<zstruct(i)-(hstruct(i)+abs(cosd(astruct(i))*wstruct(i)))/2,:,i) = 0;
    indstruct(:,x>xstruct(i)+(wstruct(i)+abs(sind(astruct(i))*hstruct(i)))/2 | x<xstruct(i)-(wstruct(i)+abs(sind(astruct(i))*hstruct(i)))/2,i) = 0;
    indstruct([1 end],:,i) = indstruct([2 end-1],:,i);
    indstruct(:,[1 end],i) = indstruct(:,[2 end-1],i);
end

% Smoothing function applied to structure indicator to minimise sharp interfaces
for i=1:smth/2
    indstruct(2:end-1,2:end-1,:) = indstruct(2:end-1,2:end-1,:) ...
                            + diff(indstruct(:,2:end-1,:),2,1)./8 ...
                            + diff(indstruct(2:end-1,:,:),2,2)./8;
    indstruct([1 end],:,:) = indstruct([2 end-1],:,:);
    indstruct(:,[1 end],:) = indstruct(:,[end-1 2],:);
end

% update initial condition within structures
for i = 1:length(zstruct)
    if ~isnan(fstruct(i)); f = indstruct(:,:,i).*fstruct(i) + (1-indstruct(:,:,i)).*f; end
    if ~isnan(Tstruct(i)); T = indstruct(:,:,i).*Tstruct(i) + (1-indstruct(:,:,i)).*T; end
    if ~isnan(Cstruct(i)); C = indstruct(:,:,i).*Cstruct(i) + (1-indstruct(:,:,i)).*C; end
end

% add smooth random perturbations
f = f + df.*rn;
T = T + dT.*rn;
C = C + dC.*rn;

% enforce bounds on porosity
f = max(1e-3,min(1-1e-3,f));

%
T = T + (Ttop-T).*exp(-max(0,Z)/h);
C = C + (Ctop-C).*exp(-max(0,Z)/h);

% store initial fields
fin = f;
Tin = T;
Cin = C;

% get density difference
Drho  = - rhol0.*(- aT.*(T-mean(T,2)) + gC.*(C-mean(C,2)));
Drhoz = (Drho(1:end-1,:)+Drho(2:end,:))./2;

% get permeability [m2]
k = k0 * f.^n;  % Kozeny-Carman relationship

% get Darcy coefficient [m2/Pas]
K  = k/mu;
Kz = (K(1:end-1,:)+K(2:end,:))./2;
Kx = (K(:,1:end-1)+K(:,2:end))./2;

% get iterative step size
dtau = (h/2)^2./K;

% prepare for plotting
load ocean
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
imagesc(x,z,f); axis equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Initial Porosity [vol]',TX{:},FS{:})
set(fh1, 'CurrentAxes', ax(2))
imagesc(x,z,log10(k)); axis equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Initial Permeability [log$_{10}$ m$^2$]',TX{:},FS{:})
set(fh1, 'CurrentAxes', ax(3))
imagesc(x,z,T); axis equal tight;  box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:});title('Initial Temperature [C]',TX{:},FS{:})
set(fh1, 'CurrentAxes', ax(4))
imagesc(x,z,C); axis equal tight;  box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:});title('Initial Salinity [wt]',TX{:},FS{:})
drawnow

% prepare solution & residual arrays for VP solver
w = zeros(N+1,N+2);  % vertical Darcy speed
u = zeros(N+2,N+1);  % horizontal Darcy speed
p = zeros(N+2,N+2);  % pore fluid pressure
F = ones(N+2,N+2)./dtau;  % residual for pressure equation
F([1,end],:) = 0;  F(:,[1,end]) = 0;


%% TIME STEPPING LOOP

% initialise timing parameters
dTdt = 0.*T(2:end-1,2:end-1);
dCdt = 0.*C(2:end-1,2:end-1);
dt   = 0;
step = 0;
time = 0;

while time <= tend

    fprintf(1,'\n\n*****  step = %d,  time = %1.3e [yr],  step = %1.3e [hr] \n\n',step,time/3600/24/365.25,dt/3600);

    % store previous rates of change
    To    = T;
    Co    = C;
    dTdto = dTdt;
    dCdto = dCdt;


    Fnorm = 1e6;
    pi    = p;
    it    = 0;


    %% NONLINEAR SOLVER LOOP

    while Fnorm >= tol && it <= maxit || it <= 100 
        
        if ~mod(it,nup)

            % UPDATE TEMPERATURE SOLUTION (SEMI-IMPLICIT UPDATE)
            
            advn_T = - advect(T(2:end-1,2:end-1),u(2:end-1,:),w(:,2:end-1),h,{ADVN,'vdf'},[1,2],BC_T);

            diff_T = kT.* (diff(T(:,2:end-1),2,1)./h^2 + diff(T(2:end-1,:),2,2)./h^2);

            dTdt = advn_T + diff_T;

            T(2:end-1,2:end-1) = To(2:end-1,2:end-1) + (dTdt + dTdto)/2 .* dt;
                        
            % apply temperature boundary conditions
            T(:,1  ) = T(:,end-1);  % left boundary: insulating
            T(:,end) = T(:,2    );  % right boundary: insulating
            T(1  ,:) = Ttop;        % top boundary: isothermal
            T(end,:) = Tbot;        % bottom boundary: constant flux
            

            % UPDATE CONCENTRATION SOLUTION (SEMI-IMPLICIT UPDATE)
            
            % calculate salinity advection
            advn_C = - advect(C(2:end-1,2:end-1),u(2:end-1,:),w(:,2:end-1),h,{ADVN,'vdf'},[1,2],BC_C);

            diff_C = kC.* (diff(C(:,2:end-1),2,1)./h^2 + diff(C(2:end-1,:),2,2)./h^2);

            dCdt = advn_C + diff_C;
            
            C(2:end-1,2:end-1) = Co(2:end-1,2:end-1) + (dCdt + dCdto)/2 .* dt;
            
            C = max(0,min(1,C));  % saveguard min/max bounds

            % apply salinity boundary conditions
            C(:,1  ) = C(:,end-1);  % left boundary: closed
            C(:,end) = C(:,2    );  % right boundary: closed
            C(1  ,:) = Ctop;        % top boundary: isochemical
            C(end,:) = Cbot;        % bottom boundary: isochemical
            
            % update density difference
            Drho  = - rhol0.*(- aT.*(T-mean(T,2)) + gC.*(C-mean(C,2)));
            Drhoz = (Drho(1:end-1,:)+Drho(2:end,:))./2;
        end

        % UPDATE VELOCITY-PRESSURE SOLUTION (PSEUDO-TRANSIENT SOLVER)

        % store previous iterative solution guesses
        pii = pi; pi = p;
        
        % calculate pressure gradient [Pa/m]
        gradPz = diff(p,1,1)./h;  % vertical gradient
        gradPx = diff(p,1,2)./h;  % horizontal gradient
        
        % calculate Darcy segregation speed [m/s]
        w = - Kz .* (gradPz + Drhoz.*grav);
        u = - Kx .* (gradPx + 0          );
        
        % calculate residual of pressure equation
        F(2:end-1,2:end-1) = diff(w(:,2:end-1),1,1)./h + diff(u(2:end-1,:),1,2)./h;

        % update pressure solution
        p = pi - alpha.*F.*dtau + beta.*(pi-pii);
        
        % apply pressure boundary conditions
        if strcmp(BC_VP{1},'closed')
            p([1 end],:) = p([2 end-1],:) + [1;-1].*(Drho([1 end],:)+Drho([2 end-1],:))./2.*grav.*h;
        else        
            p([1 end],:) = 0;
        end
        if strcmp(BC_VP{2},'closed')
            p(:,[1 end]) = p(:,[2 end-1]);
        elseif strcmp(BC_VP{2},'periodic')
            p(:,[1 end]) = p(:,[end-1 2]);
        else
            p(:,[1 end]) = 0;
        end

        % get physical time step
        dt = CFL * min([(h/2)/max(abs(w(:))) , (h/2)/max(abs(u(:))) , (h/2)^2./kT]);

        if ~mod(it,nup)
            % get preconditioned residual norm to monitor convergence
            Fnorm = norm(F(2:end-1,2:end-1).*dtau(2:end-1,2:end-1),2)./norm(p(2:end-1,2:end-1)+1,2);

            % report convergence
            fprintf(1,'---  %d,  %e\n',it,Fnorm);
        end
        
        % increment iteration count
        it = it+1;
    end

    % plot solution
    if ~mod(step,nout)
        if lvplt
            VIS = {'Visible','on'};
        else
            VIS = {'Visible','off'};
        end

        if ~exist('fh2','var'); fh2 = figure(VIS{:});
        else; set(0, 'CurrentFigure', fh2); clf;
        end
        
        colormap(ocean);
         
        axh = 6.00; axw = 7.50; %   Height and width of axis
        ahs = 1.50; avs = 1.00; %   Horzontal and vertial distance between axis
        axb = 1.50; axt = 2.00; %   Bottom and top;Size of page relative to axis
        axl = 1.75; axr = 0.90; %   Right and left; spacing of axis to page
        
        fh = axb + 2*axh + 1*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh2,UN{:},'Position',[9 9 fw fh]);
        set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh2,'Color','w','InvertHardcopy','off');
        set(fh2,'Resize','off');
        ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
        ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

        sgtitle(sprintf('Time elapsed %.3f yr', time/3600/24/365.25),TX{:},FS{:})
         
        set(fh2, 'CurrentAxes', ax(1))
        imagesc(x,z,-w(2:end-1,:).*3600); axis equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Segregation z-speed [m/hr]',TX{:},FS{:})
        
        set(fh2, 'CurrentAxes', ax(2))
        imagesc(x,z,u(:,2:end-1).*3600); axis equal tight;  box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Segregation x-speed [m/hr]',TX{:},FS{:})
        
        set(fh2, 'CurrentAxes', ax(3))
        imagesc(x,z,p(2:end-1,2:end-1)); axis equal tight;  box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Dynamic fluid pressure [Pa]',TX{:},FS{:})
      
        set(fh2, 'CurrentAxes', ax(4))
        imagesc(x,z,T(2:end-1,2:end-1)); axis equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Temperature [$^\circ$C]',TX{:},FS{:})
        
        set(fh2, 'CurrentAxes', ax(5))
        imagesc(x,z,C(2:end-1,2:end-1)); axis equal tight;  box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Salinity [wt]',TX{:},FS{:})
        
        set(fh2, 'CurrentAxes', ax(6))
        imagesc(x,z,Drho(2:end-1,2:end-1)); axis equal tight;  box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Density contrast [kg/m$^3$]',TX{:},FS{:})
        drawnow      
                            
        % print figure to file
        if svout
            print(fh2,[outdir,'/',runID,'/',runID,'_',int2str(step/nout)],'-dpng','-r200')
            save([outdir,'/',runID,'/',runID,'_',int2str(step/nout)],'u','w','p','T','C','dTdt','dCdt','K','Drho');
        end
    end
    
    % update time and step count
    step = step + 1;
    time = time + dt;  

end 
