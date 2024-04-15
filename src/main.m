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

% load colorbar
load ocean
TINY = 1e-16;

% initialise grid coordinates
h = D./N;                        % grid spacing [m]
x = linspace(h/2,D-h/2,N);       % x-coordinate vector
z = linspace(h/2,D-h/2,N);       % z-coordinate vector
[X,Z] = meshgrid(x,z);           % coordinate arrays

% set ghosted index arrays
if strcmp(BC_VP{2},'periodic')
    icx = [N,1:N,1];
    icz = [1,1:N,N];
else
    icx = [1,1:N,N];
    icz = [1,1:N,N];
end

if bnchm; mms; else  % construct manufactured solution if running benchmark [AP: Where is mms defined?]

% initialise smooth random noise
rng(15); 
smth = smth*N^2*1e-4;
rp   = randn(N,N);
for i = 1:round(smth)
    rp = rp + diffus(rp,ones(size(rp))/8,1,[1,2],{'',BC_VP{2}});
    rp = rp - mean(mean(rp));
end
rp = rp./max(abs(rp(:)));

% set basic initial conditions
switch finit  % initial porosity
    case 'linear'
        f = f0 + (f1-f0) .* Z/D;
    case 'layer'
        f = f0 + (f1-f0) .* (1+erf((Z/D-zlay)/wlay))/2;
    case 'array'
        f = fArray;
end
switch Tinit  % initial temperature
    case 'linear'
        T = T0 + (T1-T0) .* Z/D;
    case 'layer'
        T = T0 + (T1-T0) .* (1+erf((Z/D-zlay)/wlay))/2;
    case 'array'
        T = TArray;
end
switch Cinit  % initial salinity
    case 'linear'
        C = C0 + (C1-C0) .* Z/D;
    case 'layer'
        C = C0 + (C1-C0) .* (1+erf((Z/D-zlay)/wlay))/2;
end


% % % add linear structures (faults, aquifers, etc.)
% % % get indicator functions
% % % indstruct = zeros([size(f),length(fstruct)]);       % % Change zstruct to number of binary images (number of units) k in image prep stuff
% % % for i = 1:length(zstruct)                               
% % %     indstruct(:,:,i) = abs(Z-D/2)<=hstruct(i)/2 & abs(X-D/2)<=wstruct(i)/2;
% % %     indstruct(:,:,i) = imrotate(indstruct(:,:,i),astruct(i),'crop');
% % %     indstruct(:,:,i) = circshift(indstruct(:,:,i),-round((D/2-zstruct(i))/D*N),1);
% % %     indstruct(:,:,i) = circshift(indstruct(:,:,i),-round((D/2-xstruct(i))/D*N),2);
% % %     indstruct(  z>zstruct(i)+(hstruct(i)+abs(cosd(astruct(i))*wstruct(i)))/2 | z<zstruct(i)-(hstruct(i)+abs(cosd(astruct(i))*wstruct(i)))/2,:,i) = 0;
% % %     indstruct(:,x>xstruct(i)+(wstruct(i)+abs(sind(astruct(i))*hstruct(i)))/2 | x<xstruct(i)-(wstruct(i)+abs(sind(astruct(i))*hstruct(i)))/2,  i) = 0;
% % %     indstruct([1 end],:,i) = indstruct([2 end-1],:,i);
% % %     indstruct(:,[1 end],i) = indstruct(:,[2 end-1],i);
% % % end
% % 
% % % % add linear structures (faults, aquifers, etc.) ****************NEW
% % % % get indicator functions
% % % indstruct = zeros([size(f),numColors]);       % % Change zstruct to number of binary images (number of units) k in image prep stuff
% % % 
% % % % Smoothing function applied to structure indicator to minimise sharp interfaces
% % % for i=1:smth/2
% % %     indstruct(2:end-1,2:end-1,:) = indstruct(2:end-1,2:end-1,:) ...
% % %                             + diff(indstruct(:,2:end-1,:),2,1)./8 ...
% % %                             + diff(indstruct(2:end-1,:,:),2,2)./8;
% % %     indstruct([1 end],:,:) = indstruct([2 end-1],:,:);
% % %     indstruct(:,[1 end],:) = indstruct(:,[end-1 2],:);
% % % end

% update initial condition within structures
for i = 1:size(indstruct,3)
    if ~isnan(fstruct(i)); f = indstruct(:,:,i).*fstruct(i) + (1-indstruct(:,:,i)).*f; end
%     if ~isnan(Tstruct(i)); T = indstruct(:,:,i).*Tstruct(i) + (1-indstruct(:,:,i)).*T; end
    if ~isnan(Cstruct(i)); C = indstruct(:,:,i).*Cstruct(i) + (1-indstruct(:,:,i)).*C; end
end

% Smoothing function applied to structure indicator to minimise sharp
for i=1:smth/4
    f = f + diffus(f,ones(size(f))/8,1,[1,2],{'',''});
    T = T + diffus(T,ones(size(T))/8,1,[1,2],{'',''});
    C = C + diffus(C,ones(size(C))/8,1,[1,2],{'',''});
end

if exist('wat','var'); rp(wat==1) = 0; end

% add smooth random perturbations
f = f + df.*rp;
T = T + dT.*rp;
C = C + dC.*rp;

% enforce bounds on porosity
f = max(1e-3,min(1-1e-3,f));

% adjust boundary layer to top boundary conditions
T = T + (Ttop-T).*exp(-max(0,Z)/h);
C = C + (Ctop-C).*exp(-max(0,Z)/h);

end

% store initial fields
fin = f;
Tin = T;
Cin = C;

% get permeability [m2]
k = k0 * f(icz,icx).^n;  % Kozeny-Carman relationship

% get Darcy coefficient [m2/Pas]
K  = k/mu;
Kz = (K(1:end-1,:)+K(2:end,:))./2;
Kx = (K(:,1:end-1)+K(:,2:end))./2;

% get iterative step size for p-solution
dtau = (h/2)^2./K(2:end-1,2:end-1);

% prepare solution & residual arrays for VP solver
w = zeros(N+1,N+2);   % vertical Darcy speed
u = zeros(N+2,N+1);   % horizontal Darcy speed
p = zeros(N+2,N+2);   % pore fluid pressure
res_p = zeros(N,N)./dtau;  % residual for pressure equation


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



%% TIME STEPPING LOOP      


% initialise timing parameters
dTdt = 0.*T;
dCdt = 0.*C;
dt   = 0;
step = 0;
time = 0;



while time <= tend && step <= Nt      

        
    tic;

    fprintf(1,'\n\n*****  step = %d,  time = %1.3e [yr],  step = %1.3e [hr] \n\n',step,time/3600/24/365.25,dt/3600);

    % store previous rates of change
    To    = T;
    Co    = C;
    dTdto = dTdt;
    dCdto = dCdt;

    resnorm = 1e6;
    pi    = p;
    it    = 0;


    %% NONLINEAR SOLVER LOOP

    while resnorm >= tol && it <= maxit || it <= 100 
        
        if ~mod(it,nup)

            % UPDATE TEMPERATURE SOLUTION (SEMI-IMPLICIT UPDATE)
            
            advn_T = - advect(T,u(2:end-1,:),w(:,2:end-1),h,{ADVN,'vdf'},[1,2],BC_T);

            diff_T = diffus(T,kT,h,[1,2],BC_T);

            dTdt = advn_T + diff_T;

            if bnchm; dTdt = dTdt + src_T_mms(:,:,step+1); end

            res_T = (T-To)/(dt+TINY) - (dTdt + dTdto)/2;
            
            % T(wat==0) = T(wat==0) - res_T(wat==0)*dt;
            T = T - res_T*dt;
            T(wat==1) = mean(T(wat==1),'all');

            % UPDATE CONCENTRATION SOLUTION (SEMI-IMPLICIT UPDATE)
            
            % calculate salinity advection
            advn_C = - advect(C,u(2:end-1,:),w(:,2:end-1),h,{ADVN,'vdf'},[1,2],BC_C);

            diff_C = diffus(C,kC,h,[1,2],BC_C);

            dCdt = advn_C + diff_C;
            
            if bnchm; dCdt = dCdt + src_C_mms(:,:,step+1); end

            res_C = (C-Co)/(dt+TINY) - (dCdt + dCdto)/2;

            % C(wat==0) = C(wat==0)- res_C(wat==0)*dt;
            C = C - res_C*dt;
            C(wat==1) = mean(C(wat==1),'all');

            C = max(0,min(1,C));  % saveguard min/max bounds
            
            % update density difference
            Drho  = rhol0.*(- aT.*(T(icz,icx)-mean(T(icz,:),2)) + gC.*(C(icz,icx)-mean(C(icz,:),2)));
            Drhoz = (Drho(1:end-1,:)+Drho(2:end,:))./2;
            if bnchm; Drhoz = Drho_mms(:,:,step+1); end
        end

        % UPDATE VELOCITY-PRESSURE SOLUTION (PSEUDO-TRANSIENT SOLVER)

        % store previous iterative solution guesses
        pii = pi; pi = p;
        
        % calculate pressure gradient [Pa/m]
        gradPz = diff(p,1,1)./h;  % vertical gradient
        gradPx = diff(p,1,2)./h;  % horizontal gradient
        
        watfz = floor((wat([1 1:end],[1,1:end,end]) + wat([1:end end],[1,1:end,end]))/2);
        watfx = floor((wat([1,1:end,end],[1 1:end]) + wat([1,1:end,end],[1:end end]))/2);

        % calculate Darcy segregation speed [m/s]
        w = - Kz .* (gradPz - Drhoz.*grav);% .* (watfz==0);
        u = - Kx .* (gradPx              );% .* (watfx==0);
        
        % calculate residual of pressure equation
        res_p = diff(w(:,2:end-1),1,1)./h + diff(u(2:end-1,:),1,2)./h;
        if bnchm; res_p = res_p - src_p_mms(:,:,step+1); end
        % res_p(wat==1) = 0;

        % update pressure solution
        p(2:end-1,2:end-1) = pi(2:end-1,2:end-1) - alpha.*res_p.*dtau + beta.*(pi(2:end-1,2:end-1)-pii(2:end-1,2:end-1));
        
        % apply pressure boundary conditions
        if strcmp(BC_VP{1},'closed')
            p([1 end],:) = p([2 end-1],:) + [-1;1].*Drhoz([1 end],:).*grav.*h;
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
        if step>0
            dt = CFL * min([(h/2)/max(abs(w(:))) , (h/2)/max(abs(u(:))) , (h/2)^2./kT]);
            if bnchm; dt = dt_mms; end
        end

        if ~mod(it,nup)
            % get preconditioned residual norm to monitor convergence
            resnorm = norm(res_p.*dtau,2)./norm(p(2:end-1,2:end-1)+1e-6,2);

            % report convergence
            fprintf(1,'---  %d,  %e\n',it,resnorm);
        end
        
        % increment iteration count
        it = it+1;

    end  % non-linear iteration loop

    soltime = toc;  % record time to solution
    fprintf(1,'\n\n      time to solution = %1.3e [s] \n\n',soltime);


%%  PLOT AND SAVE SOLUTION
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
            print(fh2,[outdir,'/',runID,'/',runID,'_sol_',int2str(step/nout)],'-dpng','-r200')
            save([outdir,'/',runID,'/',runID,'_',int2str(step/nout)],'u','w','p','T','C','dTdt','dCdt','K','Drho');
        end
    end
    
    % update time and step count
    step = step + 1;
    time = time + dt;  

end  % time-stepping loop
