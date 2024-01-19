clear; close all; clc;
    
%addpath(genpath('/home/gary/Documents/Simulations/'))
%% SET MODEL PARAMETERS

runID   = 'bnchm_p'; % run identifier tag
outdir  = '../out'; % output directory 
nout    = 1;        % print output every 'nop' steps
lvplt   = 1;        % plot figures live (1) or in background (0)     
svout   = 1;        % save figures and data to file (1)
bnchm   = 1;        % run benchmark

% set domain parameters
D       = 1e3;      % phys. domain depth [m]

% set physical parameters
mu      = 1e-4;     % pore fluid viscosity (water) [Pa s]
k0      = 1e-9;     % background permeability [m2]
n       = 3;        % permeability powerlaw [1]
rhol0   = 1000;     % fluid density [kg/m3]
grav    = 10;       % gravity [m/s2]
kC      = 5e-8;     % chemical diffusivity [m2/s]  
kT      = 5e-7;     % thermal diffusivity [m2/s]
aT      = 1e-4;     % thermal expansivity [1/K]
gC      = 1.1;      % chemical expansivity [1/wt]

% set initial condition parameters
finit   = 'layer';  % initial condition: 'linear' or 'layer'
f0      = 0.10;     % top/background initial porosity [vol]
f1      = 0.10;     % base porosity [vol]  
df      = 0.001;    % perturbation amplitude [vol]

Tinit   = 'layer';  % initial condition: 'linear' or 'layer'
Ttop    = 50;       % top boundary temperature
Tbot    =  5;       % base boundary temperature
T0      = 50;       % top/background initial temperature [C]
T1      =  5;       % base initial temperature [C]
dT      = -0.1;     % perturbation amplitude [C]

Cinit   = 'layer';  % initial condition: 'linear' or 'layer'
Ctop    = 0.010;    % top boundary concentration [wt]
Cbot    = 0.001;    % base boundary concentration [wt]
C0      = 0.010;    % top/background concentration  [wt]
C1      = 0.001;    % base concentration [wt]
dC      = 1e-4;     % perturbation amplitude [wt]

zlay    = 0.5;      % relative depth of layer boundary
wlay    = 0.02;     % relative width of layer boundary

xstruct = [];       % midpoint x-position of structures
zstruct = [];       % midpoint z-position of structures
hstruct = [];       % height of structures
wstruct = [];       % width of structures
astruct = [];       % angle of structures to horizontal (counter-clockwise)
fstruct = [];       % porosity of structures (nan = do not set)
Tstruct = [];       % temperature of structures (nan = do not set)
Cstruct = [];       % salinity of structures (nan = do not set)

smth    = 10;       % smoothness of initial fields

% set boundary conditions
BC_T    = {'periodic','periodic'};
BC_C    = {'periodic','periodic'};
BC_VP   = {'closed','periodic'};

% set model timing parameters
tend    = 1e11;      % model stopping time [s]

% set numerical solver parameters
CFL     = 0.75;      % Courant-Friedrich-Lewy number to limit time step size
ADVN    = 'weno5';   % advection scheme
nup     = 100;       % update TC-solution and check residuals every nup iter
tol     = 1e-12;     % residual tolerance for iterative solver
maxit   = 1e4;       % maximum number of iterations
alpha   = 0.99;      % step size for iterative solver
beta    = 0.97;      % damping parameter for iterative solver

Nt      = 0;
dt_mms  = 1;
NN      = [20,40,80];

for N = NN

    %*****  RUN NAKHLA MODEL  *************************************************
    run('../src/main')
    %**************************************************************************

    Err_p = p - p_mms;
    Err_T = T - T_mms;
    Err_C = C - C_mms;

    figure(100); clf;
    colormap(ocean);
    subplot(2,3,1); imagesc(x,z, p(2:end-1,2:end-1)); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $p$ [Pa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,2); imagesc(x,z, T); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $T$ [degC]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,3); imagesc(x,z, C); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $C$ [wt]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,4); imagesc(x,z, p(2:end-1,2:end-1)-p_mms(2:end-1,2:end-1,step)); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $p$ [Pa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,5); imagesc(x,z, T-T_mms(:,:,step)); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $T$ [degC]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,6); imagesc(x,z, C-C_mms(:,:,step)); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $C$ [wt]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    drawnow;

    % get solution error
    Ep = norm(p(2:end-1,2:end-1)-p_mms(2:end-1,2:end-1,step),'fro')./norm(p_mms(2:end-1,2:end-1,step),'fro');
    ET = norm(T-T_mms(:,:,step),'fro')./norm(T_mms(:,:,step),'fro');
    EC = norm(C-C_mms(:,:,step),'fro')./norm(C_mms(:,:,step),'fro');

    % plot error convergence
    fh101 = figure(101);
    p1 = loglog(h,Ep,'rs','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('grid step [m]','Interpreter','latex')
    ylabel('rel. numerical error [1]','Interpreter','latex')
    set(gca,'TicklabelInterpreter','latex')
    title('Numerical convergence in space','Interpreter','latex','FontSize',20)

    if N == NN(1)
        p4 = loglog(D./NN,Ep.*(NN(1)./NN).^2,'k-','LineWidth',2);  % plot linear trend for comparison
    end
    if N == NN(end)
        legend([p1,p4],{'error p','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end

    % plot time to solution
    fh102 = figure(102);
    DOFS = (NN+2).*(NN+2);% + 2.*(NN).*(NN);
    dofs = (N +2).*(N +2);% + 2.*(N ).*(N );
    p5 = loglog(dofs,soltime,'r+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('\# dofs [1]','Interpreter','latex','FontSize',16)
    ylabel('time to solution [s]','Interpreter','latex','FontSize',16)
    title('Scaling of direct solver','Interpreter','latex','FontSize',20)

    if N == NN(1)
        p6 = loglog(DOFS,soltime*(DOFS./DOFS(1)).^1,'k-','LineWidth',2);  % plot linear trend for comparison
    end
    if N == NN(end)
        legend([p5,p6],{'time to solution','linear'},'Interpreter','latex','box','on','location','southeast')
    end

    drawnow
end

name = [outdir,'/',runID,'/',runID,'_bnchm'];
print(fh101,name,'-dpng','-r300','-vector');

name = [outdir,'/',runID,'/',runID,'_scale'];
print(fh102,name,'-dpng','-r300','-vector');
