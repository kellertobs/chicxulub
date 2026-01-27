% Fixed-point iterative update with Anderson acceleration

% res      : preconditioned residual of governing equation
% x        : current iterate
% g        : fixed-point updated iterate
% f        : current fixed-point update f = g - x
% x_acc    : fully accelerated iterate
% x_new    : new updated iterate with (some/no) acceleration applied
% rho.est  : estimated spectral radius at each grid point
% rho.mean : global mean of estimated spectral radius
% FHST     : history of previous fixed-point updates
% itpar... : iterative parameter structure 
%      .fp.damp  : damping coefficient for fixed-point coefficients
%      .aa.m     : depth for Anderson acceleration
%      .aa.damp  : damping coefficient for Anderson acceleration
%      .aa.reg   : regularisation coefficient for Anderson acceleration
% count    : iteration count to track history of run

function [x,GHST,FHST,rho] = iterate(x,res,rho,GHST,FHST,itpar,count)

% allocate arrays of correct shape
alpha = 0.*x;
beta  = 0.*x;
x_new = 0.*x;
x_acc = 0.*x;
f     = 0.*x;


% 1) Inertial fixed-point iterative update

% Per-DOF spectral radius rho estimates from ratio of consecutive updates
if count>2
    ratio   = abs(FHST(:,end))./abs(FHST(:,end-1) + eps);  % form ratio of two most recent updates
    rho_new = min(0.9, max(0.1, ratio));                   % clamp values to desired range
else
    rho_new = rho.mean;
end

% Moving average for stability
rho.est  = 0.7*rho.est  + 0.3*rho_new;           % moving average
rho.mean = 0.9*rho.mean + 0.1*mean(rho.est);     % moving average
rho.est  = max(rho.est, 0.5*rho.mean);           % avoid outliers
rho.est  = max(rho.est, 0.3);                    % hard lower bound

% Fixed-point update coefficients
alpha(:) = 4./(2 + rho.est).^2;
% beta (:) = alpha(:)./2;

% New fixed-point update and iterate
f(:) = itpar.fp.damp*(-alpha(:).*res(:) + beta(:).*FHST(:,end));
g    = x + f;


% 2) Anderson acceleration (Walker & Ni, 2011)

% Shift histories and store current g, f
FHST = [FHST(:,2:end) f(:)];   % previous solution updates
GHST = [GHST(:,2:end) g(:)];   % previous solution updates

if count>2 && itpar.aa.damp>eps  % only if enough history and damp>0

    n = max(0,itpar.aa.m-count+1);

    % Take differences of updates and solutions
    DF  = FHST(:,2+n:end) - FHST(:,1+n:end-1);   % history of fixed-point updates
    DG  = GHST(:,2+n:end) - GHST(:,1+n:end-1);   % history of fixed-point iterates

    % Solve min_delta || f - DF * delta + reg*I ||_2  (global regularised least squares)
    reg   = itpar.aa.reg.*rms(DF.'*DF,'all');      % get scaled regularisation level
    gamma = (DF.'*DF + reg*eye(itpar.aa.m-n)) \ (DF.'*f(:));  % solve least squares problem

    % Anderson update for fixed-point:
    x_acc(:) = g(:) - DG * gamma;

    % Damped Anderson step
    x_new(:) = g(:) + itpar.aa.damp * (x_acc(:) - g(:));

else
    % No acceleration applied
    x_new = g;
end


% 3) Update the current iterate with the new fixed-point/accelerated value
x = x_new;

end