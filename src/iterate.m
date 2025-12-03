% Chebyshev-like fixed-point iterative update with Anderson acceleration

% res      : preconditioned residual of governing equation
% x        : current iterate
% g        : fixed-point updated iterate
% r        : current residual (fixed-point update) r = g - x
% x_acc    : fully accelerated iterate
% x_new    : new updated iterate with (some/no) acceleration applied
% XHST     : history of (unaccelerated) previous iterates
% RHST     : history of previous residuals
% itpar... : iterative parameter structure 
%      .cheb.alpha : damping for first Chebychev-like coefficient
%      .cheb.beta  : damping for second Chebychev-like coefficient
%      .cheb.gamma : damping for third Chebychev-like coefficient 
%      .anda.m     : depth for Anderson acceleration
%      .anda.mix   : mixing coefficient for Anderson acceleration
%      .anda.reg   : regularisation coefficient for Anderson acceleration
% count    : iter*step count to track history of run

function [x,XHST,RHST,rho_est,rho_mean] = iterate(x,res,rho_est,rho_mean,XHST,RHST,itpar,count)

% allocate arrays of correct shape
alpha = 0.*x;
beta  = 0.*x;
gamma = 0.*x;
x_new = 0.*x;
x_acc = 0.*x;
r     = 0.*x;
g     = 0.*x;

% 1) Chebyshev-like fixed-point iterative update

% Update per-DOF rho estimates from ratio of consecutive updates
denom   = max(abs(RHST(:,end-1)), 1e-12);         % avoid blow-up
ratio   = abs(RHST(:,end)) ./ denom;            % n×1
rho_new = min(0.99, max(0.02, ratio));           % clamp to [0.02, 0.99]

% Moving average for stability
rho_est  = 0.7*rho_est  + 0.3*rho_new;           % moving average
rho_mean = 0.9*rho_mean + 0.1*mean(rho_est);     % moving average
rho_est  = max(rho_est, 0.7*rho_mean);           % avoid crazy outliers
rho_est  = max(rho_est, 0.1);                    % global lower bound

% Chebyshev-like coefficients
alpha(:) =  4./(2 + rho_est).^2              .*itpar.cheb.alpha;
beta (:) =  2.*(2 - rho_est) ./ (2 + rho_est).*itpar.cheb.beta;
gamma(:) = -1./(2 + rho_est).^2              .*itpar.cheb.gamma;

% New residual/update and fixed-point iterate
r(:) = -alpha(:).*res(:) + beta(:).*RHST(:,end) + gamma(:).*RHST(:,end-1);   % n×1
g    = x + r;

% 3) Anderson acceleration (residual-based)

% Shift histories and store current g, r
XHST = [XHST(:,2:end) g(:)];
RHST = [RHST(:,2:end) r(:)];

if (count>itpar.anda.m || ~count) && itpar.anda.mix>eps  % only if enough history and mix>0

    % Differences of iterates and residuals
    DX  = XHST(:,2:end) - XHST(:,1:end-1);   % n×m
    DR  = RHST(:,2:end) - RHST(:,1:end-1);   % n×m
    reg = itpar.anda.reg.*rms(DR.'*DR,'all');

    % Solve min_gamma || f - DF * gamma + reg*I ||_2  (global regularised least squares)
    gamma = (DR.'*DR + reg*eye(itpar.anda.m)) \ (DR.'*r(:));

    % Standard Anderson Type-I update for fixed-point:
    x_acc(:) = g(:) - DX * gamma;

    % Damped Anderson step
    x_new(:) = g(:) + itpar.anda.mix * (x_acc(:) - g(:));

else
    % No acceleration applied
    x_new = g;
end

% Update the current iterate with the new fixed-point/accelerated value
x = x_new;

end