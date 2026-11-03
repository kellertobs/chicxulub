% update density difference
rho0  = (1-V).*rhol0 + V.*rhov0;
rho   = rho0.*(1 - aT.*T - aC.*C);
Drho  = rho - mean(rho,2);
Drhoz = (Drho(icz(1:end-1),:)+Drho(icz(2:end),:))./2;
if strcmp(BC_VP{1},'closed'); Drhoz([1 end],:) = 0; end
if bnchm; Drhoz = Drho_mms(:,:,step+1); end

% update permeability
kB = k0 * f.^n;                                                            % Kozeny-Carman relation
k  = 10.^(log10(kD) + (log10(kB)-log10(kD)) .* 1./(1+exp(-(BDT-T)./50)));  % reduce permeability in ductile region

% update fluid viscosity
mu = mul.^(1-V) .* muv.^V;

% update pore velocity magnitude
vel = sqrt(((w(1:end-1,:)+w(2:end,:))/2).^2 ...
         + ((u(:,1:end-1)+u(:,2:end))/2).^2) ./ f;

% update hydrodynamic dispersivity
kd = vel.*ld .* (1-air-wat);

% update Darcy coefficient
K  = k./mu;

switch Kintp
    case 'geometric'
        Kz = (K(icz(1:end-1),:).*K(icz(2:end),:)).^0.5;
        Kx = (K(:,icx(1:end-1)).*K(:,icx(2:end))).^0.5;
    case 'arithmetic'
        Kz = (K(icz(1:end-1),:)+K(icz(2:end),:)).*0.5;
        Kx = (K(:,icx(1:end-1))+K(:,icx(2:end))).*0.5;
end

% get iterative step size for p-solution
dtau = (h/2)^2./K;

% get physical time step
dt = CFL * min([(h/2)/max(abs(w(:))) , (h/2)/max(abs(u(:))) , (h/2)^2./max(kT+kd(:))]);
if bnchm; dt = dt_mms; end
