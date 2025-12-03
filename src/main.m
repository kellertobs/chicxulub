%% CHIXCULUB: IMPACT HYDROTHERMAL CIRCULATION MODEL

fprintf('\n\n')
fprintf('*************************************************************\n');
fprintf('*****  RUN CHICXULUB MODEL | %s  **********\n',datetime('now'));
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

init;



%% TIME STEPPING LOOP      


while time <= tend && step <= Nt || max(Ra(:))<100

        
    tic;

    fprintf(1,'\n\n*****  step = %d,  time = %1.3e [yr],  step = %1.3e [hr] \n\n',step,time/3600/24/365.25,dt/3600);

    % store previous rates of change
    To    = T;
    Co    = C;
    Vo    = V;
    mqo   = mq;
    dTdto = dTdt;
    dCdto = dCdt;
    dVdto = dVdt;

    resnorm = 1e6;
    pi    = p;
    it    = 0;


    %% NONLINEAR SOLVER LOOP

    while resnorm >= tol && it <= maxit || it <= 100 
        
        if ~mod(it,nup) && step>0

            % UPDATE TEMPERATURE SOLUTION (SEMI-IMPLICIT UPDATE)
            
            advn_T = - advect(T,u,w,h,{ADVN,'vdf'},[1,2],BC_T);

            diff_T = diffus(T,kT,h,[1,2],BC_T);

            eqlb_T = -(T_wat-T_air)./tau_eqlb.*wat;

            mq     = max(0,min(1,((T-Tsol)./(Tliq-Tsol)))).^pTm;
            phsr_m = -(mq - mqo)/5/dt;

            dTdt   = advn_T + diff_T + eqlb_T - phsr_V*LHv - phsr_m*LHm;

            if bnchm; dTdt = dTdt + src_T_mms(:,:,step+1); end

            if wat_evolve
                dTdt_wat = sum(wat(:).*dTdt(:))./sum(wat(:));
            else
                dTdt_wat = 0;
            end
            dTdt = (1-wat-air).*dTdt + wat.*dTdt_wat + air*0;

            res_T = (T-To)/dt - (dTdt + dTdto)/2;

            [T,XHST.T,RHST.T,rho_est.T,rho_mean.T] = iterate(T,res_T*dt,rho_est.T,rho_mean.T,XHST.T,RHST.T,itpar,step*it);


            % set water to evolving reservoir
            if wat_evolve
                T_wat = mean(T(wat==1),'all');
            end

            % UPDATE CONCENTRATION SOLUTION (SEMI-IMPLICIT UPDATE)
            
            % calculate salinity advection
            advn_C = - advect(C,u,w,h,{ADVN,'vdf'},[1,2],BC_C);

            diff_C = diffus(C,kC,h,[1,2],BC_C);

            dCdt   = advn_C + diff_C;
            
            if bnchm; dCdt = dCdt + src_C_mms(:,:,step+1); end

            if wat_evolve
                dCdt_wat = sum(wat(:).*dCdt(:))./sum(wat(:));
            else
                dCdt_wat = 0;
            end
            dCdt = (1-wat-air).*dCdt + wat.*dCdt_wat + air*0;

            res_C = (C-Co)/(dt+TINY) - (dCdt + dCdto)/2;

            [C,XHST.C,RHST.C,rho_est.C,rho_mean.C] = iterate(C,res_C*dt,rho_est.C,rho_mean.C,XHST.C,RHST.C,itpar,step*it);

            C = max(0,min(1,C));  % saveguard min/max bounds

            % set water to evolving reservoir
            if wat_evolve
                C_wat = mean(C(wat==1),'all');
            end

            % UPDATE VAPOUR SOLUTION (SEMI-IMPLICIT UPDATE)
            
            % calculate vapour advection
            advn_V = - advect(V,u,w,h,{ADVN,'vdf'},[1,2],BC_V);

            diff_V = diffus(V,kV,h,[1,2],BC_V);

            Vq = vapour(T,C,Plith);

            phsr_V = -(V - Vq)/5/dt;

            dVdt = advn_V + diff_V + phsr_V;

            if bnchm; dVdt = dVdt + src_V_mms(:,:,step+1); end

            if wat_evolve
                dVdt_wat = sum(wat(:).*dVdt(:))./sum(wat(:));
            else
                dVdt_wat = 0;
            end
            dVdt = (1-wat-air).*dVdt + wat.*dVdt_wat + air*0;

            res_V = (V-Vo)/(dt+TINY) - (dVdt + dVdto)/2;

            [V,XHST.V,RHST.V,rho_est.V,rho_mean.V] = iterate(V,res_V*dt,rho_est.V,rho_mean.V,XHST.V,RHST.V,itpar,step*it);

            V = max(0,min(1,V));  % saveguard min/max bounds
 
            % set water to evolving reservoir
            if wat_evolve
                V_wat = mean(V(wat==1),'all');
            end

            % update density difference
            rho   = rhol0.*(1 - aT.*T - aC.*C - aV.*V );
            Drho  = rho - mean(rho,2);
            Drhoz = (Drho(icz(1:end-1),:)+Drho(icz(2:end),:))./2;
            if bnchm; Drhoz = Drho_mms(:,:,step+1); end

            % update permeability
            kB = k0 * f.^n;                                                            % Kozeny-Carman relation
            k  = 10.^(log10(kD) + (log10(kB)-log10(kD)) .* 1./(1+exp(-(BDT-T)./100))); % reduce permeability in ductile region

            % update Darcy coefficient
            K  = k/mu;
            Kz = (K(icz(1:end-1),:)+K(icz(2:end),:))./2;
            Kx = (K(:,icx(1:end-1))+K(:,icx(2:end)))./2;

            % get iterative step size for p-solution
            dtau = (h/2)^2./K;
        end

        % UPDATE VELOCITY-PRESSURE SOLUTION (PSEUDO-TRANSIENT SOLVER)

        % store previous iterative solution guesses
        pii = pi; pi = p;
        
        % calculate pressure gradient [Pa/m]
        gradPz = diff(p(:,2:end-1),1,1)./h;  % vertical gradient
        gradPx = diff(p(2:end-1,:),1,2)./h;  % horizontal gradient

        % calculate Darcy segregation speed [m/s]
        w = - Kz .* (gradPz - Drhoz.*grav);% .* (watfz==0);
        u = - Kx .* (gradPx              );% .* (watfx==0);
        
        w(2:end,:) = w(2:end,:).*(1-air);  % set air to zero flow

        % calculate residual of pressure equation
        res_p = diff(w,1,1)./h + diff(u,1,2)./h;
        if bnchm; res_p = res_p - src_p_mms(:,:,step+1); end

        res_p(air>=0.5) = 0;  % set air to zero

        % update pressure solution
        [p(2:end-1,2:end-1),XHST.p,RHST.p,rho_est.p,rho_mean.p] = iterate(p(2:end-1,2:end-1),res_p*dt,rho_est.p,rho_mean.p,XHST.p,RHST.p,itpar,step*it);
        
        % apply pressure boundary conditions
        if strcmp(BC_VP{1},'closed')
            p([1 end],:) = p([2 end-1],:) + [-1;1].*Drhoz([1 end],icx).*grav.*h;
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
            resnorm = norm(RHST.p(:,end))./norm(p+eps);% ... 
                    % + norm(upd_T,2)./norm(T+1,2) ...
                    % + norm(upd_C,2)./norm(C+1,2) ...
                    % + norm(upd_V,2)./norm(V+1,2);

            % report convergence
            fprintf(1,'---  %d,  %e\n',it,resnorm);
        end
        
        % increment iteration count
        it = it+1;

    end  % non-linear iteration loop

    fprintf(1,'\n\n      mean water T = %2.4f [C] \n',T_wat);
    fprintf(1,'      mean water C = %2.4f [wt] \n',C_wat);

    soltime = toc;  % record time to solution
    fprintf(1,'\n\n      time to solution = %1.3e [s] \n\n',soltime);

    % update dimensionless numbers
    Vel = sqrt(((w(1:end-1,:)+w(2:end,:))/2).^2 ...
             + ((u(:,1:end-1)+u(:,2:end))/2).^2);
    Ra = Vel.*D./(kT+kC+kV); 


%%  PLOT AND SAVE SOLUTION
    output;
    
    % update time and step count
    step = step + 1;
    time = time + dt;  

    
end  % time-stepping loop
