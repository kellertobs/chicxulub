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
            
            advn_T = - advect(T,u(2:end-1,:),w(:,2:end-1),h,{ADVN,'vdf'},[1,2],BC_T);

            diff_T = diffus(T,kT,h,[1,2],BC_T);

            eqlb_T = -(T_wat-T_air)./tau_eqlb.*wat;

            dTdt = advn_T + diff_T + eqlb_T - phsr_V*LH;

            if bnchm; dTdt = dTdt + src_T_mms(:,:,step+1); end

            if wat_evolve
                dTdt_wat = sum(wat(:).*dTdt(:))./sum(wat(:));
            else
                dTdt_wat = 0;
            end
            dTdt = (1-wat-air).*dTdt + wat.*dTdt_wat + air*0;

            res_T = (T-To)/dt - (dTdt + dTdto)/2;

            upd_T = - gamma*res_T*dt + delta*upd_T;
            T = T + upd_T;

            % set water to evolving reservoir
            if wat_evolve
                T_wat = mean(T(wat==1),'all');
            end

            % UPDATE CONCENTRATION SOLUTION (SEMI-IMPLICIT UPDATE)
            
            % calculate salinity advection
            advn_C = - advect(C,u(2:end-1,:),w(:,2:end-1),h,{ADVN,'vdf'},[1,2],BC_C);

            diff_C = diffus(C,kC,h,[1,2],BC_C);

            dCdt = advn_C + diff_C;
            
            if bnchm; dCdt = dCdt + src_C_mms(:,:,step+1); end

            if wat_evolve
                dCdt_wat = sum(wat(:).*dCdt(:))./sum(wat(:));
            else
                dCdt_wat = 0;
            end
            dCdt = (1-wat-air).*dCdt + wat.*dCdt_wat + air*0;

            res_C = (C-Co)/(dt+TINY) - (dCdt + dCdto)/2;

            upd_C = - gamma*res_C*dt + delta*upd_C;
            C = C + upd_C;

            C = max(0,min(1,C));  % saveguard min/max bounds

            % set water to evolving reservoir
            if wat_evolve
                C_wat = mean(C(wat==1),'all');
            end

            % UPDATE VAPOUR SOLUTION (SEMI-IMPLICIT UPDATE)
            
            % calculate vapour advection
            advn_V = - advect(V,u(2:end-1,:),w(:,2:end-1),h,{ADVN,'vdf'},[1,2],BC_V);

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

            upd_V = - gamma*res_V*dt + delta*upd_V;
            V = V + upd_V;

            V = max(0,min(1,V));  % saveguard min/max bounds
 
            % set water to evolving reservoir
            if wat_evolve
                V_wat = mean(V(wat==1),'all');
            end

            % update density difference
            rho  = rhol0.*(1 - aT.*T(icz,icx) ...
                             - aC.*C(icz,icx) ...
                             - aV.*V(icz,icx) );% .* (1-air(icz,icx));
            Drho = rho - mean(rho,2);
            % Drho(air(icz,icx)+wat(icz,icx)>=0.5) = 0;  % set air and water to zero
            Drhoz = (Drho(1:end-1,:)+Drho(2:end,:))./2;
            if bnchm; Drhoz = Drho_mms(:,:,step+1); end
        end

        % UPDATE VELOCITY-PRESSURE SOLUTION (PSEUDO-TRANSIENT SOLVER)

        % store previous iterative solution guesses
        pii = pi; pi = p;
        
        % calculate pressure gradient [Pa/m]
        gradPz = diff(p,1,1)./h;  % vertical gradient
        gradPx = diff(p,1,2)./h;  % horizontal gradient

        % calculate Darcy segregation speed [m/s]
        w = - Kz .* (gradPz - Drhoz.*grav);% .* (watfz==0);
        u = - Kx .* (gradPx              );% .* (watfx==0);
        
        w(2:end,:) = w(2:end,:).*(1-air(:,icx));  % set air to zero flow

        % calculate residual of pressure equation
        res_p = diff(w(:,2:end-1),1,1)./h + diff(u(2:end-1,:),1,2)./h;
        if bnchm; res_p = res_p - src_p_mms(:,:,step+1); end

        res_p(air>=0.5) = 0;  % set air to zero

        % update pressure solution
        upd_p = - alpha*res_p.*dtau + beta*upd_p;
        p(2:end-1,2:end-1) = p(2:end-1,2:end-1) + upd_p;
        
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
            resnorm = norm(res_p.*dtau.*(1-air-wat),2)./norm(p+1,2) ... 
                    + norm(res_T.*dt*gamma.*(1-air-wat),2)./norm(T+1,2) ...
                    + norm(res_C.*dt*gamma.*(1-air-wat),2)./norm(C+1,2) ...
                    + norm(res_V.*dt*gamma.*(1-air-wat),2)./norm(V+1,2);

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
    Vel = sqrt(((w(1:end-1,2:end-1)+w(2:end,2:end-1))/2).^2 ...
             + ((u(2:end-1,1:end-1)+u(2:end-1,2:end))/2).^2);
    Ra = Vel.*D./(kT+kC+kV); 


%%  PLOT AND SAVE SOLUTION
    output;
    
    % update time and step count
    step = step + 1;
    time = time + dt;  

    
end  % time-stepping loop
