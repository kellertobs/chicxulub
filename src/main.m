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

    % reset update history
    specrad.T.est = specrad.T.mean + 0.*specrad.T.est;
    specrad.C.est = specrad.C.mean + 0.*specrad.C.est;
    specrad.V.est = specrad.V.mean + 0.*specrad.V.est;
    specrad.p.est = specrad.p.mean + 0.*specrad.p.est;

    GHST.T        = 0.*GHST.T;
    GHST.C        = 0.*GHST.C;
    GHST.V        = 0.*GHST.V;
    GHST.p        = 0.*GHST.p;

    FHST.T        = 0.*FHST.T;
    FHST.C        = 0.*FHST.C;
    FHST.V        = 0.*FHST.V;
    FHST.p        = 0.*FHST.p;

    resnorm = 1e6;
    pi    = p;
    it    = 0;


    %% NONLINEAR SOLVER LOOP

    while resnorm >= tol && it <= maxit
        
        if ~mod(it,nup) && step>0

            % update thermo-chemical evolution
            thermochem;

            % update model parameter
            update;

        end

        % UPDATE VELOCITY-PRESSURE SOLUTION (PSEUDO-TRANSIENT SOLVER)
        switch solver

            case 'direct'

                darcy_direct;

            case 'pseudotransient'

                darcy_pseudotrans;

        end

        if ~mod(it,nup)
            % get preconditioned residual norm to monitor convergence
            resnorm = norm(FHST.p(:,end))./norm(p+1) ... 
                    + norm(FHST.T(:,end))./norm(T+1) ...
                    + norm(FHST.C(:,end))./norm(C+1) ...
                    + norm(FHST.V(:,end))./norm(V+1);

            % report convergence
            fprintf(1,'---  %d,  %e\n',it,resnorm);
        end
        
        % increment iteration count
        it = it+1;

    end  % non-linear iteration loop

    fprintf(1,'\n\n      mean water T = %2.4f [C] \n',T_wat);
    fprintf(1,'      mean water C = %2.4f [wt] \n',C_wat);
    if any(mlt(:))
    fprintf(1,'      mean magma T = %2.4f [C] \n',T_mlt);
    fprintf(1,'      mean magma C = %2.4f [wt] \n',C_mlt);
    end

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
