% UPDATE TEMPERATURE SOLUTION (SEMI-IMPLICIT UPDATE)

advn_T = - advect(T,u,w,h,{ADVN,'vdf'},[1,2],BC_T);

diff_T = diffus(T,kT+kd,h,[1,2],BC_T);

eqlb_T = -(T_wat-T_air)./tau_eqlb.*wat;

mq     = max(0,min(1,((T-Tsol)./(Tliq-Tsol)))).^pTm;
mlt    = mq>0.6;
phsr_m = -(mq - mqo)/5/dt;

dTdt   = advn_T + diff_T + eqlb_T - phsr_V*LHv - phsr_m*LHm;

if bnchm; dTdt = dTdt + src_T_mms(:,:,step+1); end

if any(mlt(:))
    dTdt_mlt = sum(mlt(:).*dTdt(:))./sum(mlt(:)+eps);
else
    dTdt_mlt = 0;
end
if wat_evolve
    dTdt_wat = sum(wat(:).*dTdt(:))./sum(wat(:));
else
    dTdt_wat = 0;
end

dTdt = (1-mlt-wat-air).*dTdt + mlt.*dTdt_mlt + wat.*dTdt_wat + air*0;

res_T = (T-To)/dt - (dTdt + dTdto)/2;

[T,GHST.T,FHST.T,specrad.T] = iterate(T,res_T*dt,specrad.T,GHST.T,FHST.T,itpar,it/nup);


% set water to evolving reservoir
if any(mlt(:))
    T_mlt = mean(T(mlt==1),'all');
end
if wat_evolve
    T_wat = mean(T(wat==1),'all');
end

% UPDATE CONCENTRATION SOLUTION (SEMI-IMPLICIT UPDATE)

% calculate salinity advection
advn_C = - advect(C,u,w,h,{ADVN,'vdf'},[1,2],BC_C);

diff_C = diffus(C,kC+kd,h,[1,2],BC_C);

dCdt   = advn_C + diff_C;

if bnchm; dCdt = dCdt + src_C_mms(:,:,step+1); end

if any(mlt(:))
    dCdt_mlt = sum(mlt(:).*dCdt(:))./sum(mlt(:)+eps);
else
    dCdt_mlt = 0;
end
if wat_evolve
    dCdt_wat = sum(wat(:).*dCdt(:))./sum(wat(:));
else
    dCdt_wat = 0;
end

dCdt = (1-mlt-wat-air).*dCdt + mlt.*dCdt_mlt+ wat.*dCdt_wat + air*0;

res_C = (C-Co)/(dt+TINY) - (dCdt + dCdto)/2;

[C,GHST.C,FHST.C,specrad.C] = iterate(C,res_C*dt,specrad.C,GHST.C,FHST.C,itpar,it/nup);

C = max(0,min(1,C));  % saveguard min/max bounds

if any(mlt(:))
    C_mlt = mean(C(mlt==1),'all');
end
if wat_evolve
    C_wat = mean(C(wat==1),'all');
end

% UPDATE VAPOUR SOLUTION (SEMI-IMPLICIT UPDATE)

% calculate vapour advection
advn_V = - advect(V,u,w,h,{ADVN,'vdf'},[1,2],BC_V);

diff_V = diffus(V,kV+kd,h,[1,2],BC_V);

Vq = vapour(T,C,Plith);

phsr_V = -(V - Vq)/5/dt;

dVdt = advn_V + diff_V + phsr_V;

if bnchm; dVdt = dVdt + src_V_mms(:,:,step+1); end

if any(mlt(:))
    dVdt_mlt = sum(mlt(:).*dVdt(:))./sum(mlt(:)+eps);
else
    dVdt_mlt = 0;
end
if wat_evolve
    dVdt_wat = sum(wat(:).*dVdt(:))./sum(wat(:));
else
    dVdt_wat = 0;
end

dVdt = (1-mlt-wat-air).*dVdt + mlt.*dVdt_mlt + wat.*dVdt_wat + air*0;

res_V = (V-Vo)/(dt+TINY) - (dVdt + dVdto)/2;

[V,GHST.V,FHST.V,specrad.V] = iterate(V,res_V*dt,specrad.V,GHST.V,FHST.V,itpar,it/nup);

V = max(0,min(1,V));  % saveguard min/max bounds

% set water to evolving reservoir
if wat_evolve
    V_wat = mean(V(wat==1),'all');
end