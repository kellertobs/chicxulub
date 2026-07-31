
% store previous iterative solution guesses
pii = pi; pi = p;

% calculate pressure gradient [Pa/m]
gradPz = diff(p(:,2:end-1),1,1)./h;  % vertical gradient
gradPx = diff(p(2:end-1,:),1,2)./h;  % horizontal gradient

% calculate Darcy segregation speed [m/s]
w = - Kz .* (gradPz - Drhoz.*grav);% .* (watfz==0);
u = - Kx .* (gradPx              );% .* (watfx==0);

w = w.*(1-airw).*(1-watw);  % set air, wat to zero flow
u = u.*(1-airu).*(1-watu);  % set air, wat to zero flow

% calculate residual of pressure equation
res_p = diff(w,1,1)./h + diff(u,1,2)./h;
if bnchm; res_p = res_p - src_p_mms(:,:,step+1); end

res_p(air>=0.5) = 0;  % set air to zero

% update pressure solution
[p(2:end-1,2:end-1),GHST.p,FHST.p,specrad.p] = iterate(p(2:end-1,2:end-1),res_p.*dtau,specrad.p,GHST.p,FHST.p,itpar,it);

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