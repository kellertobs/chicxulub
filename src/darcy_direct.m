IIL = [];       % equation indeces into A
JJL = [];       % variable indeces into A
AAL = [];       % coefficients for A

% top boundary
ii  = MapP(1,:);
jj1 = ii; jj2 = MapP(2,:);

aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+topbot];

% bottom boundary
ii  = MapP(end  ,:);
jj1 = ii; jj2 = MapP(end-1,:);

aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+topbot];

% left & right boundaries
ii  = [MapP(2:end-1,1    ); MapP(2:end-1,end)]; % left & right
jj1 = ii;
if periodic
    jj2 = [MapP(2:end-1,end-1); MapP(2:end-1,2    )];
else
    jj2 = [MapP(2:end-1,2    ); MapP(2:end-1,end-1)];
end

aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sides];

% internal points
ii  = MapP(2:end-1,2:end-1);
jj1 = MapP(1:end-2,2:end-1);
jj2 = MapP(3:end-0,2:end-1);
jj3 = MapP(2:end-1,1:end-2);
jj4 = MapP(2:end-1,3:end-0);

% coefficients multiplying pressure
switch Kintp
    case 'geometric'
        K1 = (K(icz(1:end-2),:).*K(icz(2:end-1),:)).^0.5;   K2 = (K(icz(2:end-1),:).*K(icz(3:end-0),:)).^0.5;
        K3 = (K(:,icx(1:end-2)).*K(:,icx(2:end-1))).^0.5;   K4 = (K(:,icx(2:end-1)).*K(:,icx(3:end-0))).^0.5;
    case 'arithmetic'
        K1 = (K(icz(1:end-2),:)+K(icz(2:end-1),:)).*0.5;   K2 = (K(icz(2:end-1),:)+K(icz(3:end-0),:)).*0.5;
        K3 = (K(:,icx(1:end-2))+K(:,icx(2:end-1))).*0.5;   K4 = (K(:,icx(2:end-1))+K(:,icx(3:end-0))).*0.5;
end

aa  = (K1+K2+K3+K4)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL;-aa(:)    ];      % P on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; K1(:)/h^2];      % P one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; K2(:)/h^2];      % P one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; K3(:)/h^2];      % P one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; K4(:)/h^2];      % P one below

% assemble coefficient matrix
LL  = sparse(IIL,JJL,AAL,NP,NP);

% Right hand side
IIR = [];
JJR = ones(Nz*Nx,1);
AAR = [];

ii  = MapP(2:end-1,2:end-1);

rr  = (K2.*Drhoz(2:end,:)-K1.*Drhoz(1:end-1,:))./h.*grav;

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% assemble right-hand side vector
RR  = sparse(IIR,JJR,AAR,NP,1);

if strcmp(BC_VP{1},'closed')
    nzp = round(Nz/2);
    nxp = round(Nx/2);
    nnp = MapP(nzp,nxp);
    LL(nnp,:)   = 0;
    LL(nnp,nnp) = K(nzp,nxp)./h^2;
    RR(nnp,:)   = 0;
end

SCL = (abs(diag(LL))).^0.5;
SCL = diag(sparse( 1./SCL ));

LL  = SCL*LL*SCL;
RR  = SCL*RR;

if ~exist('pcol','var'); pcol = colamd(LL); end  % update column permutation for sparsity pattern once per time step
dLL          = decomposition(LL(:,pcol), 'lu');  % get LU-decomposition for consistent performance of LL \ RR
SS(pcol,1)   = dLL \ RR;                         % solve permuted decomposed system

SOL = SCL * (LL\RR);

p = full(reshape(SOL,Nz+2,Nx+2));

% calculate pressure gradient [Pa/m]
gradPz = diff(p(:,2:end-1),1,1)./h;  % vertical gradient
gradPx = diff(p(2:end-1,:),1,2)./h;  % horizontal gradient

% calculate Darcy segregation speed [m/s]
w = - Kz .* (gradPz - Drhoz.*grav);% .* (watfz==0);
u = - Kx .* (gradPx              );% .* (watfx==0);

% set air and water to zero
w = w.*(1-airw).*(1-watw);
u = u.*(1-airu).*(1-watu);
p = p.*(1-air([1,1:end,end],[1,1:end,end])).*(1-wat([1,1:end,end],[1,1:end,end]));  % set air to zero flow
