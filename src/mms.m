% prepare workspace
clear x z t f T C p w u f_mms T_mms C_mms p_mms pi Drho_mms
TINY = 1e-16;
syms f_mms(x,z,t) T_mms(x,z,t) C_mms(x,z,t) p_mms(x,z,t)

fprintf(1,'\n\n  ***  compose manufactured solution\n\n');

% compose manufactured solution variables
f_mms(x,z,t) = 1e-1 .* (1 + (sin(2*(x)*pi/D).*cos(2*(z)*pi/D))/4);
T_mms(x,z,t) = 1e+2 .* (1 + (cos(2*(x)*pi/D).*sin(2*(z)*pi/D)).*exp(-t/1e7)/2);
C_mms(x,z,t) = 5e-2 .* (1 + (sin(2*(x)*pi/D).*sin(2*(z)*pi/D)).*exp(-t/1e7)/2);
p_mms(x,z,t) = 1e+4 .* (  -  cos(2*(x)*pi/D).*cos(2*(z)*pi/D));

fprintf(1,'       f   = %s \n',char(f_mms));
fprintf(1,'       T   = %s \n',char(T_mms));
fprintf(1,'       C   = %s \n',char(C_mms));
fprintf(1,'       P   = %s \n',char(p_mms));
fprintf(1,'       . ');

% update coefficients
k_mms    = k0 * f_mms.^n;
K_mms    = k_mms/mu;
Drho_mms = rhol0.*(- aT.*(T_mms-100) + gC.*(C_mms-0.05));
fprintf(1,' . ');

% update porous flow velocities
w_mms = - K_mms .* (diff(p_mms,z) - Drho_mms.*grav);
u_mms = - K_mms .* (diff(p_mms,x)                 );
fprintf(1,' . ');

% manufactured solution residuals
res_p_mms = (diff(w_mms,z) + diff(u_mms,x));
res_T_mms = diff(T_mms,t) + u_mms * diff(T_mms,x) + w_mms * diff(T_mms,z) - kT * (diff(T_mms,x,2) + diff(T_mms,z,2));
res_C_mms = diff(C_mms,t) + u_mms * diff(C_mms,x) + w_mms * diff(C_mms,z) - kC * (diff(C_mms,x,2) + diff(C_mms,z,2));
fprintf(1,' . ');

% plot manufactured solution
figure(15);
colormap(ocean);
subplot(2,3,1); fcontour( -w_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufcat. $w$ [m/s]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,2); fcontour(  u_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $u$ [m/s]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,3); fcontour(  p_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $p$ [Pa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
fprintf(1,' . ');
subplot(2,3,4); fcontour( f_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $f$ [vol]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,5); fcontour( T_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $T$ [degC]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,6); fcontour( C_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $C$ [wt]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
drawnow;
fprintf(1,' . \n');

% plot manufactured residuals
figure(16);
colormap(ocean);
subplot(2,3,1); fcontour( res_p_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $p$-res','Interpreter','latex');
subplot(2,3,2); fcontour( res_T_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $T$-res','Interpreter','latex');
subplot(2,3,3); fcontour( res_C_mms(x,z,0) ,[0,D],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $C$-res','Interpreter','latex');
drawnow;

% evaluate mms source terms on appropriate grids
fprintf(1,'\n  ***  evaluate manufactured solution\n\n');
t_mms  =  0:dt_mms:Nt*dt_mms;
x_mms  = -h/2:h:D+h/2;
z_mms  = -h/2:h:D+h/2;
xu_mms = (x_mms(1:end-1)+x_mms(2:end))./2;
zw_mms = (z_mms(1:end-1)+z_mms(2:end))./2;

fprintf(1,'       Patience, my young Padawan!\n');
fprintf(1,'       . ');

[x,z,t] = meshgrid(x_mms(2:end-1),z_mms(2:end-1),t_mms);
src_p_mms = double(subs(res_p_mms)); fprintf(1,' . ');
src_T_mms = double(subs(res_T_mms)); fprintf(1,' . ');
src_C_mms = double(subs(res_C_mms)); fprintf(1,' . ');

% plot evaluated source terms
figure(16);
subplot(2,3,4); imagesc(x_mms(2:end-1),z_mms(2:end-1),src_p_mms(:,:,1)); axis ij equal tight; colorbar; box on; title('evaluated $p$-res','Interpreter','latex');
subplot(2,3,5); imagesc(x_mms(2:end-1),z_mms(2:end-1),src_T_mms(:,:,1)); axis ij equal tight; colorbar; box on; title('evaluated $T$-res','Interpreter','latex');
subplot(2,3,6); imagesc(x_mms(2:end-1),z_mms(2:end-1),src_C_mms(:,:,1)); axis ij equal tight; colorbar; box on; title('evaluated $C$-res','Interpreter','latex');
drawnow;

% evaluate analytical solution on appropriate grids
[x,z,t]  = meshgrid(xu_mms,z_mms,t_mms);
u_mms    = double(subs(u_mms)); fprintf(1,' . ');
[x,z,t]  = meshgrid(x_mms,zw_mms,t_mms);
w_mms    = double(subs(w_mms));
Drho_mms = double(subs(Drho_mms)); fprintf(1,' . ');
[x,z,t]  = meshgrid(x_mms,z_mms,t_mms);
p_mms    = double(subs(p_mms)); fprintf(1,' . ');
[x,z,t]  = meshgrid(x_mms(2:end-1),z_mms(2:end-1),t_mms);
f_mms    = double(subs(f_mms));
T_mms    = double(subs(T_mms));
C_mms    = double(subs(C_mms)); fprintf(1,' . ');

% set initial condition
f = f_mms(:,:,1);
T = T_mms(:,:,1);
C = C_mms(:,:,1);

% reset coordinate vectors for numerical model
x = linspace(h/2,D-h/2,N);
z = linspace(h/2,D-h/2,N);
t = 0;

fprintf(1,' . \n');
