%% MAE 185 Midterm 
clear; clc; close all;

%defining variables
L=10^-5;H=8*10^-6;nx=75;ny=80;
cp=1005;cv=718;R=cp-cv;
gamma=cp/cv;Pr=0.71;Tinf=288.15;pinf=101300;rho0=1.225;M=4;
a0=sqrt(gamma*R*Tinf); % speed of sound
uinf=M*a0;
t=0; 
dt=2.35e-11;

%creating grid
[x,y]=ndgrid(0:L/(nx-1):L,0:H/(ny-1):H);
dx=L/(nx-1);
dy=H/(ny-1);

%ICs
u=zeros(size(x));
v=zeros(size(x));
T=zeros(size(x));
p=zeros(size(x));

T(:,:)=Tinf;
p(:,:)=pinf;
u(:,:)=uinf;

mu=sutherland(T);
k=cp/Pr.*mu;
rho=p./R./T;

%Initialize U, Et, E, and F
U=prim2cons(rho,u,v,T,cv);
Et=U(:,:,4);
E=zeros(nx,ny,4);
F=zeros(nx,ny,4);

% while loop
count=0;

while count<1500
    %vp=max(4/3.*mu.*gamma.*mu./Pr./rho,[],'all');
    %dt=(abs(u)/dx + abs(v)/dy + sqrt(gamma.*R.*T.*(dx^-2 + dy^-2))+ 2*vp*(dx^-2 + dy^-2)).^-1;
    %dt=min(dt,[],'all');
    
    %predictor step, fwd FD for x and y

    %For E, tau gets backward in x, central in y
    
    tauxx=2*mu.*(ddx_bwd(u,dx)-(1/3*(ddx_bwd(u,dx)+ddy_central(v,dy))));
    tauxy=mu.*(ddy_central(u,dy)+ddx_bwd(v,dx));

    %For E, qdotx gets backward in x

    qdotx=-k.*ddx_bwd(T,dx);

    E(:,:,1)=rho.*u;
    E(:,:,2)=rho.*u.^2+p-tauxx;
    E(:,:,3)=rho.*u.*v-tauxy;
    E(:,:,4)=(Et+p).*u-u.*tauxx-v.*tauxy+qdotx;

    %For F, tau gets central in x, backwards in y

    tauyy=2*mu.*(ddy_bwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_bwd(v,dy))));
    tauxy=mu.*(ddy_bwd(u,dy)+ddx_central(v,dx));

    %For F, qdoty gets backward in y

    qdoty=-k.*ddy_bwd(T,dy);

    F(:,:,1)=rho.*v;
    F(:,:,2)=rho.*u.*v-tauxy;
    F(:,:,3)=rho.*v.^2+p-tauyy;
    F(:,:,4)=(Et+p).*v-v.*tauyy-u.*tauxy+qdoty;

    %the predictor step calculation, fwd and bwd difference functions account for 3D arrays E and F

    Ubar=U-dt*(ddx_fwd(E,dx)+ddy_fwd(F,dy));

    %now get primitive variables from Ubar so we can plug into Ebar, Fbar

    [~,u,v,T,p,~,~] = cons2prim(Ubar,R,cv);

    %ENFORCE BCs

    %at the outflow
    u(nx,:)=2*u(nx-1,:)-u(nx-2,:);
    v(nx,:)=2*v(nx-1,:)-v(nx-2,:);
    p(nx,:)=2*p(nx-1,:)-p(nx-2,:);
    T(nx,:)=2*T(nx-1,:)-T(nx-2,:);

    %at the inflow and far-field
    u(1,:)=uinf;
    u(:,ny)=uinf;
    p(1,:)=pinf;
    p(:,ny)=pinf;
    T(1,:)=Tinf;
    T(:,ny)=Tinf;
    v(1,:)=0;
    v(:,ny)=0;

    %at the wall
    u(:,1)=0;
    v(:,1)=0;
    T(:,1)=Tinf;
    p(:,1)=2*p(:,3)-p(:,2);

    %at the leading edge
    u(1,1)=0;
    p(1,1)=pinf;
    T(1,1)=Tinf;

    rho=p./R./T;
    mu=sutherland(T);
    k=cp/Pr.*mu;

    %now converting these boundary conditions to reflect in Ubar
    Ubar=prim2cons(rho,u,v,T,cv);
    Et=Ubar(:,:,4);

    %then the corrector step, backward fd in x and y

    %For Ebar, tau gets forward in x, central in y
    
    tauxx=2*mu.*(ddx_fwd(u,dx)-(1/3*(ddx_fwd(u,dx)+ddy_central(v,dy))));
    tauxy=mu.*(ddy_central(u,dy)+ddx_fwd(v,dx));

    %For Ebar, qdotx gets forward in x

    qdotx=-k.*ddx_fwd(T,dx);

    E(:,:,1)=rho.*u;
    E(:,:,2)=(rho.*(u.^2))+p-tauxx;
    E(:,:,3)=(rho.*u.*v)-tauxy;
    E(:,:,4)=(Et+p).*u-(u.*tauxx)-v.*tauxy+qdotx;

    %For Fbar, tau gets central in x, forward in y

    tauyy=2*mu.*(ddy_fwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_fwd(v,dy))));
    tauxy=mu.*(ddy_fwd(u,dy)+ddx_central(v,dx));

    %For Fbar, qdoty gets forward in y

    qdoty=-k.*ddy_fwd(T,dy);

    F(:,:,1)=rho.*v;
    F(:,:,2)=rho.*u.*v-tauxy;
    F(:,:,3)=rho.*v.^2+p-tauyy;
    F(:,:,4)=(Et+p).*v-v.*tauyy-u.*tauxy+qdoty;

    %now the corrector step calculation, fwd and bwd difference functions account for 3D arrays E and F
   
    U=0.5*(U+Ubar-dt*(ddx_bwd(E,dx)+ddy_bwd(F,dy)));

    %update primitive variables
    [~,u,v,T,p,~,~] = cons2prim(U,R,cv);

    %ENFORCE BCs

    %at the outflow
    u(nx,:)=2*u(nx-1,:)-u(nx-2,:);
    v(nx,:)=2*v(nx-1,:)-v(nx-2,:);
    p(nx,:)=2*p(nx-1,:)-p(nx-2,:);
    T(nx,:)=2*T(nx-1,:)-T(nx-2,:);

    %at the inflow and far-field
    u(1,:)=uinf;
    u(:,ny)=uinf;
    p(1,:)=pinf;
    p(:,ny)=pinf;
    T(1,:)=Tinf;
    T(:,ny)=Tinf;
    v(1,:)=0;
    v(:,ny)=0;

    %at the wall
    u(:,1)=0;
    v(:,1)=0;
    T(:,1)=Tinf; %constant temperature wall
    p(:,1)=2*p(:,3)-p(:,2);

    %at the leading edge
    u(1,1)=0;
    p(1,1)=pinf;
    T(1,1)=Tinf;
    
    rho=p./R./T;
    mu=sutherland(T);
    k=cp/Pr.*mu;

    %update conservative variables with the new primitive variables
    U=prim2cons(rho,u,v,T,cv);
    Et=U(:,:,4);

    t=t+dt;
    count=count+1;
end

%calculate Mach Angle
M_theta = asin(1/M);
M_theta1 = asin(a0/max(u(1,:)));
    
[gx,gy] = gradient(rho ,dx,dy);
G       = hypot(gx,gy);

% choose two x-rows where the shock is still straight
iA =  round(0.05*nx);            % 5% of plate length
iB =  round(0.2*nx);            % 20% of plate length

% in each row pick the strongest gradient ABOVE the boundary layer
jA = find(G(iA,5:end)==max(G(iA,5:end)),1,'first') + 4;   % skip first 4 y-cells
jB = find(G(iB,5:end)==max(G(iB,5:end)),1,'first') + 4;

% physical coordinates of those two shock points
xA = x(iA,jA);   yA = y(iA,jA);
xB = x(iB,jB);   yB = y(iB,jB);

% Mach angle from slope
M_theta2 = atan( (yB - yA) / (xB - xA) );



fprintf('Theoretical θ = %.2f°\n', rad2deg(M_theta));
fprintf('Top-edge estimate θ = %.2f°\n', rad2deg(M_theta1));
fprintf('Numerical slope-based θ = %.2f°\n', rad2deg(M_theta2));