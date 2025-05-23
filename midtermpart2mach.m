%% MAE 185 Midterm 
clear; clc; close all;

%defining variables
L=10^-5;
H=8*10^-6;
nx=75;
ny=80;
cp=1005;
cv=718;
R=cp-cv;
gamma=cp/cv;
Pr=0.71;
Tinf=288.15;
pinf=101300;
rho0=1.225;
M=4;
a0=sqrt(gamma*R*Tinf);
uinf=M*a0;
<<<<<<< HEAD
t=0;
=======
t=0; 
>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed
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

<<<<<<< HEAD
%for convergence plot
T_old = Tinf;
u_old = u;
convergence_T=zeros(1,1500);
convergence_u=zeros(1,1500);

=======
>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed
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
<<<<<<< HEAD

count=0;
figure
=======
count=0;

>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed
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
    %qdotx(:,1)=0; %adiabatic wall condition

    E(:,:,1)=rho.*u;
    E(:,:,2)=rho.*u.^2+p-tauxx;
    E(:,:,3)=rho.*u.*v-tauxy;
    E(:,:,4)=(Et+p).*u-u.*tauxx-v.*tauxy+qdotx;

    %For F, tau gets central in x, backwards in y

    tauyy=2*mu.*(ddy_bwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_bwd(v,dy))));
    tauxy=mu.*(ddy_bwd(u,dy)+ddx_central(v,dx));

    %For F, qdoty gets backward in y

    qdoty=-k.*ddy_bwd(T,dy);
    %qdoty(:,1)=0;  %adiabatic wall condition

    F(:,:,1)=rho.*v;
    F(:,:,2)=rho.*u.*v-tauxy;
    F(:,:,3)=rho.*v.^2+p-tauyy;
    F(:,:,4)=(Et+p).*v-v.*tauyy-u.*tauxy+qdoty;
<<<<<<< HEAD
  
=======

>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed
    %the predictor step calculation, fwd and bwd difference functions account for 3D arrays E and F

    Ubar=U-dt*(ddx_fwd(E,dx)+ddy_fwd(F,dy));

    %now get primitive variables from Ubar so we can plug into Ebar, Fbar

    [~,u,v,T,p,~,~] = cons2prim(Ubar,R,cv);

    %ENFORCE BCs

<<<<<<< HEAD
=======
    %at the outflow
    u(nx,:)=2*u(nx-1,:)-u(nx-2,:);
    v(nx,:)=2*v(nx-1,:)-v(nx-2,:);
    p(nx,:)=2*p(nx-1,:)-p(nx-2,:);
    T(nx,:)=2*T(nx-1,:)-T(nx-2,:);

>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed
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
    T(:,1)=T(:,2);
    p(:,1)=2*p(:,2)-p(:,3);

    %at the leading edge
    u(1,1)=0;
    p(1,1)=pinf;
<<<<<<< HEAD
   
    %at the outflow
    u(nx,:)=2*u(nx-1,:)-u(nx-2,:);
    v(nx,:)=2*v(nx-1,:)-v(nx-2,:);
    p(nx,:)=2*p(nx-1,:)-p(nx-2,:);
    T(nx,:)=2*T(nx-1,:)-T(nx-2,:);
=======
    T(1,1)=Tinf;
>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed

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
    %qdotx(:,1)=0;   %adiabatic wall condition

    E(:,:,1)=rho.*u;
    E(:,:,2)=(rho.*(u.^2))+p-tauxx;
    E(:,:,3)=(rho.*u.*v)-tauxy;
    E(:,:,4)=(Et+p).*u-(u.*tauxx)-v.*tauxy+qdotx;

    %For Fbar, tau gets central in x, forward in y

    tauyy=2*mu.*(ddy_fwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_fwd(v,dy))));
    tauxy=mu.*(ddy_fwd(u,dy)+ddx_central(v,dx));

    %For Fbar, qdoty gets forward in y

    qdoty=-k.*ddy_fwd(T,dy);
    %qdoty(:,1)=0; %adiabatic wall condition

    F(:,:,1)=rho.*v;
    F(:,:,2)=rho.*u.*v-tauxy;
    F(:,:,3)=rho.*v.^2+p-tauyy;
    F(:,:,4)=(Et+p).*v-v.*tauyy-u.*tauxy+qdoty;

    %now the corrector step calculation, fwd and bwd difference functions account for 3D arrays E and F
   
    U=0.5*(U+Ubar-dt*(ddx_bwd(E,dx)+ddy_bwd(F,dy)));

    %update primitive variables
    [~,u,v,T,p,~,~] = cons2prim(U,R,cv);

    %ENFORCE BCs

<<<<<<< HEAD
=======
    %at the outflow
    u(nx,:)=2*u(nx-1,:)-u(nx-2,:);
    v(nx,:)=2*v(nx-1,:)-v(nx-2,:);
    p(nx,:)=2*p(nx-1,:)-p(nx-2,:);
    T(nx,:)=2*T(nx-1,:)-T(nx-2,:);

>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed
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
    T(:,1)=T(:,2);
    p(:,1)=2*p(:,2)-p(:,3);

    %at the outflow
    u(nx,:)=2*u(nx-1,:)-u(nx-2,:);
    v(nx,:)=2*v(nx-1,:)-v(nx-2,:);
    p(nx,:)=2*p(nx-1,:)-p(nx-2,:);
    T(nx,:)=2*T(nx-1,:)-T(nx-2,:);

    %at the leading edge
    u(1,1)=0;
    p(1,1)=pinf;
    
    rho=p./R./T;
    e=cv*T;
    mu=sutherland(T);
    k=cp/Pr.*mu;

    %update conservative variables with the new primitive variables
    U=prim2cons(rho,u,v,T,cv);
    Et=U(:,:,4);

<<<<<<< HEAD
    %convergence plot
    convergence_T(1,count+1) = norm(T-T_old,2)/norm(T);
    T_old = T; %update prev step

    t=t+dt;
    if mod(count,100)==0 || count==0
        tiledlayout(2,3)
        ax1=nexttile;
        pcolor(x,y,rho), shading interp, axis equal tight;
        grid off
        title('\rho [kg/m^3]')
        colormap(ax1,'jet');
        colorbar
        caxis([1.1 3.5])
        ax2=nexttile;
        pcolor(x,y,u), shading interp, axis equal tight;
        grid off
        title('u [m/s]')
        colormap(ax2,'jet');
        colorbar
        caxis([0 1350])
        ax3=nexttile;
        pcolor(x,y,v), shading interp, axis equal tight;
        grid off
        title('v [m/s]')
        colormap(ax3,'jet');
        colorbar
        caxis([0 175])
        ax4=nexttile;
        pcolor(x,y,e), shading interp, axis equal tight;
        grid off
        title('e [J/kg]')
        colormap(ax4,'jet');
        colorbar
        caxis([210000 350000])
        ax5=nexttile;
        pcolor(x,y,p), shading interp, axis equal tight;
        grid off
        title('P [Pa]')
        colormap(ax5,'jet');
        colorbar
        caxis([110000 280000])
        ax6=nexttile;
        pcolor(x,y,T), shading interp, axis equal tight;
        grid off
        title('T [K]')
        colormap(ax6,'hot');
        colorbar
        caxis([290 490])
        drawnow
    end
    count=count+1;
end

figure
plot(dt*(1:1500),convergence_T)
title('Convergence for T')
xlabel('Time [s]')
ylabel('Residual')

=======
    t=t+dt;
    count=count+1;
end

>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed
%calculate Mach Angle
M_theta = asin(1/M);
M_theta1 = asin(a0/max(u(1,:)));
    
[gx,gy] = gradient(rho ,dx,dy);
G       = hypot(gx,gy);

% choose two x-rows where the shock is still straight
<<<<<<< HEAD
iA =  round(0.03 *nx);            % 3 % of plate length
iB =  round(0.15*nx);            % 15 % of plate length
=======
iA =  round(0.05*nx);            % 5% of plate length
iB =  round(0.2*nx);            % 20% of plate length
>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed

% in each row pick the strongest gradient ABOVE the boundary layer
jA = find(G(iA,6:end)==max(G(iA,5:end)),1,'first') + 5;   % skip first 5 y-cells
jB = find(G(iB,6:end)==max(G(iB,5:end)),1,'first') + 5;

% physical coordinates of those two shock points
xA = x(iA,jA);   yA = y(iA,jA);
xB = x(iB,jB);   yB = y(iB,jB);

% Mach angle from slope
M_theta2 = atan( (yB - yA) / (xB - xA) );


fprintf('Theoretical θ = %.2f°\n', rad2deg(M_theta));
fprintf('Top-edge estimate θ = %.2f°\n', rad2deg(M_theta1));
<<<<<<< HEAD
fprintf('Numerical slope-based θ = %.2f°\n', rad2deg(M_theta2));



% FUNCTIONS

function U = prim2cons(rho,u,v,T,cv)

[nx,ny]=size(u);
e=cv.*T;
U=zeros(nx,ny,4);

U(:,:,1)=rho;
U(:,:,2)=rho.*u;
U(:,:,3)=rho.*v;
U(:,:,4)=rho.*(e+(u.^2+v.^2)/2);

end


function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)

rho=squeeze(U(:,:,1));
u=squeeze(U(:,:,2))./rho;
v=squeeze(U(:,:,3))./rho;
T=(squeeze(U(:,:,4))./rho-(u.^2+v.^2)/2)/cv;
e=cv*T;
p=rho.*R.*T;
Et=squeeze(U(:,:,4));

end



function mu = sutherland(T)

mu0=1.735e-5;
S1=110.4;
T0=288.15;

mu=mu0*(T/T0).^(3/2).*(T0+S1)./(T+S1);

end

function sol = ddx_bwd(f,dx)
    if ndims(f)==3
        [x,y,n]=size(f);
        sol=zeros(x,y,n);
        for k=1:n
            for i=1:x
                for j=1:y
                    if i>1
                        sol(i,j,k)=(f(i,j,k)-f(i-1,j,k))/dx;
                    else
                        sol(i,j,k)=(f(i+1,j,k)-f(i,j,k))/dx;
                    end
                end
            end
        end
    else      
        [x,y]=size(f);
        sol=zeros(x,y);
        for i=1:x
            for j=1:y 
                if i>1                       
                    sol(i,j)=(f(i,j)-f(i-1,j))/dx;                    
                else                       
                    sol(i,j)=(f(i+1,j)-f(i,j))/dx;
                end
            end
        end
    end
end







function sol = ddx_central(f,dx)
    [x,y]=size(f);
    sol=zeros(x,y);
    for i=1:x
        for j=1:y
            if i==1
                sol(i,j)=(-3*f(i,j)+4*f(i+1,j)-f(i+2,j))/2/dx;
            elseif i==x
                sol(i,j)=(3*f(i,j)-4*f(i-1,j)+f(i-2,j))/2/dx;
            else
                sol(i,j)=(f(i+1,j)-f(i-1,j))/2/dx;
            end
        end
    end
end







function sol = ddx_fwd(f,dx)
    if ndims(f)==3
        [x,y,n]=size(f);
        sol=zeros(x,y,n);
        for k=1:n
            for i=1:x
                for j=1:y
                    if i<x
                        sol(i,j,k)=(f(i+1,j,k)-f(i,j,k))/dx;
                    else
                        sol(i,j,k)=sol(i-1,j,k);
                    end
                end
            end
        end
    else      
        [x,y]=size(f);
        sol=zeros(x,y); 
        for i=1:x    
            for j=1:y                
                if i<x                        
                    sol(i,j)=(f(i+1,j)-f(i,j))/dx;                    
                else                        
                    sol(i,j)=sol(i-1,j);                   
                end                
            end            
        end
    end
end








function sol = ddy_central(f,dy)
    [x,y]=size(f);
    sol=zeros(x,y);
    for i=1:x
        for j=1:y
            if j==1
                sol(i,j)=(-3*f(i,j)+4*f(i,j+1)-f(i,j+2))/2/dy;
            elseif j==y
                sol(i,j)=(3*f(i,j)-4*f(i,j-1)+f(i,j-2))/2/dy;
            else
                sol(i,j)=(f(i,j+1)-f(i,j-1))/2/dy;
            end
        end
    end
end









function sol = ddy_fwd(f,dy)
    if ndims(f)==3
        [x,y,n]=size(f);
        sol=zeros(x,y,n);
        for k=1:n
            for i=1:x
                for j=1:y
                    if j<y
                        sol(i,j,k)=(f(i,j+1,k)-f(i,j,k))/dy;
                    else
                        sol(i,j,k)=sol(i,j-1,k);
                    end
                end
            end
        end
    else      
        [x,y]=size(f);
        sol=zeros(x,y); 
        for i=1:x    
            for j=1:y                
                if j<y                        
                    sol(i,j)=(f(i,j+1)-f(i,j))/dy;                    
                else                        
                    sol(i,j)=sol(i,j-1);                   
                end                
            end            
        end
    end
end










function sol = ddy_bwd(f,dy)
    if ndims(f)==3
        [x,y,n]=size(f);
        sol=zeros(x,y,n);
        for k=1:n
            for i=1:x
                for j=1:y
                    if j>1
                        sol(i,j,k)=(f(i,j,k)-f(i,j-1,k))/dy;
                    else
                        sol(i,j,k)=(f(i,j+1,k)-f(i,j,k))/dy;
                    end
                end
            end
        end
    else      
        [x,y]=size(f);
        sol=zeros(x,y);
        for i=1:x
            for j=1:y 
                if j>1                       
                    sol(i,j)=(f(i,j)-f(i,j-1))/dy;                    
                else                       
                    sol(i,j)=(f(i,j+1)-f(i,j))/dy;
                end
            end
        end
    end
end
=======
fprintf('Numerical slope-based θ = %.2f°\n', rad2deg(M_theta2));
>>>>>>> be7848fafb9be178045d9cf0e42aa34e5b9b94ed
