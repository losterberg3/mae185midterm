%% MAE 185 Midterm 
clear; clc; close all;

%defining variables
L=10^-5;H=8*10^-6;nx=75;ny=80;
cp=1005;cv=718;R=cp-cv;
gamma=cp/cv;Pr=0.71;Tinf=288.15;pinf=101300;rho0=1.225;M=4;
a0=sqrt(gamma*R*Tinf); % speed of sound
uinf=M*a0;
t=0; 
%add variables for midterm part 2.1
beta = 0.8; kappa = 10;

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
%for convergence plot
T_old = Tinf;
u_old = u;
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

%% while loop

count=0;
i=1;

figure(1)
while count<1500
    vp=max(4/3.*mu.*gamma.*mu./Pr./rho,[],'all');
    dt=(abs(u)/dx + abs(v)/dy + sqrt(gamma.*R.*T.*(dx^-2 + dy^-2))+ 2*vp*(dx^-2 + dy^-2)).^-1;
    dt=min(dt,[],'all');
    
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
    E=real(E);

    %For F, tau gets central in x, backwards in y

    tauyy=2*mu.*(ddy_bwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_bwd(v,dy))));
    tauxy=mu.*(ddy_bwd(u,dy)+ddx_central(v,dx));

    %For F, qdoty gets backward in y

    qdoty=-k.*ddy_bwd(T,dy);

    F(:,:,1)=rho.*v;
    F(:,:,2)=rho.*u.*v-tauxy;
    F(:,:,3)=rho.*v.^2+p-tauyy;
    F(:,:,4)=(Et+p).*v-v.*tauyy-u.*tauxy+qdoty;
    F=real(F);

    %the predictor step calculation, fwd and bwd difference functions account for 3D arrays E and F

    Ubar=U-dt*(ddx_fwd(E,dx)+ddy_fwd(F,dy));

    %now get primitive variables from Ubar so we can plug into Ebar, Fbar

    [~,u,v,T,p,~,~] = cons2prim(Ubar,R,cv);

    %ENFORCE BCs

    %at the outflow
    u(nx,:)=2*u(nx-2,:)-u(nx-1,:);
    v(nx,:)=2*v(nx-2,:)-v(nx-1,:);
    p(nx,:)=2*p(nx-2,:)-p(nx-1,:);
    T(nx,:)=2*T(nx-2,:)-T(nx-1,:);

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

    %clamping these values
    T = max(T, 1e-3);
    p = max(p, 1e-3);

    rho=p./R./T;
    mu=sutherland(T);
    k=cp/Pr.*mu;

    if any(T(:) <= 0)
        error('Negative or zero temperature encountered!');
    end

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
    E=real(E);

    %For Fbar, tau gets central in x, forward in y

    tauyy=2*mu.*(ddy_fwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_fwd(v,dy))));
    tauxy=mu.*(ddy_fwd(u,dy)+ddx_central(v,dx));

    %For Fbar, qdoty gets forward in y

    qdoty=-k.*ddy_fwd(T,dy);

    F(:,:,1)=rho.*v;
    F(:,:,2)=rho.*u.*v-tauxy;
    F(:,:,3)=rho.*v.^2+p-tauyy;
    F(:,:,4)=(Et+p).*v-v.*tauyy-u.*tauxy+qdoty;
    F=real(F);

    %now the corrector step calculation, fwd and bwd difference functions account for 3D arrays E and F
   
    U=0.5*(U+Ubar-dt*(ddx_bwd(E,dx)+ddy_bwd(F,dy)));

    %update primitive variables
    [~,u,v,T,p,~,~] = cons2prim(U,R,cv);

    %ENFORCE BCs

    %at the outflow
    u(nx,:)=2*u(nx-2,:)-u(nx-1,:);
    v(nx,:)=2*v(nx-2,:)-v(nx-1,:);
    p(nx,:)=2*p(nx-2,:)-p(nx-1,:);
    T(nx,:)=2*T(nx-2,:)-T(nx-1,:);

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
    
    %clamping these values
    T = max(T, 1e-3);
    p = max(p, 1e-3);

    rho=p./R./T;
    mu=sutherland(T);
    k=cp/Pr.*mu;

    if any(T(:) <= 0)
        error('Negative or zero temperature encountered!');
    end

    %update conservative variables with the new primitive variables
    U=prim2cons(rho,u,v,T,cv);
    Et=U(:,:,4);

    u=real(u);
    v=real(v);
    T=real(T);
    convergence_T(i) = norm(T-T_old,2)/norm(T);
    convergence_u(i) = norm(u-u_old,2)/norm(u); %ELEPHANT
    u_old = u;
    T_old = T; %update prev step
    p=real(p);
    k=real(k);
    mu=real(mu);
    rho=real(rho);
    %schlieren phtography
    drho_dx = ddx_bwd(rho,dx);  % compute gradients
    drho_dy = ddy_bwd(rho,dy);

    curl_rho = sqrt(drho_dx.^2 + drho_dy.^2); %calculate curl of rho
    S = beta*exp(-kappa*(abs(curl_rho)/abs(max(max(curl_rho))))); %calculate schlieren photography



    Et=real(Et);
    t=t+dt;
    i=i+1;
    if mod(count,50)==0 || count==0
        figure(1); grid off
        tiledlayout(2,4)
        nexttile %plot rho
        pcolor(x,y,rho), shading interp, axis equal tight;
        title('\rho')
        colorbar
        lim = caxis;
        caxis(lim);

        nexttile %plot Schlieren Photography [S(x,y)]
        pcolor(x,y,S), colormap(gray),shading interp, axis equal tight;
        title('Numerical Schlieren Image');
        colorbar
        caxis([0.1 0.8]);
        colormap(gray);
        lim = caxis;
        caxis(lim);

        nexttile %plot u
        pcolor(x,y,u), colormap(parula), shading interp, axis equal tight;
        title('u')
        colorbar
        lim = caxis;
        caxis(lim);

        nexttile %plot v
        pcolor(x,y,v), shading interp, axis equal tight;
        title('v')
        colorbar
        lim = caxis;
        caxis(lim);

        nexttile %plot Et
        pcolor(x,y,Et), shading interp, axis equal tight;
        title('Et')
        colorbar
        lim = caxis;
        caxis(lim);

        nexttile %plot p
        pcolor(x,y,p), shading interp, axis equal tight;
        title('p')
        colorbar
        lim = caxis;
        caxis(lim);

        nexttile %plot T
        pcolor(x,y,T), shading interp, axis equal tight;
        title('T')
        colorbar
        lim = caxis;
        caxis(lim);

        drawnow   

    end
    
    count=count+1;
% convergence plot
    % if mod(count, 50) == 0
    %     figure(2);tiledlayout(2,2);
    %     nexttile;
    %     semilogy(1:count, convergence_u, 'r-');
    %     nexttile;
    %     plot(1:count, convergence_u, 'r-');
    %     title('Convergence of velocity');
    %     xlabel('Iteration'); ylabel('Residual');
    %     nexttile;
    %     semilogy(1:count, convergence_T, 'r-');
    %     nexttile;
    %     plot(1:count, convergence_T, 'r-');
    %     title('Convergence of temperature');
    %     xlabel('Iteration'); ylabel('Residual');
    %     drawnow;
    % end

end

    %calculate Mach Angle
    M_theta = asin(1/M);
    M_theta1 = asin(a0/max(max(u(1,:))));
    
[gx,gy] = gradient(rho ,dx,dy);
G       = hypot(gx,gy);

% choose two x-rows where the shock is still straight
iA =  round(0.05*nx);            % 5 % of plate length
iB =  round(0.2*nx);            % 20 % of plate length

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

%% FUNCTIONS

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
