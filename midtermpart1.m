%defining variables
L=10^-5;
H=8*10^-5;
nx=75;
ny=80;
dt=2.35*10^-11;
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
t=0;

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

%U, E, and F
U=prim2cons(rho0,u,v,T,cv);
E=zeros(4,nx,ny);
F=zeros(4,nx,ny);
Ebar=zeros(4,nx,ny);
Fbar=zeros(4,nx,ny);

%%


while t<1
    mu=sutherland(T);
    k=cp/Pr.*mu;
    rho=p./R./T;
    e=cv*T;
    Et=rho.*(e+(u.^2+v.^2)/2);

    %predictor step, fwd FD for x and y

    %getting U

    U=prim2cons(rho,u,v,T,cv);

    %For E, tau gets backward in x, central in y
    
    tauxx=2*mu.*(ddx_bwd(u,dx)-(1/3*(ddx_bwd(u,dx)+ddy_central(v,dy))));
    tauxy=mu.*(ddy_central(u,dy)+ddx_bwd(v,dx));

    %For E, qdotx gets backward in x

    qdotx=-k.*ddx_bwd(T,dx);

    E(1,:,:)=rho.*u;
    E(2,:,:)=rho.*u.^2+p-tauxx;
    E(3,:,:)=rho.*u.*v-tauxy;
    E(4,:,:)=(Et+p).*u-u.*tauxx-v.*tauxy+qdotx;

    %For F, tau gets central in x, backwards in y

    tauyy=2*mu.*(ddy_bwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_bwd(v,dy))));
    tauxy=mu.*(ddy_bwd(u,dy)+ddx_central(v,dx));

    %For F, qdoty gets backward in y

    qdoty=-k.*ddy_bwd(T,dy);

    F(1,:,:)=rho.*v;
    F(2,:,:)=rho.*u.*v-tauxy;
    F(3,:,:)=rho.*v.^2+p-tauyy;
    F(4,:,:)=(Et+p).*v-v.*tauyy-u.*tauxy+qdoty;

    %the predictor step

    Ubar=U-dt*(ddx_fwd(E,dx)+ddy_fwd(F,dy));

    %now get primitive variables from Ubar so we can plug into Ebar, Fbar

    [rho,u,v,T,p,e,Et] = cons2prim(Ubar,R,cv);

    %then the corrector step, backward fd in x and y

    %For Ebar, tau gets forward in x, central in y
    
    tauxx=2*mu.*(ddx_fwd(u,dx)-(1/3*(ddx_fwd(u,dx)+ddy_central(v,dy))));
    tauxy=mu.*(ddy_central(u,dy)+ddx_fwd(v,dx));

    %For Ebar, qdotx gets forward in x

    qdotx=-k.*ddx_fwd(T,dx);

    Ebar(1,:,:)=rho.*u;
    Ebar(2,:,:)=rho.*u.^2+p-tauxx;
    Ebar(3,:,:)=rho.*u.*v-tauxy;
    Ebar(4,:,:)=(Et+p).*u-u.*tauxx-v.*tauxy+qdotx;

    %For Fbar, tau gets central in x, forward in y

    tauyy=2*mu.*(ddy_fwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_fwd(v,dy))));
    tauxy=mu.*(ddy_fwd(u,dy)+ddx_central(v,dx));

    %For Fbar, qdoty gets forward in y

    qdoty=-k.*ddy_fwd(T,dy);

    Fbar(1,:,:)=rho.*v;
    Fbar(2,:,:)=rho.*u.*v-tauxy;
    Fbar(3,:,:)=rho.*v.^2+p-tauyy;
    Fbar(4,:,:)=(Et+p).*v-v.*tauyy-u.*tauxy+qdoty;

    %ENFORCE BCs, not done yet

    %
    %
    %
    %
    
    %now the corrector step calculation
    U=0.5*(U+Ubar-dt*(ddx_bwd(E,dx)+ddy_bwd(F,dy)));







    t=t+dt;
end



a=sqrt(gamma*R*T);
M=uinf/a;
rho=p./R./T;
mu=sutherland(T);
k=cp/Pr*mu;
e=cv*T;








%forward differences in x and central differences in y
tauxx1=2*mu.*(ddx_fwd(u,dx)-(1/3*(ddx_fwd(u,dx)+ddy_central(v,dy))));
tauyy1=2*mu.*(ddy_central(v,dy)-(1/3*(ddx_fwd(u,dx)+ddy_central(v,dy))));
tauxy1=mu.*(ddy_central(u,dy)+ddx_fwd(v,dx));

%part2
%central differences in x and backward differences in y
tauxx2=2*mu.*(ddx_central(u,dx)-(1/3*(ddx_central(u,dx)+ddy_bwd(v,dy))));
tauyy2=2*mu.*(ddy_bwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_bwd(v,dy))));
tauxy2=mu.*(ddy_bwd(u,dy)+ddx_central(v,dx));













% FUNCTIONS

function U = prim2cons(rho,u,v,T,cv)

[nx,ny]=size(u);
e=cv.*T;
U=zeros(4,nx,ny);

U(1,:,:)=rho;
U(2,:,:)=rho.*u;
U(3,:,:)=rho.*v;
U(4,:,:)=rho.*(e+(u.^2+v.^2)/2);

end


function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)

rho=squeeze(U(1,:,:));
u=squeeze(U(2,:,:))./rho;
v=squeeze(U(3,:,:))./rho;
T=(squeeze(U(4,:,:))./rho-(u.^2+v.^2)/2)/cv;
e=cv*T;
p=rho.*R.*T;
Et=squeeze(U(4,:,:));

end



function mu = sutherland(T)

mu0=1.735e-5;
S1=110.4;
T0=288;

mu=mu0*(T/T0).^(3/2).*(T0+S1)./(T+S1);

end

function sol = ddx_bwd(f,dx)
    if ndims(f)==3
        [n,x,y]=size(f);
        sol=zeros(n,x,y);
        for k=1:n
            f=flip(f,2);
            for i=1:x
                for j=1:y
                    if i<x
                        sol(k,i,j)=(f(k,i+1,j)-f(k,i,j))/dx;
                    else
                        sol(k,i,j)=sol(k,i-1,j);
                    end
                end
            end
            SOL=squeeze(sol(k,:,:));
            sol(k,:,:)=flip(SOL,2);
        end
    else      
        [x,y]=size(f);
        sol=zeros(x,y);
        f=flip(f,2);
        for i=1:x
            for j=1:y
                if i<x
                    sol(i,j)=(f(i+1,j)-f(i,j))/dx;
                else
                    sol(i,j)=sol(i-1,j);
                end
            end
        end
        sol=flip(sol,2);
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
        [n,x,y]=size(f);
        sol=zeros(n,x,y);
        for k=1:n
            for i=1:x
                for j=1:y
                    if i<x
                        sol(k,i,j)=(f(k,i+1,j)-f(k,i,j))/dx;
                    else
                        sol(k,i,j)=sol(k,i-1,j);
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
        [n,x,y]=size(f);
        sol=zeros(n,x,y);
        for k=1:n
            for i=1:x
                for j=1:y
                    if j<y
                        sol(k,i,j)=(f(k,i,j+1)-f(k,i,j))/dy;
                    else
                        sol(k,i,j)=sol(k,i,j-1);
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
        [n,x,y]=size(f);
        sol=zeros(n,x,y);
        for k=1:n
            f=flip(f,2);
            for i=1:x
                for j=1:y
                    if i<x
                        sol(k,i,j)=(f(k,i,j+1)-f(k,i,j))/dy;
                    else
                        sol(k,i,j)=sol(k,i,j-1);
                    end
                end
            end
            SOL=squeeze(sol(k,:,:));
            sol(k,:,:)=flip(SOL,1);
        end
    else      
        [x,y]=size(f);
        sol=zeros(x,y);
        f=flip(f,2);
        for i=1:x
            for j=1:y
                if i<x
                    sol(i,j)=(f(i,j+1)-f(i,j))/dy;
                else
                    sol(i,j)=sol(i,j-1);
                end
            end
        end
        sol=flip(sol,1);
    end
end
