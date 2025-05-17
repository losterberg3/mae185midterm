%defining variables

L=10^-5;
H=8*10^-5;
nx=75;
ny=80;
dt=2.35*10^-11;
M=4;

%creating grid
[x,y]=ndgrid(0:L/nx:L,0:H/ny:H);






%part 1
%forward differences in x and central differences in y
tauxx1=2*sutherland(T).*(ddx_fwd(u,dx)-(1/3*(ddx_fwd(u,dx)+ddy_central(v,dy))));
tauyy1=2*sutherland(T).*(ddy_central(v,dy)-(1/3*(ddx_fwd(u,dx)+ddy_central(v,dy))));
tauxy1=sutherland(T).*(ddy_central(u,dy)+ddx_fwd(v,dx));

%part2
%central differences in x and backward differences in y
tauxx2=2*sutherland(T).*(ddx_central(u,dx)-(1/3*(ddx_central(u,dx)+ddy_bwd(v,dy))));
tauyy2=2*sutherland(T).*(ddy_bwd(v,dy)-(1/3*(ddx_central(u,dx)+ddy_bwd(v,dy))));
tauxy2=sutherland(T).*(ddy_bwd(u,dy)+ddx_central(v,dx));


















function mu = sutherland(T)

mu0=1.735e-5;
S1=110.4;
T0=288;

mu=mu0*(T/T0).^(3/2).*(T0+S1)./(T+S1);

end

function sol = ddx_bwd(f,dx)
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

function sol = ddy_bwd(f,dy)
    [x,y]=size(f);
    sol=zeros(x,y);
    f=flip(f,1);
    for i=1:x
        for j=1:y
            if j<y
                sol(i,j)=(f(i,j+1)-f(i,j))/dy;
            else
                sol(i,j)=sol(i,j-1);
            end
        end
    end
    sol=flip(sol,1);
end

