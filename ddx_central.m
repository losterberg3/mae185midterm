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