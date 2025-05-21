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