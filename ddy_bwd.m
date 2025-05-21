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