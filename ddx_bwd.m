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
