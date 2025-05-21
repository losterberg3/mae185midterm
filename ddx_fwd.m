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