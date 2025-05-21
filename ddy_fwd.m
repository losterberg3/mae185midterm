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