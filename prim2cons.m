function U = prim2cons(rho,u,v,T,cv)

    [nx,ny]=size(u);
    e=cv.*T;
    U=zeros(nx,ny,4);

    U(:,:,1)=rho;
    U(:,:,2)=rho.*u;
    U(:,:,3)=rho.*v;
    U(:,:,4)=rho.*(e+(u.^2+v.^2)/2);

end