function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)

rho=squeeze(U(:,:,1));
u=squeeze(U(:,:,2))./rho;
v=squeeze(U(:,:,3))./rho;
T=(squeeze(U(:,:,4))./rho-(u.^2+v.^2)/2)/cv;
e=cv*T;
p=rho.*R.*T;
Et=squeeze(U(:,:,4));

end