function mu = sutherland(T)

mu0=1.735e-5;
S1=110.4;
T0=288.15;

mu=mu0*(T/T0).^(3/2).*(T0+S1)./(T+S1);

end