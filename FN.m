Eb=0.5;

T=300;
kb = 1.3806503e-23; % m^2kgs^{-2}K^{-1}
q = 1.60217646e-19; % C
hbar = 1.05457148e-34; %m^2kg/s
m0 = 9.10938188e-31; % kg;
kbT = kb*T; 
meff=1*m0;
d=10E-9; %m
factor=1E4;

for ii=1:30
V(ii)=ii*0.1;    
if V(ii)<=Eb
J(ii)=q^3*meff/(4*hbar*meff*((q*Eb)^0.5-(q*Eb-q*V(ii))^0.5)^2)*(V(ii)/d)^2*exp(-4*(2*meff)^0.5/(3*hbar*q*V(ii))*d*((q*Eb)^1.5-(q*Eb-q*V(ii))^1.5));
else
J(ii)=q^3*meff/(4*hbar*meff*q*Eb)*(V(ii)/d)^2*exp(-4*(2*meff*(q*Eb)^3)^0.5/(3*hbar*q*V(ii))*d);
end
end
J=J/factor;

semilogy(V,J,'ro--');