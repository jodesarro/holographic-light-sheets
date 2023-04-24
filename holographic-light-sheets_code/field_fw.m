function Psi = field_fw(An,in,n,L,K,Q,rho0,phi0)
syms rho phi z

rhog = sqrt(rho.^2+rho0^2-2*rho*rho0.*cos(phi-phi0));

betan(in) = Q+(2*pi/L).*n;
krhon(in) = sqrt(K^2-betan.^2);

aux2(in) = besselj(0,krhon.*rhog).*exp(-1i.*betan.*z);

Psi(rho,phi,z) = sum(An.*aux2);

end