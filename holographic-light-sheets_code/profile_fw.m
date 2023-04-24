function [An,in,n] = profile_fw (I,L,N)
syms z

F(z) = I;

zz = linspace(0,L,4000);

in = 1:2*N+1;
n = -N-1+in;

fF = matlabFunction(F);

An = zeros(size(n));
for j=1:length(n)
An(j) = (1/L)*trapz(zz,fF(zz).*exp(1i*((2*pi/L)*n(j)).*zz));
end


end