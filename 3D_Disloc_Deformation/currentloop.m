function [Bx, Br] = currentloop(x,r,i,a)

al = r./a;
be = x./a;
ga = x./r;
Q = ((1 + al).^2 + be.^2);
k = 2*sqrt(al./Q);

mu0 = 4*pi*10^-7; % magnetic permeability of free space

B0 = i*mu0./(2*a); % field at the middle of the wire

[K,E] = ellipke(k.^2,eps*100);

Bx = B0./(pi*sqrt(Q)).*(E.*(1-al.^2-be.^2)./(Q-4.*al) + K);
Br = B0.*ga./(pi*sqrt(Q)).*(E.*(1+al.^2+be.^2)./(Q-4.*al) - K);
end
