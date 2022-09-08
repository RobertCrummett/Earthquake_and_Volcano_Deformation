% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 3, pg.63
% Edge dislocation fault
clear; clf; close all

s = 10;
d = 10;
xi = linspace(-200,200,1000); % points to realize displacements at

u = edge_diloc(s,d,xi);
u_1 = u{1};
u_2 = u{2};

figure(1)
plot(xi,u_1/s,'r-',LineWidth=2)
hold on
plot(xi,u_2/s,'b-',LineWidth=2)
hold off
legend("$u_{1}$/s","$u_{2}$/s","Interpreter","Latex")
xlabel("\xi / d (m)")
ylabel("u / s (m)")
grid on
set(gca,"GridLineStyle","--")

function u = edge_diloc(s, d, xi)
% vertical fault slipping uniformly below d
zeta = xi/d;
u1 = -s/pi*(1./(1+zeta.^2));
u2 = s/pi*(zeta./(1+zeta.^2) + atan(zeta));
u = {u1, u2};
end
