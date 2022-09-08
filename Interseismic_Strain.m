% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 2, pg.38-39
% Interseismic deformation
% Let d1 -> inf, while d2 is fixed.
% This model was preposed by Savage and Burford (1970) as a first order
% model of interseismic deformation at places like the San Andreas fault
clear; clf; close all

s = 10;
dsdt = 0.05; % fault slip rate
d2 = 10;
x1 = linspace(-d2*10, d2*10);

% Surface displacement
u3 = s/pi*atan(x1/d2);
% Surface veclocity
nu3 = dsdt/pi*atan(x1/d2);
% Surface strain rate
deps13dt = dsdt/(2*pi*d2)*(1./(1+(x1/d2).^2));

figure(1)
plot(x1/d2, u3/s, 'k-', LineWidth=2)
title("Interseismic Displacment")
xlabel("Distance from fault (x_{1}/d_{2})")
ylabel("Displacement (u_{3}/s)")
grid on
set(gca,'GridLineStyle','--')

figure(2)
plot(x1/d2, nu3, 'k-', LineWidth=2)
title("Interseismic Ground Velocity")
xlabel("Distance from fault (x_{1}/d_{2})")
ylabel("Velocity (\nu_{3})")
grid on
set(gca,'GridLineStyle','--')

figure(3)
plot(x1/d2, deps13dt, 'k-', LineWidth=2)
title("Interseismic Strain Rate")
xlabel("Distance from fault (x_{1}/d_{2})")
ylabel("Strain Rate (d\epsilon_{13}/dt)")
grid on
set(gca,'GridLineStyle','--')

x = [1,2; 3, 4];
xi = [0,5];

[g1_1, g1_2, g2_1, g2_2] = g_displacements(x,xi);