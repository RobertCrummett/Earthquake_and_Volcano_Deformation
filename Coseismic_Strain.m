% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 2, pg.38-39
% Coseismic Faulting
% Good model for certain parts of the San Andreas where d2 = 0
clear; clf; close all

s = 10;
d1 = 10;
x1 = linspace(-d1*10, d1*10);

% Displacement across fault
u3 = s*atan(d1./x1)/pi; % do not use atan2()

% Strain across fault
eps13 = -s/(2*pi*d1)*(1./(1+(x1/d1).^2));

% strain is everywhere negative because earthquakes release shear strain

figure(1)
plot(x1/d1, u3/s, 'k-', LineWidth=2)
title("Coseismic Displacment")
xlabel("Distance from fault (x_{1}/d_{1})")
ylabel("Displacement (u_{3}/s)")
grid on
set(gca,'GridLineStyle','--')

figure(2)
plot(x1/d1, eps13*d1/s, 'k-', 'LineWidth', 2)
title("Coseismic Strain")
xlabel("Distance from fault (x_{1}/d_{1})")
ylabel("Strain (\epsilon_{13}d_{1}/s)")
grid on
set(gca,'GridLineStyle','--')