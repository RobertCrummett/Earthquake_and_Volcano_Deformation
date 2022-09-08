% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 2, pg.35-36
% Slip occurs over a finite interval 2*d in the u2 direction in a fullspace
% Note that the stress function gives the stress increase/decrease relative
% to a preexisting field. Therefore, upon slipping, we expect a decrease of
% stress in he slipping region and an increase outside. This is observed
% clearly in Fig. 2
clear; clf; close all;

s = 10; % m of displacement
d = 10; % halfwidth of fault (m)
u1 = linspace(-20,20,40); % profile distance across fault
u2 = linspace(-20, 20, 40); % depth

[U1,U2] = meshgrid(u1,u2);

disp = finite_fault_disp(s,U1,U2,d);
sigma13 = finite_fault_stress(s,U1,U2,d);

figure(1)
contourf(U1,U2,disp,10)
colorbar()
title("Displacement due to finite vertical fault")
ylabel("Depth (m)")
xlabel("Profile Distance (m)")
set(gca(),"YDir","reverse")

figure(2)
contourf(U1,U2,sigma13,10)
colorbar()
title("Stress (1,3) due to finite vertical fault")
ylabel("Depth (m)")
xlabel("Profile Distance (m)")
set(gca(),"YDir","reverse")

function displacement = finite_fault_disp(s,U1,U2,d)
displacement = -s/(2*pi)*(atan2(U1,U2+d) - atan2(U1,U2-d));
end

function sigma13 = finite_fault_stress(s,U1,U2,d,shear_mod)
% for oppositely signed dipole dislocations
if ~exist("shear_mod","var")
    disp("Default shear modulus used: shear_mod = 35e9")
    shear_mod = 35e9;
end
R1 = sqrt(U1.^2 + (U2 + d).^2);
R2 = sqrt(U1.^2 + (U2 - d).^2);
sigma13 = -s*shear_mod/(2*pi)*((U2 + d)./(R1.^2) - (U2 - d)./(R2.^2));
end