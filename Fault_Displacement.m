% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 2, pg.32-34
% Slip occurs over an infinite interval in the u2 direction in a fullspace
% Notice, stresses decay as 1/r from the dislocation line, not the 
% dislocation plane - it is the the dislocation line that is the source of
% strain
clear; clf; close all;

s = 10; % m of displacement
u1 = linspace(-20,20,40); % profile distance across fault
u2 = linspace(-20,20,40); % depth

[U1,U2] = meshgrid(u1,u2);

disp = fault_disp(s,U1,U2);

[sigma23, sigma13] = fault_stress(s,U1,U2);

figure(1)
contourf(U1,U2,disp,10)
colorbar()
title("Displacement due to vertical fault")
ylabel("Depth (m)")
xlabel("Profile Distance (m)")
set(gca(),"YDir","reverse")

figure(2)
subplot(1,2,1)
contourf(U1,U2,sigma23,10)
colorbar()
title("Stress (2,3) due to vertical fault")
ylabel("Depth (m)")
xlabel("Profile Distance (m)")
set(gca(),"YDir","reverse")

subplot(1,2,2)
contourf(U1,U2,sigma13,10)
colorbar()
title("Stress (1,3) due to vertical fault")
ylabel("Depth (m)")
xlabel("Profile Distance (m)")
set(gca(),"YDir","reverse")

function displacement = fault_disp(s,U1,U2)
displacement = - s/(2*pi)*atan2(U1,U2);
end

function [sigma23, sigma13] = fault_stress(s,U1,U2,shear_mod)
if ~exist("shear_mod","var")
    disp("Default shear modulus used: shear_mod = 35e9")
    shear_mod = 35e9;
end
R = sqrt(U1.^2 + U2.^2);
THETA = acos(U1./R);
sigma23 = s*shear_mod/(2*pi)*sin(THETA)./R;
sigma13 = -s*shear_mod/(2*pi)*cos(THETA)./R;
end