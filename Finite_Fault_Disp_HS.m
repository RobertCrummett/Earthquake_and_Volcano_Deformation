% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 2, pg.37-38
% Displacement due to finite fault in a half-space
clear; clf; close all

s = 10; % slip
d1 = 10; % depth (greater)
d2 = 5; % depth (less)
x1 = linspace(-50,50,100); % profile distance across fault
x2 = linspace(-25,0,50); % depth

[X1,X2] = meshgrid(x1,x2);

u3 = finite_fault_disp(s,X1,X2,d1,d2);
surfu = finite_fault_surfd(s,x1,d1,d2);

figure(1)
contourf(X1,X2,u3,15)
hold on
plot(x1,zeros(length(x1)),"k-") % free surface boundary
hold off
ylabel("Depth (m)")
ylim([min(x2),5])
xlabel("Profile Distance (m)")
title("Cross-cut displacement across finite fault in halfspace")

% *Useful because this is where measurements are made*
figure(2)
plot(x1,surfu,"k-","LineWidth",2)
hold on; 
plot(x1,zeros(length(x1)),"k--")
hold off
grid on
xlabel("Profile distance across fault (m)")
ylabel("Displacement across fault (m)")
title("Displacements on the free-surface")

function u3 = finite_fault_disp(s,X1,X2,d1,d2)
u3 = -s/(2*pi)*(atan2(X1,X2+d1)-atan2(X1,X2-d1)-atan2(X1,X2+d2)+atan2(X1,X2-d2));
end

function surf_disp = finite_fault_surfd(s,x1,d1,d2)
% surface displacements, ie x2 = 0
% this function has no surface expression
surf_disp = -s/pi*(atan2(x1,d1) - atan2(x1,d2));
end