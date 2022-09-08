% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 3, pg.67-68
% Dip-Slip Faults
% This code recreates Figure 3.10, the displacements associated with thrust
% faulting
clear; clf; close all
s = -2;
dip = 20;
d = 3.8;
x = (-10:0.1:30);

dip = deg2rad(dip);
zeta = (x - cot(dip)*d)/d;

u1_rigid = s/pi*(cos(dip)*(atan2(x,d) - dip + pi/2) - x.*(d*cos(dip) + x*sin(dip))./(x.^2 + d^2));
u2_rigid = -s/pi*(sin(dip)*(atan2(x,d) - dip + pi/2) + d.*(d*cos(dip) + x*sin(dip))./(x.^2 + d^2));

u1 = s/pi*(cos(dip)*atan2(x - d*cot(dip), d) + (sin(dip) - zeta*cos(dip))./(1 + zeta.^2));
u2 = -s/pi*(sin(dip)*atan2(x - d*cot(dip), d) + (cos(dip) + zeta*sin(dip))./(1 + zeta.^2));

figure(1)
plot(x,u1_rigid/abs(s),'--k','LineWidth',2); hold on
plot(x,u2_rigid/abs(s),'--b','LineWidth',2); 
plot(x,u1/abs(s),'-k','LineWidth',2);
plot(x,u2/abs(s),'-r','LineWidth',2);hold off
legend('Horizontal', 'Vertical', 'Horizontal Shifted','Vertical Shifted','Location', 'best')
grid on
set(gca,'GridLineStyle','--')

u1_s = -s/pi*(cos(dip)*(atan2(x - d*cot(dip),d) - pi/2*sign(x)) + (sin(dip) - zeta*cos(dip))./(1 + zeta.^2));
u2_s = s/pi*(sin(dip)*(atan2(x - d*cot(dip),d) - pi/2*sign(x)) + (cos(dip) + zeta*sin(dip))./(1 + zeta.^2));

figure(2)
subplot(2,1,1)
plot(x,u1/abs(s),'--k','LineWidth',2); hold on
plot(x,u2/abs(s),'--b','LineWidth',2); hold off
title('Equations 3.70')
legend('Horizontal Corrected', 'Vertical Corrected','Location','best')
grid on
set(gca,'GridLineStyle','--')
subplot(2,1,2)
plot(x,u1_s/abs(s),'-k','LineWidth',2); hold on
plot(x,u2_s/abs(s),'-r','LineWidth',2);hold off
legend('Horizontal Displacement','Vertical Displacment','Location', 'best')
title('Figure 3.10, Thrust Faulting')
xlabel('Distance from trench')
ylabel('u / s')
grid on
set(gca,'YLim',[-1,0.5])
set(gca,'GridLineStyle','--')