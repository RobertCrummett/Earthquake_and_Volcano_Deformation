% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 4, pg 109 = 110
% Recreation of Figure 4.12 and Figure 4.13
clear; clf; close all

% input params
z = linspace(0,2,1000);
x = linspace(-20,20,1000);
s = 1;
ad = 1;
ac = 1.25;
as = 1.6;
[slip,displ] = uniform_disloc(s,x,z,ad);
[cslip,cdispl] = cstressdrop_crack(s,x,z,ac);
[tslip,tdispl] = tapered_crack(s,x,z,as);

% plots
figure(1)
subplot(1,2,1)
plot(slip,z,"k-","LineWidth",2)
hold on
plot(cslip,z,"b-","LineWidth",2)
plot(tslip,z,"r-","LineWidth",2)
hold off
set(gca,"YDir","reverse")
grid on
set(gca,"GridLineStyle",":")
title("Cross-Section")
ylim([0,2])
xlim([-0.1,1.6])
xlabel("Slip")
ylabel("Depth from Surface")
legend("Uniform-Slip Disloc.","Constant Stress Drop Crack",...
    "Tapered Crack","Location","se")

subplot(1,2,2)
plot(x/ad,displ,"k-","LineWidth",1.25)
hold on
plot(x/ac,cdispl,"b-","LineWidth",1.25)
plot(x/as,tdispl,"r-","LineWidth",1.25)
hold off
grid on
title("Map View")
set(gca,"GridLineStyle",":")
xlabel("Distance from fault (x/d)")
ylabel("Displacement")
xlim([-10,10])
legend("Uniform-Slip Disloc.","Constant Stress Drop Crack",...
    "Tapered Crack","Location","se")

sgtitle("Figure 4.12")

figure(2)
plot(x/ad,cdispl - displ, 'LineWidth',2)
hold on
plot(x/ad,tdispl - displ, 'LineWidth',2)
hold off
xlim([0,5])
legend("Crack - Dislocation","Tapered Crack - Dislocation",...
    "Location","e")
grid on
set(gca,"GridLineStyle",":")
xlabel("Distance from fault (x/d)")
ylabel("Displacement differences")
title("Figure 4.13")

% dlip and displacement scenario functions
function [slip, displ] = uniform_disloc(s,x,z,ad)
slip = zeros(size(z));
slip(z < ad) = s;
displ = s/pi*atan(ad./x); % not atan2(), use atan
end

function [slip, displ]  = cstressdrop_crack(s,x,z,ac)
slip = s*sqrt(1 - (z/ac).^2);
displ = s/2*(sign(x).*sqrt(1 + (x./ac).^2) - x/ac);
end

function [slip,displ] = tapered_crack(s,x,z,as)
slip = s*(1 - (z/as).^2).^(3/2);
displ = -3/4*s*((x/as) + 2/3*(x/as).^3 - 2/3*sign(x).*(1 + (x/as).^2).^(3/2));
end

