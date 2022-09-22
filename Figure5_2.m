% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 5, pg 119 - 120
% Displacements across and infinitely long strike-slip fault in an elastic
% half-space with different shear modulus on either sife of the fault. mu1
% refers to the shear modulus on the right side of the fault. See figure
% 5.1 in the text for more information
clear; clf; close all

s = 1;
x1_1 = linspace(0, 100, 1000);
x1_2 = linspace(-100, 0, 1000);
d = 10;

ratio = [1,2,5,10];
colors = ["#d3d3d3","#aeaeae","#757575","#000000"];

figure(1)
for idx = 1:length(ratio)
    mu1 = 35e9;
    mu2 = mu1/ratio(idx);
    A = s/pi*mu2/(mu1 + mu2);
    B = s/pi*mu1/(mu1 + mu2);
    
    u1 = 2*A*atan2(x1_1,d);
    u2 = 2*B*atan2(x1_2,d);

    hold on
    plot(x1_1/d,u1,'Color',colors(5-idx),'LineStyle','-','LineWidth',1.75)
    plot(x1_2/d,u2,'Color',colors(5-idx),'LineStyle','-','LineWidth',1.75)
    hold off
end
grid on
set(gca,"GridLineStyle",":")
xlabel("Distance from fault (x_{1}/d)")
ylabel("Displacement")
legend("$\mu_{1}/\mu_{2} = 1$","$\mu_{1}/\mu_{2} = 2$",...
    "$\mu_{1}/\mu_{2} = 5$","$\mu_{1}/\mu_{2} = 10$",...
    "Interpreter","Latex","Location","se")