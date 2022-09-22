% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 5, pg 124 - 125
% Recreation of Figure 5.6
% I multiplied all of Equation 5.23 by a factor of 2 to recreate the
% figures. I am not sure why the Equation 5.23 as given in the text does
% not reproduce the figure exactly.
clear; clf; close all

s = 1;
d = 10;
x = linspace(-d,d,100);
h = 3;

mu = 35e9;
mu_star = mu*[1, 0.8, 0.5];
kappa = (mu - mu_star)./(mu + mu_star);
colors = ['r','b','m'];
ls = ['-','-','-'];

x1 = linspace(-d,-h,100);
x2 = 0;
nupp = 100000;
last = zeros(1, length(x1));
tol = 1e-20;

figure(1)
subplot(2,1,1)
% u^(1) calculation
for i = 1:length(mu_star)
    k = kappa(i);
    c = colors(i);
    l = ls(i);
    hold on
    for n = 0:nupp
        u = k^n*atan2(x1 - 2*n*h, x2 + d);
        u_S = last + u;
        if max(abs(u_S - last)) < tol
            u_S = 2*(1 - k)*s/(2*pi)*u_S;
            break
        end
        last = u_S;
    end
    plot(x1/d,u_S/s,'Color',c,'LineStyle',l,'LineWidth',1.5)
    clear u_S
    last = zeros(1, length(x1));
    hold off
end

x1 = linspace(-h,h,100);
% u^(2) calculation
for i = 1:length(mu_star)
    k = kappa(i);
    c = colors(i);
    l = ls(i);
    hold on
    for n = 1:nupp
        u = k^n*(atan((x1 - 2*n*h)/(x2 + d))...
            + atan((x1 + 2*n*h)/(x2 + d)));
        u_S = last + u;
        if max(abs(u_S - last)) < tol
            u_S = 2*u_S*s/(2*pi) + 2*s/(2*pi)*atan(x1/(x2 + d));
            break
        end
        last = u_S;
    end
    plot(x1/d,u_S/s,'Color',c,'LineStyle',l,'LineWidth',1.5)
    clear u_S
    last = zeros(1, length(x1));
    hold off
end

x1 = linspace(h,d,100);
% u^(3) calculation
for i = 1:length(mu_star)
    k = kappa(i);
    c = colors(i);
    l = ls(i);
    hold on
    for n = 0:nupp
        u = k^n*atan2(x1 + 2*n*h, x2 + d);
        u_S = last + u;
        if max(abs(u_S - last)) < tol
            u_S = 2*(1 - k)*s/(2*pi)*u_S;
            break
        end
        last = u_S;
    end
    plot(x1/d,u_S/s,'Color',c,'LineStyle',l,'LineWidth',1.5)
    clear u_S
    last = zeros(1, length(x1));
    hold off
end

xlabel("Distance (x/d)")
ylabel("Surface displacement (u_{3}/s)")
legend("$\mu*/\mu = 1$","$\mu*/\mu = 0.8$","$\mu*/\mu = 0.5$",...
    "Location","se","Interpreter","Latex")
box on

subplot(2,1,2)
x1 = linspace(-d,-h,100);
% strain^(1) calculation
for i = 1:length(mu_star)
    k = kappa(i);
    c = colors(i);
    l = ls(i);
    hold on
    for n = 0:nupp
        strain = k^n*(x2 + d)./((x2 + d)^2 + (x1 - 2*n*h).^2);
        strain_S = last + strain;
        if max(abs(strain_S - last)) < tol
            strain_S = 2*(1 - k)*s/(2*pi)*strain_S;
            break
        end
        last = strain_S;
    end
    plot(x1/d,strain_S/s*d,'Color',c,'LineStyle',l,'LineWidth',1.5)
    clear strain_S
    last = zeros(1, length(x1));
    hold off
end

x1 = linspace(-h,h,100);
% strain^(2) calculation
for i = 1:length(mu_star)
    k = kappa(i);
    c = colors(i);
    l = ls(i);
    hold on
    for n = 1:nupp
        strain = k^n*((x2 + d)./((x2 + d)^2 + (x1 - 2*n*h).^2)...
            + (x2 + d)./((x2 + d)^2 + (x1 + 2*n*h).^2));
        strain_S = last + strain;
        if max(abs(strain_S - last)) < tol
            strain_S = 2*strain_S*s/(2*pi) + ...
                2*s/(2*pi)*(x2 + d)./((x2 + d)^2 + (x1).^2);
            break
        end
        last = strain_S;
    end
    plot(x1/d,strain_S/s*d,'Color',c,'LineStyle',l,'LineWidth',1.5)
    clear strain_S
    last = zeros(1, length(x1));
    hold off
end

x1 = linspace(h,d,100);
% strain^(3) calculation
for i = 1:length(mu_star)
    k = kappa(i);
    c = colors(i);
    l = ls(i);
    hold on
    for n = 0:nupp
        strain = k^n*(x2 + d)./((x2 + d)^2 + (x1 + 2*n*h).^2);
        strain_S = last + strain;
        if max(abs(strain_S - last)) < tol
            strain_S = 2*(1 - k)*s/(2*pi)*strain_S;
            break
        end
        last = strain_S;
    end
    plot(x1/d,strain_S/s*d,'Color',c,'LineStyle',l,'LineWidth',1.5)
    clear strain_S
    last = zeros(1, length(x1));
    hold off
end

xlabel("Distance (x/d)")
ylabel("Surface strain (\epsilon_{13}*d/s)")
legend("$\mu*/\mu = 1$","$\mu*/\mu = 0.8$","$\mu*/\mu = 0.5$",...
    "Location","ne","Interpreter","Latex")
box on
sgtitle("Figure 5.6")