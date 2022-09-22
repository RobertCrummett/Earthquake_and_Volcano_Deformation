% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 4, pg 111 - 113
% Recreation of Figure 4.15
%
% This would not have been possible without the help of Dr. Paul Segall,
% who took the time to email me a script to produce the correct solution.
% Thank you very much for your time Dr. Segall.

clear; clf; close all
N = 2^10;
lambda = 200;
sd = 1e-4; % stress drop
a = 20; % crack half-length
x0 = 50; % crack center
G = 1; % shear modulus

x = linspace(-lambda/2, lambda/2, N+1);
x = x(1:end-1);

u = 2*sd*sqrt(a^2 - (x - x0).^2);
u = real(u); % slip

n = -0.5*N:(0.5*N-1);
STIFF = -pi*G*abs(n)/lambda; % stiffness in the transform domain

U = fft(u,N);
Uhat = fftshift(U); % shifting to align with STIFF

Shat = STIFF.*Uhat; % convolution of stiffness and slip
S = ifftshift(Shat); % unshifting to align with x
s = ifft(S);

figure(1)
plot(x,u/(2*sd*a),"m-","LineWidth",1.75)
hold on
plot(x,s/sd,"k-","LineWidth",1.75)
text(50,2,["Normalized","slip"],"HorizontalAlignment","center",...
    "FontSize",11)
text(28,9,["Normalized","stress"],"HorizontalAlignment","right",...
    "FontSize",11)
hold off
box on
grid on
set(gca,"GridLineStyle",":")
xlim([0,100])
ylim([-2,12])
xlabel("Distance","FontSize",11)
title("Figure 4.15","FontSize",13,"FontWeight","normal")

