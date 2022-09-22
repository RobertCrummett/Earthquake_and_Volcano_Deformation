% Code from Dr. Paul Segall, personal correspondence (9.21.2022)
% This is the solution to reproducing Figure 4.15 (very clever shifting)

% Test of fft stress computation for anti-plane strain crack
%
% Note that for the stress to accurately match the analytical solution (in
% this case stress inside crack = -1.00) the half-length must be small
% compared to computational domain. Otherwise cyclical replications of the
% crack contribute to the stress.


clear Uhat S

N = 2^10;
% set displacement to be zero outside crack and elliptical inside
Lambda = 200;  % size of domain
u = zeros(N,1);
sd = 1.0e-4;    % stress drop
a  = 20;        % crack half-length
x0 = 0;         % crack center
G = 1;          % shear modulus

%x-coordinates
x = linspace(-Lambda/2,Lambda/2,N+1);  
x = x(1:end-1);                     

for k=1:N
	if abs( x(k) - x0 ) < a
		u(k) = 2*sd*sqrt(  a^2 - ( x(k) - x0)^2);
	end
end

% Compute the stiffness in the transform domain
% Note: this  could be done outside a subroutine and passed in
	n = -0.5*N:1:0.5*N-1;
	Stiff = -pi*G*abs(n)/Lambda;

    t0 = clock;
% Given u(x) FFT the displacement
	Uhat = fft(u, N);

% rearrange the indicies so that n ranges from -N/2 to +N/2
    U = fftshift(Uhat)';

% Convolve stiffness with slip distribution
	S = Stiff .* U;

% Rearrange indicies so that n ranges from 0 to N-1
    Shat = ifftshift(S);
    
% Inverse FFT to compute stress
	stress = ifft(Shat, N);
 	etime(clock,t0)
    
% Note: could be written more compactly as 
%    stress = ifft(ifftshift(Stiff .* fftshift(fft(u, N))'),N);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(x, u/(2*sd*a), 'b', x, stress/sd, 'r', 'LineWidth',2)
xlabel('Distance','FontSize',14)
title('Slip and Stress for Crack via FFT','FontSize',16)
grid
h = legend('Normalized Slip', 'Normalized Stress');
h.FontSize = 12;