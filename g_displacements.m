function [output_struct] = g_displacements(x, xi, mu, nu)
% Calculates displacements at (x_1, x_2) due to a line force acting at
% (0, xi_2)
%
% Inputs:
% x - A two component vector in the one-two space
% xi - A two component vector in the one-two space
%      If the first element is not zero (ie, xi_1 ~= 0), then the x_1
%      coordinate is shifted by xi_1
% mu - Shear modulus (optional)
% nu - Lame constant (optional)
%
% Outputs:
% g1_1 - Displacements in the one direction of a line force acting in the 
%        one direction.
% g1_2 - Displacements in the two direction of a line force acting in the 
%        one direction.
% g2_1 - Displacements in the one direction of a line force acting in the 
%        two direction.
% g2_2 - Displacements in the two direction of a line force acting in the 
%        two direction.
%
% Source:
% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 3, pg.60-63

% Default elastic parameter values
if ~exist("mu","var")
    mu = 35e9;
    disp("Default mu = 35e9 used")
end
if ~exist("nu","var")
    nu = 0.25;
    disp("Default nu = 0.25 used")
end

if numel(x) ~= numel(xi) % same size check
    xi = repmat(xi, size(x,1), 1);
    if numel(x) ~= numel(xi)
        error("Both x and xi must be full of inputs")
    end
end
if size(x,2) ~= 2 || size(xi,2) ~= 2 % correct length check
    error("x and xi must be two component vectors with real inputs")
end
if x ~= real(x) || xi ~= real(xi) || mu ~= real(mu) || nu ~= real(nu)
    error("All inputs must be real") % real check
end

x_1 = x(:,1);
x_2 = x(:,2);
xi_2 = xi(:,2);

for idx = 1:length(x_1)
    if xi(idx,1) ~= 0 % if the line force is not at x_1 = 0
        x_1(idx) = x_1(idx) - xi(idx,1); % shift x
    end
end

r_1 = sqrt(x_1.^2 + (x_2 - xi_2).^2);
r_2 = sqrt(x_1.^2 + (x_2 + xi_2).^2);
theta_2 = atan2(x_1, x_2 + xi_2);

% displacements
g1_1 = -1/(2*pi*mu*(1-nu))*((3 - 4*nu)*log(r_1)/4 + (8*nu^2 - 12*nu + ...
    5)*log(r_2)/4 + (x_2 - xi_2).^2./(4*r_1.^2) + ((3 - 4*nu)*...
    (x_2 + xi_2).^2 + 2*xi_2.*(x_2 + xi_2) - 2*xi_2.^2)./(4*r_2.^2) - ...
    xi_2.*x_2.*(x_2 + xi_2).^2./(r_2.^4));
g1_2 = 1/(2*pi*mu*(1-nu))*((1-2*nu)*(1-nu)*theta_2 + (x_2 - xi_2).*x_1./...
    (4*r_1.^2) + (3-4*nu)*(x_2-xi_2).*x_1./(4*r_2.^2) + xi_2.*x_2.*x_1.*...
    (x_2 + xi_2)./(r_2.^4));
g2_1 = 1/(2*pi*mu*(1-nu))*(-(1-2*nu)*(1-nu)*theta_2 + (x_2 - xi_2).*x_1...
    ./(4*r_1.^2) + (3-4*nu)*(x_2-xi_2).*x_1./(4*r_2.^2) + xi_2.*x_2.*...
    x_1.*(x_2 + xi_2)./(r_2.^4));
g2_2 = 1/(2*pi*mu*(1-nu))*(-(3 - 4*nu)*log(r_1)/4 - (8*nu^2 - 12*nu +...
    5)*log(r_2)/4 - (x_1).^2./(4*r_1.^2) + (2.*xi_2.*x_2 - (3-4*nu)*x_1...
    .^2)./(4*r_2.^2) - xi_2.*x_2.*x_1.^2./(r_2.^4));

% stresses
s1_11 = -1/(2*pi*(1-nu))*(x_1.^3/(r_1.^4) + x_1.*(x_1.^2 - 4*xi_2.*x_2 ...
    - 2*xi_2.^2)./(r_2.^4) + 8*xi_2.*x_2.*x_1.*(x_2 + xi_2).^2./(r_2.^6)...
    + (1-2*nu)/2*(x_1./(r_1.^2) + 3*x_1./(r_2.^2) - 4*x_2.*x_1.*(x_2 + ...
    xi_2)./(r_2.^4)));
s1_12 = -1/(2*pi*(1-nu))*((x_2 - xi_2).*x_1.^2./(r_1.^4) +(x_2 + xi_2).*...
    (2*xi_2.*x_2 + x_1.^2)./(r_2.^4) - 8*xi_2.*x_2.*x_1.^2.*(x_2 + xi_2)...
    ./(r_2.^6) + (1-2*nu)/2*((x_2 - xi_2)./(r_1.^2) + (3*x_2 + xi_2)./...
    (r_2.^2) - 4*x_2.*(x_2 + xi_2)./(r_2.^4)));
s1_22 = -1/(2*pi*(1-nu))*((x_2 - xi_2).^2.*x_1./(r_2.^4) -x_1.*(xi_2.^2 ...
    - x_2.^2 + 6*xi_2.*x_2)./(r_2.^4) + 8*xi_2.*x_2.*x_1.^3./(r_2.^6) -...
    (1-2*nu)/2*(x_1./(r_1.^2) - x_1./(r_2.^2) - 4*x_2.*x_1.*(x_2 + xi_2)...
    ./(r_2.^4)));
s2_11 = -1/(2*pi*(1-nu))*((x_2 - xi_2).*x_1.^2./(r_1.^4) + ((x_2 + xi_2)...
    .*(2*xi_2.^2 + x_1.^2) - 2*xi_2.*x_1.^2)./(r_2.^4) + 8*xi_2.*x_2.*...
    (x_2 + xi_2).*x_1.^2./(r_2.^6) + (1-2*nu)/2*(-(x_2 - xi_2)./(r_1.^2)...
    + (3*xi_2 + x_2)./(r_2.^2) + 4*x_2.*x_1.^2./(r_2.^4)));
s2_12 = -x_1/(2*pi*(1-nu)).*((x_2 - xi_2).^2./(r_1.^4) + (x_2.^2 - ...
    2*xi_2.*x_2 - xi_2.^2)./(r_2.^4) + 8*xi_2.*x_2.*(x_2 + xi_2).^2./...
    (r_2.^6) + (1-2*nu)/2*(1./(r_1.^2) - 1./(r_2.^2) + 4*x_2.*(x_2 + ...
    xi_2)./(r_2.^4)));
s2_22 = -1/(2*pi*(1-nu))*((x_2 - xi_2).^3./(r_1.^4) + (x_2 + xi_2).*...
    ((x_2 +xi_2).^2 + 2*xi_2.*x_2)./(r_2.^4) + 8*xi_2.*x_2.*(x_2 + xi_2)...
    .^2./(r_2.^6) + (1-2*nu)/2*((x_2 - xi_2)./(r_1.^2) + (3*x_2 + xi_2)...
    ./(r_2.^2) - 4*x_2.*x_1.^2./(r_2.^4)));

% output structure
output_struct.g1 = {g1_1, g1_2};
output_struct.g2 = {g2_1, g2_2};
output_struct.s1 = {s1_11, s1_12, s1_22};
output_struct.s2 = {s2_11, s2_12, s2_22};
end