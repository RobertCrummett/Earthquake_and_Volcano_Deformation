function [strike_slip,dip_slip,opening] = FiniteRectDisloc(slip,x1,x2,d,dip,...
    W,L,mu,lambda)
% FiniteRectDisloc() calculates the displacements in 3D due to a strike-
% slip fault, a dip-slip fault, and an opening
% See Figure 3.20, pg 78
%
% Inputs: 
% s - uniform amount fo slip or opening along faults surface
% x1 - x range
% x2 - y range
% d - The depth to the base of the fault (Figure 3.20, pg 78)
% dip - Dip of fault in degrees from the positive horizontal
% W - Down-dip width of the fault
% L - Horizontal span of the fault
% mu - Shear modulus, default is 35e9
% lambda - Lamé constant, default is 35e9
%
% Outputs:
% strike_slip - Structure contianing dislocations corresponding to slip in 
%               the {1} direction
% dip_slip - Structure contianing dislocations corresponding to slip in the
%            {2} direction
% opening - Structure contianing dislocations corresponding to opening 
%           in the {3} direction
%
% Source:
% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 3, pg 80-82

if ~exist("mu","var")
    mu = 35e9;
    disp("Default mu (shear modulus) = 35e9 used")
end
if ~exist("lambda","var")
    lambda = 35e9;
    disp("Default lambda (Lamé constant) = 35e9 used")
end
con = mu/(lambda+mu); clear lambda mu
dip = pi*dip/180; % degrees to radians

p = x2*cos(dip) + d*sin(dip);
q = x2*sin(dip) - d*cos(dip);

y_tilde = @(eta) eta*cos(dip) + q*sin(dip);
d_tilde = @(eta) eta*sin(dip) - q*cos(dip);
R = @(zeta,eta) sqrt(zeta.^2 + y_tilde(eta).^2 + d_tilde(eta).^2);
X = @(zeta) sqrt(zeta.^2 + q.^2);

% check if dip = pi/2 or -pi/2
if cos(dip) < 1e-14 && cos(dip) > -1e-14 % numeric zero
    I_5 = @(zeta,eta) -con*zeta.*sin(dip)./(R(zeta,eta) + d_tilde(eta));
    I_4 = @(zeta,eta) -con*q./(R(zeta,eta) + d_tilde(eta));
    I_3 = @(zeta,eta) con/2*(eta./(R(zeta,eta) + d_tilde(eta)) + y_tilde(eta)...
        .*q./((R(zeta,eta) + d_tilde(eta)).^2) - log(R(zeta,eta) + eta));
    I_2 = @(zeta,eta) con*(-log(R(zeta,eta) + eta)) - I_3(zeta,eta);
    I_1 = @(zeta,eta) -con/2*(zeta.*q./((R(zeta,eta) + d_tilde(eta)).^2));
else
    I_5 = @(zeta,eta) con*2/cos(dip)*atan2(eta*(X(zeta) + q*cos(dip)) + ...
        X(zeta).*(X(zeta) + R(zeta,eta))*sin(dip), zeta*(R(zeta,eta) + ...
        X(zeta))*cos(dip));
    I_4 = @(zeta,eta) con/cos(dip)*(log(R(zeta,eta) + d_tilde(eta)) - ...
        sin(dip)*log(R(zeta,eta) + eta));
    I_3 = @(zeta,eta) con*(y_tilde(eta)./(cos(dip)*(R(zeta,eta) + ...
        d_tilde(eta))) - log(R(zeta,eta) + eta)) + sin(dip)*I_4(zeta,eta)/...
        cos(dip);
    I_2 = @(zeta,eta) con*(-log(R(zeta,eta) + eta)) - I_3(zeta,eta);
    I_1 = @(zeta,eta) con*(-zeta/(cos(dip)*(R(zeta,eta)+d_tilde(eta)))) - ...
        sin(dip)*I_5(zeta,eta)/cos(dip);
end

ss_u1 = @(zeta,eta) -slip/(2*pi)*(zeta.*q./(R(zeta,eta).*(R(zeta,eta) + eta))...
    + atan2(zeta.*eta, q.*R(zeta,eta)) + I_1(zeta,eta)*sin(dip));
ss_u2 = @(zeta,eta) -slip/(2*pi)*(y_tilde(eta).*q./(R(zeta,eta).*(R(zeta,eta) + eta))...
    + q*cos(dip)./(R(zeta,eta) + eta) + I_2(zeta,eta)*sin(dip));
ss_u3 = @(zeta,eta) -slip/(2*pi)*(d_tilde(eta).*q./(R(zeta,eta).*(R(zeta,eta) + eta))...
    + q*sin(dip)./(R(zeta,eta) + eta) + I_4(zeta,eta)*sin(dip));

ds_u1 = @(zeta,eta) -slip/(2*pi)*(q./R(zeta,eta) - I_3(zeta,eta)*sin(dip)*cos(dip));
ds_u2 = @(zeta,eta) -slip/(2*pi)*(y_tilde(eta).*q./(R(zeta,eta).*...
    (R(zeta,eta) + zeta)) + cos(dip)*atan(zeta.*eta./(q.*R(zeta,eta))) -...
    I_1(zeta,eta)*sin(dip)*cos(dip));
ds_u3 = @(zeta,eta) -slip/(2*pi)*(d_tilde(eta).*q./(R(zeta,eta).*...
    (R(zeta,eta) + zeta)) + sin(dip)*atan(zeta.*eta./(q.*R(zeta,eta))) -...
    I_5(zeta,eta)*sin(dip)*cos(dip));

op_u1 = @(zeta,eta) slip/(2*pi)*(q.^2./(R(zeta,eta).*(R(zeta,eta) + eta))...
    - I_3(zeta,eta)*sin(dip)^2);
op_u2 = @(zeta,eta) slip/(2*pi)*(-q.*d_tilde(eta)./(R(zeta,eta).*...
    (R(zeta,eta) + zeta)) - sin(dip)*(zeta.*q./(R(zeta,eta).*...
    (R(zeta,eta) + eta)) - atan(zeta.*eta./(q.*R(zeta,eta)))) - ...
    I_1(zeta,eta)*sin(dip)^2);
op_u3 = @(zeta,eta) slip/(2*pi)*(q.*y_tilde(eta)./(R(zeta,eta).*...
    (R(zeta,eta) + zeta)) + cos(dip)*(zeta.*q./(R(zeta,eta).*...
    (R(zeta,eta) + eta)) - atan(zeta.*eta./(q.*R(zeta,eta)))) - ...
    I_5(zeta,eta)*sin(dip)^2);

% strike slip displacements
strike_slip.u1 = Chinnery(ss_u1,x1,p,W,L);
strike_slip.u2 = Chinnery(ss_u2,x1,p,W,L);
strike_slip.u3 = Chinnery(ss_u3,x1,p,W,L);
% dip slip displacements
dip_slip.u1 = Chinnery(ds_u1,x1,p,W,L);
dip_slip.u2 = Chinnery(ds_u2,x1,p,W,L);
dip_slip.u3 = Chinnery(ds_u3,x1,p,W,L);
% opening displacements
opening.u1 = Chinnery(op_u1,x1,p,W,L);
opening.u2 = Chinnery(op_u2,x1,p,W,L);
opening.u3 = Chinnery(op_u3,x1,p,W,L);
end

function f = Chinnery(f, x, p, W, L)
% Chinnery (1961) notation function
f = f(x,p) - f(x,p-W) - f(x-L,p) + f(x-L,p-W);
end



