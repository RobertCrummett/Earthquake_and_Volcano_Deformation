function [strike_slip, dip_slip, opening] = Point_Source(potency, x1, x2,...
    depth, dip, nu)
% Point_Source() calculates the displacement due to a point source in a 
% half space.
%
% Inputs: 
% potency - (slip_{i} * dArea), where {i} is the direction of choice
%           In the Okada coordinate system (Figure 3.20 pg 78), {1} is the
%           strike slip direction, {2} is the dip-slip direction, and {3}
%           corresponds to the opening direction
% x1 - x_{1} range
% x2 - x_{2} range
% depth - The depth to the base of the fault (Figure 3.20, pg 78)
% dip - Dip of fault in degrees from the positive horizontal
% nu - Poissons constant, default is 0.25
%
% Outputs:
% strike_slip - Structure contianing dislocations corresponding to slip in 
%               the {1} direction
% dip_slip - Structure contianing dislocations corresponding to slip in the
%            {2} direction
% opening - Structure contianing dislocations corresponding to opening 
%               in the {3} direction
%
% Source:
% Earthquake and Volcano Deformation, Paul Segall (2010)
% Chapter 3, pg 78-80

if ~exist("nu","var")
    nu = 0.25;
    disp('Default nu = 0.25 used')
end

dip = deg2rad(dip);
d = depth; clear depth

p = x2*cos(dip) + d*sin(dip);
q = x2*sin(dip) - d*cos(dip);
r = sqrt(x1.^2 + p.^2 + q.^2);

I0_1 = (1-2*nu)*x2.*((1./(r.*(r + d).^2)) - x1.^2.*(3*r+d)./(r.^3.*(r+d).^3));
I0_2 = (1-2*nu)*x1.*((1./(r.*(r + d).^2)) - x2.^2.*(3*r+d)./(r.^3.*(r+d).^3));
I0_3 = (1-2*nu)*(x1./(r.^3)) - I0_2;
I0_4 = -(1-2*nu)*(x1.*x2.*(2*r+d)./(r.^3.*(r+d).^2));
I0_5 = (1-2*nu)*((1./(r.*(r + d))) - x1.^2.*(2*r+d)./(r.^3.*(r+d).^2));

strike_slip.u1 = -potency/(2*pi)*((3*x1.^2.*q)./(r.^5) + I0_1*sin(dip));
strike_slip.u2 = -potency/(2*pi)*((3*x1.*x2.*q)./(r.^5) + I0_2*sin(dip));
strike_slip.u3 = -potency/(2*pi)*((3*x1.*d.*q)./(r.^5) + I0_4*sin(dip));

dip_slip.u1 = -potency/(2*pi)*((3*x1.*p.*q)./(r.^5) - I0_3*sin(dip)*cos(dip));
dip_slip.u2 = -potency/(2*pi)*((3*x2.*p.*q)./(r.^5) - I0_1*sin(dip)*cos(dip));
dip_slip.u3 = -potency/(2*pi)*((3*d.*p.*q)./(r.^5) - I0_5*sin(dip)*cos(dip));

opening.u1 = potency/(2*pi)*((3*x1.*q.^2)./(r.^5) - I0_3*sin(dip)^2);
opening.u2 = potency/(2*pi)*((3*x2.*q.^2)./(r.^5) - I0_1*sin(dip)^2);
opening.u3 = potency/(2*pi)*((3*d*q.^2)./(r.^5) - I0_5*sin(dip)^2);

end