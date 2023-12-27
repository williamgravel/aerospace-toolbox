function [llh] = ecef2llh(ecef,opts)
arguments
    ecef (:,3) double
    opts.ReferenceFrame (1,1) string = 'gps'
    opts.AngleUnits (1,1) string {mustBeMember(opts.AngleUnits,{'rad','deg'})} = 'deg'
end

x = ecef(:,1);
y = ecef(:,2);
z = ecef(:,3);

params = TRF(opts.ReferenceFrame);

w = vecnorm([x,y],2,2);

l = params.e^2/2;
m = (w./params.a).^2;
n = ((1 - params.e^2)*z/params.b).^2;
p = (m + n - 4*l^2)/6;

G = m.*n*l^2;
H = 2*p.^3 + G;

% if H < H_min then abort

C = ((H + G + 2*sqrt(H.*G))/2).^(1/3);

i = -(2*l^2 + m + n)/2;
P = p.^2;

beta = i/3 - C - P./C;

k = l^2*(l^2 - m - n);
t = sqrt(sqrt(beta.^2 - k) - (beta + i)/2) - sign(m - n).*sqrt(abs((beta - i)/2));

F = t.^4 + 2*i.*t.^2 + 2*l*(m - n).*t + k;
dFdt = 4*t.^3 + 4*i.*t + 2*l*(m - n);
Delta_t = -F./(dFdt);

u = t + Delta_t + l;
v = t + Delta_t - l;

phi = atan2d(z.*u,w.*v);

Delta_w = w.*(1 - 1./u);
Delta_z = z.*(1 - (1 - params.e^2)./v);

h = sign(u - l).*vecnorm([Delta_w,Delta_z],2,2);
lambda = atan2d(y,x);

llh = [phi, lambda, h];

if strcmp(opts.AngleUnits,'rad')
    llh(1:2) = deg2rad(llh(1:2));
end

end
