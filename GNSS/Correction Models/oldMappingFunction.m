function [outputArg1,outputArg2] = oldMappingFunction(inputArg1,inputArg2)
b_h = 0.0029; % mean b coefficient
c_0_h = 0.062; % mean c coefficient

if lat >= 0 % northern hemisphere
    c_10_h = 0.001;
    c_11_h = 0.005;
    psi = 0;
else % southern hemisphere
    c_10_h = 0.002;
    c_11_h = 0.007;
    psi = pi;
end

% note that lat is in radians here
c_h = c_0_h + ((cos((doy - 28)/365*2*pi + psi) + 1)*c_11_h/2 + c_10_h)*(1 - cos(lat));

a_0_h = 0;
A_h = 0;
i = 0;
for n = 0:9
    for m = 0:n
        i = i + 1;
        a_0_h = a_0_h + GMF.a_0_h(i)*P(n+1,m+1) + GMF.b_0_h(i)*Q(n+1,m+1);
        A_h = A_h + GMF.A_h(i)*P(n+1,m+1) + GMF.B_h(i)*Q(n+1,m+1);
    end
end
a_h = a_0_h + A_h*cos(doy/365.25*2*pi)*1e-5;

b_w = 0.00146;
c_w = 0.04391;

a_0_w = 0;
A_w = 0;
i = 0;
for n = 0:9
    for m = 0:n
        i = i + 1;
        a_0_w = a_0_w + GMF.a_0_w(i)*P(n+1,m+1) + GMF.b_0_w(i)*Q(n+1,m+1);
        A_w = A_w + GMF.A_w(i)*P(n+1,m+1) + GMF.B_w(i)*Q(n+1,m+1);
    end
end
a_w = a_0_w + A_w*cos(doy/365.25*2*pi)*1e-5;

end

