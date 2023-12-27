function [m_h,m_w] = mappingFunction(lat,lon,h,doy,el,a_h,a_w)

if nargin == 7
    [b_h,b_w,c_h,c_w] = mappingCoefficients(lat,lon,doy);
else
    NMF = load('Lookup Tables\NMF.mat');

    doy_p = (doy - 28)/365.25*2*pi; % day of year angle (annual amplitude) [rad]

    a_h = interp1(NMF.phi,NMF.a_0.a_h,lat) + interp1(NMF.phi,NMF.A.a_h,lat)*cos(doy_p);
    b_h = interp1(NMF.phi,NMF.a_0.b_h,lat) + interp1(NMF.phi,NMF.A.b_h,lat)*cos(doy_p);
    c_h = interp1(NMF.phi,NMF.a_0.c_h,lat) + interp1(NMF.phi,NMF.A.c_h,lat)*cos(doy_p);

    a_w = interp1(NMF.phi,NMF.a_0.a_w,lat);
    b_w = interp1(NMF.phi,NMF.a_0.b_w,lat);
    c_w = interp1(NMF.phi,NMF.a_0.c_w,lat);
end

% continued fraction definition
f_t = @(a,b,c) 1 + a/(1 + b/(1 + c));
f_b = @(el,a,b,c) sind(el) + a/(sind(el) + b/(sind(el) + c));
f = @(el,a,b,c) f_t(a,b,c)/f_b(el,a,b,c);

m_0_h = f(el,a_h,b_h,c_h);

% height correction from Niell [1996]
a_ht = 2.53e-5;
b_ht = 5.49e-3;
c_ht = 1.14e-3;

dm_dh = 1/sind(el) - f(el,a_ht,b_ht,c_ht);
m_h = m_0_h + dm_dh*h;

m_w = f(el,a_w,b_w,c_w);

end
