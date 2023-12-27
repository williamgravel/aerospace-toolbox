function [ecef] = llh2ecef(llh,opts)
arguments
    llh (:,3) double
    opts.ReferenceFrame (1,1) string = 'GPS'
end

phi = llh(1);
lambda = llh(2);
h = llh(3);

params = TRF(opts.ReferenceFrame);

N = params.a/sqrt(1 - params.e^2*sind(phi)^2);

x = (N + h)*cosd(phi)*cosd(lambda);
y = (N + h)*cosd(phi)*sind(lambda);
z = (N*(1 - params.e^2) + h)*sind(phi);

ecef = [x, y, z];

end
