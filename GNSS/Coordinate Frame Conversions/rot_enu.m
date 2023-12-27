function [R_L] = rot_enu(ecef,opts)
arguments
    ecef (:,3) double
    opts.ReferenceFrame (1,1) string = 'GPS'
end

llh = ecef2llh(ecef,'ReferenceFrame',opts.ReferenceFrame,'AngleUnits','deg');
R_L = rot(1,90 - llh(1),'Unit','deg')*rot(3,llh(2) + 90,'Unit','deg');

end
