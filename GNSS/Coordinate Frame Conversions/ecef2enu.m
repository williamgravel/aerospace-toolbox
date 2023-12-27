function [enu] = ecef2enu(ecef)
arguments
    ecef (:,3) double
end

enu = rot_enu(ecef)*ecef;

end
