function [T] = tropoModel(lat,lon,h,doy,el,opts)
arguments
    lat (1,1) double
    lon (1,1) double
    h (1,1) double
    doy (1,1) double
    el (1,1) double
    opts.MeteoModel (1,1) string = ""
    opts.DelayModel (1,1) string = ""
    opts.MappingFunction (1,1) string = ""
    opts.RefinedDecimals (1,1) logical = true
end

[~,ZHD,ZWD,a_h,a_w] = zenithModel(lat,lon,h,doy,"MeteoModel",opts.MeteoModel,"DelayModel",opts.DelayModel,"RefinedDecimals",true);

% switch opts.DelayModel
%     case "Saastamoinen"
%         if opts.RefinedDecimals % courtesy of Davis et al. [1985]
%             ZHD = 0.0022768*(1 + 0.00266*cosd(2*lat) + 0.00028e-3*h)*P_0;
%             ZWD = 0.0022768*(1255/T_0 + 0.05)*e_0;
%         else % original constants from Saastamoinen [1972]
%             ZHD = 0.002277*(1 + 0.0026*cosd(2*lat) + 0.00028e-3*h)*P_0;
%             ZWD = 0.002277*(1255/T_0 + 0.05)*e_0;
%         end
%         
%     case "Hopfield"
%         % orthometric heights H (above geoid)
%         h_d = 43;
%         h_w = 12;
% 
%         N_d = 77.6*P_0/T_0;
%         N_w = 77.6*4810*e_0/T_0^2;
% 
%         ZHD = 1e-6/5*N_d*h_d;
%         ZWD = 1e-6/5*N_w*h_w;
% end

switch opts.MappingFunction
    case "NMF"
        [m_h,m_w] = mappingFunction(lat,lon,h,doy,el);
    case "VMF3"
        [m_h,m_w] = mappingFunction(lat,lon,h,doy,el,a_h,a_w);
end

T = ZHD*m_h + ZWD*m_w;

% References
% Davis et al. [1985]: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/RS020i006p01593

end