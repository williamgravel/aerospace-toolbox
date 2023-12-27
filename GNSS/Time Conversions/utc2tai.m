function [tai] = utc2tai(utc,ls)
ls.DateUTC = datetime(ls.Date,"TimeZone","UTC") + days(1);

utc_ct = datetime(utc,"TimeZone","UTC");

ls_initial = seconds(10);
ls_issued = ls.CumulativeAdjustment(find(ls.DateUTC <= utc_ct,1,"last"));
ls_current = seconds(second(utc) == 60);

if isempty(ls_issued)
    ls_issued = seconds(0);
end

tai = utc_ct + ls_initial + ls_issued + ls_current;

end
