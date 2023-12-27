function [utc] = tai2utc(tai,ls)
ls.DateTAI = datetime(ls.Date,"TimeZone","UTC") + days(1) + ls.CumulativeAdjustment + seconds(9);

[Y,M,D] = ymd(tai);
[H,MI,S] = hms(tai);
tai_ls = datetime([Y,M,D,H,MI,S],"TimeZone","UTCLeapSeconds");

ls_initial = seconds(10);
ls_issued = ls.CumulativeAdjustment(find(ls.DateTAI <= tai,1,'last'));
ls_current = seconds(any(ls.DateTAI - seconds(1) == tai));

if isempty(ls_issued)
    ls_issued = seconds(0);
end

utc = tai_ls - ls_initial - ls_issued - ls_current;

end
