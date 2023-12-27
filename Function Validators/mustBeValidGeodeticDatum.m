function mustBeValidGeodeticDatum(value)
    if ~any(strcmpi(value,{'WGS84','PZ-90',''}))
        eidType = 'mustBeValidGeodeticDatum:notValidGeodeticDatum';
        msgType = 'Input must be a valid geodetic datum.';
        throwAsCaller(MException(eidType,msgType))
    end
end
