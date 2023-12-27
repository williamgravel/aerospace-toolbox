function mustBeValidTimeStandard(value)
    if ~any(strcmpi(value,{'gps','glonass','galileo','beidou','qzss','irnss','utc','tai'}))
        eidType = 'mustBeValidTimeScale:notValidTimeScale';
        msgType = 'Input must be a valid time scale.';
        throwAsCaller(MException(eidType,msgType))
    end
end
