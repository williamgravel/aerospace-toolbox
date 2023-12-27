function mustBeValidSatelliteSystem(value)
    if ~any(strcmpi(value,{'gps','glonass','galileo','beidou','qzss','irnss'}))
        eidType = 'mustBeValidSatelliteSystem:notValidSatelliteSystem';
        msgType = 'Input must be a valid satellite system.';
        throwAsCaller(MException(eidType,msgType))
    end
end
