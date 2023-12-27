%
function [T P E TM]=UNB3MM(LATRAD,HEIGHTM,DAYOYEAR)
%
%---|---1|0---|---2|0---|---3|0---|---4|0---|---5|0---|---6|0---|---7|0-
%=======================================================================
%
%  Subroutine UNB3mm
%  =================
%
%     This subroutine uses UNB3m model to calculate the surface 
%     temperature, surface pressure, surface relative humidity and
%     mean temperature of water vapor for a given latitude, height 
%     and day of year (DOY).
%
%   How this function works:
%   ------------------------
%
%     1. A look-up table is used to calculate mean sea level (MSL)
%        values for pressure, temperature, relative humidity (RH) 
%        and lapse rate parameters (for given latitude and DOY).
%
%     2. MSL water vapor pressure (WVP) is computed according to
%        IERS Conventions 2003.
%
%     3. Pressure, temperature and WVP values are computed for 
%        station height.
%
%     4. Mean temperature of water vapor is computed.
%
%   General Comments:
%   -----------------
%
%     1. This subroutine was adapted from subroutine UNB3M.
%
%     2. The internal variables for RH are called E0, EAVG and EAMP.
%
%   References: 
%   ----------
%
%     Leandro R.F., Santos, M.C. and Langley R.B., (2006). 
%     UNB Neutral Atmosphere Models: Development and Performance. 
%     Proceedings of ION NTM 2006, pp. 564-573, Monterey, California, 
%     January, 2006.
%
%     Davis, J.L., T.A. Herring, I.I. Shapiro, A.E.E. Rogers, 
%     and G. Elgered (1985). Geodesy by Radio Interfereometry: 
%     Effects of Atmospheric Modeling Errors on Estimates of 
%     Baseline Length. Radio Science, Vol. 20, No. 6, pp. 1593-1607.
%
%     IERS (2004). IERS Conventions (2003), eds. D.D. McCarthy and 
%     G. Petit, IERS Technical Note No. 32, International Earth Rotation 
%     and Reference Systems Service, Verlag des Bundesmates für Kartographie 
%     und Geodasie, Frankfurt am Main.
%
%   Created/Modified
%   ----------------
%
%	DATE        WHO                WHAT
%	----        ---                ----	
%   2006/01        Rodrigo Leandro          SUBROUTINE UNB3M.m created
%   2006/06/22     Rodrigo Leandro          SUBROUTINE UNB3MM.m created 
%                                           (addapted from UNB3M.m)
%   2006/11/23     Rodrigo Leandro          Code re-formatted. Capital case
%                                           should be always used in code 
%                                           from now on.
%
%   INPUT/OUTPUT
%   ------------
%
%     PARAMETER   I/O   TYPE       DESCRIPTION
%     ---------   ---   ----       -----------
%     LATRAD      IN    DOUBLE     Station geodetic latitude (radians)
%     HEIGHTM     IN    DOUBLE     Station orthometric height (m)
%     DAYOYEAR    IN    DOUBLE     Day of year
%     T           OUT   DOUBLE     Surface temperature (k)
%     P           OUT   DOUBLE     Surface pressure (mbar)
%     E           OUT   DOUBLE     Surface water vapor pressure (mbar)
%     TM          OUT   DOUBLE     Mean temperature of water vapor (k)
%
%=======================================================================
%
%
%=======================================================================
% Initialize UNB3m look-up table
%-----------------------------------------------------------------------
AVG=[  15.0  1013.25  299.65  75.00  6.30  2.77
       30.0  1017.25  294.15  80.00  6.05  3.15
       45.0  1015.75  283.15  76.00  5.58  2.57
       60.0  1011.75  272.15  77.50  5.39  1.81
       75.0  1013.00  263.65  82.50  4.53  1.55];
%      lat      P       T      RH    beta lambda
%
AMP=[  15.0   0.00   0.00   0.00  0.00  0.00
       30.0  -3.75   7.00   0.00  0.25  0.33
       45.0  -2.25  11.00  -1.00  0.32  0.46
       60.0  -1.75  15.00  -2.50  0.81  0.74
       75.0  -0.50  14.50   2.50  0.62  0.30];
%      lat     P      T      RH    beta lambda
%=======================================================================
%
EXCEN2 = 6.6943799901413e-03;
MD     = 28.9644;
MW     = 18.0152;
K1     = 77.604;
K2     = 64.79;
K3     = 3.776e5;
R      = 8314.34;
C1     = 2.2768e-03;
K2PRIM = K2 - K1*(MW/MD);
RD     = R / MD;
DTR = 1.745329251994329e-02;
DOY2RAD=(0.31415926535897935601e01)*2/365.25;

%---|---1|0---|---2|0---|---3|0---|---4|0---|---5|0---|---6|0---|---7|0-

%
%
%=======================================================================
% Transform latitude from radians to decimal degrees
%-----------------------------------------------------------------------
LATDEG = LATRAD * 180.0 / pi;
%=======================================================================
%
%
%=======================================================================
% Deal with southern hemisphere and yearly variation
%-----------------------------------------------------------------------
TD_O_Y = DAYOYEAR;
if LATDEG<0
    TD_O_Y = TD_O_Y + 182.625;
end
COSPHS = cos((TD_O_Y - 28) * DOY2RAD );
%=======================================================================
%
%
%=======================================================================
% Initialize pointers to lookup table
%-----------------------------------------------------------------------
LAT = abs( LATDEG );

if LAT>=75
    P1 = 5;
    P2 = 5;
    M = 0;
elseif LAT<=15
    P1 = 1;
    P2 = 1;
    M = 0;
else
    P1 = fix((LAT - 15)/15.D0) + 1;
    P2 = P1 + 1;
    M = (LAT - AVG(P1,1) ) / ( AVG(P2,1) - AVG(P1,1) );
end
%=======================================================================
%
%
%=======================================================================
% Compute average surface tropo values by interpolation
%-----------------------------------------------------------------------
PAVG = M * ( AVG(P2,2) - AVG(P1,2) ) + AVG(P1,2);
TAVG = M * ( AVG(P2,3) - AVG(P1,3) ) + AVG(P1,3);
EAVG = M * ( AVG(P2,4) - AVG(P1,4) ) + AVG(P1,4);
BETAAVG   = M * ( AVG(P2,5) - AVG(P1,5) ) + AVG(P1,5);
LAMBDAAVG = M * ( AVG(P2,6) - AVG(P1,6) ) + AVG(P1,6);
%=======================================================================
%
%
%=======================================================================
% Compute variation of average surface tropo values
%-----------------------------------------------------------------------
PAMP = M * ( AMP(P2,2) - AMP(P1,2) ) + AMP(P1,2);
TAMP = M * ( AMP(P2,3) - AMP(P1,3) ) + AMP(P1,3);
EAMP = M * ( AMP(P2,4) - AMP(P1,4) ) + AMP(P1,4);
BETAAMP   = M * ( AMP(P2,5) - AMP(P1,5) ) + AMP(P1,5);
LAMBDAAMP = M * ( AMP(P2,6) - AMP(P1,6) ) + AMP(P1,6);
%=======================================================================
%
%
%=======================================================================
% Compute surface tropo values
%-----------------------------------------------------------------------
P0 = PAVG - PAMP * COSPHS;
T0 = TAVG - TAMP * COSPHS;
E0 = EAVG - EAMP * COSPHS;
BETA = BETAAVG - BETAAMP * COSPHS;
BETA   = BETA / 1000.D0;
LAMBDA = LAMBDAAVG - LAMBDAAMP * COSPHS;
%=======================================================================
%
%
%=======================================================================
% Transform from relative humidity to WVP (IERS Conventions 2003)
%-----------------------------------------------------------------------
ES = 0.01 * exp(1.2378847e-5 * (T0 ^ 2) - 1.9121316e-2 * ...
    T0 + 3.393711047e1 - 6.3431645e3 * (T0 ^ -1));
FW = 1.00062 + 3.14e-6 * P0 + 5.6e-7 * ((T0 - 273.15) ^ 2);
E0 = (E0 / 1.00e2) * ES * FW;
%=======================================================================
%
%
%=======================================================================
% Compute power value for pressure & water vapour
%-----------------------------------------------------------------------
EP = 9.80665 / 287.054 / BETA;
%=======================================================================
%
%
%=======================================================================
% Scale surface values to required height
%-----------------------------------------------------------------------
T = T0 - BETA * HEIGHTM;
P = P0 * ( T / T0 ) ^ EP;
E = E0 * ( T / T0 ) ^ ( EP * (LAMBDA+1) );
%=======================================================================
%
%
%=======================================================================
% Compute the acceleration at the mass center
% of a vertical column of the atmosphere
%-----------------------------------------------------------------------
GEOLAT = atan((1.0-EXCEN2)*tan(LATRAD));
DGREF = 1.0 - 2.66e-03*cos(2.0*GEOLAT) - 2.8e-07*HEIGHTM;
GM    = 9.784 * DGREF;
DEN   = ( LAMBDA + 1.0 ) * GM;
%=======================================================================
%
%
%=======================================================================
% Compute mean temperature of the water vapor
%-----------------------------------------------------------------------
TM  = T * (1 - BETA * RD / DEN);
%=======================================================================
%
%
%=======================================================================
% End of function
%=======================================================================