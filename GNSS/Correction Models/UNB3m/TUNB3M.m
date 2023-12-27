function TUNB3M()
%
%---|---1|0---|---2|0---|---3|0---|---4|0---|---5|0---|---6|0---|---7|0-
%=======================================================================
%
%  Program TUNB3m
%  ==============
%
%     This program calls UNB3m subroutines to compute neutral atmosphere
%     delays, delay rates and meteorlogical parameters for several Latitudes, 
%     Heights, days of year and elevation angles. Created for test purpouses. 
%     T stands for 'Test'.
%
%   Created/Modified
%   ----------------
%
%      DATE          WHO                 WHAT
%      ----          ---                 ----	
%      2006/01       Rodrigo Leandro     PROGRAM created
%      2006/12/08    Rodrigo Leandro     Program modified to be used for
%                                        unb3m, unb3mr and unb3mm simultaneously
%
%   INPUT/OUTPUT
%   ------------
%
%	This program does not have any input, and it outputs a file called 
%      'tunb3m.txt' containing all computed values.
%
%=======================================================================
%---|---1|0---|---2|0---|---3|0---|---4|0---|---5|0---|---6|0---|---7|0-
%
%
%
%=======================================================================
% Open tunb3m.txt file for output
%-----------------------------------------------------------------------
OUTF=fopen('tunb3m.txt','w');
%=======================================================================
%
%
%=======================================================================
% Run UNB3m and write results
%-----------------------------------------------------------------------
%
%
% =================
% Write file header
% -----------------
fprintf(OUTF, ...
    'LATITUDE    HEIGHT       DOY  E. ANGLE      TEMP ');
fprintf(OUTF, ...
    '    PRESS         WVP        TM     DELAY   D. RATE\n');
fprintf(OUTF, ...
    '   (DEG)       (M)               (DEG)       (K)');
fprintf(OUTF, ...
    '    (MBAR)      (MBAR)       (K)       (M)   (M/RAD)\n\n');
% =================
%
% =======================================================
% %ompute values for several latitudes, heights, doy, and 
% elevation angles and print then in the output file.
% -------------------------------------------------------
      for LAT = 0:30:90
       for HEI = 0:500:1000
        for DOY = 1:100:351
         for ELV = 30:30:90
          LATRAD = LAT * pi/180;
          ELVRAD= ELV * pi/180;
          [RTROP HZD HMF WZD WMF]=UNB3M(LATRAD,HEI,DOY,ELVRAD);
          [DRATE HZD DHMFDEL WZD DWMFDEL]=UNB3MR(LATRAD,HEI,DOY,ELVRAD);
          [T P E TM]=UNB3MM(LATRAD,HEI,DOY);
          fprintf(OUTF, ...
              '%8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f    %8.3f  %8.3f  %8.3f  %8.3f  \n', ...
              LAT,HEI,DOY,ELV,T,P,E,TM,RTROP,DRATE);
         end
        end
       end
      end
% =======================================================
%
%=======================================================================
%
%
%=======================================================================
% End of Program
%=======================================================================
fclose('all');