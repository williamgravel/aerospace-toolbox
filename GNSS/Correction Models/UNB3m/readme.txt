----|---1|0---|---2|0---|---3|0---|---4|0---|---5|0---|---6|0---|---7|0-
========================================================================

UNB3m_pack.zip
==============

Contents
--------

readme.txt              -> This file
UNB3M.f                 -> Subroutine for unb3m - FORTRAN
UNB3M.m                 -> Function for unb3m - MatLab
UNB3MR.f                -> Subroutine for unb3mr - FORTRAN
UNB3MR.m                -> Function for unb3mr - MatLab
UNB3MM.f                -> Subroutine for unb3mm - FORTRAN
UNB3MM.m                -> Function for unb3mm - MatLab
TUNB3M.f                -> UNB3m_pack test program - FORTRAN
TUNB3M.m                -> UNB3m_pack test function - MatLab 
tunb3m_.txt             -> Test program output file
Leandroetal2006.pdf     -> Reference paper for unb3m


Test procedure
--------------

FORTRAN users: Compile TUNB3M program with the other subroutines (UNB3M, 
 UNB3MM, UNB3MR). Run the executable. Compare the file created 
 (tunb3m.txt) with the one provided with the package (tunb3m_.txt).

MatLab users: Put all *.m files in the same directory. Run TUNB3M 
 function. Compare the file created (tunb3m.txt) with the one provided
 with the package (tunb3m_.txt).
 
========================================================================
----|---1|0---|---2|0---|---3|0---|---4|0---|---5|0---|---6|0---|---7|0-

