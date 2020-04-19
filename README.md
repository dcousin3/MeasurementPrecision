# MeasurementPrecision
Measurement Precision toolkit for R

This repository contains source code and a package suitable for R 3.0 and above. It computes 
a few basic statistics and round them according to the measurement instrument precision (delta_x).

See :
Cousineau, D. (2020) How many decimals? Rounding descriptive and inferential statistics based on 
measurement precision, Journal of Mathematical Psychology. doi: 10.1016/j.jmp.2020.102362
for the formal mathematical derivations of the results and an explanation of the four scenarios considered.

You can install this library on you computer if the library devtools is installed with:

devtools::install_github("dcousin3/MeasurementPrecision")<br>
library(MeasurementPrecision)

Check the UserManual.pdf document for more on the functions.
