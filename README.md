# SITU
Codes and data files to run the Sea Ice Tracking Utility (SITU)

Usage of these files and algorithms is permitted without license or fee.

Users should cite: 

DeRepentigny, P., L.B. Tremblay, R. Newton, and S. Pfirman, (2016), 
Patterns of sea ice retreat in the transition to a seasonally ice-free Arctic. 
Journal of Climate, DOI: 10.1175/JCLI-D-15-0733.1. For the SITU system.

The codes have been written and are maintained by Patricia DeRepentigny, 
L. Bruno Tremblay and Robert Newton, of the University of Colorado,
McGill University and Columbia Univesity, respectively.  

Inquiries should be directed to Robert Newton: bnewton@ldeo.columbia.edu.

This repository contains the codes used to run simulations reported on in 
in the manuscript, Defining the 'ice shed' of the Arctic Ocean�s Last Ice 
Area and its future evolution. 
Robert Newton, Stephanie Pfirman, L. Bruno Tremblay, Patricia Derepentigny

Users interested in downloading and using the codes are encouraged to 
reach out to the corresponding author (Newton) to see whether there are
more recent versions that might be more suitable for their project(s). 

The 'main' script for ice tracking is "LagTracks...".  We've found that 
it works well to clone LagTracks and edit for each project. 

The LagTracks scripts require:

FAITPaths... : A script to set path statements for your environment.

north_x_y_lat_lon: an ascii data file with two columns that establishes
the EASE-25 grid. 

CESMmask.mat: a mask for use with CESM ouptut; Ocean = 1; Land = 0.

XY_lia.mat: a matlab data file; two columns with the coordinates of the
Last Ice Area in the EASE-25 grid.  

The data inputs are CESM output interpolated onto the EASE-25 grid. 
We used monthly average ice drift in the CESM 'u' and 'v' directions:
uice and vice.   

The SITU directory includes several other utilities that are not
strictly needed to get Lagrangian tracks with SITU, but which may
prove useful as you customize the scripts for your project.  

LIAMovies is a basic script to animate ice tracks.
YearMonth2JMonth, YearWeek2JWeek, and date2jd are utilities 
to translate into Julian month, week and date.  
BathyEASE.mat is the Arctic bathymetry translated into the EASE-25
grid.  It is useful for plotting.  
