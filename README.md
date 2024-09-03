# SunRISE-LWA
Part I. Solar Eclipse data analysis
Use the code to plot light curves for the desired frequencies to be analyzed.
-Plot_light_curve_filtering_butterworth.py

Part II. Image processing
(a)Use standardization & normalization and then manual threshold filter to remove noise patches and vertical lines, however, this might also remove parts of the signal.
-Multi_FIT_Newkrik_processed_format.py
(b)In the cases where the SRB signature is split into two or more fit files, use this code to concatenate them and do analysis.
-Multi_FIT_Newkrik_concatenate.py

Part III. Type II SRBs and CMEs analysis
(a)Use this code to do CME values estimation. The code uses the frequencies to get the electron densities, and then use the Newkirk model to derive the shock height, and then use some differential equations to estimate parameters such as drift rates and speeds. To use this code, plot your spectrogram and then left-click on the band you want to analyze to place dots, where the curve connecting the dots should follow the general shape of the band of the type II signal in order to reflect the trend of the drift rate variation. After placing the dots, right click to plot the characteristic curves.
-Multi_FIT_Newkrik_processed_format.py
(b)Use this Matlab code for CME characteristics calculation, including drift rates, shock speeds, Alfvén velocity, ambient magnetic field etc. To get the values in the equations, still use the dynamic spectrogram to place dots, and this time remember to place dots across the entire fundamental band, and calculate the average drift rate. Take the first point of shock height, the starting frequency of the fundamental band and the average drift rate to calculate the shock speed. For other parameters and their calculation, refer to the paper in the next link.
-CME_chars_calculation.m
Refer to this paper for the equations used:
ANGEO - Low-frequency solar radio type II bursts and their association with space weather events during the ascending phase of solar cycle 25 (copernicus.org)
