# beacon-modeling
Mathematical model for complex molecular beacon thermodynamics based on fluorescent data

This is a mathematical model for the binding of targets to molecular beacons. The model uses fluorescent melting data to calculate the dissociation constants for the beacon (B) to the targets (Y). Three code files for different types of data along with three example data files are provided.

Code Files:

OneSite_Fit_KDBY
	This MATLAB code file fits fluorescent traces to a one-site model to find the dissociation constant for 	one target binding (KDBY).
  
TwoSite_Fit_KDBYY
	This MATLAB code file fits fluorescent traces to a two-site model to find the dissociation constant for 	a second target binding (KDBYY). KDBY parameters can be estimated using the OneSite model.
  
TwoSite_Fit_KDBY_KDYY
	This MATLAB code file fits fluorescent traces to a two-site model to find the dissociation constant for 	both targets binding (KDBY and KDBYY). 


Data File Format:

xlsx data files require the following sheets:

x = number of target concentrations
y = temperature points

1. NormFluoro - raw data sheet with columns Temperature, NTC, [Y] = 0, data1, data2,...,datax
2. weights - matrix of size (x,y)
3. conc - row with target concentrations 
4. KD - initial dissociation constants in the following order:
	One-site code: K_DB, K_DBY, K_DYY
	Two-site code: K_DB, K_DBY, K_DBYY, K_DYY

	*Ensure that all KDs are in over the same temperature range. 	

5. Y - free target inputs of size (x,y)
6. B - free beacon inputs of size (x,y)

Example data files are provided for clarity.


