# CGMC-XYZmodel

This code can be used to simulate quantum trajectories for the driven-dissipative XYZ Heisenberg model using a cluster Gutzwiller ansatz for the wave function. 
More information on this model can be found in arXiv:1812.00643.

An example code can be found in TESTFILE.m for a specific set of system parameters.

The code can be adapted to work for other models. The most important files to be changed would then be H1.m, H2.m and the used wave function ansatz (which in the present form of the code can be found in the file CalculateSssxx.m on line 21). Depending on the model new operators may have to be defined, this can be done in GetAllOperators.m.
