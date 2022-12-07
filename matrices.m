function m=matrices(p)
%---------------------------------------------------------------------------------------------------------------------------
% This function houses vectors and matrices of the model in the structure "m" 
%---------------------------------------------------------------------------------------------------------------------------
m.a_grid=gridspec(p.Amin,p.Amax,p.naa,1); % assets grid points (log-spaced)
m.z_grid=[0.1 1.0]; % idiosyncratic productivity grid points 
m.Pi=[0.9 0.1; 0.1 0.9]; % idiosyncratic productivity transition matrix
m.mu=stationarydistZ(p,m.Pi); % stationary distribution of idiosyncratic shocks
m.Z_grid=[0.99 1.01]; % aggregate productivity grid points 
m.PZ=[0.5 0.5; 0.1 0.9]; % aggregate productivity transition matrix
m.K_grid=gridspec(p.Kmin,p.Kmax,p.nKK,0); % aggregate capital grid points (linearly-spaced)
m.Zsim=simulshocks(p,m); % p.simulT simulated aggregagte shocks