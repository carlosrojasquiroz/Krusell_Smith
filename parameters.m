function p=parameters()
%---------------------------------------------------------------------------------------------------------------------------
% This function houses parameters of the model in the structure "p" 
% (you can call it whatever you want, but be careful of being coherent in the whole code)
%---------------------------------------------------------------------------------------------------------------------------
p.naa=100; % # of asset grid points
p.nKK=5; % # of aggregate capital grid points
p.Amin=0; % borrowing limit
p.Amax=30; % max level of assets
p.A_ss=5; % assets steady state 
p.nzz=2; % # of idiosyincratic productivity grid points
p.nZZ=2;% # of aggregate productivity grid points
p.sigma=1; % intertemporal elasticity of substitution
p.beta=0.96; % discount factor
p.alpha=0.33; % participation of capital
p.delta=0.05; % depreciation rate
p.r=0.03; % interest rate (initial)
p.w=1.3; % wage (initial)
p.N=1; % mass of households
p.Kmin=(p.alpha*0.99./(0.0220 + p.delta))^(1/(1-p.alpha))*0.55; % low level for K
p.Kmax=(p.alpha*1.01./(0.0220 + p.delta))^(1/(1-p.alpha))*0.55; % max level for K
%---------------------------------------------------------------------------------------------------------------------------
% Moreover, this script also includes parameters to configurate the algorithm,
% the display of information from the code, and plotting some figures
%---------------------------------------------------------------------------------------------------------------------------
p.tol=1e-6; % tolerance criterion
p.intp=500; % # of asset grid points for interpolation
p.maxiter=1000; % maximum # iterations in the algorithms
p.simulT=10000; % # of time periods for aggregate shocks simulation
p.algo=1; % algo=0, EGM iterates over the VF / algo=1, EGM iterates over the PF
p.disp1=1; % disp=1, displays the tolerance error and the time to execute the EGM
p.disp2=1; % disp=1, displays the tolerance error and the time to solve the model
p.fig=0; % fig=1, plots figures, 0 otherwise