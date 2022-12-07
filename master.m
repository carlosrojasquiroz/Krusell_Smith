clc,close,clear;
%---------------------------------------------------------------------------------------------------------------------------
%% Aiyagari model with aggregate shocks
%---------------------------------------------------------------------------------------------------------------------------
% This code solves the Aiyagari model with aggregate shocks by the K-S method.
% (c) Carlos Rojas Quiroz - December 2022
%---------------------------------------------------------------------------------------------------------------------------
%% 1. Parameters
%---------------------------------------------------------------------------------------------------------------------------
p=parameters();
%---------------------------------------------------------------------------------------------------------------------------
%% 2. Matrices (grid specifications, stationary distribution, etc.)
%---------------------------------------------------------------------------------------------------------------------------
m=matrices(p);
%---------------------------------------------------------------------------------------------------------------------------
%% 3. Demand of capital and labor (from firms)
%---------------------------------------------------------------------------------------------------------------------------
[p.K,p.L]=firms(p,m);
%---------------------------------------------------------------------------------------------------------------------------
%% 4. Partial equilibrium
%---------------------------------------------------------------------------------------------------------------------------
% 4.1. Policy and value functions
%---------------------------------------------------------------------------------------------------------------------------
[pe.V,pe.g_a,pe.g_c]=egm(p,m);
%---------------------------------------------------------------------------------------------------------------------------
% 4.2 Euler errors
%---------------------------------------------------------------------------------------------------------------------------
pe.EE=eulererrors(pe.g_a,p,m);
%---------------------------------------------------------------------------------------------------------------------------
% 4.3 Transition matrix
%---------------------------------------------------------------------------------------------------------------------------
pe.tm=transitionmatrix(p,m,pe.g_a,p.intp);
%---------------------------------------------------------------------------------------------------------------------------
% 4.4 Assets distribution
%---------------------------------------------------------------------------------------------------------------------------
pe.dist=stationarydist(p,pe.tm,p.intp);
%---------------------------------------------------------------------------------------------------------------------------
%% 5. General equilibrium
%---------------------------------------------------------------------------------------------------------------------------
ge=solveAiyagari(p,m);
%---------------------------------------------------------------------------------------------------------------------------
% 5.1 Market for capital
%---------------------------------------------------------------------------------------------------------------------------
[ge.Kd_curve,ge.Ks_curve]=capitalmarket(p,m,ge.r,p.naa/2,1);