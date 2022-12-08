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
%% 4. Individual Policy and value functions
%---------------------------------------------------------------------------------------------------------------------------
[sol.V,sol.g_a,sol.g_c]=vfi(p,m);