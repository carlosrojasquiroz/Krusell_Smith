function A=gridspec(Amin,Amax,naa,gridform)
%---------------------------------------------------------------------------------------------------------------------------
% This function creates two different grid specifications for assets, given
% the borrowing constraint and the maximum level calibrated in parameters
% object. Note that Index is an indicator function such that:
% gridform=0, linearly-spaced grid
% gridform=1, log-spaced grid
%---------------------------------------------------------------------------------------------------------------------------
if gridform==0
    A=linspace(Amin,Amax,naa);
elseif gridform==1
    A=exp(linspace(log(Amin+1),log(Amax+1),naa))-1;
else
    disp('You should choose between p.gridform= 0 or 1 only')
    disp('p.gridform=0, linearly-spaced grid')
    disp('p.gridform=1, log-spaced grid')
    A=NaN(1,naa);
end