function mu0=stationarydistZ(p,matrix)
%---------------------------------------------------------------------------------------------------------------------------
% This function obtains the stationary distribution of productivity shocks. 
% It is used for computing the aggregate labor and the grid for capital.
%---------------------------------------------------------------------------------------------------------------------------
mu0=(1./p.nzz)*ones(1,p.nzz);
err=1;
iter=0;
while err>p.tol && iter<p.maxiter                
    mu1=mu0*matrix;
    err=max(abs(mu1-mu0)); 
    mu0=mu1;
    iter=iter+1;
end