function mu=stationarydist(p,nn,Pi,intp)
%---------------------------------------------------------------------------------------------------------------------------
% This function obtains the stationary distribution of assets
%---------------------------------------------------------------------------------------------------------------------------
mu0=(1/(nn*intp))*ones(1,nn*intp); 
err=1;
iter=0;
while err> p.tol && iter<p.maxiter                
   mu1=mu0*Pi;
   err=max(abs(mu1-mu0));
   mu0=mu1;
   iter=iter+1;
end
mu=reshape(mu0,nn,intp);