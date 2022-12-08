function aggregateNS(p,m)

Km_ts=zeros(p.simulT,1); % time series of the mean of capital distribution 
L=m.z_grid*m.mu; % constant labor supply
Z=m.Zsim; % aggregate productivity
r=zeros(p.simulT-1,1); % interest rate
w=zeros(p.simulT-1,1);  % wage
% Beginning-of-period capital distributions
dist_a_zLow=zeros(p.simulT,p.intp); % low idio. shock
dist_a_zHigh=zeros(p.simulT,p.intp); % high idio. shock

a_grid=gridspec(p.Amin,p.Amax,p.intp,1);


