function m=matrices(p)
%---------------------------------------------------------------------------------------------------------------------------
% This function houses vectors and matrices of the model in the structure "m" 
%---------------------------------------------------------------------------------------------------------------------------
m.a_grid=gridspec(p.Amin,p.Amax,p.naa,1); % assets grid points (log-spaced)
m.z_grid=[0.1 1.0]; % idiosyncratic productivity grid points 
m.Pi=[0.9 0.1; 0.1 0.9]; % idiosyncratic productivity transition matrix
m.mu=stationarydist(p,p.nzz,m.Pi,1); % stationary distribution of idiosyncratic shocks
m.Z_grid=[0.99 1.01]; % aggregate productivity grid points 
m.PZ=[0.5 0.5; 0.1 0.9]; % aggregate productivity transition matrix
m.mU=stationarydist(p,p.nZZ,m.PZ,1); % stationary distribution of aggregate shocks
m.P=kron(m.PZ,m.Pi);  % matrix of transition probabilities
m.K_grid=gridspec(3,7,p.nKK,0); % aggregate capital grid points (linearly-spaced)
[m.Zsim,m.agshock]=simulshocks(p,m); % p.simulT # of simulated aggregate shocks
m.B=[0 1 0 1]; % Initial vector of coefficients B of the ALM
%---------------------------------------------------------------------------------------------------------------------------
% Initial assets function (a')
%---------------------------------------------------------------------------------------------------------------------------
m.g_a=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
for d_1=1:p.nzz
   for d_3=1:p.nZZ
      for d_4=1:p.nKK
         m.g_a(d_1,:,d_3,d_4)=m.a_grid;
      end
   end
end
%---------------------------------------------------------------------------------------------------------------------------
% Initial density on the grid for the non-stochastic simulation
%---------------------------------------------------------------------------------------------------------------------------
m.a_cross=zeros(p.nzz,p.naa); 
m.a_cross(1:p.nzz,round(p.A_ss/((p.Amax-p.Amin)/p.naa)))=1; 
%---------------------------------------------------------------------------------------------------------------------------
% Auxilary matrices of transition probabilities (for the next period)
%---------------------------------------------------------------------------------------------------------------------------
m.P_Ll=zeros(p.nzz,p.naa,p.nZZ,p.nKK); % for a low agg. + low idios. shocks
m.P_Lh=zeros(p.nzz,p.naa,p.nZZ,p.nKK); % for a low agg. + high idios. shocks
m.P_Hl=zeros(p.nzz,p.naa,p.nZZ,p.nKK); % for a high agg. + low idios. shocks
m.P_Hh=zeros(p.nzz,p.naa,p.nZZ,p.nKK); % for a high agg. + high idios. shocks
%---------------------------------------------------------------------------------------------------------------------------
m.P_Ll(1,:,1,:)=m.P(1,1)*ones(p.naa,p.nKK);
m.P_Ll(2,:,1,:)=m.P(2,1)*ones(p.naa,p.nKK);
m.P_Ll(1,:,2,:)=m.P(3,1)*ones(p.naa,p.nKK);
m.P_Ll(2,:,2,:)=m.P(4,1)*ones(p.naa,p.nKK);
%---------------------------------------------------------------------------------------------------------------------------
m.P_Lh(1,:,1,:)=m.P(1,2)*ones(p.naa,p.nKK);
m.P_Lh(2,:,1,:)=m.P(2,2)*ones(p.naa,p.nKK);
m.P_Lh(1,:,2,:)=m.P(3,2)*ones(p.naa,p.nKK);
m.P_Lh(2,:,2,:)=m.P(4,2)*ones(p.naa,p.nKK);
%---------------------------------------------------------------------------------------------------------------------------
m.P_Hl(1,:,1,:)=m.P(1,3)*ones(p.naa,p.nKK);
m.P_Hl(2,:,1,:)=m.P(2,3)*ones(p.naa,p.nKK);
m.P_Hl(1,:,2,:)=m.P(3,3)*ones(p.naa,p.nKK);
m.P_Hl(2,:,2,:)=m.P(4,3)*ones(p.naa,p.nKK);
%---------------------------------------------------------------------------------------------------------------------------
m.P_Hh(1,:,1,:)=m.P(1,4)*ones(p.naa,p.nKK);
m.P_Hh(2,:,1,:)=m.P(2,4)*ones(p.naa,p.nKK);
m.P_Hh(1,:,2,:)=m.P(3,4)*ones(p.naa,p.nKK);
m.P_Hh(2,:,2,:)=m.P(4,4)*ones(p.naa,p.nKK);
%---------------------------------------------------------------------------------------------------------------------------