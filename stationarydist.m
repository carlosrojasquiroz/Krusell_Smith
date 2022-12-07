function mu=stationarydist(p,Pi,intp)
%---------------------------------------------------------------------------------------------------------------------------
% This function obtains the stationary distribution of assets
%---------------------------------------------------------------------------------------------------------------------------
mu0=(1/(p.nzz*intp))*ones(1,p.nzz*intp); 
err=1;
iter=0;
while err> p.tol && iter<p.maxiter                
   mu1=mu0*Pi;
   err=max(abs(mu1-mu0));
   mu0=mu1;
   iter=iter+1;
end
mu=reshape(mu0,p.nzz,intp);

if p.fig==1
    LW=1.5;
    FS=16;
    a=gridspecA(p,intp);
    disT=sum(mu);
%---------------------------------------------------------------------------------------------------------------------------
% Assets distribution
%---------------------------------------------------------------------------------------------------------------------------
    figure()
    plot(a,disT, 'LineWidth',LW+1)
    xlabel('Assets, $a$','fontsize',FS,'interpreter','latex')
    ylabel('Fraction of households','fontsize',FS,'interpreter','latex')
    grid on
    ax=gca;
    ax.FontSize =FS;
    saveas(gcf,'Histogram','epsc')
end