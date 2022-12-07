function [Kd,Ks]=capitalmarket(p,m,eqR,nnr,Index)
%---------------------------------------------------------------------------------------------------------------------------
% This function graphs supply and demand curves in the market for capital
% iff Index=1
%---------------------------------------------------------------------------------------------------------------------------
lowhalfintv=linspace(0.001,eqR,nnr/2);
uphalfintv=linspace(eqR+lowhalfintv(2)-lowhalfintv(1),0.1,nnr/2);
rvals=[lowhalfintv uphalfintv];
Kd=(p.alpha*p.Z./(rvals + p.delta)).^(1/(1-p.alpha))*p.L;
Ks=zeros(1,nnr);
m.a_grid=gridspecA(p,p.naa);
p.disp1=0;
p.fig=0;
for ind_i=1:nnr
    p.r=rvals(ind_i);
    p.w=p.Z*(1-p.alpha)*(p.Z*p.alpha/(p.r+p.delta))^(p.alpha/(1-p.alpha));
    [~,g_a,~]=egm(p,m);
    tm=transitionmatrix(p,m,g_a,p.naa);
    dist=stationarydist(p,tm,p.naa);
    Ks(1,ind_i)=sum(dist*m.a_grid');
end
%---------------------------------------------------------------------------------------------------------------------------
% Assets policy function
%---------------------------------------------------------------------------------------------------------------------------
if Index==1
    figure()
    LW=1.5;
    FS=16;
    plot(Ks,rvals,'LineWidth',LW);
    hold on;
    grid on;
    plot(Kd,rvals,'LineWidth',LW);
    xlim([2 10])
    ylim([0 0.1])
    xlabel('$K$','fontsize',FS,'interpreter','latex')
    ylabel('$r$','fontsize',FS,'interpreter','latex')
    title('Equilibrium in the market for capital')
    legend('$Supply$','$Demand$','fontsize',FS,'interpreter','latex'...
        ,'Location','best')
    ax=gca;
    ax.FontSize=FS; 
    [xi,yi] = polyxpoly(Ks,rvals,Kd,rvals);
    scatter(xi,yi,'filled','MarkerFaceColor','k','HandleVisibility','off')
    saveas(gcf,'MarketforK','epsc')  
end