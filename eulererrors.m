function EE=eulererrors(g_a,p,m)
%---------------------------------------------------------------------------------------------------------------------------
% This function computes Euler errors from the interpolation of the assets policy function. It
% interpolates the actual object m.g_a using p.intp grid points.
%---------------------------------------------------------------------------------------------------------------------------
a_intp=gridspecA(p,p.intp);
EE=zeros(p.nzz,p.intp);
for ind_j=1:p.nzz                   
    for ind_i=1:p.intp
    ap=interp1(m.a_grid,g_a(ind_j,:),a_intp(ind_i),'linear');
    c=m.z_grid(ind_j)*p.w+(1+p.r)*a_intp(ind_i)-ap;
    uc=c.^(-p.sigma);
    Euc=0;
    for ind_j1=1:p.nzz
        app=interp1(m.a_grid,g_a(ind_j1,:),ap,'linear');
        cp=m.z_grid(ind_j1)*p.w+(1+p.r)*ap-app;
        Euc=m.Pi(ind_j,ind_j1)*cp.^(-p.sigma)+Euc;
    end
    EE(ind_j,ind_i)=(1+p.r)*p.beta*Euc/uc-1;
    end
end

if p.fig==1
    LW=1.5;
    FS=16;
    zlow=1;
    zhigh=2;
%---------------------------------------------------------------------------------------------------------------------------
% Euler errors
%---------------------------------------------------------------------------------------------------------------------------
    figure()
    plot(a_intp, log(abs(EE(zlow,:))), 'LineWidth',LW)
    hold on;
    grid on;
    xlabel('$a$','fontsize',FS,'interpreter','latex')
    ylabel('$\log(|EE|)$','fontsize',FS,'interpreter','latex')
    title('Euler errors')
    plot(a_intp, log(abs(EE(zhigh,:))), 'LineWidth',LW)
    legend('$z_{low}$','$z_{high}$', 'fontsize',FS,'interpreter','latex'...
        ,'Location','best')
    yline(0,'LineStyle',':', 'LineWidth',LW,'HandleVisibility','off')
    ax=gca;
    ax.FontSize=FS;
    saveas(gcf,'Eulererrors','epsc')
end