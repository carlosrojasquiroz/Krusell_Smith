function [V,g_a,g_c]=egm(p,m)
%---------------------------------------------------------------------------------------------------------------------------
% This function computes the policy functions for assets and consumption
% using the endogenous grid method.
%---------------------------------------------------------------------------------------------------------------------------
g_a=zeros(p.nzz,p.naa);
g_c=zeros(p.nzz,p.naa);
err=1;
iter=0;
c_aux=p.w*m.z_grid'+p.r*m.a_grid;
V=utility(c_aux,p.sigma)./(1-p.beta);

if p.disp1==1   
        disp('I start the EGM algorithm')
    if p.algo==0
        disp('iterating over the value function...')
    elseif p.algo==1
        disp('iterating over the policy function...')
    end
end

tic;
while err>p.tol && iter<p.maxiter
    iter=iter+1;
    Vaux=V;
    EUc= m.Pi*c_aux.^(-p.sigma);
    Ucp= ((1+p.r)*p.beta*EUc).^(-1/p.sigma);
    a=(Ucp + m.a_grid - m.z_grid'*p.w)/(1+p.r);
    for ind_j=1:p.nzz
        g_a(ind_j,:)=interp1(a(ind_j,:),m.a_grid,m.a_grid,'linear','extrap');
    end
    g_a(g_a<m.a_grid(1)) = m.a_grid(1);
    g_a(g_a>m.a_grid(end)) = m.a_grid(end);
    g_c=m.z_grid'*p.w + (1+p.r)*m.a_grid - g_a;
    Vo=utility(g_c,p.sigma)+p.beta*m.Pi*Vaux;
    for ind_j=1:p.nzz
        V(ind_j,:)=interp1(a(ind_j,:),Vo(ind_j,:),m.a_grid,'linear','extrap');
    end
    if p.algo==0
        err=max(max(abs(V-Vaux)));
    elseif p.algo==1
        err=max(max(abs(c_aux-g_c)));
    else
        disp('You should choose between options 0 or 1 only')
        disp('p.algo=0, iterate over the value function')
        disp('p.algo=1, iterate over the consumption policy function')
        break
    end
    c_aux=g_c;
    if p.disp1==1   
        disp(['EGM iteration number ',num2str(iter)])
        disp(['error = ',num2str(err)])
    end
end
Endtime=toc;

if p.disp1==1   
    disp(['Time to execute the EGM loop ',num2str(Endtime)])
end

if p.fig==1
    LW=1.5;
    FS=16;
    zlow=1;
    zhigh=2;
%---------------------------------------------------------------------------------------------------------------------------
% Assets policy function
%---------------------------------------------------------------------------------------------------------------------------
    figure()
    plot(m.a_grid,g_a(zlow,:),'LineWidth',LW)
    hold on;
    grid on;
    xlabel('$a$','fontsize',FS,'interpreter','latex')
    ylabel('$a^{\prime}=g_a(z,a)$','fontsize',FS,'interpreter','latex')
    title('Assets policy function')
    plot(m.a_grid,g_a(zhigh,:), 'LineWidth',LW)
    plot(m.a_grid,m.a_grid, 'LineWidth',LW,'LineStyle','--')
    legend('$z_{low}$','$z_{high}$','$45^o \ line$','fontsize',FS,'interpreter','latex'...
        ,'Location','best')
    yline(0,'LineStyle',':', 'LineWidth',LW,'HandleVisibility','off')
    ax=gca;
    ax.FontSize=FS;
    saveas(gcf,'Assets','epsc')
%---------------------------------------------------------------------------------------------------------------------------
% Consumption policy function
%---------------------------------------------------------------------------------------------------------------------------
    figure()
    plot(m.a_grid,g_c(zlow,:),'LineWidth',LW)
    hold on;
    grid on;
    xlabel('$a$','fontsize',FS,'interpreter','latex')
    ylabel('$c=g_c(z,a)$','fontsize',FS,'interpreter','latex')
    title('Consumption policy function')
    plot(m.a_grid,g_c(zhigh,:), 'LineWidth',LW)
    legend('$z_{low}$','$z_{high}$', 'fontsize',FS,'interpreter','latex'...
        ,'Location','best')
    yline(0,'LineStyle',':', 'LineWidth',LW,'HandleVisibility','off')
    ax=gca;
    ax.FontSize=FS;
    saveas(gcf,'Consumption','epsc')
%---------------------------------------------------------------------------------------------------------------------------
% Value function
%---------------------------------------------------------------------------------------------------------------------------
    figure()
    plot(m.a_grid,V(zlow,:),'LineWidth',LW)
    hold on;
    grid on;
    xlabel('$a$','fontsize',FS,'interpreter','latex')
    ylabel('$V(z,a)$','fontsize',FS,'interpreter','latex')
    title('Value function')
    plot(m.a_grid,V(zhigh,:), 'LineWidth',LW)
    legend('$z_{low}$','$z_{high}$', 'fontsize',FS,'interpreter','latex'...
        ,'Location','best')
    yline(0,'LineStyle',':', 'LineWidth',LW,'HandleVisibility','off')
    ax=gca;
    ax.FontSize=FS;
    saveas(gcf,'Valuefunction','epsc')    
end