function [V,g_a,g_c]=egm(p,m)
%---------------------------------------------------------------------------------------------------------------------------
% This function computes the policy functions for assets and consumption
% using the endogenous grid method.
%---------------------------------------------------------------------------------------------------------------------------
g_a=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
g_c=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
r=zeros(p.nZZ,p.nKK);
w=zeros(p.nZZ,p.nKK);
c_aux=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
V=zeros(p.nzz,p.naa,p.nZZ,p.nKK);

err=1;
iter=0;
L=m.z_grid*m.mu';

for d_3=1:p.nZZ
    for d_4=1:p.nKK
        r(d_3,d_4)=m.Z_grid(d_3)*p.alpha*(L/m.K_grid(d_4))^(1-p.alpha)-p.delta;
        w(d_3,d_4)=m.Z_grid(d_3)*(1-p.alpha)*(m.Z_grid(d_3)*p.alpha/(r(d_3,d_4)+p.delta))^(p.alpha/(1-p.alpha));
        c_aux(:,:,d_3,d_4)=w(d_3,d_4)*m.z_grid'+r(d_3,d_4)*m.a_grid;
        V(:,:,d_3,d_4)=utility(c_aux(:,:,d_3,d_4),p.sigma)./(1-p.beta);
    end
end

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
    V_aux=V;
    for d_1=1:p.nzz
        for d_3=1:p.nZZ
            for d_4=1:p.nKK
                c_oo=reshape(c_aux(:,:,:,d_4),p.nzz*p.nZZ,p.naa);
                EUc= (1+r(d_3,d_4))*m.P(d_1*d_3+(d_1-d_3)*(d_1>d_3),:)*c_oo.^(-p.sigma);
                Ucp= (p.beta*EUc).^(-1/p.sigma);
                a=(Ucp + m.a_grid - m.z_grid(d_1)*w(d_3,d_4))/(1+r(d_3,d_4));
                g_a(d_1,:,d_3,d_4)=interp1(a,m.a_grid,m.a_grid,'linear','extrap');
                for d_2=1:p.naa
                    if g_a(d_1,d_2,d_3,d_4)<m.a_grid(1)
                        g_a(d_1,d_2,d_3,d_4)=m.a_grid(1);
                    elseif g_a(d_1,d_2,d_3,d_4)>m.a_grid(end)
                        g_a(d_1,d_2,d_3,d_4)=m.a_grid(end);
                    end
                end
                g_c(d_1,:,d_3,d_4)=m.z_grid(d_1)*w(d_3,d_4) + (1+r(d_3,d_4))*m.a_grid - g_a(d_1,:,d_3,d_4);
                V_oo=utility(g_c(d_1,:,d_3,d_4),p.sigma)+p.beta*m.P(d_1*d_3+(d_1-d_3)*(d_1>d_3),:)*reshape(V_aux(:,:,:,d_4),p.nzz*p.nZZ,p.naa);
                V(d_1,:,d_3,d_4)=interp1(a,V_oo,m.a_grid,'linear','extrap');                
             end
        end
    end
    if p.algo==0
        err=max(max(max(max(abs(V-V_aux)))));
    elseif p.algo==1
        err=max(max(max(max(abs(c_aux-g_c)))));
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