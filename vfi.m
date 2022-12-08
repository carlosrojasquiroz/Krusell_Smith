function [V,g_a,p_a,g_c]=vfi(p,m)
V=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
g_a=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
g_c=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
distV=1;
iter=0;
r=zeros(p.nZZ,p.nKK);
w=zeros(p.nZZ,p.nKK);
L=m.z_grid*m.mu;

for d_3=1:p.nZZ
    for d_4=1:p.nKK
        r(d_3,d_4)=m.Z_grid(d_3)*p.alpha*(L/m.K_grid(d_4))^(1-p.alpha)-p.delta;
        w(d_3,d_4)=m.Z_grid(d_3)*(1-p.alpha)*(m.Z_grid(d_3)*p.alpha/(r(d_3,d_4)+p.delta))^(p.alpha/(1-p.alpha));
    end
end

while distV>p.tol && iter < p.maxiter
    [Vnew,p_a]=vfi_step(p,m,V);
    distV=max(abs(V(:)-Vnew(:)));
    V=Vnew;
    iter=iter+1;
    if p.disp1==1   
        disp(['VFI iteration number ',num2str(iter)])
        disp(['error = ',num2str(distV)])
    end
end

for d_1=1:p.nzz
    for d_2=1:p.naa
        for d_3=1:p.nZZ
            for d_4=1:p.nKK
                g_a(d_1,d_2,d_3,d_4)=m.a_grid(p_a(d_1,d_2,d_3,d_4));
                g_c(d_1,d_2,d_3,d_4)=w(d_3,d_4)*m.z_grid(d_1)+(1+r(d_3,d_4))*m.a_grid(d_2)-g_a(d_1,d_2,d_3,d_4);
            end
        end
    end
end