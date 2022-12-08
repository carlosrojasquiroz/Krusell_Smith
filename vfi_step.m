function [Vnew,p_a]=vfi_step(p,m,V)
Vnew=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
p_a=zeros(p.nzz,p.naa,p.nZZ,p.nKK);
for d_1=1:p.nzz
    for d_3=1:p.nZZ
        for d_4=1:p.nKK
            EV=m.P(d_1*d_3+(d_1-d_3)*(d_1>d_3),:)*reshape(V(:,:,:,d_4),p.nzz*p.nZZ,p.naa);
            for d_2=1:p.naa
                c=p.w*m.z_grid(d_1)+(1+p.r)*m.a_grid(d_2)-m.a_grid;
                Vaux=utility(c,p.sigma)+p.beta*EV;
                [Vnew(d_1,d_2,d_3,d_4),p_a(d_1,d_2,d_3,d_4)]=max(Vaux(:));
            end
        end
    end
end
