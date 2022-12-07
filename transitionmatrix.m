function TM=transitionmatrix(p,m,g_a,intp)
%---------------------------------------------------------------------------------------------------------------------------
% This function computes the transition matrix for a grid specification
% with "intp" grid points
%---------------------------------------------------------------------------------------------------------------------------
a=gridspecA(p,intp);    
TM=zeros(p.nzz,intp,p.nzz,intp);
for ind_i = 1:intp     
    for ind_j = 1:p.nzz
        a_oo=interp1(m.a_grid,g_a(ind_j,:),a(ind_i),'linear');
        a_low=min(intp-1,find(a<=a_oo,1,'last'));
        a_low=max(a_low,1);
        a_sup=a_low+1;
        w_s=(a_oo-a(a_low))/(a(a_sup)-a(a_low));
        w_l=1-w_s;
        for ind_j1=1:p.nzz
            TM(ind_j,ind_i,ind_j1,a_low)=m.Pi(ind_j,ind_j1)*w_l;
            TM(ind_j,ind_i,ind_j1,a_sup)=m.Pi(ind_j,ind_j1)*w_s;
        end
    end
end
TM=reshape(TM,intp*p.nzz,intp*p.nzz);