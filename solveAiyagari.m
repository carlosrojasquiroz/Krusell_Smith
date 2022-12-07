function R=solveAiyagari(p,m)
%---------------------------------------------------------------------------------------------------------------------------
% This function solves the general equilibrium of the Aiyagari model 
%---------------------------------------------------------------------------------------------------------------------------
p.disp1=0;
p.fig=0;
r_min=-p.delta;
r_max=(1-p.beta)/p.beta;
err_R=1;
iter=0;

if p.disp2==1   
    disp('I start to solve the Aiyagari model...')
end

tic;
while abs(err_R)>p.tol && iter<p.maxiter
    iter = iter + 1;
    p.r=0.5*(r_max+r_min);
    [p.K0,p.L]=firms(p,m);
    p.w=(1-p.alpha)*p.Z*(p.K0/p.L)^p.alpha;
    if p.r<0
        p.Amin=0;
    else
        p.Amin=max(0,-p.w*m.z_grid(1)/p.r);
    end
    p.naa=p.intp;
    m.a_grid=gridspecA(p,p.naa);
    [V,g_a,g_c]=egm(p,m);
    tm=transitionmatrix(p,m,g_a,p.naa);
    dist=stationarydist(p,tm,p.naa);
    Ks=sum(dist*m.a_grid');
    [Kd,~]=firms(p,m);
    r1=p.alpha*p.Z*(p.L./max(Ks,0.0001)).^(1-p.alpha)-p.delta;
    w1=(1-p.alpha)*p.Z*(max(Ks,0.0001)./p.L).^p.alpha;
    err_R=r1-p.r;
    err_W=w1-p.w;
    if err_R < 0
        r_max = p.r;
    else
        r_min = p.r;
    end
    if p.disp2==1   
    disp(['Iteration number ',num2str(iter)])
    disp(['error(r) = ',num2str(abs(err_R)),', error(w) = ',num2str(abs(err_W))])
    end
end
Endtime=toc;

if p.disp2==1   
    disp(['Time to solve the model ',num2str(Endtime)])
end

R.r=p.r;
R.w=p.w;
R.Ks=Ks;
R.Kd=Kd;
R.L=p.L;
R.g_c=g_c;
R.g_a=g_a;
R.V=V;
R.dist=dist;
R.hist=sum(dist);
R.a_grid=m.a_grid;
R.Y=p.Z*R.Ks^p.alpha*R.L^(1-p.alpha);