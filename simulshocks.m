function [Z,agshock]=simulshocks(p,m)
agshock=zeros(1,p.simulT);
agshock(1)=1; 

for t=2:p.simulT
   raux=rand; 
   if raux<=m.PZ(1,agshock(t-1)) 
      agshock(t)=1; 
   else
      agshock(t)=2;
   end
end

Z=m.Z_grid(agshock);