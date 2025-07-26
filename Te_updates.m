function [pop,flag]=Te_updates(weight,new_p,new_m,zzold,zz,pop,k)
global minreference_point maxreference_point;
flag=0;
part1=abs(zz-minreference_point);
part2=abs(maxreference_point-minreference_point);
part=part1./part2;
newobj=max(weight.*part);
part1=abs(zzold-minreference_point);
part=part1./part2;
oldobj=max(weight.*part); 
if newobj<oldobj
    pop(k).OS_chrom=new_p;
    pop(k).MV_chrom=new_m;
    pop(k).Fitness=zz(1,:);
    flag=1;
end
end