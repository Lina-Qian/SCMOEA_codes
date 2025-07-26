function [pop,flag]= LS1(pop,i)
fit1=zeros(1,2);
[schedule,pop(i).OS_chrom,fit1(1,1),fit1(1,2)]=Decoding(pop(i).OS_chrom,pop(i).MV_chrom);
cpath = findCriticalPath(schedule);
r1=ceil(rand*size(cpath,1));r2=ceil(rand*size(cpath,1));
while r1==r2
    r2=ceil(rand*size(cpath,1));
end
xx=pop(i).OS_chrom;
x=pop(i).OS_chrom;
tmp = find(x==cpath(r1,1));
index1 = tmp(cpath(r1,2));
tmp = find(x==cpath(r2,1));
index2 = tmp(cpath(r2,2));
tmpj = x(index1);
x(index1) = x(index2);
x(index2) = tmpj;
newx = x;
fit2=zeros(1,2);
[~,newx,fit2(1,1),fit2(1,2)]=Decoding(newx,pop(i).MV_chrom);
result = NDS(fit2, fit1);
if result==1
    pop(i).OS_chrom=newx;
    pop(i).Fitness=fit2(1,:);
    flag=1;
elseif result==2||result==0
    pop(i).OS_chrom=xx;
    pop(i).Fitness=fit1(1,:);
    flag=0;
end
end
