function [pop,flag]= LS2(pop,i)
global H
fit1=zeros(1,2);
[schedule,pop(i).OS_chrom,fit1(1,1),fit1(1,2)]=Decoding(pop(i).OS_chrom,pop(i).MV_chrom);
cpath = findCriticalPath(schedule);
r1=ceil(rand*size(cpath,1));
xx=pop(i).MV_chrom;
x=pop(i).MV_chrom;
a = cpath(r1,1);
j = cpath(r1,2);
idx = sum(H(1,1:a-1))+j;
if x(idx)~=1
    x(idx)=1;
end
newx = x;
fit2=zeros(1,2);
[~,~,fit2(1,1),fit2(1,2)]=Decoding(pop(i).OS_chrom,newx);
result = NDS(fit2, fit1);
if result==1
    pop(i).MV_chrom=newx;
    pop(i).Fitness=fit2(1,:);
    flag=1;
elseif result==2||result==0
    pop(i).MV_chrom=xx;
    pop(i).Fitness=fit1(1,:);
    flag=0;
end
end