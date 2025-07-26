function [pop,flag]= LS4(pop,i)
global message_matrix H;
fit1=zeros(1,2);
os=pop(i).OS_chrom;
mvv=pop(i).MV_chrom;
mv=pop(i).MV_chrom;
[schedule,os,fit1(1,1),fit1(1,2)]=Decoding(os,mv);
schedule=schedule(schedule(:,2)~=0,:);
num=size(schedule,1);
set=[];
for j=1:num
    start_hour = mod(schedule(j,5),24); 
    end_hour = mod(schedule(j,6),24); 
    mid_hour=schedule(j,6)-schedule(j,5);
    if mid_hour>19
        set=[set;schedule(j,1),schedule(j,2),schedule(j,3)];
    elseif 10>start_hour&&15<=end_hour
        set=[set;schedule(j,1),schedule(j,2),schedule(j,3)];
    elseif isElectricityPriceHigh(start_hour)||isElectricityPriceHigh(end_hour)
        set=[set;schedule(j,1),schedule(j,2),schedule(j,3)];
    end
end
set=unique(set,'rows');
selectM=set(:, end);
selectM=unique(selectM,'rows');
num_selectM=size(selectM);
for kk=1:num_selectM
    bb=set(set(:,3)==selectM(kk),:);
    sizebb=size(bb,1);
    select_operation=ceil(rand*sizebb);
    message=message_matrix{bb(select_operation,1),bb(select_operation,2)};
    logicalIndex = message(3, :) == 1;
    midmessage = message(:, logicalIndex);
    selectmv=midmessage(1,:);
    randomIndex = randi(length(selectmv));
    randomValue = selectmv(randomIndex);
    mv(1,sum(H(1,1:bb(select_operation,1)-1))+bb(select_operation,2))=randomValue;
end
newmv=mv;
fit2=zeros(1,2);
[~,~,fit2(1,1),fit2(1,2)]=Decoding(pop(i).OS_chrom,newmv);
result = NDS(fit2, fit1);
if result==1
    pop(i).MV_chrom=newmv;
    pop(i).Fitness=fit2(1,:);
    flag=1;
elseif result==2||result==0
    pop(i).MV_chrom=mvv;
    pop(i).Fitness=fit1(1,:);
    flag=0;
end
end
