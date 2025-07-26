function [pop,flag] = LS3(pop,i)
global   message_matrix H
OS_chrom=pop(i).OS_chrom;
MV_chrom=pop(i).MV_chrom;
[schedule,New_OS_chrom,Ee] = Decode_Tec(OS_chrom,MV_chrom);
[~,machine_index] = max(Ee);
schedule = schedule(schedule(:,2)~=0,:);
set= [schedule(schedule(:,3)==machine_index,1),schedule(schedule(:,3)==machine_index,2)];
for idx = 1:size(set,1)
    a = set(idx,1);
    b = set(idx,2);
    set(idx,3) = size(message_matrix{a,b},2);
end
idx2 = find(set(:,3)>=6);
if isempty(idx2)
    new_OS_chrom=OS_chrom;
    new_MV_chrom=MV_chrom;
else
    idx3 = randperm(length(idx2),1);
    c = set(idx2(idx3),1);
    j = set(idx2(idx3),2);
    ms_idx = sum(H(1,1:c-1))+j;
    my_machine_num = set(idx2(idx3),3);
    idx4 = randperm(my_machine_num,1);
    while ceil(idx4/3)==ceil(MV_chrom(ms_idx)/3)
        idx4 = randperm(my_machine_num,1);
    end
    newMV_chrom=MV_chrom;
    newMV_chrom(ms_idx) = idx4;
    new_MV_chrom = newMV_chrom;
    new_OS_chrom=New_OS_chrom;
end
[~,OS_chrom,fit1(1,1),fit1(1,2)]=Decoding(OS_chrom,MV_chrom);
[~,new_OS_chrom,fit2(1,1),fit2(1,2)]=Decoding(new_OS_chrom,new_MV_chrom);
result = NDS(fit2, fit1);
if result==1
    pop(i).OS_chrom=new_OS_chrom;
    pop(i).MV_chrom=new_MV_chrom;
    pop(i).Fitness=fit2(1,:);
    flag=1;
elseif result==2||result==0
    pop(i).OS_chrom=OS_chrom;
    pop(i).MV_chrom=MV_chrom;
    pop(i).Fitness=fit1(1,:);
    flag=0;
end
end

