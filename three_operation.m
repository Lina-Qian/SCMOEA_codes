function [pop,pfit3,normfit3,ns3,nf3]=three_operation(pfit3,normfit3,ns3,nf3,pop,weights,indices_3,LP,neighbour,R)
global PS minreference_point maxreference_point;
k=size(indices_3,1);
numst=3;
learngen = LP;
rr = rand;
spacing = 1/k;
randnums = sort(mod(rr : spacing : 1 + rr - 0.5 * spacing, 1));
sumfit=0;
for j=1:numst
    sumfit=pfit3(j)+sumfit;
end
for j=1:numst
    normfit3(j)=pfit3(j)/sumfit;
end
partsum = 0;
count(1) = 0;
stpool = [];
for j = 1 : length(pfit3)
    partsum = partsum + normfit3(j);
    count(j + 1) = length(find(randnums < partsum));
    select(j, 1) = count(j + 1) - count(j);
    stpool = [stpool; ones(select(j, 1), 1) * j];
end
stpool = stpool(randperm(k));
scount=zeros(1,numst);
lcount=zeros(1,numst);
for i=1:k
    chooseA=stpool(i);
    indices = find([pop.Category] == chooseA);
    num=size(indices,2);
    if num==1
        aa=indices;
    elseif num>1
        aa=randsample(indices, 1);
        while aa==indices_3(i)
            aa=randsample(indices, 1);
        end
    else
        aa=ceil(rand*PS);
        while aa==indices_3(i)
            aa=ceil(rand*PS);
        end
    end
    [new_pchrom1,new_mchrom1,new_pchrom2,new_mchrom2] = crossover(pop(aa).OS_chrom,pop(aa).MV_chrom,pop(indices_3(i)).OS_chrom,pop(indices_3(i)).MV_chrom);
    if rand<0.8
        [new_pchrom1,new_mchrom1] = mutation(new_pchrom1,new_mchrom1);
        [new_pchrom2,new_mchrom2] = mutation(new_pchrom2,new_mchrom2);
    end
    f1=zeros(1,2);
    f2=zeros(1,2);
    [~,new_pchrom1,f1(1,1),f1(1,2)]=Decoding(new_pchrom1,new_mchrom1);
    [~,new_pchrom2,f2(1,1),f2(1,2)]=Decoding(new_pchrom2,new_mchrom2);
    minreference_point=min(minreference_point,f1);
    minreference_point=min(minreference_point,f2);
    maxreference_point=max(maxreference_point,f1);
    maxreference_point=max(maxreference_point,f2);
    nei=neighbour(indices_3(i),:);
    for a=1:R
        kk=nei(a);
        tmp(1,:)=pop(kk).Fitness;
        [pop,flag]=Te_updates(weights(kk,:),new_pchrom1,new_mchrom1,tmp,f1,pop,kk);
        if flag==1
            scount(stpool(i))=scount(stpool(i))+1;
        elseif(flag==0)
            lcount(stpool(i))=lcount(stpool(i))+1;
        end
    end
    for a=1:R
        kk=nei(a);
        tmp(1,:)=pop(kk).Fitness;
        [pop,flag]=Te_updates(weights(kk,:),new_pchrom2,new_mchrom2,tmp,f2,pop,kk);
        if flag==1
            scount(stpool(i))=scount(stpool(i))+1;
        elseif(flag==0)
            lcount(stpool(i))=lcount(stpool(i))+1;
        end
    end
end
ns3=[ns3;scount];
nf3=[nf3;lcount];
[i,~]=size(ns3);

if i >= learngen
    [ns_row,~]=size(ns3);
    [nf_row,~]=size(ns3);
    for j = 1 : numst
        sum_ns=0;sum_nf=0;
        for k=1:ns_row
            sum_ns=sum_ns+ns3(k,j);
        end
        for k=1:nf_row
            sum_nf=sum_nf+nf3(k,j);
        end
        if (sum_ns + sum_nf) == 0
            pfit3(j) = 0.01;
        else
            pfit3(j) = sum_ns / (sum_ns + sum_nf) + 0.01;
        end
    end
    if ~isempty(ns3), ns3(1, :) = [];  end
    if ~isempty(nf3), nf3(1, :) = [];  end
end
end