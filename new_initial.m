function [Total_OS,Total_MV] = new_initial()
global  N H SH PS ;
X_OS=zeros(PS*0.9,SH);
Total_MV=[];
X0_OS=zeros(1,SH);
for i=1:N
    for j=1:H(i)
        k=sum(H(1,1:i))-H(i)+j;
        X0_OS(k)=i;
    end
end
tmp=X0_OS;
X_OS(1,:)=tmp(randperm(length(tmp)));
for i=2:PS*0.9
    tmp=X_OS(i-1,:);
    X_OS(i,:)=tmp(randperm(length(tmp)));
end
[Pchrom,Mchrom]=H3();
Total_OS=[X_OS;Pchrom];
Total_MV(1,:) = 1*ones(1,SH);
Total_MV(2,:) = LMS(X_OS(i,:));
X_MS1 = H1(); 
X_MS2 = H2();
Total_MV = [Total_MV;X_MS1;X_MS2];
for i = 0.7*PS+3:PS*0.9
    Total_MV(i,:) = RMS(); 
end
Total_MV = [Total_MV;Mchrom];
end

function [ms,pt,pe] = RMS()
global N SH H message_matrix PB
ms = zeros(1,SH);
pt=0;
pe=0;
for i = 1:N
    for j=1:H(i)
        myselections = size(message_matrix{i,j},2);
        t=ceil(rand*myselections);
        t1=sum(H(1,1:i-1))+j;
        ms(t1)=t;
        pt = pt+message_matrix{i,j}(4,t);
    end
end
pe = pe+pt*PB(message_matrix{i,j}(3,t));
end

function ms = LMS(os)
global SH H message_matrix
ms = zeros(1,SH);
for index = 1:length(os)
    i = os(index); 
    for j = 1:H(i) 
        myselections = message_matrix{i,j}(1,message_matrix{i,j}(3,:)==1);
        myselection = myselections(1);
        ms(sum(H(1,1:i-1))+j) = myselection;
    end
end
end

function MS = H1()
global SH N H message_matrix PS
MS = zeros(0.35*PS,SH);
for n = 1:2000
    ms(n,:) = randi([1, 2], SH, 1);
    pt(n) = 0;
    for i = 1:N
        for j=1:H(i)
            t = sum(H(1,1:i-1))+j;
            pt(n) = pt(n)+message_matrix{i,j}(4,ms(n,t));
        end
    end
end
[~,index] = sort(pt,2);
for n = 1:0.35*PS
    MS(n,:) = ms(index(n),:);
end
end

function MS = H2()
global SH N H message_matrix PS PB
MS = zeros(0.35*PS,SH);
for n = 1:2000
    pe(n) = 0;
    for i = 1:N
        for j=1:H(i)
            t = sum(H(1,1:i-1))+j;
            index = find(message_matrix{i,j}(3,:)==1);
            myselections = size(index,2);
            idx=ceil(rand*myselections);
            ms(n,t) = index(idx);
            pe(n) = pe(n)+message_matrix{i,j}(4,ms(n,t))*PB(message_matrix{i,j}(3,ms(n,t)));
        end
    end
end
[~,index2] = sort(pe,2);
for n = 1:0.35*PS
    MS(n,:) = ms(index2(n),:);
end
end

function [Pchrom,Mchrom]=H3()
global N H NM time message_matrix PS SH TM M
subpopsize=PS*0.1;
Pchrom=[];
Mchrom=[];
m_chrom=zeros(subpopsize,SH);
p_chrom=zeros(subpopsize,SH);
chrom=zeros(1,SH);
for i=1:N
    for j=1:H(i)
        k=sum(H(1,1:i))-H(i)+j;
        chrom(k)=i;
    end
end
s1=chrom;
s2=zeros(1,SH);
p=zeros(1,N);
for i=1:SH
    p(s1(i))=p(s1(i))+1;
    s2(i)=p(s1(i));
end
OPM=zeros(1,SH);
for i=1:SH
    OPM(i)=NM{s1(i),s2(i)};
end
[~,index]=sort(OPM,'ascend');
chrom=chrom(index);
for i=1:subpopsize
    p_chrom(i,:)=chrom;
end
e=[0];
for k=1:subpopsize
    mt=zeros(1,TM);
    for i=1:TM
        mt(i)=e;
    end
    mm=zeros(1,SH);
    s1=p_chrom(k,:);
    s2=zeros(1,SH);
    p=zeros(1,N);
    for i=1:SH
        p(s1(i))=p(s1(i))+1;
        s2(i)=p(s1(i));
    end
    for i=1:SH      
        n=NM{s1(i),s2(i)};
        if n==1
            mm(i)=M{s1(i),s2(i),1};
            mt(mm(i))=mt(mm(i))+time{s1(i),s2(i),mm(i)};
            continue;
        else
            avalible_m=zeros(1,n);
            avalible_m_load=zeros(1,n);
            for j=1:n
                avalible_m(j)=M{s1(i),s2(i),j};
                avalible_m_load(j)=mt(avalible_m(j));
            end
            [minvalue,~]=min(avalible_m_load);
            candidateM=find(avalible_m_load==minvalue);
            sizeM=size(candidateM,2);
            if(sizeM==1)
                mm(i)=avalible_m(candidateM);
                mt(mm(i))=mt(mm(i))+time{s1(i),s2(i),mm(i)};
            else
                t=time{s1(i),s2(i),avalible_m(candidateM(1))};
                mm(i)=avalible_m(candidateM(1));
                for kk=2:sizeM
                    tmp=time{s1(i),s2(i),avalible_m(candidateM(kk))};
                    if t>tmp
                        mm(i)=avalible_m(candidateM(kk));
                    end
                end
                mt(mm(i))=mt(mm(i))+time{s1(i),s2(i),mm(i)};
            end
        end
    end  
    for i=1:SH
        t1=s1(i);
        t2=s2(i);
        m_chrom(k,sum(H(1,1:t1-1))+t2)=mm(i);
    end
    for i=1:N
        for j=1:H(i)
            midmessage=message_matrix{i,j};
            indices = find( midmessage(2, :) == m_chrom(k,sum(H(1,1:i-1))+j) & midmessage(3, :) == 3);
            m_chrom(k,sum(H(1,1:i-1))+j)=indices;
        end
    end
end
Pchrom=[Pchrom;p_chrom];
Mchrom=[Mchrom;m_chrom];
end