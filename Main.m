%SCMOEA
clc;
clear all;
dbstop if error;
global N H SH NM TM PS M time;
global message_matrix;
global setuptime PA PB Pidle Pstandby Eshift Eturnonfirst Estandby2process 
global minreference_point maxreference_point;
T=[1,2,3];
R=20;
[~,T_size]=size(T);
numst=T_size;
PS=100;
Maxiter=200;
o=2;
pMutation = 0.8;
LP=40;

empty_individual.OS_chrom = [];
empty_individual.MV_chrom = [];
empty_individual.Fitness = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.CrowdingDistance = [];
empty_individual.Category = [];
pop = repmat(empty_individual, PS, 1);


path = 'C:\Users\Desktop\main\Mk\';
data = cell(1,10);
for i = 1:100
    newData = upper(['Mk', num2str(i)]);
    data{i} = newData; 
end
txt_tail='.txt';

%%main
for file=1:100
    numbers = [];
    for i = 1:length(data)
        match = regexp(data{i}, '\d+', 'match');
        if ~isempty(match)
            numbers(i) = str2double(match{1});
        end
    end
    MK_name = [path,'MK',num2str(numbers(file)),'.mat'];
    E_name = [path,'Energy',num2str(numbers(file)),'.mat'];
    respath='result\';
    tmp6='\';
    respath=[respath,data{file},tmp6];
    clear tmp6
    a=['mkdir' respath];
    system(a);
    fprintf('%s %s\r\n','caculating',data{file});
    L_array=zeros(1,10);
    totalPF=[];
    for round=1:10
        load(MK_name);
        load(E_name);
        [weights,neighbour] = init_weight(PS+1,o,R);
        [os_chrom,mv_chrom] = new_initial();
        fitness=zeros(PS,2);
        for i=1:PS
            [~,os_chrom(i,:),fitness(i,1),fitness(i,2)]=Decoding(os_chrom(i,:),mv_chrom(i,:));
        end
        [minreference_point,~]=min(fitness);
        [maxreference_point,~]=max(fitness);
        category= classifyIndividuals(fitness);
        for i=1:PS
            pop(i).OS_chrom= os_chrom(i,:);
            pop(i).MV_chrom= mv_chrom(i,:);
            pop(i).Fitness=fitness(i,:);
            pop(i).Category=category(i);
        end     
        ns1 =[];
        nf1 = [];
        pfit1 = ones(1, numst);
        normfit1=zeros(1, numst);
        ns2 =[];
        nf2 = [];
        pfit2 = ones(1, numst);
        normfit2=zeros(1, numst);
        ns3 =[];
        nf3 = [];
        pfit3 = ones(1, numst);
        normfit3=zeros(1, numst);
        PF_elite=[];
        for iter=1:Maxiter
            fprintf('%s %s %d %s %d\r\n',data{file},'round',round,'iter',iter);

            indices_1 = find([arrayfun(@(x) x.Category, pop)] == 1);
            if indices_1~=0
                [pop,pfit1,normfit1,ns1,nf1]=one_operation(pfit1,normfit1,ns1,nf1,pop,weights,indices_1,LP,neighbour,R);
            end

            indices_2 = find([arrayfun(@(x) x.Category, pop)] == 2);
            if indices_2~=0
                [pop,pfit2,normfit2,ns2,nf2]=two_operation(pfit2,normfit2,ns2,nf2,pop,indices_2,weights,LP,neighbour,R);
            end

            indices_3 = find([arrayfun(@(x) x.Category, pop)] == 3);
            if indices_3~=0
                [pop,pfit3,normfit3,ns3,nf3]=three_operation(pfit3,normfit3,ns3,nf3,pop,weights,indices_3,LP,neighbour,R);
            end

            for i=1:PS
                [pop,flag]= LS1(pop,i);
                if flag~=1
                    [pop,flag]= LS2(pop,i);
                    if flag~=1
                        [pop,flag]= LS3(pop,i);
                        if flag~=1
                            [pop,flag]= LS4(pop,i);
                        end
                    end
                end
            end

            newfitness=zeros(PS,2);
            for i=1:PS
                newfitness(i,:)=pop(i).Fitness;
            end

            category=classifyIndividuals(newfitness);
            for i=1:PS
                pop(i).Category=category(i);
            end
         
            [pop, F] = NonDominatedSorting(pop);
            F1 = pop(F{1});
            PF_elite=[PF_elite;F1];        
        end

        unique_pop = PF_elite(1);
        for i = 2:length(PF_elite)
            is_duplicate = false;
            for j = 1:length(unique_pop)
                if isequal(PF_elite(i).Fitness, unique_pop(j).Fitness)
                    is_duplicate = true;
                    break;
                end
            end
            if ~is_duplicate
                unique_pop(end+1) =PF_elite(i); 
            end
        end
        PF_elite = unique_pop; 

        [PF_elite, F] = NonDominatedSorting(PF_elite);
        F1 = PF_elite(F{1});
        L=numel(F1);
        for i=1:L
            obj(i,:)=F1(i).Fitness;
        end
        obj=unique(obj,'rows');  
        [L,~]=size(obj);
        L_array(round)=L;
        for cc=1:L
            totalPF=[totalPF;obj(cc,:)];
        end
    end
    current_index=1;
    for round=1:10
        obj=[];
        endindex=current_index+L_array(round)-1;
        for i=current_index:endindex
            obj=[obj;totalPF(i,:)];
        end
        current_index=current_index+L_array(round);
        tmp5=obj';
        tmp1='res';
        tmp2=num2str(round);
        tmp3='.txt';
        resPATH=[respath tmp1 tmp2 tmp3];
        fout=fopen(resPATH,'w');
        fprintf(fout,'%5.2f %6.3f\r\n',tmp5);
        fclose(fout);
    end
    fprintf('%s %s\r\n','Finish ',data{file})
end