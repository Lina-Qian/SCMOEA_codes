clear
clc
fclose all;
num = 6;

al{1} = 'C:\Users\Desktop\DDEPSO\';
al{2} = 'C:\Users\Desktop\QPSO\';
al{3} = 'C:\Users\Desktop\MOPSO\';
al{4} = 'C:\Users\Desktop\NSGAII\';
al{5} = 'C:\UsersDesktop\HGA\';
al{6} = 'C:\Users\Desktop\RLMOMAD\';

for idx = 1:100 
   
    if isfolder([al{1}, 'Mk', num2str(idx), '/'])
        ds = sprintf('Mk%d/', idx);
    elseif isfolder([al{1}, 'MK', num2str(idx), '/'])
        ds = sprintf('MK%d/', idx);
    else
        continue; 
    end
    
    name = 'res';
    txt='.txt';
    x0 = cell(1,num);
    y0 = cell(1,num);
    X0 = [];
    Y0 = [];
    x = {};
    y = {};
    k=1;
    for file = 1:10  
        for n = 1:num
            RealPath{n,k}=[al{n},ds,name,num2str(file),txt];
            [a,b]=DataRead(RealPath{n,k});
            x0{n} = [x0{n},a];
            y0{n} = [y0{n},b];
            X0 = [X0 a];
            Y0 = [Y0 b];
        end
        k = k+1;
    end
    all0 = [X0;Y0]';
    minXf1 = min(X0);
    maxXf1 = max(X0);
    minXf2 = min(Y0);
    maxXf2 = max(Y0);
    for i = 1:length(x0)
        x0{i} = (x0{i} - minXf1)/(maxXf1-minXf1);
        y0{i} = (y0{i} - minXf2)/(maxXf2-minXf2);
    end
    all0(:,1) = (all0(:,1) - minXf1)/(maxXf1-minXf1);
    all0(:,2) = (all0(:,2) - minXf2)/(maxXf2-minXf2);
    PF_all=pareto1(all0);
    PF_true = all0(PF_all,:);
    for n = 1:num
        tmp = [x0{n};y0{n}]';
        PF{n}=pareto1(tmp);
        x{n}(:,1) = tmp(PF{n},1);
        y{n}(:,2) = tmp(PF{n},2);
        PF_approx{n} = [x{n}(:,1),y{n}(:,2)];
        igd_value = calculate_igd(PF_true, PF_approx{n});
        igd(idx,n) = igd_value;
    end
    respath='IGD.txt';
    fout=fopen(respath,'w');
    fprintf(fout,'%4f %4f %4f %4f %4f %4f\r\n',igd');
    fclose(fout);
end


function igd = calculate_igd(PF_true, PF_approx)
igd = 0;
for i = 1:size(PF_true, 1)
    distances = sqrt(sum((PF_true(i,:) - PF_approx).^2, 2));
    [~, minIndex] = min(distances);
    minDistance = distances(minIndex);
    igd = igd + minDistance;
end
igd = igd / size(PF_true, 1);
end

function [x,y]=DataRead(RealPath)
fin=fopen(RealPath,'r'); 
sizeA = [2 Inf];
A=fscanf(fin,'%f %f',sizeA);
x = A(1,:);
y = A(2,:);
fclose(fin);
end

function PF=pareto1(obj)
PF=[];
M=2;
[obj_size,~]=size(obj);
pn=zeros(1,obj_size);
S=0;
for i=1:obj_size
    for j=1:obj_size
        dom_less=0;
        dom_equal=0;
        dom_more=0;
        if (obj(i,1)>obj(j,1))
            dom_more = dom_more + 1;
        elseif (obj(i,1)==obj(j,1))
            dom_equal = dom_equal + 1;
        else
            dom_less = dom_less + 1;
        end
        if (obj(i,2)>obj(j,2))
            dom_more = dom_more + 1;
        elseif (obj(i,2)==obj(j,2))
            dom_equal = dom_equal + 1;
        else
            dom_less = dom_less + 1;
        end
        if dom_less == 0 && dom_equal ~= M 
            pn(i) = pn(i)+ 1;
        end
    end
    if pn(i)== 0 
        PF=[PF i];
        S=S+1;
    end
end
end