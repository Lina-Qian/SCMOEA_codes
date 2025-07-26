function [weights,neighbour] = init_weight(ps,M,R)   
Popsize=ps-1;                                         
weights=zeros(Popsize,M);                           
count=1;
for i=1:Popsize                                    
    weights(count,1)=i/ps;
    weights(count,2)=1-i/ps;
    count=count+1;
end
distance = zeros(Popsize,Popsize);               
neighbour=zeros(Popsize,R);                        
for i=1:Popsize
    for j=i+1:Popsize
        A=weights(i,:);B=weights(j,:);
        distance(i,j)=(A-B)*(A-B)';           
        distance(j,i)=distance(i,j);            
    end
    [~,sindex]=sort(distance(i,:));           
    neighbour(i,:)=sindex(1:R);                
end
end