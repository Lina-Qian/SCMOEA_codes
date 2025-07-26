function [p_parent_1,m_parent_1,p_parent_2,m_parent_2]=crossover(pchrom1,mchrom1,pchrom2,mchrom2)
    global N SH;
    p_parent_1=pchrom1;
    m_parent_1=mchrom1;
    p_parent_2=pchrom2;
    m_parent_2=mchrom2;
    J1=[];
    c1_p=zeros(1,SH);
    c2_p=zeros(1,SH);
    while size(J1,1)==0 && size(J1,2)==0
        J1=find(round(rand(1,N))==1);
    end
    for j=1:SH
        if ismember(p_parent_1(j),J1) 
            c1_p(j)=p_parent_1(j);
        end    
        if ~ismember(p_parent_2(j),J1) 
            c2_p(j)=p_parent_2(j);
        end
    end    
    index_1_1=find(c1_p==0);
    index_1_2=find(c2_p~=0);
    index_2_1=find(c2_p==0);
    index_2_2=find(c1_p~=0);
    for j=1:size(index_1_1,2)
        c1_p(index_1_1(j))=p_parent_2(index_1_2(j));
    end
    for j=1:size(index_2_1,2)
        c2_p(index_2_1(j))=p_parent_1(index_2_2(j));
    end  
    p_parent_1=c1_p;
    p_parent_2=c2_p; 
    s=round(rand(1,SH))==1;
    for i=1:SH
      if(s(i)==1)
         t=m_parent_1(i);
         m_parent_1(i)=m_parent_2(i);
         m_parent_2(i)=t;
       end
    end    
end
