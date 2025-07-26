function [pchrom,mchrom]=mutation(pchrom,mchrom)
global SH N H NM ;
        p1=ceil(rand*SH);   
        p2=ceil(rand*SH);   
        while (p1==p2)||(pchrom(p1)==pchrom(p2))  
           p2=ceil(rand*SH);
        end
        t=pchrom(p1);             
        pchrom(p1)=pchrom(p2);           
        pchrom(p2)=t;            
         s1=pchrom;
         s2=zeros(1,SH);
         p=zeros(1,N);
        for i=1:SH
             p(s1(i))=p(s1(i))+1;
             s2(i)=p(s1(i));
        end     
        s3=mchrom;
        p1=ceil(rand*SH);      
        p2=ceil(rand*SH);      
        while(p1==p2)
           p2=ceil(rand*SH);
        end        
        n=NM{s1(p1),s2(p1)};
        m=ceil(rand*n*3);
     
        if n>1
        while(s3(p1)==m)
               m=ceil(rand*n*3);                
        end
        end
        mchrom(1,sum(H(1,1:s1(p1)-1))+s2(p1))=m;        
        n=NM{s1(p2),s2(p2)};
        m=ceil(rand*n*3);      
        if n>1
        while(s3(p1)==m)
               m=ceil(rand*n*3);              
        end
        end
        mchrom(1,sum(H(1,1:s1(p2)-1))+s2(p2))=m;
end