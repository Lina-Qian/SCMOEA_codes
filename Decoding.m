function [schedule,os_chrom,Cmax,TEC]=Decoding(os_chrom,mv_chrom)
global N H message_matrix setuptime Eshift TM PA PB Pstandby Pidle Estandby2process Eturnonfirst
bigM=1000000;
EC=0;
ED=0;
schedule=zeros(0,6);
alls=[];
job_operation=zeros(1,N);
job_pre_endtime=zeros(1,N);
for idx=1:length(os_chrom)
    i=os_chrom(idx);
    job_operation(i)=job_operation(i)+1;
    mindex=sum(H(1,1:i-1))+job_operation(i);
    my_machineid=message_matrix{i,job_operation(i)}(2,mv_chrom(mindex)); 
    my_machinev=message_matrix{i,job_operation(i)}(3,mv_chrom(mindex)); 
    my_ptime=message_matrix{i,job_operation(i)}(4,mv_chrom(mindex)); 
    my_can_time=job_pre_endtime(i);
    midschedule=schedule(schedule(:,3)==my_machineid,:);
    midschedule=sortrows(midschedule,5);
    emptysolt=zeros(0,2);
    num_schedules=size(midschedule,1);
    if num_schedules==0
        midsolt=[0 bigM];
        emptysolt=[emptysolt;midsolt];
    else
        pre_endtime=0;
        for kk=1:num_schedules
            if (midschedule(kk,5)>pre_endtime)
                midsolt=[pre_endtime midschedule(kk,5)];
                emptysolt=[emptysolt;midsolt];
            end
            pre_endtime=midschedule(kk,6);
        end
        midsolt=[pre_endtime bigM];
        emptysolt=[emptysolt;midsolt];
    end
    emptysolt=sortrows(emptysolt,1);
    num_empty=size(emptysolt,1);
    for kk=1:num_empty
        mid_mystime=max(my_can_time,emptysolt(kk,1));
        if isempty(midschedule(midschedule(:,6)==emptysolt(kk,1),1))==1||midschedule(midschedule(:,6)==emptysolt(kk,1),1)~=i
            if emptysolt(kk,2)-mid_mystime>=my_ptime+setuptime(i)
                mid_schedule=[i,0,my_machineid,my_machinev,mid_mystime,mid_mystime+setuptime(i)];
                schedule=[schedule;mid_schedule];
                my_stime=mid_mystime+setuptime(i);
                my_ctime=my_stime+my_ptime;
                mid_schedule=[i,job_operation(i),my_machineid,my_machinev,my_stime,my_ctime];
                schedule=[schedule;mid_schedule];
                job_pre_endtime(i)=my_ctime;
                break;
            end
        else
            if emptysolt(kk,2)-mid_mystime>=my_ptime
                v1=midschedule(midschedule(:,6)==emptysolt(kk,1),4);
                switch v1/my_machinev
                    case{1}
                        v=v1;
                        ED=ED+0;
                    case{2,1/2}
                        v = 1;
                        ED = ED+Eshift(v);
                    case{2/3,3/2}
                        v=2;
                        ED = ED+Eshift(v);
                    case{3,1/3}
                        v=3;
                        ED = ED+Eshift(v);
                end
                my_stime=mid_mystime;
                my_ctime=my_stime+my_ptime;
                mid_schedule=[i,job_operation(i),my_machineid,my_machinev,my_stime,my_ctime];
                schedule=[schedule;mid_schedule];
                job_pre_endtime(i)=my_ctime;
                break;
            end
        end
    end
end
addmschedule = [];
mschedule = {};
EA = 0; 
EB = 0; 
EC = 0;
ED = 0;
EE = 0; 
for m = 1:TM
    midschedule = schedule(schedule(:,3)==m,:);
    midschedule = sortrows(midschedule,5); 
    pre_endtime = 0;
    for index = 1:size(midschedule,1)
        start_hour = mod(midschedule(index,5),24); 
        end_hour = mod(midschedule(index,6),24); 
        mid_hour=midschedule(index,6)-midschedule(index,5);
        current_price = getElectricityPrice(start_hour); 
        if midschedule(index,2)==0
            for h = start_hour:start_hour+setuptime(midschedule(index,1))
                if h>24
                    h=mod(h,24);
                end
                period_price = getElectricityPrice(h);
                EA = EA +  PA * period_price ; 
            end
        else
            for h = start_hour:start_hour+mid_hour
                if h>24
                    h=mod(h,24);
                end
                period_price = getElectricityPrice(h);
                EB = EB + PB(midschedule(index,4)) * period_price ; 
            end
        end
        if index>1
            if midschedule(index,5)> pre_endtime
                midslot = [pre_endtime midschedule(index,5)];       
                idlev = 0;
                Eidle=0;
                Estandby=0;
                mid_pre_endtime=mod(pre_endtime,24);
                mid=midschedule(index,5)-pre_endtime;
                if midschedule(index-1,4)/midschedule(index,4)==1
                    idlev = midschedule(index-1,4);
                    for h=mid_pre_endtime:mid_pre_endtime+mid
                        if h>24
                            h=mod(h,24);
                        end
                        current_price= getElectricityPrice(h);
                        Eidle = Eidle+Pidle(idlev)*current_price;
                    end
                elseif midschedule(index-1,4)/midschedule(index,4)==2 || midschedule(index-1,4)/midschedule(index,4)==0.5
                    idlev = 1;
                    for h=mid_pre_endtime:mid_pre_endtime+mid
                        if h>24
                            h=mod(h,24);
                        end
                        current_price= getElectricityPrice(h);
                        Eidle = Eidle+Pidle(1)*current_price + Eshift(idlev);
                    end
                elseif midschedule(index-1,4)/midschedule(index,4)==2/3 || midschedule(index-1,4)/midschedule(index,4)==3/2
                    idlev = 2;
                    for h=mid_pre_endtime:mid_pre_endtime+mid
                        if h>24
                            h=mod(h,24);
                        end
                        current_price= getElectricityPrice(h);
                        Eidle = Eidle+Pidle(2)*current_price + Eshift(idlev);
                    end
                else
                    idlev = 3;
                    for h=mid_pre_endtime:mid_pre_endtime+mid
                        if h>24
                            h=mod(h,24);
                        end
                        current_price= getElectricityPrice(h);
                        Eidle =Eidle+Pidle(3)*current_price + Eshift(idlev);
                    end
                end
                for h=mid_pre_endtime:mid_pre_endtime+mid
                    if h>24
                        h=mod(h,24);
                    end
                    current_price= getElectricityPrice(h);
                    Estandby = Estandby+Pstandby*current_price + Estandby2process(midschedule(index,4));
                end
                if Eidle <= Estandby
                    EC = EC + Eidle;
                    a = [0 0 m idlev pre_endtime midschedule(index,5)];
                    addmschedule = [addmschedule;a];
                else
                    EC = EC + Estandby;
                    a = [0 0 m 0 pre_endtime midschedule(index,5)];
                    addmschedule = [addmschedule;a];
                end
            end
        end
        pre_endtime = midschedule(index,6);
    end
    mschedule{1,m} = [midschedule;addmschedule];
    mschedule{1,m} = sortrows(mschedule{1,m},5);
    alls = [alls;mschedule{1,m}];
end
EE = 0;
for m = 1:TM
    index = find(mschedule{1,m}(:,4)~=0,1);
    if isempty(index)
        EE = EE + 0;
    else
        EE = EE + Eturnonfirst(mschedule{1,m}(index,4));
    end
end
E = EA + EB + EC + ED + EE;
schedule2=schedule(schedule(:,2)~=0,:);
[~,index]= sort(schedule2(:,5));
os_newchrom= schedule2(index,1)';
os_chrom = os_newchrom ;
Cmax = max(schedule(:,6));
TEC = E;
end