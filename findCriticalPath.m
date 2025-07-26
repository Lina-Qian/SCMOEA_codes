function cpath = findCriticalPath(schedule)
 
    for i = 1:size(schedule,1)
        if schedule(i,2)==0
            schedule(i+1,5) = schedule(i,5);
        end
    end
    schedule(schedule(:,2)==0,:) = [];
   
    [~,im]=max(schedule(:,6));
    starti = schedule(im,5);
    cpath(1,1) = schedule(im,1);
    cpath(1,2) = schedule(im,2);
    k = 1;
    while(starti~=0)
       
        if schedule(im,2)~=1 
            im1 = find(schedule(:,1)==schedule(im,1)&schedule(:,2)==(schedule(im,2)-1));
            cpj = schedule(im1,6);
        else
            cpj = 0;
        end
     
        if find(schedule(schedule(:,3)==schedule(im,3),5)<schedule(im,5))~=-1
            im2 = find(schedule(:,3)==schedule(im,3)&schedule(:,6)<=schedule(im,5),1, 'last');
            cpm = schedule(im2,6);
        else
            cpm = 0;
        end
        if cpj == 0 & cpm == 0
            starti = 0;
        elseif cpj>cpm
            im = im1;
            cpath(k+1,1) = schedule(im,1);
            cpath(k+1,2) = schedule(im,2);
            starti = schedule(im,5);
        else
            im = im2;
            cpath(k+1,1) = schedule(im,1);
            cpath(k+1,2) = schedule(im,2);
            starti = schedule(im,5);
        end
        k = k+1;
    end
end