function Fitness= classifyIndividuals(fitness)
global PS 
         min1 = min(fitness(:, 1));
         max1 = max(fitness(:, 1));
         min2 = min(fitness(:, 2));
         max2 = max(fitness(:, 2));
 
         fitness_norm = zeros(size(fitness));
         fitness_norm(:, 1) = (fitness(:, 1) - min1) / (max1 - min1);
         fitness_norm(:, 2) = (fitness(:, 2) - min2) / (max2 - min2);
     
         angles = atan2(fitness_norm(:, 2), fitness_norm(:, 1));
    
         angles = mod(angles, 2*pi);
     
         angles_deg = angles * (180 / pi);
       
         mid=size(fitness,1);
         Fitness = zeros(mid, 1);
       
         for i = 1:PS
             if angles_deg(i) >= 0 && angles_deg(i) < 30
                 Fitness(i) = 1; 
             elseif angles_deg(i) >= 30 && angles_deg(i) < 60
                 Fitness(i) = 2; 
             elseif angles_deg(i) >= 60 && angles_deg(i) <= 90
                 Fitness(i) = 3; 
             end
         end
end