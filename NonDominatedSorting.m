function [pop, F] = NonDominatedSorting(pop)
nPop = numel(pop);
for i = 1:nPop
    pop(i).DominationSet = [];
    pop(i).DominatedCount = 0;
end
F = {};
F{1} = [];
for i = 1:nPop
    for j = i+1:nPop
        if Dominates(pop(i).Fitness, pop(j).Fitness)
            pop(i).DominationSet = [pop(i).DominationSet, j];
            pop(j).DominatedCount = pop(j).DominatedCount + 1;
        elseif Dominates(pop(j).Fitness, pop(i).Fitness)
            pop(j).DominationSet = [pop(j).DominationSet, i];
            pop(i).DominatedCount = pop(i).DominatedCount + 1;
        end
    end
    if pop(i).DominatedCount == 0
        F{1} = [F{1}, i];
        pop(i).Rank = 1;
    end
end

k = 1;
while true
    Q = [];
    for i = F{k}
        p = pop(i);
        for j = p.DominationSet
            q = pop(j);
            q.DominatedCount = q.DominatedCount - 1;
            if q.DominatedCount == 0
                Q = [Q, j]; 
                q.Rank = k + 1;
            end
            pop(j) = q;
        end
    end
    if isempty(Q)
        break;
    end
    F{k+1} = Q;
    k = k + 1;
end

for i = 1:nPop
    if ~isfield(pop(i), 'Rank')
        pop(i).Rank = k;
    end
end
end