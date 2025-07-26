function isDominate = Dominates(fit1, fit2)
    isBetter = all(fit1 <= fit2) && any(fit1 < fit2);
    isDominate = isBetter;
end