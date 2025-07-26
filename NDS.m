function result = NDS(fit1, fit2)
    if ~isequal(size(fit1), [1, 2]) || ~isequal(size(fit2), [1, 2])
        error('输入必须是1x2向量');
    end         
    if all(fit1 <= fit2) && any(fit1 < fit2)
        result = 1;
    elseif all(fit2 <= fit1) && any(fit2 < fit1)
        result = 2; 
    else
        result = 0; 
    end
end