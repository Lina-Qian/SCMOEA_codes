function price = getElectricityPrice(hour)
    if hour >= 0 && hour < 7
        price = 0.4;
    elseif hour >= 7 && hour < 10
        price = 0.8;
    elseif hour >= 10 && hour < 15
        price = 1.3;
    elseif hour >= 15 && hour < 18
        price = 0.8;
    elseif hour >= 18 && hour < 21
        price = 1.3;
    elseif hour >= 21 && hour < 23
        price = 0.8;
    else % hour >= 23 || hour < 0
        price = 0.4;
    end
end