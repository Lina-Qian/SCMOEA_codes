function isHighPrice = isElectricityPriceHigh(hour)
isHighPrice = false;
if hour > 10 && hour < 15
    isHighPrice = true;
end
end