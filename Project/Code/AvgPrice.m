%Calculate the Average Price for a firm given the within firm shares and
%sales and the respective fares.
function sol = AvgPrice(p,shr)
    sol = sum(p.*shr);
end
