%Calculate Profits for all parents in a year using marginal cost
%coefficients, price and marketshares
function sol = calcprofit(p,shr,marketsize,mc)
    sol = sum((p - mc).*(shr.*marketsize));
    if sol < 0
        disp('hi');
    end
end