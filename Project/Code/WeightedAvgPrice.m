%Computes the Weighted Average Price
function sol = WeightedAvgPrice(x,w)
    sol = sum(w.*x)./sum(w);
end