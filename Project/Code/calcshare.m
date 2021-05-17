%Estimate the Share using Demand Coefficients
function sol = calcshare(beta, X)
    num = exp(X*beta);
    denom = 1 + sum(exp(X*beta));
    sol = num./denom;
end