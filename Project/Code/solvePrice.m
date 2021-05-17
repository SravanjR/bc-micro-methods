%Price Estimate using Firm price FOCs. Uses the estimated ownership matrix for a
%given year, the estimated market share, the marginal cost coefficient and
%the marginal cost residual
function sol = solvePrice(mrkup,mc)
    sol = mc + mrkup;
end