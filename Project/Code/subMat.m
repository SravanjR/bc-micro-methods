%Creates estimated Ownership matrix for a given year using demand
%coefficents, the year, the parents and the estimated marketshare
function DELTA = subMat(alpha,sigma,share,withinGroupShare,Year,CarrID,Nest,nObs)
ownPriceDeriv = alpha*share.*((1/(1-sigma))-(sigma/(1-sigma))*withinGroupShare - share);
DELTA = zeros(nObs,nObs,'single');
for j = 1:nObs 
    tempnest = Nest(j);
    tempID = CarrID(j);
    tempyear = Year(j);
    DELTA(:,j) = DELTA(:,j) + (tempnest == Nest).* cast(-alpha*share.*((sigma/(1-sigma))*withinGroupShare(j) + share(j)),'single');
    DELTA(:,j) = DELTA(:,j) + (tempnest ~= Nest).* cast(-alpha.*share(j).*share,'single');
    DELTA(:,j) = DELTA(:,j).*(tempID == CarrID).*(tempyear == Year);
    DELTA(j,j) = cast(ownPriceDeriv(j),'single');
end
end