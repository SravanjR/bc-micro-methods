function [objfun, newWeightMat, varcovar] = project_GMM2(y,mc,p,d,x1,x2,z1,z2,params, wMat, isConstrained)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Methods
%  GMM function file with Supply and Demand
% Output Variables:
%   objfun is the value of the GMM objective function
%   newWeightMat is self explanatory
%   varcovar is the estimated variance covariance matrix
%        note that this is calculated AFTER updating the weighting matrix,
%        based on the updated weighting matrix
% Input Variables:
%  y is the quality vector
%   p is price col vector
%   mc is the marginal cost
%   d is dept col vector
%   x1 is covariates matrix for demand
%   x2 is the covariates matrix for supply
%   z1 is the matrix of instruments for demand
%   z2 is the matrix of instruments for supply
%   params are the parameters
%   wMat is the weighting matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigmaInd = 15; %subscript index of sigma coefficient
priceInd = 16; %subscript index of price coefficient
deptInd = 17; %subscript index of deptature coefficient

if isConstrained == 1 %constrain the nesting parameter to be in (0,1)
    params(sigmaInd) = exp(params(sigmaInd)) / (1 + exp(params(sigmaInd)));
end

nObs = size(p,1);
nMoms = size(wMat,1);

% First set of moments (demand)
xi = y - (x1*params(1:sigmaInd) + params(priceInd)*p + params(deptInd)*d);
mom1 = z1'*xi/nObs;

% Second set of moments (supply)
g = params(18:30);
omega = mc-(x2*g);
mom2 = z2'*omega/nObs;

% Combine the moment conditions
moms = [mom1 ; mom2];
if nargout == 1
    objfun = moms'*wMat*moms;
elseif nargout == 2
    objfun = moms'*wMat*moms;
    temp = [(z1.*repmat(xi,1,size(z1,2))) (z2.*repmat(omega,1,size(z2,2)))];
    newWeightMat = inv(temp'*temp/nObs);    
else
    objfun = moms'*wMat*moms;
    temp = [(z1.*repmat(xi,1,size(z1,2))) (z2.*repmat(omega,1,size(z2,2)))];
    newWeightMat = inv(temp'*temp/nObs);    
    % Need to compute derivative matrix
    gamma = [(z1'*[-x1 -p -d]) zeros(nMoms/2,13); (zeros(nMoms/2,17)) z2'*[-x2]]/nObs;
    varcovar = inv(gamma'*newWeightMat*gamma)/(nObs);
end



