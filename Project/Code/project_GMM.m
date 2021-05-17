function [objfun, newWeightMat, varcovar] = project_GMM(y,p,d,x,z,params, wMat, isConstrained)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Methods
% Nested Logit function file 
% Output Variables:
%   objfun is the value of the GMM objective function
%   newWeightMat is self explanatory
%   varcovar is the estimated variance covariance matrix
%        note that this is calculated AFTER updating the weighting matrix,
%        based on the updated weighting matrix
% Input Variables:
%   y is the nested logit LHS variable 
%   p is price col vector
%   d is the departures col vector
%   x is covariates matrix
%   z is the matrix of covariates and excluded instruments
%   params are the parameters
%   wMat is the weighting matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigmaInd = 15; %subscript index of sigma coefficient
priceInd = 16; %subscript index of price coefficient
deptInd = 17; %subscript index of deptature coefficient

if isConstrained == 1 %constrain the nesting parameter to be in (0,1)
    params(sigmaInd) = exp(params(sigmaInd)) / (1 + exp(params(sigmaInd)));
end

nObs = size(p,1);
nMoms = size(wMat,1);

xi = y - (x*params(1:sigmaInd) + params(priceInd)*p + params(deptInd)*d);
moms = z'*xi/nObs;

if nargout == 1
    objfun = moms'*wMat*moms;
elseif nargout == 2
    objfun = moms'*wMat*moms;
    temp = [z.*repmat(xi,1,nMoms)];
    newWeightMat = inv(temp'*temp/nObs);    
else
    objfun = moms'*wMat*moms;
    temp = [z.*repmat(xi,1,nMoms)];
    newWeightMat = inv(temp'*temp/nObs);    
    % Need to compute derivative matrix
    gamma = z'*[-x -p -d]/nObs;
    varcovar = inv(gamma'*newWeightMat*gamma)/nObs;
        %best to avoid using inv when possible.
end