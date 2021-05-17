%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Methods
% Project - Main File
% Code Outline: The code consists of nested logit. The code allows for various
%   normalizations of the covariates and instruments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
disp(' '); % to separate results from previous executions
diary Work.txt
format short g

% Set input values 

NORM_STYLE = 1; % use 0 for not normalizing at all
                % use 1 for dividing by mean
                % use 2 for dividing by standard deviation
                % use 3 for substracting mean and then dividing by standard
                % deviation
NORM_EXTENT = 1; % use 1 for normalizing before creating instruments and not
                 % normalizing anything afterward
                 % use 2 for normalizing the continuous covariates (not 
                 % quantity, of course, and not logWithinGroupShare)
                 % after creating the instruments, without normalizing the 
                 % instruments
                 % use 3 for normalizing the continuous covariates (not 
                 % quantity, of course, and not logWithinGroupShare) and the
                 % instruments, after creating the instruments
NORM_NESTED = 1; % use 1 to make no changes to logWithinGroupShare
                 % use 2 to normalize logWithinGroupShare
                 
% Label output clearly
disp('Empirical Methods Project');                 

% Load and separate data
load("Premerger_db1b")
load("Premerger_MkDist")
Premerger_db1b = sortrows(Premerger_db1b,30);
Premerger_MkDist = sortrows(Premerger_MkDist,1);
Premerger_db1b.Distsq = Premerger_db1b.Dist.^2;
NUM_GROUPS = size(unique(Premerger_db1b.Nest),1); % controls the number of groups used in nested logit
MKT_SIZE= Premerger_db1b.AvgPop;  % set market size to the geometric mean of all endpoint cities' MSA Populations


% Store parameter values
nObs=size(Premerger_db1b.Year,1);
T=size(unique(Premerger_db1b.Year),1);

% Calculate market shares
share =Premerger_db1b.ProdSh;

% Calculate outside good shares for each year
outsideShare=zeros(nObs,1);
yearList=unique(Premerger_db1b.Year);
for j=1:T
    for i=1:nObs
    if Premerger_db1b.Year(i)==yearList(j)
        outsideShare(i)=1-sum(share(Premerger_db1b.Year==yearList(j)),1);
    end
    end
end

if NORM_EXTENT == 1
    switch NORM_STYLE
        case 0
            normPrice = Premerger_db1b.BinFare;
            normDept = Premerger_db1b.dept;
            normDist = Premerger_db1b.Dist;
            normDistSq = Premerger_db1b.Distsq;
        case 1
            normPrice = Premerger_db1b.BinFare/mean(Premerger_db1b.BinFare);
            normDept = Premerger_db1b.dept/mean(Premerger_db1b.dept);
            normDist = Premerger_db1b.Dist/mean(Premerger_db1b.Dist);
            normDistSq = Premerger_db1b.Distsq/mean(Premerger_db1b.Distsq);
        case 2
            normPrice = Premerger_db1b.BinFare / std(Premerger_db1b.BinFare);
            normDept = Premerger_db1b.dept/ std(Premerger_db1b.dept);
            normDist = Premerger_db1b.Dist/std(Premerger_db1b.Dist);
            normDistSq = Premerger_db1b.Distsq/std(Premerger_db1b.Distsq);
        case 3
            normPrice = (Premerger_db1b.BinFare - mean(Premerger_db1b.BinFare)) / std(Premerger_db1b.BinFare);
            normDept = (Premerger_db1b.dept - mean(Premerger_db1b.dept)) / std(Premerger_db1b.dept);
            normDist = (Premerger_db1b.Dist- mean(Premerger_db1b.Dist)) / std(Premerger_db1b.Dist);
            normDistSq = (Premerger_db1b.Distsq- mean(Premerger_db1b.Distsq)) / std(Premerger_db1b.Distsq);
    end

     %% Create instruments
firmList = unique(Premerger_db1b.CarrID);
yearList=unique(Premerger_db1b.Year);
nFirms = size(firmList,1);
insts = zeros(nObs,7); 
    for iFirm=1:nFirms 
        for iYear=1:T
        isThisFirm=(Premerger_db1b.CarrID==firmList(iFirm)); % flag characteristics of the same firm
        isThisYear=(Premerger_db1b.Year==yearList(iYear)); % flag characteristics in the relevant year
        isThisObs=logical(isThisFirm.*isThisYear); % flag current obs
        isValidObs=logical((1-isThisFirm).*isThisYear); % flag characteristics of competing products in a given year
          % create ivs using flagged observations
        insts(isThisObs,1)=(sum(isValidObs.*Premerger_db1b.NoCarr))./sum(isValidObs); 
        insts(isThisObs,2)=(sum(isValidObs.*Premerger_db1b.Dist))./sum(isValidObs);
        end
    end
        insts(:,3) =  Premerger_db1b.Hubdest;
        insts(:,4) =  Premerger_db1b.NdestEnd;
        insts(:,5) =  Premerger_db1b.P25;
        insts(:,6) =  Premerger_db1b.P75;
        insts(:,7) =  Premerger_db1b.fit_dept;
        normInsts=insts;
        
% covars = [ones(nObs,1) Premerger_db1b.NumDest Premerger_db1b.Dist Premerger_db1b.Distsq Premerger_db1b.FLLAS Premerger_db1b.Slot Premerger_db1b.OtherCarr Premerger_db1b.American Premerger_db1b.Continental ...
%     Premerger_db1b.Delta Premerger_db1b.United Premerger_db1b.US_Airways Premerger_db1b.JetBlue Premerger_db1b.Southwest]; % covariates
% Z = [covars normInsts];
else %NORM_EXTENT ==  2 | 3
     % Create instruments
firmList = unique(Premerger_db1b.CarrID);
yearList=unique(Premerger_db1b.Year);
nFirms = size(firmList,1);
insts = zeros(nObs,7); %instruments
    for iFirm=1:nFirms 
        for iYear=1:T
        isThisFirm=(Premerger_db1b.CarrID==firmList(iFirm)); % flag characteristics of the same firm
        isThisYear=(Premerger_db1b.Year==yearList(iYear)); % flag characteristics in the relevant year
        isThisObs=logical(isThisFirm.*isThisYear);  % flag current obs
        isValidObs=logical((1-isThisFirm).*isThisYear); % flag characteristics of competing products in a given year
          % create ivs using flagged observations
        insts(isThisObs,1)=(sum(isValidObs.*Premerger_db1b.NoCarr))./sum(isValidObs); 
        insts(isThisObs,2)=(sum(isValidObs.*Premerger_db1b.Dist))./sum(isValidObs);
        end
    end
        insts(:,3) =  Premerger_db1b.Hubdest;
        insts(:,4) =  Premerger_db1b.NdestEnd;
        insts(:,5) =  Premerger_db1b.P25;
        insts(:,6) =  Premerger_db1b.P75;
        insts(:,7) =  Premerger_db1b.fit_dept;
 
    switch NORM_STYLE
        case 0
            normPrice = Premerger_db1b.BinFare;
            normDept = Premerger_db1b.dept;
            if NORM_EXTENT == 3 % normalize the instruments
                normInsts = insts ./ (ones(nObs,1) * mean(insts,1));
            else
                normInsts = insts;
            end
            
        case 1
            normPrice = Premerger_db1b.BinFare/mean(Premerger_db1b.BinFare);
            normDept = Premerger_db1b.dept/mean(Premerger_db1b.dept);
            normDist = Premerger_db1b.Dist/mean(Premerger_db1b.Dist);
            normDistSq = Premerger_db1b.Distsq/mean(Premerger_db1b.Distsq);
            if NORM_EXTENT == 3 % normalize the instruments
                normInsts = insts ./ (ones(nObs,1) * mean(insts,1));
            else
                normInsts = insts;
            end
        case 2
            normPrice = Premerger_db1b.BinFare / std(Premerger_db1b.BinFare);
            normDept = Premerger_db1b.dept/std(Premerger_db1b.dept);
            if NORM_EXTENT == 3 % normalize the instruments
                normInsts = insts ./ (ones(nObs,1) * std(insts,1,1));
            else
                normInsts = insts;
            end
        case 3
            normPrice = (Premerger_db1b.BinFare - mean(Premerger_db1b.BinFare)) / std(Premerger_db1b.BinFare);
            normDept = (Premerger_db1b.dept - mean(Premerger_db1b.dept)) / std(Premerger_db1b.dept);
            if NORM_EXTENT == 3 % normalize the instruments
                normInsts = (insts - (ones(nObs,1) * mean(insts,1))) ...
                    ./ (ones(nObs,1) * std(insts,1,1));
            else
                normInsts = insts;
            end
    end

%    covars = [ones(nObs,1) Premerger_db1b.NumDest normDist normDistSq Premerger_db1b.FLLAS Premerger_db1b.Slot Premerger_db1b.OtherCarr Premerger_db1b.American Premerger_db1b.Continental ...
%     Premerger_db1b.Delta Premerger_db1b.United Premerger_db1b.US_Airways Premerger_db1b.JetBlue Premerger_db1b.Southwest]; % covariates
%  % covariates
%     Z = [covars normInsts];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NESTED LOGIT MODEL   %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let nesting variable vary based on the number of nests
groupId=Premerger_db1b.Nest;

% Calculate within group share 
groupQ = zeros(size(Premerger_db1b.CarrID));
% Construct 'time invariant' nests
for iGroup = 1:NUM_GROUPS
    groupIndVec = find(ismember(groupId,iGroup)); % identify time invariant group members
    groupQ(groupIndVec) = sum(Premerger_db1b.total_passengers_by_product(groupIndVec)); % calculate group quantity and record value for each group member
end

withinGroupShare =  Premerger_db1b.total_passengers_by_product ./ groupQ;
logWithinGroupShare = log(withinGroupShare);
clear groupQ;

% Identify group members by creating a binary matrix 
% such that (i,j) = 1 if i and j are in the same group
% Use this matrix for iv construction
groupIdMat = zeros(nObs,nObs,'int16');
for j = 1:nObs
    tempgroupID = groupId(j);
    tempYear = Premerger_db1b.Year(j);
    groupIdMat(:,j) = (tempgroupID == groupId & tempYear == Premerger_db1b.Year);
end

instGroupIdMat = groupIdMat - eye(nObs,'int16'); % remove product itself from the 
    % group that will be used to construct an instrument for 
    % its within group share. The matrix is symmetric, but I think about
    % row (i) as containing 1 in column j whenever j should be used to
    % construct the instrument for i.

clear tempgroupID isGroupMember;
 
% priceInsts = normInsts;

% Choose to normalize log within group share
if NORM_NESTED == 1
    normLogWithinGroupShare = logWithinGroupShare;
else
    switch NORM_STYLE
        case 0
            normLogWithinGroupShare = logWithinGroupShare; 
        case 1
            normLogWithinGroupShare = logWithinGroupShare / mean(logWithinGroupShare);
        case 2
            normLogWithinGroupShare = logWithinGroupShare / std(logWithinGroupShare);
        case 3
            normLogWithinGroupShare = (logWithinGroupShare - mean(logWithinGroupShare)) ...
                / std(logWithinGroupShare);
    end
end
    
if NORM_EXTENT == 1 
    % create within group share instruments from normalized demand covariates
    wInGroupShareInsts(:,1) = sum(repmat(cast(Premerger_db1b.NumDest,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,2) = sum(repmat(cast(Premerger_db1b.Dist,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,3) = sum(repmat(cast(Premerger_db1b.Distsq,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,4) = sum(repmat(cast(Premerger_db1b.FLLAS,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,5) = sum(repmat(cast(Premerger_db1b.Slot,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,6) = sum(repmat(cast(Premerger_db1b.OtherCarr,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,7) = sum(repmat(cast(Premerger_db1b.American,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,8) = sum(repmat(cast(Premerger_db1b.Continental,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,9) = sum(repmat(cast(Premerger_db1b.Delta,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,10) = sum(repmat(cast(Premerger_db1b.United,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,11) = sum(repmat(cast(Premerger_db1b.US_Airways,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,12) = sum(repmat(cast(Premerger_db1b.JetBlue,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,13) = sum(repmat(cast(Premerger_db1b.Southwest,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
else %NORM_EXTENT == 2 | 3
    % create within group share instruments from non-normalized demand covariates
    % create within group share instruments from normalized demand covariates
    wInGroupShareInsts(:,1) = sum(repmat(Premerger_db1b.NumDest',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,2) = sum(repmat(Premerger_db1b.Dist',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,3) = sum(repmat(Premerger_db1b.Distsq',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,4) = sum(repmat(Premerger_db1b.FLLAS',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,5) = sum(repmat(Premerger_db1b.Slot',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,6) = sum(repmat(Premerger_db1b.OtherCarr',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,7) = sum(repmat(Premerger_db1b.American',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,8) = sum(repmat(Premerger_db1b.Continental',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,9) = sum(repmat(Premerger_db1b.Delta',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,10) = sum(repmat(Premerger_db1b.United',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,11) = sum(repmat(Premerger_db1b.US_Airways',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,12) = sum(repmat(Premerger_db1b.JetBlue',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts(:,13) = sum(repmat(Premerger_db1b.Southwest',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    if NORM_EXTENT == 3
        switch NORM_STYLE
            case 1 % divide by mean
                wInGroupShareInsts = wInGroupShareInsts ./ ...
                    (ones(nObs,1) * mean(wInGroupShareInsts,1));
            case 2 % divide by std
                wInGroupShareInsts = wInGroupShareInsts ./ ...
                    (ones(nObs,1) * std(wInGroupShareInsts,1));
            case 3 % substract mean, divide by std
                wInGroupShareInsts = (wInGroupShareInsts - ...
                    (ones(nObs,1) * mean(wInGroupShareInsts))) ...
                    ./ (ones(nObs,1) * std(wInGroupShareInsts,1));
        end
    end
end
%%
FUN_TOL = 1e-17;
% create variables for GMM function
covars = [ones(nObs,1) Premerger_db1b.NumDest normDist normDistSq Premerger_db1b.FLLAS Premerger_db1b.Slot Premerger_db1b.OtherCarr Premerger_db1b.American Premerger_db1b.Continental ...
    Premerger_db1b.Delta Premerger_db1b.United Premerger_db1b.US_Airways Premerger_db1b.JetBlue Premerger_db1b.Southwest logWithinGroupShare]; % covariates
sigmaInd = 15;
priceInd = 16;
deptInd = 17;
z = [covars wInGroupShareInsts normInsts];
y = log(share) - log(outsideShare);

% 2SLS Estimation 
x = [covars normPrice normDept];
theta2sls = real((z'*x)\(z'*y));

% GMM Estimation (constraints on sigma)
% Find an initial guess for the parameters, using the identity matrix as
% the weight.
isConstrained = 1;
weightMat = eye(size(z,2));
thetaStart = theta2sls;
thetaStart(sigmaInd) = 0.5;
options = optimset('MaxFunEvals', 1000000, 'Display', 'none', 'TolFun', FUN_TOL);

thetaNested = fminsearch(@(theta) project_GMM(y,normPrice,normDept,covars,z,theta, weightMat, isConstrained),thetaStart, options);

% Optimal weight matrix
[temp, optw] = project_GMM(y,normPrice,normDept,covars,z,thetaNested, weightMat, isConstrained);

% Go again with optimal matrix
thetaNested = fminsearch(@(theta) project_GMM(y,normPrice,normDept,covars,z, theta, optw, isConstrained), thetaNested, options);

[obj1, ~, v1] = project_GMM(y,normPrice,normDept,covars,z, thetaNested, optw, isConstrained);

se = real(abs(diag(v1)));
demandbetasnomc = thetaNested;
demandsenomc = sqrt(real(se));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now to set up GMM for the  Multi-Product Firms in a Bertrand Nash Equilibrium

%%
%Price Derivative Matrix
sigma = thetaNested(sigmaInd);
alpha = thetaNested(priceInd);
%DELTA = subMat(alpha, sigma, share,withinGroupShare, Premerger_db1b.Year, Premerger_db1b.CarrID, Premerger_db1b.Nest,nObs);
markup = -(subMat(alpha, sigma, share,withinGroupShare, Premerger_db1b.Year, Premerger_db1b.CarrID, Premerger_db1b.Nest,nObs)\cast(share,'single'));
clear DELTA se
markup = cast(markup,'double');
mc = Premerger_db1b.BinFare - markup;

%%
if NORM_EXTENT == 1
    switch NORM_STYLE
        case 0
            normmc = mc;
        case 1
            normmc = mc/mean(mc);
        case 2
            normmc = mc / std(mc);
        case 3
            normmc = (mc - mean(mc)) / std(mc);
    end
end
%%
%Within group share supply instruments
 if NORM_EXTENT == 1 
    % create within group share instruments from covariates for supply
    wInGroupShareInsts2(:,1) = sum(repmat(cast(Premerger_MkDist.SmDist,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts2(:,2) = sum(repmat(cast(Premerger_MkDist.LgDist,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts2(:,3) = sum(repmat(cast(Premerger_db1b.HubMC,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
else %NORM_EXTENT == 2 | 3
    % create within group share instruments from non-normalized covariates
    % create within group share instruments from normalized covariates
    wInGroupShareInsts2(:,1) = sum(repmat(cast(Premerger_MkDist.SmDist,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts2(:,2) = sum(repmat(cast(Premerger_MkDist.LgDist,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    wInGroupShareInsts2(:,3) = sum(repmat(cast(Premerger_db1b.HubMC,'int16')',nObs,1) .*instGroupIdMat,2) ...
        ./ sum(instGroupIdMat,2);
    if NORM_EXTENT == 3
        switch NORM_STYLE
            case 1 % divide by mean
                wInGroupShareInsts2 = wInGroupShareInsts2 ./ ...
                    (ones(nObs,1) * mean(wInGroupShareInsts2,1));
            case 2 % divide by std
                wInGroupShareInsts2 = wInGroupShareInsts2 ./ ...
                    (ones(nObs,1) * std(wInGroupShareInsts2,1));
            case 3 % substract mean, divide by std
                wInGroupShareInsts2 = (wInGroupShareInsts2 - ...
                    (ones(nObs,1) * mean(wInGroupShareInsts2))) ...
                    ./ (ones(nObs,1) * std(wInGroupShareInsts2,1));
        end
    end
end
 %%
% create variables for GMM function 
 %Supply Dependant Variable
 mc = cast(mc,'double');
 %Demand Covariates
 covars = [ones(nObs,1) Premerger_db1b.NumDest normDist normDistSq Premerger_db1b.FLLAS Premerger_db1b.Slot Premerger_db1b.OtherCarr Premerger_db1b.American Premerger_db1b.Continental ...
    Premerger_db1b.Delta Premerger_db1b.United Premerger_db1b.US_Airways Premerger_db1b.JetBlue Premerger_db1b.Southwest logWithinGroupShare]; % covariates
%Supply Covariates 
covars2 = [ones(nObs,1) Premerger_MkDist.SmDist Premerger_MkDist.LgDist Premerger_db1b.HubMC Premerger_db1b.Slot ...
     Premerger_db1b.OtherCarr Premerger_db1b.American Premerger_db1b.Continental Premerger_db1b.Delta Premerger_db1b.United Premerger_db1b.US_Airways Premerger_db1b.JetBlue Premerger_db1b.Southwest];
%Covariate Set
 covars3 = [ones(nObs,1) Premerger_db1b.NumDest normDist normDistSq Premerger_db1b.FLLAS Premerger_db1b.Slot ...
     Premerger_MkDist.SmDist Premerger_MkDist.LgDist Premerger_db1b.HubMC ...
     Premerger_db1b.OtherCarr Premerger_db1b.American Premerger_db1b.Continental Premerger_db1b.Delta Premerger_db1b.United Premerger_db1b.US_Airways Premerger_db1b.JetBlue Premerger_db1b.Southwest];
z = [covars3 wInGroupShareInsts wInGroupShareInsts2 normInsts];
thetaStart2 = cast((z'*covars2)\(z'*normmc),'double');
%%
FUN_TOL = 1e-17;
% GMM Estimation (constraints on sigma)
% Find an initial guess for the parameters, using the identity matrix as
% the weight.
isConstrained = 1;
weightMat = eye(size(z,2)*2);
thetaStart3 = [thetaNested;thetaStart2];
options = optimset('MaxFunEvals', 1000000, 'Display', 'none', 'TolFun', FUN_TOL);

thetaNested2 = fminsearch(@(theta) project_GMM2(y,normmc,normPrice,normDept,covars,covars2,z,z,theta, weightMat, isConstrained),thetaStart3, options);

% Optimal weight matrix
[temp, optw] = project_GMM2(y,normmc,normPrice,normDept,covars,covars2,z,z, thetaNested2, weightMat, isConstrained);

% Go again with optimal matrix
thetaNested3 = fminsearch(@(theta) project_GMM2(y,normmc,normPrice,normDept,covars,covars2,z,z, theta, optw, isConstrained), thetaNested2, options);

%Variance-Covariance Matrix
[obj2, ~, v1] = project_GMM2(y,normmc,normPrice,normDept,covars,covars2,z,z, thetaNested3, optw, isConstrained);
se = real(abs(diag(v1)));
demandbetas = thetaNested3(1:17);
demandse = sqrt(real(se(1:17)));
supplybetas = thetaNested3(18:30);
supplyse = sqrt(real(se(18:30)));
writematrix(demandbetas, "demandbetas.txt")
writematrix(demandse, "demandse.txt")
writematrix(supplybetas, "supplybetas.txt")
writematrix(supplyse, "supplyse.txt")
demandbetas = load("demandbetas.txt");
demandse = load("demandse.txt");
supplybetas = load("supplybetas.txt");
supplyse = load("supplyse.txt");
%%
% Demand and Supply Estimates with Real Standard Errors
disp(['Nested Logit Model with ', num2str(NUM_GROUPS), ' groups (Constrained), demand only, 2-step GMM']);
disp(['Premerger']);
disp('       Const       NumDest        Dist        Distsq        FLLAS       Slot  OtherCarr     American   Continental     Delta     United     US Airways  JetBlue      Southwest     Sigma     Fare     Departures');
disp(['Coef:', num2str(demandbetas')]);
disp(['SE:   ', num2str(demandse')]);
disp('-------------------');

% Estimates    Standard Errors
disp(['Nested Logit Model with ', num2str(NUM_GROUPS), ' groups (Constrained), Supply, 2-step GMM']);
disp(['Premerger']);
disp('       Const       SmDist        LgDist        HubMC        Slot           OtherCarr     American      Continental      Delta      United       US Airways      JetBlue       Southwest');
disp(['Coef:', num2str(supplybetas')]);
disp(['SE:   ', num2str(supplyse')]);
disp('-------------------');

clear alpha covars covars2 DELTA isConstrained insts instGroupIdMat logWithinGroupShare groupId
clear groupIndVec markup normDept normInsts normLogWithinGroupShare normPrice optw theta2sls thetaNested
clear thetaNested2 z Z x y weightMat wInGroupShareInsts thetaStart thetaStart2
clear thetaStart3 v1 priceInsts options mc groupIdMat alpha sigma temp wInGroupShareInsts2 se
clear covars3 deptInd isThisFirm isThisGroup isThisObs isThisYear isValidObs iYear j
clear NORM_EXTENT NORM_NESTED NORM_STYLE outsideShare priceInd sigmaInd T yearList nFirms i iFirm firmList iGroup
clear iRow thetaNested3


%% Clear off Unneeded Values
diary off
















