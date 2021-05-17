%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Methods
% Project - Main File 2
% Code Outline: The code consists of Computing the Model Counterfactuals
% Run main.m first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
disp(' '); % to separate results from previous executions
diary Work.txt
format short g
rng(1);
% Load and separate data
load("Premerger_db1b")
load("Premerger_MkDist")
Premerger_db1b = sortrows(Premerger_db1b,30);
Premerger_MkDist = sortrows(Premerger_MkDist,1);
Premerger_db1b.Distsq = Premerger_db1b.Dist.^2;
NUM_GROUPS = size(unique(Premerger_db1b.Nest),1); % controls the number of groups used in nested logit
MKT_SIZE= Premerger_db1b.AvgPop;  % set market size to the geometric mean of  endpoint cities' MSA Populations

Premerger_db1b.normPrice = Premerger_db1b.BinFare/mean(Premerger_db1b.BinFare);
Premerger_db1b.normDept = Premerger_db1b.dept/mean(Premerger_db1b.dept);
Premerger_db1b.normDist = Premerger_db1b.Dist/mean(Premerger_db1b.Dist);
Premerger_db1b.normDistSq = Premerger_db1b.Distsq/mean(Premerger_db1b.Distsq);


% Store parameter values
nObs=size(Premerger_db1b.Year,1);
T=size(unique(Premerger_db1b.Year),1);

% Calculate market shares
share =Premerger_db1b.ProdSh;

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

demandbetas = load("demandbetas.txt");
demandse = load("demandse.txt");
supplybetas = load("supplybetas.txt");
supplyse = load("supplyse.txt");


%% Recalculate Substitution Matrix and Marginal Cost
currSubmat = subMat(demandbetas(16), demandbetas(15), share,withinGroupShare, Premerger_db1b.Year, Premerger_db1b.CarrID, Premerger_db1b.Nest,nObs);
Currmarkup = -(currSubmat\cast(share,'single'));
clear currSubmat
Currmc = Premerger_db1b.BinFare - Currmarkup;
Currmc = cast(Currmc,'double');
Premerger_db1b.mc = Currmc;

%% Average Price by Market By Firm and Firm Profitability
uniqueCarr = unique(Premerger_db1b.TICKET_CARRIER);
uniqueYears = unique(Premerger_db1b.Year);
nYears = size(uniqueYears,1);
nFirms = size(uniqueCarr,1);
year = [0];
firm = ["AA"];
profit = [0];
price = [0];
revenue = [0];
%aggregate the quantity for each time the year appears and record the
%market size
for Count = 1:nYears
    for count2 = 1:nFirms
        tempYear = uniqueYears(Count);
        tempFirm = uniqueCarr(count2);
         check = (tempYear == Premerger_db1b.Year & strcmp(tempFirm ,Premerger_db1b.TICKET_CARRIER));
        tempPrice = Premerger_db1b.BinFare(check);
        tempshr = Premerger_db1b.ProdSh(check);
        tempmktsize = MKT_SIZE(check);
        tempmc = Currmc(check);
        tempfirmsales = sum(Premerger_db1b.total_passengers_by_product(check));
        temprodsales = Premerger_db1b.total_passengers_by_product(check);
        tempwithinfirmshr = temprodsales./tempfirmsales;
        year = [year;tempYear];
        firm = [firm; tempFirm];
        profit = [profit ;calcprofit(tempPrice,tempshr,tempmktsize,tempmc)];
       price = [price ;AvgPrice(tempPrice, tempwithinfirmshr)];
       revenue = [revenue; sum(tempPrice.*tempshr.*tempmktsize)];
    end
end
avgPriceAndProfitByFirmYear = table(year, firm, profit, price, revenue);
avgPriceAndProfitByFirmYear = avgPriceAndProfitByFirmYear(2:end,:);
clear uniqueCarr uniqueYears nYears nFirms start tempYear tempFirm tempPrice tempshr tempmc tempfirmsales
clear tempprodsales tempwithinfirmshr Count count2 check year firm profit price

%% Counterfactual: AA and US Merger: Change all US Products to AA Products
counter_premerger_db1b = Premerger_db1b;
counter_premerger_mktdist = Premerger_MkDist;
counter_premerger_db1b.CarrID(strcmp("US" ,Premerger_db1b.TICKET_CARRIER)) = 1;
counter_premerger_db1b.US_Airways(strcmp("US" ,Premerger_db1b.TICKET_CARRIER)) = 0;
counter_premerger_db1b.American(strcmp("US" ,Premerger_db1b.TICKET_CARRIER)) = 1;
counter_premerger_db1b.TICKET_CARRIER(strcmp("US" ,Premerger_db1b.TICKET_CARRIER)) = {'AA'};

%Give US Airway Hub Airports to American Airlines
USHub = unique(Premerger_db1b.MkID(find(strcmp("US" ,Premerger_db1b.TICKET_CARRIER)&Premerger_db1b.HubMC == 1)));
for i = 1:size(USHub,1)
    tempHub = USHub(i);
    counter_premerger_db1b.HubMC(strcmp("AA" ,counter_premerger_db1b.TICKET_CARRIER)&counter_premerger_db1b.MkID == tempHub) = 1;
end

%Remove Nonunique products by USAirways and AA and add departures
AAMkt = unique(counter_premerger_db1b.MkID(counter_premerger_db1b.CarrID == 1));
AABin = unique(counter_premerger_db1b.BinFare(counter_premerger_db1b.CarrID == 1));
for i = 1:size(AAMkt,1)
    tempmkt = AAMkt(i);
    for j = 1:size(AABin,1)
        tempbin = AABin(j);
        check = (counter_premerger_db1b.CarrID == 1)&(counter_premerger_db1b.MkID == tempmkt)&(counter_premerger_db1b.BinFare == tempbin);
        counter_premerger_db1b.dept(check) = sum(counter_premerger_db1b.dept(check));
        counter_premerger_db1b.total_passengers_by_product(check) = sum(counter_premerger_db1b.total_passengers_by_product(check));
    end
end
counter_premerger_db1b.ProdSh = counter_premerger_db1b.total_passengers_by_product./counter_premerger_db1b.AvgPop;
c = [counter_premerger_db1b.CarrID counter_premerger_db1b.MkID counter_premerger_db1b.BinFare];
[C,ia] = unique(c(:,1:3),'rows','first');
counter_premerger_mktdist = counter_premerger_mktdist(ia,:);
counter_premerger_db1b = counter_premerger_db1b(ia,:);
counter_withinGroupShare = calcwithinshare(counter_premerger_db1b.ProdSh, counter_premerger_db1b.CarrID, NUM_GROUPS, counter_premerger_db1b.Nest, counter_premerger_db1b.AvgPop);


%% Use Pricing FOCs to reestimate pricing equilibrium and shares
Currmc2 = Currmc(ia);
tempwithinShare = counter_withinGroupShare;
tempmarkup = Currmarkup(ia);
clear Currmarkup Currmc 
clear c C ia AAMkt AABin
%%
nObs2=size(counter_premerger_db1b.Year,1);
for l = 1:2
      if(l == 1)
           tempprice = counter_premerger_db1b.BinFare;
           tempmc = Currmc2;
           normtempprice = tempprice/mean(tempprice);
      else
           tempprice = counter_premerger_db1b.BinFare;
           normtempprice = tempprice/mean(tempprice);
           tempmc = tempprice - tempmarkup;
      end
      %Using Prices reestimate shares
      X = [ones(nObs2,1) counter_premerger_db1b.NumDest counter_premerger_db1b.normDist counter_premerger_db1b.normDistSq counter_premerger_db1b.FLLAS counter_premerger_db1b.Slot counter_premerger_db1b.OtherCarr counter_premerger_db1b.American counter_premerger_db1b.Continental ...
    counter_premerger_db1b.Delta counter_premerger_db1b.United counter_premerger_db1b.US_Airways counter_premerger_db1b.JetBlue counter_premerger_db1b.Southwest log(tempwithinShare) normtempprice counter_premerger_db1b.normDept]; % covariates
        counter_premerger_db1b.ProdSh = calcshare(demandbetas, X);
        tempshare = counter_premerger_db1b.ProdSh;
        %Using shares, reestimate within Group Shares
        tempwithinShare = calcwithinshare(tempshare, counter_premerger_db1b.CarrID, NUM_GROUPS, counter_premerger_db1b.Nest, counter_premerger_db1b.AvgPop);
        %Using shares and within Group shares, reestimate ownership matrix
        tempsub = subMat(demandbetas(16), demandbetas(15), tempshare,tempwithinShare, counter_premerger_db1b.Year, counter_premerger_db1b.CarrID, counter_premerger_db1b.Nest,nObs2);
        %Using reestimated shares, reestimated ownership  matrix,
        %within group shares and pricing FOCs, reestimate markup and prices
        tempmarkup = -(tempsub\cast(tempshare,'single'));
        clear tempsub
        tempmarkup = cast(tempmarkup,'double');
        counter_premerger_db1b.BinFare = solvePrice(tempmarkup,tempmc);
        abs(tempprice - counter_premerger_db1b.BinFare)
        if(abs(tempprice - counter_premerger_db1b.BinFare) < 10e-6)
%             Currmc2 = counter_premerger_db1b.BinFare - tempmarkup;
           break;
%         elseif(sum(isnan(abs(tempprice - counter_premerger_db1b.BinFare))) > 0)
%             counter_premerger_db1b.BinFare = tempprice;
%             Currmc2 = tempmc;
%             break;
         end      
end
Currmc2 = tempmc;
counter_premerger_db1b.mc = Currmc2;
counter_withinGroupShare = tempwithinShare;
clear tempprice tempwithinShare tempshare X tempmc tempmarkup
%% Average Price by Market By Firm and Firm Profitability
%Set the minimum price of flights to be 10.
counter_premerger_mktdist = counter_premerger_mktdist(counter_premerger_db1b.BinFare >= 10,:);
counter_premerger_db1b = counter_premerger_db1b(counter_premerger_db1b.BinFare >= 10,:);
counter_withinGroupShare = calcwithinshare(counter_premerger_db1b.ProdSh, counter_premerger_db1b.CarrID, NUM_GROUPS, counter_premerger_db1b.Nest, counter_premerger_db1b.AvgPop);

uniqueCarr = unique(counter_premerger_db1b.TICKET_CARRIER);
uniqueYears = unique(counter_premerger_db1b.Year);
nYears = size(uniqueYears,1);
nFirms = size(uniqueCarr,1);
year = [0];
firm = ["AA"];
profit = [0];
price = [0];
revenue = [0];
%aggregate the profit and note the average price for each firm and year in
%the counterfactual
for Count = 1:nYears
    for count2 = 1:nFirms
        tempYear = uniqueYears(Count);
        tempFirm = uniqueCarr(count2);
        check = (tempYear == counter_premerger_db1b.Year & strcmp(tempFirm ,counter_premerger_db1b.TICKET_CARRIER));
        tempPrice = counter_premerger_db1b.BinFare(check);
        tempshr = counter_premerger_db1b.ProdSh(check);
        tempmc = Currmc2(check);
        tempfirmsales = sum(counter_premerger_db1b.total_passengers_by_product(check));
        temprodsales = counter_premerger_db1b.total_passengers_by_product(check);
        tempwithinfirmshr = temprodsales./tempfirmsales;
        tempmktsize = counter_premerger_db1b.AvgPop(check);
        year = [year;tempYear];
        firm = [firm; tempFirm];
        profit = [profit ;calcprofit(tempPrice,tempshr,tempmktsize,tempmc)];
        price = [price ;AvgPrice(tempPrice, tempwithinfirmshr)];
        revenue = [revenue; sum(tempPrice.*tempshr.*tempmktsize)];
    end
end
counter_avgPriceAndProfitByFirmYear = table(year, firm, profit, price, revenue);
counter_avgPriceAndProfitByFirmYear = counter_avgPriceAndProfitByFirmYear(2:end,:);

clear uniqueCarr uniqueYears nYears nFirms start tempYear tempFirm tempPrice tempshr tempmc tempfirmsales
clear Count count2 tempprodsales tempwithinfirmshr check year firm profit price

%% Compare Average Change in Prices to Naive Diff-in-Diff Result
counter_avgPriceAndProfitByFirmYear2 = counter_avgPriceAndProfitByFirmYear(counter_avgPriceAndProfitByFirmYear.price ~= 0 & counter_avgPriceAndProfitByFirmYear.profit ~= 0,:);
avgPriceAndProfitByFirmYear2 = avgPriceAndProfitByFirmYear(avgPriceAndProfitByFirmYear.price ~= 0 & avgPriceAndProfitByFirmYear.profit ~= 0,:);

year = [0];
firm = ["AA"];
profitchange = [0];
pricechange = [0];
revchange = [0];

avgChange = zeros(size(counter_avgPriceAndProfitByFirmYear2,1),3);
for k = 1:size(counter_avgPriceAndProfitByFirmYear2,1)
    if strcmp("AA",counter_avgPriceAndProfitByFirmYear2.firm(k))
        firm = [firm; "AA/US"];
        currYear = counter_avgPriceAndProfitByFirmYear2.year(k);
        year = [year;currYear];
        pricechange = [pricechange; mean(counter_avgPriceAndProfitByFirmYear2.price(k) - avgPriceAndProfitByFirmYear2.price(find(currYear == avgPriceAndProfitByFirmYear2.year & (strcmp("AA",avgPriceAndProfitByFirmYear2.firm) | strcmp("US",avgPriceAndProfitByFirmYear2.firm)))))];
        profitchange = [profitchange; counter_avgPriceAndProfitByFirmYear2.profit(k) - sum(avgPriceAndProfitByFirmYear2.profit(find(currYear == avgPriceAndProfitByFirmYear2.year & (strcmp("AA",avgPriceAndProfitByFirmYear2.firm) | strcmp("US",avgPriceAndProfitByFirmYear2.firm)))))];
        revchange = [revchange;counter_avgPriceAndProfitByFirmYear2.revenue(k) - sum(avgPriceAndProfitByFirmYear2.revenue(find(currYear == avgPriceAndProfitByFirmYear2.year & (strcmp("AA",avgPriceAndProfitByFirmYear2.firm) | strcmp("US",avgPriceAndProfitByFirmYear2.firm)))))];
    else
        currYear = counter_avgPriceAndProfitByFirmYear2.year(k);
        year = [year;currYear];
        currFirm = counter_avgPriceAndProfitByFirmYear2.firm(k);
        firm = [firm; currFirm];
        prevMeanPrice = avgPriceAndProfitByFirmYear2.price(find(currYear == avgPriceAndProfitByFirmYear2.year & strcmp(currFirm,avgPriceAndProfitByFirmYear2.firm)));
        prevMeanProfit = avgPriceAndProfitByFirmYear2.profit(find(currYear == avgPriceAndProfitByFirmYear2.year & strcmp(currFirm,avgPriceAndProfitByFirmYear2.firm)));
        prevmeanRev = avgPriceAndProfitByFirmYear2.revenue(find(currYear == avgPriceAndProfitByFirmYear2.year & strcmp(currFirm,avgPriceAndProfitByFirmYear2.firm)));
        pricechange = [pricechange;counter_avgPriceAndProfitByFirmYear2.price(k) - prevMeanPrice];
        profitchange = [profitchange;counter_avgPriceAndProfitByFirmYear2.profit(k) - prevMeanProfit];
        revchange = [revchange; counter_avgPriceAndProfitByFirmYear2.revenue(k) - prevmeanRev];
    end
end
avgChange = table(year,firm, pricechange,profitchange,revchange);
avgChange = avgChange(2:end,:);
initAvgPrice = WeightedAvgPrice(Premerger_db1b.BinFare,Premerger_db1b.ProdSh);
mergerAvgPrice = WeightedAvgPrice(counter_premerger_db1b.BinFare,counter_premerger_db1b.ProdSh);
avgProfitChange = mean(avgChange.profitchange);
avgRevChange = mean(avgChange.revchange);
PriceElas = sum(sum(subMat(demandbetas(16),demandbetas(15),Premerger_db1b.ProdSh,withinGroupShare,Premerger_db1b.Year,Premerger_db1b.CarrID,Premerger_db1b.Nest,nObs)));
nObs2=size(counter_premerger_db1b.Year,1);
counter_PriceElas = sum(sum(subMat(demandbetas(16),demandbetas(15),counter_premerger_db1b.ProdSh,counter_withinGroupShare,counter_premerger_db1b.Year,counter_premerger_db1b.CarrID,counter_premerger_db1b.Nest,nObs2)));
marginalCost = WeightedAvgPrice(Premerger_db1b.mc,Premerger_db1b.ProdSh);
counter_marginalCost = WeightedAvgPrice(counter_premerger_db1b.mc,counter_premerger_db1b.ProdSh);
lernerIndex = WeightedAvgPrice((Premerger_db1b.BinFare - Premerger_db1b.mc)./Premerger_db1b.BinFare,Premerger_db1b.ProdSh);
counter_lernerIndex = WeightedAvgPrice((counter_premerger_db1b.BinFare - counter_premerger_db1b.mc)./counter_premerger_db1b.BinFare,counter_premerger_db1b.ProdSh);


%% Save Values
writetable(avgChange, "avgChange.txt")
writetable(avgPriceAndProfitByFirmYear2, "avgPriceAndProfitByFirmYear2.txt")
writematrix(avgProfitChange, "avgProfitChange.txt")
writematrix(avgRevChange, "avgRevChange.txt")
writetable(counter_avgPriceAndProfitByFirmYear2, "counter_avgPriceAndProfitByFirmYear2.txt")
writematrix(counter_lernerIndex, "counter_lernerIndex.txt")
writematrix(counter_marginalCost, "counter_marginalCost.txt")
writematrix(counter_PriceElas, "counter_PriceElas.txt")
writematrix(counter_withinGroupShare, "counter_withinGroupShare.txt")
writematrix(withinGroupShare, "withinGroupShare.txt")
writematrix(initAvgPrice, "initAvgPrice.txt")
writematrix(lernerIndex, "lernerIndex.txt")
writematrix(marginalCost, "marginalCost.txt")
writematrix(mergerAvgPrice, "mergerAvgPrice.txt")
writematrix(nObs, "nObs.txt")
writematrix(nObs2, "nObs2.txt")
writematrix(PriceElas, "PriceElas.txt")
%% Clear off Unneeded Values
clear k l MKT_SIZE NUM_GROUPS revchange revenue temprodsales FUN_TOL
clear Currmc2 prevmeanRev share T iGroup groupId USHub groupIndVec tempmktsize i
clear currYear currFirm prevMeanPrice prevMeanProfit year firm pricechange profitchange
clear firm year profitchange pricechange counter_avgPriceAndProfitByFirmYear avgPriceAndProfitByFirmYear
diary off




