%Estimate the WithinGroupShare using Estimated Shares, MKT_SIZE and Nests
%Estimate the WithinGroupShare using Estimated Shares, MKT_SIZE and Nests
function withinGroupShare = calcwithinshare(share, CarrID, NUM_GROUPS, groupId, MKT_SIZE)
    % Calculate within group share 
groupQ = zeros(size(CarrID));
for iGroup = 1:NUM_GROUPS
    groupIndVec = find(ismember(groupId,iGroup)); % identify group members
    groupQ(groupIndVec) = sum(share(groupIndVec).*MKT_SIZE(groupIndVec)); % calculate group quantity and record value for each group member
end
withinGroupShare =  (share.*MKT_SIZE)./groupQ;
end