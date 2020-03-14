function plotMultipleBarsWithError(ax, allMean, allStd)
bar(allMean, 'grouped');
hold on
ngroups = size(allMean, 1);
nbars = size(allMean, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, allMean(:,i), allStd(:,i), 'k', 'linestyle', 'none');
end