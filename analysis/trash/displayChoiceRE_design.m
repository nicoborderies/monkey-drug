%% display_choiceRE_design


reward = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
effort = [repmat([1 2 3 4],1,4)];

f = figure;

s = scatter(reward,effort,'filled');

colMatrix = nan(16,3);
for i =1:16
    colMatrix(i,:) = ...
        [0+effort(i)/4 , 0+reward(i)/4 ,0];
end

s.CData = colMatrix;
s.SizeData = 250;

ax = gca;
ax.XLim = [0.5 5]; ax.XTick = [ 1 2 3 4]+0.5;
ax.YLim = [0.5 5]; ax.YTick = [ 1 2 3 4]+0.5;
ax.TickDir = 'out';
ax.XTickLabel = {'','','',''};
ax.YTickLabel = {'','','',''};

% ax.TickLength = [0 0];