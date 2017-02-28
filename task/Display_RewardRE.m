function [] = Display_RewardRE(WindowPtr,rewardImg,xcenter,ycenter, rewardFactor )
% Display_RewardRE draws on the screen the option of the choice (reward
% quantity, effort level) & visual feedback on exerted force
% 
% written by Nicolas Borderies - 02/2015.

heigth = 165;
lowerLevel = 410;

% draw the reward level
for iR = 1:rewardFactor
    rectIcentive=CenterRectOnPoint(Screen('Rect',rewardImg),xcenter,ycenter+lowerLevel-(iR-1)*heigth);
    Screen('DrawTexture',WindowPtr,rewardImg,[],rectIcentive,[],[],0);
end



end