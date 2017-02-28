function [] = Display_spot(WindowPtr,Img,xcenter,ycenter)
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