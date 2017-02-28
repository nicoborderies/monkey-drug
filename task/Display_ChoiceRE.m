function [] = Display_ChoiceRE(WindowPtr,rewardImg,xcenter,ycenter, level, rewardFactor , effortFactor , effortPrecision )

% Display_ChoiceRE draws on the screen the option of the choice (reward
% quantity, effort level) & visual feedback on exerted force
% 
% written by Nicolas Borderies - 02/2015.

% linear scale
    calibFB = 100;
    visualFB = level/calibFB;
    threshold = effortFactor;

% exponential scale
%     convexity = 2;
%     visualFB =  (exp(level*0.01*convexity)-1)/(exp(convexity)-1);
%     threshold = (exp(effortFactor*convexity)-1)/(exp(convexity)-1);


width = 180;
heigth = 700;
lineWidth = 10;
lowerLevel = 450;

% draw the effort rectangle
baseRect = [0, 0, width, threshold*heigth];
centeredRect = CenterRectOnPointd(baseRect, xcenter, ycenter+lowerLevel-(threshold*heigth)/2);
rectColor = [150 150 150];
Screen('FillRect', WindowPtr, rectColor, centeredRect);

% draw the orange & blue cursor
col = [255 153 0];
col2 = [0 0 255];

    Screen('FillRect',WindowPtr,col,[(xcenter-width/2) (ycenter+lowerLevel-lineWidth-visualFB*heigth) (xcenter+width/2) (ycenter+lowerLevel)]);
if (visualFB >= threshold)
    Screen('FillRect',WindowPtr,col2,[(xcenter-width/2) (ycenter+lowerLevel-lineWidth-visualFB*heigth) (xcenter+width/2) (ycenter+lowerLevel-threshold*heigth)]);
end


% draw the contours
Screen('DrawLine',WindowPtr,[0 0 0],(xcenter-width/2), (ycenter+lowerLevel-threshold*heigth), (xcenter+width/2), (ycenter+lowerLevel-threshold*heigth),lineWidth);
Screen('DrawLine',WindowPtr,[0 0 0],(xcenter-width/2), (ycenter+lowerLevel), (xcenter+width/2), (ycenter+lowerLevel),lineWidth);
Screen('DrawLine',WindowPtr,[0 0 0],(xcenter-width/2), (ycenter+lowerLevel-threshold*heigth-lineWidth/2), (xcenter-width/2), (ycenter+lowerLevel+lineWidth/2),lineWidth);
Screen('DrawLine',WindowPtr,[0 0 0],(xcenter+width/2), (ycenter+lowerLevel-threshold*heigth-lineWidth/2), (xcenter+width/2), (ycenter+lowerLevel+lineWidth/2),lineWidth);



end

