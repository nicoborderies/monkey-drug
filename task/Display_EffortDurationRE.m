function [] = Display_EffortDurationRE(WindowPtr,xcenter,ycenter, level, effortFactor, forceDuration, durationUnit  )

% Display_ChoiceRE draws on the screen the option of the choice (reward
% quantity, effort level) & visual feedback on exerted force
% 
% written by Nicolas Borderies - 02/2015.

visualFB = (forceDuration/durationUnit);

spotSize = 75;
width = spotSize*visualFB;
heigth = spotSize*visualFB;
lineWidth = 5;

% draw 
col = [0 0 255];
Screen('FillRect',WindowPtr,col,[(xcenter-width/2) (ycenter-heigth/2) (xcenter+width/2) (ycenter+heigth/2)]);

% draw the contours
Screen('DrawLine',WindowPtr,[0 0 0],(xcenter-spotSize/2), (ycenter+spotSize/2), (xcenter+spotSize/2), (ycenter+spotSize/2),lineWidth);
Screen('DrawLine',WindowPtr,[0 0 0],(xcenter-spotSize/2), (ycenter-spotSize/2), (xcenter+spotSize/2), (ycenter-spotSize/2),lineWidth);
Screen('DrawLine',WindowPtr,[0 0 0],(xcenter-spotSize/2), (ycenter-spotSize/2), (xcenter-spotSize/2), (ycenter+spotSize/2),lineWidth);
Screen('DrawLine',WindowPtr,[0 0 0],(xcenter+spotSize/2), (ycenter-spotSize/2), (xcenter+spotSize/2), (ycenter+spotSize/2),lineWidth);


end

