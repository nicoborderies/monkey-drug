%% init_computer

% screen
    % dimensions
        L=1366;
        H=768;
        
        % background position
        xSpot=L/2;ySpot=H/2;
        x=L/2; 
        y=H/2;
        
        % effort visual feedback
        xEffortScale(1)=(0.35)*L; % left option
        xEffortScale(2)=(0.85)*L; % right option
        xReward = xEffortScale - (0.20)*L; 
        yEffortScale=y;
        
    % display
        Screen('Preference', 'SkipSyncTests', 1);
        [window]=Screen('OpenWindow',0,[0 0 0],[]);
        Screen('TextSize', window, 40);
        Screen('TextFont', window, 'arial');
        HideCursor;

% keyboard
    KbName('UnifyKeyNames');
    key.left=37;
    key.right=39;
    key.space=32;
    key.escape=27;