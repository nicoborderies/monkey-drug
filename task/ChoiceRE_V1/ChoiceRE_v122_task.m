function [data,gripdata,gripbaseline] = ChoiceRE_v122_task(root,calib,reward_unit_duration,effort_unit_duration)
%% REWARD-EFFORT CHOICE TASK: v1.1 (training procedure)
% Experimental design selection: Nicolas Borderies, Sébastien Bouret, Mathias Pessiglione, Caroline Jahn, Chiara Varazzani, Simon Garret
% Data acquisition: Nicolas Borderies, Sophie Gilardeau
% Subjects (macaques rhesus): Esmeralda(f) , Bob(m)
% Task development: Nicolas Borderies
% Data analysis: Nicolas Borderies

% TASK STRUCTURE:


% PSEUDO-CODE:
%
% INITIALIZE Configurations 
% INITIALIZE Task Parametrization
% WHILE trial < ntrial
%   WAIT inter-trial interval
%   PROCESS adaptive design
%   DISPLAY preparatory signal
%   DISPLAY reward-effort option(s)
%   RECORD choice & effort
%   SEND reward
%   IF error
%       repeat trial
%   END
% END
% SAVE data
% PROCESS offline analysis



% clear all;
% close all;



%% Software Configurations     
% ------- Directory Configuration ------- %
    taskdir = [root '\task'];
    picdir= [taskdir '\display'];
    
% ------- Grip Configuration ------- %
% setting manually transducter parameters (offset~100 with a constant gain during the session)
% MatOs;
MatOs_multichannel;clear Handle;

% iGrip = 1; % left=1 / right=2 ;
    Handle=InitializeGrip('SerialMBB',[1 2]); % 2 grips


% ------- Screen Configuration ------- %
cd(picdir);

L=1366;
H=768;
xSpot=L/2;ySpot=H/2;
x=L/2; % background position
y=H/2;

% force visual fb
xEffortScale(1)=(0.35)*L; % left option
xEffortScale(2)=(0.85)*L; % right option
xReward = xEffortScale - (0.20)*L; 
yEffortScale=y;
% wh_px    = 400;   % window's heigh
% ww_px    = 400;   % window's width

Screen('Preference', 'SkipSyncTests', 1);

[window]=Screen('OpenWindow',0,[0 0 0],[]);

Screen('TextSize', window, 40);
Screen('TextFont', window, 'arial');
HideCursor;



% ------- Keyboard Configuration ------- %
KbName('UnifyKeyNames');
key.left=37;
key.right=39;
key.space=32;
key.escape=27;


% ------- Parallel Port Configuration ------- %
config_io;
outp(888,0); % put all pins to zero


% ------- Generator reset ------- %
rand('state',sum(100*clock));


% ------- Loads images and creates positions ------- %
cd(picdir);

background=Screen('MakeTexture',window,imread('bg_100.bmp'));
background_omission=Screen('MakeTexture',window,imread('bg_100.bmp'));
background_anticipation=Screen('MakeTexture',window,imread('bg_lateralization.bmp'));
background_lateralizationError=Screen('MakeTexture',window,imread('bg_lateralization.bmp'));


pic_incentive = Screen('MakeTexture',window,imread(['reward_drop2.bmp']));

pic_spot_red=Screen('MakeTexture',window,imread('spot_1_red.bmp'));
pic_spot_green=Screen('MakeTexture',window,imread('spot_1_green.bmp'));
pic_spot_blue=Screen('MakeTexture',window,imread('spot_1_blue.bmp'));

rect_spot_red=CenterRectOnPoint(Screen('Rect',pic_spot_red),xSpot,ySpot);
rect_spot_green=CenterRectOnPoint(Screen('Rect',pic_spot_green),xSpot,ySpot);
rect_spot_blue=CenterRectOnPoint(Screen('Rect',pic_spot_blue),xSpot,ySpot);
rect_bg=CenterRectOnPoint(Screen('Rect',background),x,y);


%% Task parametrization
cd(root);

% ------- Experimental Factors ------- %

rewardLvls = [0.6, 1, 2, 4]; % incitation : cardinal levels
effortLvls = [0.10 , 0.40 , 0.70 , 1.00 ];
laterality = [1 2]; % laterality: left/right


% ------- Trials ------- %

nTrialByBlock=8;
nBlock=200;

% temporal design: blocks 
    trialNumber = nTrialByBlock*nBlock;
    trials = [1:trialNumber];


% current specifications
    % only extreme effort options offered
%         effortBlock1 = [ repmat(1,1,nTrialByBlock/4) , repmat(2,1,nTrialByBlock/4) , repmat(3,1,nTrialByBlock/4), repmat(4,1,nTrialByBlock/4)]; 
%         effortBlock2 = [ repmat([repmat(1,1,nTrialByBlock/16) , repmat(2,1,nTrialByBlock/16) , repmat(3,1,nTrialByBlock/16), repmat(4,1,nTrialByBlock/16)],1,4)]; 
        effortBlock1 = [ 1,1,2,2,3,4,3,4]; 
        effortBlock2 = [ 4,3,4,3,2,2,1,1]; 
        twoLevelBlock = [2,2,2,2,3,3,3,3];

    % only currentReward offered
        currentReward = 3; currentEffort = 3; % it could be 2 or 3 or 4
        rewardBlock1 = repmat(currentReward,1,nTrialByBlock);
        rewardBlock2 = rewardBlock1;
        constantBlock1 = repmat(currentReward,1,nTrialByBlock);
        constantBlock1 = twoLevelBlock;
        constantBlock2 = repmat(currentEffort,1,nTrialByBlock);
        constantBlock2 = twoLevelBlock;


%         rewardBlock1 = [ repmat([repmat(1,1,nTrialByBlock/64) , repmat(2,1,nTrialByBlock/64) , repmat(3,1,nTrialByBlock/64), repmat(4,1,nTrialByBlock/64)],1,16)];
%         rewardBlock2 = [ repmat([1,2,3,4],1,64)];  % similar rewards for both options
        sideBlock = [repmat(1,1,nTrialByBlock/2) , repmat(2,1,nTrialByBlock/2)];

        
% alternative specifications

% pseudo-random permutation across options sequence
    e1 = []; r1 = [];e2 = []; r2 = []; side = [];
    for iBlock = 1:nBlock
        index = randperm(nTrialByBlock); 
        switch mod(iBlock,2)  % alternace between blocks of varying effort & blocks of varying rewards
            case 1
            e1 = [e1, effortBlock1(index) ];
            r1 = [r1, constantBlock1(index)];
            e2 = [e2, effortBlock2(index) ];
            r2 = [r2, constantBlock1(index)];
            case 0
            e1 = [e1, constantBlock2(index) ];
            r1 = [r1, effortBlock1(index)];
            e2 = [e2, constantBlock2(index) ];
            r2 = [r2, effortBlock2(index)];
        end
            side = [side, sideBlock(index)];
    end


% --- Performance Criterions --- %
omissionCriterion = 30; % maximal n° of successive omissions authorized
withdrawalCriterion = 0; % maximal n° of successive failed trials which are repeated
anticipationCriterion = 20; % (%) of force level exerted before the time of choice
reactionTimeCriterion = 20; % (%) of force level exerted after starting onsert
effortPrecision = 0.20;
effortReleaseCriterion = 0;
calib_increment = 10;

% ------- Temporal parameters (sec.)------- %

starting_duration = 3; % time before the experiment start running
% spotduration=1;
min_redSpotDuration=0.5;
max_redSpotDuration=2.5; 
% min_cueduration=0.5;
% max_cueduration=1.5;  
anticipationduration=0.1;
responseduration=2.5;
% intertrialduration=3;
min_intertrialduration=0.5;
max_intertrialduration=1.5;
failure_delay = 4;


redSpotDuration=rand(1,10*numel(trials))*(max_redSpotDuration-min_redSpotDuration)+min_redSpotDuration;
intertrialduration=rand(1,10*numel(trials))*(max_intertrialduration-min_intertrialduration)+min_intertrialduration;


% ------- Pre-Allocating variables to save ------- %

% rewards = [r1(1) ; r2(1)]; % initial trial's factors
% efforts = [e1(1) ; e2(1)];

responsetime=nan(1,10*numel(trials));
reactiontime=nan(1,10*numel(trials));
feedbacktime=nan(1,10*numel(trials));
leveltime   =nan(1,10*numel(trials));

force=nan(1,10*numel(trials));
sumforce=nan(1,10*numel(trials));
forceduration=zeros(2,10*numel(trials));
baseline=nan(1,10*numel(trials));
perf=nan(1,10*numel(trials));
sumperf=nan(1,10*numel(trials));
gain=zeros(1,10*numel(trials));
choice=zeros(1,10*numel(trials));
% choice=side;
error=nan(1,10*numel(trials));
repetition=zeros(1,10*numel(trials));

rewards=nan(2,10*numel(trials));
efforts=nan(2,10*numel(trials));


gripdata=struct;
gripbaseline=struct;
gripanticipation=struct;

for s=1:2
    calibration(s,:) = repmat(calib(s),1,numel(trials));
end


% ------- Initializating variables ------- %

total=0; %total gain counter
task_exit=0;  %exit flag to end task
err=0;   %error counter
omission = 0; % omission counter
withdrawal = 0; % withdrawal counter
ntrial = 0; %starting trial
iTrial = 0;
iBlock = 1;


%% TASK STRUCTURE %%

cd(taskdir);

% ready to start the task 
Screen('TextSize', window, 20);
[width,hight]=RectSize(Screen('TextBounds',window,'appuyer sur une touche pour commencer'));
Screen('DrawText',window,'appuyer sur une touche pour commencer',x-width/2,y+300,[100 100 100]);
Screen(window,'Flip');
KbWait;
Screen(window,'Flip');
WaitSecs(starting_duration);
   

while task_exit==0    
    ntrial = ntrial +1; % trial incrementation
    iTrial = iTrial +1; % trial incrementation

    % Adaptive Control of Design
    if  isnan(rewards(:,ntrial))
        % display intructions
        Screen(window,'Flip'); 
        Screen('TextSize', window, 20);
        [width,hight]=RectSize(Screen('TextBounds',window,'Sélectionner les dimensions qui varient'));
        Screen('DrawText',window,'Sélectionner les dimensions qui varient',x-width/2,y,[100 100 100]);
        [width,hight]=RectSize(Screen('TextBounds',window,'recompenses'));
        Screen('DrawText',window,'recompenses',x-width/2-x/2,y+300,[100 100 100]);
        [width,hight]=RectSize(Screen('TextBounds',window,'efforts'));
        Screen('DrawText',window,'efforts',x-width/2+x/2,y+300,[100 100 100]);
        [width,hight]=RectSize(Screen('TextBounds',window,'mixed'));
        Screen('DrawText',window,'mixed',x-width/2,y+400,[100 100 100]);
        Screen(window,'Flip');
        
        % monitor response
        [~, keyCode] = KbWait;
        Block1 = [ 1,1,2,2,3,4,3,4]; 
        Block2 = [ 4,3,4,3,2,2,1,1];
        Block3 = ones(1,nTrialByBlock);
        index = randperm(nTrialByBlock); 
        Block1 = Block1(index); Block2 = Block2(index);
        if  keyCode(key.left)==1
            currentEffort = 3;
            twoLevelBlock = [2,2,2,2,3,3,3,3];
            rewards(:,ntrial:ntrial+nTrialByBlock-1) = [Block1 ; Block2]; 
            efforts(:,ntrial:ntrial+nTrialByBlock-1) = [ twoLevelBlock ; twoLevelBlock ];
        elseif  keyCode(key.right)==1
            currentReward = 1;
            twoLevelBlock = [2,2,2,2,3,3,3,3];
            rewards(:,ntrial:ntrial+nTrialByBlock-1) = [ twoLevelBlock ; twoLevelBlock ]; 
            efforts(:,ntrial:ntrial+nTrialByBlock-1) = [Block1; Block2 ];   
        elseif keyCode(key.space)==1 
            Block1 = [ 1,2,3,4,2,2,3,3];
            Block2 = [ 4,3,2,1,2,2,3,3];
            Block3 = fliplr(Block1);
            Block4 = fliplr(Block2);
            rewards(:,ntrial:ntrial+nTrialByBlock-1) = [ Block1(index) ; Block2(index) ]; 
            efforts(:,ntrial:ntrial+nTrialByBlock-1) = [ Block3(index) ; Block4(index) ]; 
        else
            task_exit = 1;
        end
    end
    trials(iTrial) = ntrial;

        
    
    
    % Display Background
    if ntrial==1 
        Screen('DrawTexture',window,background,[],rect_bg);
    else
        switch error(ntrial-1)
            case 0
                Screen('DrawTexture',window,background,[],rect_bg);
            case 1
                Screen('DrawTexture',window,background_omission,[],rect_bg);
            case 2
                Screen('DrawTexture',window,background_anticipation,[],rect_bg);
            case 3
                Screen('DrawTexture',window,background_lateralizationError,[],rect_bg);
            otherwise
                Screen('DrawTexture',window,background,[],rect_bg);
        end
    end
    backgroundtime=Screen(window,'Flip');
    
    
    % Inter-Trial Interval:
    % Monitor Baseline & Exit before end of all trials
    i=0; 
    while GetSecs<(backgroundtime + intertrialduration(ntrial))
        i=i+1;
        for s = 1:2
            [grip,Tgrip]=ReadGripValue('SerialMBB',Handle,s);
            gripbaseline.grip{s,ntrial}(i)=grip;
            gripbaseline.time{s,ntrial}(i)=Tgrip;

            [keyisdown, secs, keycode] = KbCheck;
            if keyisdown==1 && keycode(KbName('ESCAPE'))
                task_exit = 1;
            end
        end
        WaitSecs('UntilTime',backgroundtime+0.04*i); % sample at 50Hz (for MIE device)
    end
    for s=1:2
         baseline(s,ntrial) = mean(gripbaseline.grip{s,ntrial});
    end
    
    % Display preparatory red spot
    Screen('DrawTexture',window,background,[],rect_bg);
    Screen('DrawTexture',window,pic_spot_red,[],rect_spot_red);
    prepareTime = Screen(window,'Flip');
    
    
    % Monitor anticipation error
    i=0;level = [0 0];
    while GetSecs<(prepareTime + redSpotDuration(ntrial) + anticipationduration)
        i=i+1;
        for s = 1:2
            [grip,Tgrip]=ReadGripValue('SerialMBB',Handle,s);
            gripanticipation.grip{s,ntrial}(i)=grip - baseline(s,ntrial);
            gripanticipation.time{s,ntrial}(i)=Tgrip;
            level(s)=positivepart((grip-baseline(s,ntrial))/calibration(s,ntrial)*100);
            if level(s)>=anticipationCriterion % anticipation error
                error(ntrial)= 2 ;
                choice(ntrial)=s;
            end
        end
        if error(ntrial)== 2  % interruption precoce if effort onset too early
            break
        end
            % Display go signal
            Screen('DrawTexture',window,background,[],rect_bg);
                % Reward/effort information display during preparation
                % phase
                if GetSecs>( prepareTime + redSpotDuration(ntrial) )
                     Screen('DrawTexture',window,pic_spot_green,[],rect_spot_green);
                else
                     Screen('DrawTexture',window,pic_spot_red,[],rect_spot_red);
                end

            % Display effort scale
            if choice(ntrial)==0
                side2display = [1:2];
            else
                side2display =  choice(ntrial);
            end
            for iSide = side2display
                Display_ChoiceRE(window,pic_incentive,xEffortScale(iSide),yEffortScale,level(iSide),rewards(iSide,ntrial),effortLvls(efforts(iSide,ntrial)), effortPrecision);
                Display_RewardRE(window,pic_incentive,xReward(iSide),yEffortScale,rewards(iSide,ntrial));
            end  
            Screen(window,'Flip'); 

        WaitSecs('UntilTime',backgroundtime+0.04*i); % sample at 50Hz (for MIE device)
    end
    responsetime(ntrial)=GetSecs-anticipationduration;
    
%% Monitor performance
    if error(ntrial)~=2
        i=0; level = [0 0]; stop_monitoring=0; feedback = 0;
        while GetSecs<(responsetime(ntrial)+responseduration+reward_unit_duration*gain(ntrial))
            i=i+1;
            for s = 1:2
                [grip,Tgrip]=ReadGripValue('SerialMBB',Handle,[s]);
                gripdata.grip{s,ntrial}(i)=grip - baseline(s,ntrial);
                gripdata.time{s,ntrial}(i)=Tgrip;
                level(s)=positivepart((grip-baseline(s,ntrial))/calibration(s,ntrial)*100);
                if (level(s) >= reactionTimeCriterion ) && isnan(reactiontime(ntrial))
                    reactiontime(ntrial)= GetSecs-responsetime(ntrial);
                    error(ntrial) = 4 ; % precision error (lower /upper to the range of imposed force)
                    if choice(ntrial)==0 % define choice for this only if not a repetion of the option
                        choice(ntrial) = find(level == max(level)) ; 
                    end
                    gripdata.responseOnset{s,ntrial}(i)=1;
                else
                    gripdata.responseOnset{s,ntrial}(i)=0;
                end
                % Display go signal
                Screen('DrawTexture',window,background,[],rect_bg);
                Screen('DrawTexture',window,pic_spot_green,[],rect_spot_green);
                % Display effort scale
                if choice(ntrial)==0
                    side2display = [1:2];
                else
                    side2display =  choice(ntrial);
                end
                for iSide = side2display
                   Display_ChoiceRE(window,pic_incentive,xEffortScale(iSide),yEffortScale,level(iSide),rewards(iSide,ntrial),effortLvls(efforts(iSide,ntrial)), effortPrecision);
                   Display_RewardRE(window,pic_incentive,xReward(iSide),yEffortScale,rewards(iSide,ntrial));
                end                    
            end
            
           
            % Monitor effort duration
            for iSide = 1:2
                if (level(iSide) >= 100*effortLvls(efforts(iSide,ntrial))) 
                    forceduration(iSide,ntrial) = forceduration(iSide,ntrial)+0.04;
                    leveltime(ntrial) = GetSecs - responsetime(ntrial);
                    
                elseif stop_monitoring == 0
                         forceduration(iSide,ntrial) = 0;
                end
            end
            
            if stop_monitoring == 0 && (forceduration(1,ntrial) >= effort_unit_duration || forceduration(2,ntrial) >= effort_unit_duration ) 
                if find(level == max(level)) ~= choice(ntrial) ;
                    error(ntrial)=3; % change of mind errors
                elseif error(ntrial) ~= 5
                    error(ntrial)=0; % do not consider anymore response as a precision error but as a correct response
                end 
                stop_monitoring = 1;
            end;
            
            if stop_monitoring == 1 %&& (level(choice(ntrial)) < 100*(effortLvls(efforts(choice(ntrial),ntrial)) - effortReleaseCriterion )) % exit recording criterion
%                 break
                % Reward Feedback:
                if choice(ntrial)~=0 && feedback==0                                              % if a choice is made
                    if forceduration(choice(ntrial),ntrial) >= effort_unit_duration      % && this choice is a correct exerted force (lower amplitude reached & time duration reached )
%                         Display_RewardRE(window,pic_incentive,xReward(choice(ntrial)),yEffortScale,rewards(choice(ntrial),ntrial));
                        gain(ntrial) = rewardLvls(rewards(choice(ntrial),ntrial));
                        feedbacktime(ntrial)=GetSecs;
                        outp(888,00000001);
                        feedback = 1;
                    end
                elseif feedback==1   
                    if GetSecs> feedbacktime(ntrial) + reward_unit_duration*gain(ntrial)
                        outp(888,0);
                        break;
                    end
                end      
                
            end
            
            
            Screen(window,'Flip'); 
            WaitSecs('UntilTime',responsetime(ntrial)+0.04*i); % sample at 50Hz (for MIE device)
        end
        outp(888,0);
        total=total+gain(ntrial);



        % Data to save
        if choice(ntrial)~=1 && choice(ntrial)~=2;choice(ntrial)=0;end
%         if choice(ntrial)~=side(ntrial) && choice(ntrial)~=0 ;error(ntrial)=3;end % lateralization errors
        if choice(ntrial)~=0;
            force(ntrial)=positivepart(max(gripdata.grip{choice(ntrial),ntrial}));
            perf(ntrial)= force(ntrial)/calibration(choice(ntrial),ntrial)*100;
            sumforce(ntrial)= sum(positivepart(gripdata.grip{choice(ntrial),ntrial})./responseduration); 
            sumperf(ntrial)= sumforce(ntrial)/calibration(choice(ntrial),ntrial)*100;
        end
            
% %         % Reward Feedback:
%             %    association rule & display
%             Screen('DrawTexture',window,background,[],rect_bg);
% %             if error(ntrial)==0 
%             if choice(ntrial)~=0                                                   % if a choice is made
%                 if forceduration(choice(ntrial),ntrial) >= effort_unit_duration      % && this choice is a correct exerted force (lower amplitude reached & time duration reached )
%                     Display_RewardRE(window,pic_incentive,xReward(choice(ntrial)),yEffortScale,rewards(choice(ntrial),ntrial));
%                     gain(ntrial) = rewardLvls(rewards(choice(ntrial),ntrial));
%                     
%                     outp(888,00000001);
%                     WaitSecs(reward_unit_duration*gain(ntrial));
%                     outp(888,0);
%                 end
%             end      
%             total=total+gain(ntrial);
%             Screen(window,'Flip');
            
%             %    delivery (Send TTL Trigger)
%             if choice(ntrial)~=0 
%                 if forceduration(choice(ntrial),ntrial) >= effort_unit_duration 
%                     outp(888,00000001);
%                     WaitSecs(reward_unit_duration*gain(ntrial));
%                     outp(888,0);
%     %                 WaitSecs('UntilTime',feedbacktime(ntrial) + reward_unit_duration*max(rewardLvls) + 0.010 ); 
%                 end
%             end
%             feedbacktime(ntrial)=GetSecs;
%             Screen('DrawTexture',window,background,[],rect_bg);
%             Screen(window,'Flip');
    
    end
    
    
    % Adatping design to behavior: repeating trials/ exit experiment/...
    if error(ntrial)==0; %correct
%         error(ntrial)=0;
        omission = 0;
%         withdrawal = 0;
        
    elseif error(ntrial)~=0
        if  error(ntrial)~=3 && error(ntrial)~=4  % repeat choice if omission/anticipation error
             err=err+1;  % error counter
             if repetition(ntrial)~=2
                 repetition(ntrial+1)=1;
             else
                 repetition(ntrial+1)=2;
                 choice(ntrial+1)=choice(find(choice~=0,1,'last'));
             end
             if error(ntrial)~=2
                error(ntrial)=1; 
                omission = omission + 1; 
             end
        elseif error(ntrial)==3 || error(ntrial)==4 % repeat option if change of mind/failure error 
             if  error(ntrial)==4
                 withdrawal = withdrawal + 1;
             end
             if withdrawal <= withdrawalCriterion
                  err=err+1;  % error counter
                  repetition(ntrial+1)=2;
                  choice(ntrial+1)=choice(find(choice~=0,1,'last'));
             else
                  withdrawal = 0;
             end
        end
%          if ntrial ~= trials(end)             
%             side(end+1) = NaN;
%             side((ntrial+1):end) = [ side(ntrial) side((ntrial+1):(end-1))];
            trials(end+1) = trials(end)+1;
            calib(end+1) = calib(end);
%          end
         
%         if effortLvls(efforts(choice(ntrial),ntrial))== max(effortLvls)
%             if error(ntrial)==4 
%                 withdrawal = withdrawal + 1;
%             end
%         end
    end
        
%             if omission > omissionCriterion
%                 task_exit = 1;
%             end
%             if withdrawal > withdrawalCriterion
%                 task_exit = 1;
%             end

        
    
    if ntrial ~= trials(end)
        if repetition(ntrial+1)~=0
            rewards = [rewards , nan(2,1) ];
            rewards(:,ntrial+1:end) = [ rewards(:,ntrial) , rewards(:,ntrial+1:end-1)];
            efforts = [efforts , nan(2,1) ];
            efforts(:,ntrial+1:end) = [ efforts(:,ntrial) , efforts(:,ntrial+1:end-1)];
        end
    else
        task_exit = 1;
    end
    
    
     
end
Screen(window,'Flip');

% End of task
Screen('TextSize', window, 20);
[width,hight]=RectSize(Screen('TextBounds',window,'appuyer sur une touche pour terminer'));
Screen('DrawText',window,'appuyer sur une touche pour terminer',x-width/2,y+300,[100 100 100]);
Screen(window,'Flip');
KbWait;

Screen('CloseAll');
CloseGripDevice('SerialMBB',Handle);
clc;

% compile all recorded data
n = trials(1:iTrial);
ordinalRewardLeft = rewards(1,n);
ordinalRewardRight = rewards(2,n);
ordinalEffortLeft = efforts(1,n);
ordinalEffortRight = efforts(2,n);
cardinalRewardLeft = rewardLvls(ordinalRewardLeft);
cardinalRewardRight = rewardLvls(ordinalRewardRight);
cardinalEffortLeft = effortLvls(ordinalEffortLeft);
cardinalEffortRight = effortLvls(ordinalEffortRight);

data=[(n);rewards(1,n);rewards(2,n);efforts(1,n);efforts(2,n);
      cardinalRewardLeft ; cardinalRewardRight ; cardinalEffortLeft ;cardinalEffortRight ;
      choice(n);error(n);gain(n);repetition(n);
      force(n);sumforce(n);forceduration(n);perf(n);sumperf(n);baseline(n);
      reactiontime(n);responsetime(n);feedbacktime(n);leveltime(n)]'; 
data = mat2dataset(data);
data = set(data,'VarNames', {'trialNumber','ordinalRewardLeft','ordinalRewardRight','ordinalEffortLeft','ordinalEffortRight',...
                                'cardinalRewardLeft','cardinalRewardRight','cardinalEffortLeft','cardinalEffortRight',...
                                'sideChoice','errorType','gain','isRepeatedTrial',...
                                'force','cumulativeForce','forceDuration','perf','cumulativePerf','baseline',...
                                'reactionTime','responseTimeMarker','feedbackTimeMarker','perfTimeMarker'}); 
                            



end
