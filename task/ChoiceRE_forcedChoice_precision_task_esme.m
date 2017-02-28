function [data,gripdata,gripbaseline] = ChoiceRE_forcedChoice_precision_task_esme(root,calib,reward_unit_duration,effort_unit_duration)
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



%% 1] Initialize parameters     
    % directory (every directory called should be in the root dir)
        taskdir = [root '\task'];
        picdir= [taskdir '\display'];
        
    % grip (setting the interface, manual parametrization)
        % offset = 100
        % lower bound left = 60/ lower bound right = 50 
        Screen('Preference', 'SkipSyncTests', 1);
        MatOs_multichannel;
        clear Handle;
        Handle=InitializeGrip('SerialMBB',[1 2]); % 2 grips
        fs = 50; % sampling frequency ( Hz )

    % computer devices
    
    
    % task parameters
    effortPrecision = input('intervale de précision (ex: 0.1 = 10% force de calibration)? :');


cd(picdir);

L=1440;
H=900;
% xSpot=L/2;ySpot=H/2;
x=L/2; % background position
y=H/2;
% force visual fb
xEffortScale(1)=(0.38)*L; % left option
xEffortScale(2)=(0.38+0.5)*L; % right option
xReward = xEffortScale - (0.20)*L; 
xSpot = [0.25 , 0.75]*L; ySpot=H/2;
yEffortScale=y;

% Screen('Preference', 'SkipSyncTests', 1);
% [window]=Screen('OpenWindow',0,[0 0 0],[]);
%% extended_display
     [L,H] = Screen('WindowSize',0);
     L=1440;
     H=900;
     l2=1920;
     h2=1200;
    Screen('Preference', 'SkipSyncTests', 1); %--> or see 'help SyncTrouble'
    scrAll = Screen('Screens');
    screenNum = max(scrAll);
    x0 = 0; xx = L*1;
    y0 = 0; yy = H*1;
    [window, bound] = Screen('OpenWindow', screenNum, [], [x0 y0 xx yy]);

Screen('TextSize', window, 40);
Screen('TextFont', window, 'arial');


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

% rect_spot_red=CenterRectOnPoint(Screen('Rect',pic_spot_red),xSpot,ySpot);
% rect_spot_green=CenterRectOnPoint(Screen('Rect',pic_spot_green),xSpot,ySpot);
% rect_spot_blue=CenterRectOnPoint(Screen('Rect',pic_spot_blue),xSpot,ySpot);
for i=1:2
    rect_spot{i}=CenterRectOnPoint(Screen('Rect',pic_spot_blue),xSpot(i),ySpot);
end

rect_bg=CenterRectOnPoint(Screen('Rect',background),x,y);


%% Task parametrization
cd(root);

% ------- Experimental Factors ------- %

rewardLvls = [1, 2, 3, 4]; % incitation : cardinal levels

effortLvls = [0.2 , 0.4 , 0.6 , 0.8 ];
% effortLvls = [0.2 , 0.4 , 0.6 , 0.8 ];


laterality = [1 2]; % laterality: left/right


% ------- Trials ------- %

nTrialByBlock=32;
nBlock=1000;

% temporal design: blocks 
    trialNumber = nTrialByBlock*nBlock;
    trials = [1:trialNumber];


% current specifications
    % only extreme effort options offered
        effortBlock1 = [ 1 2 3 4  1 2 3 4  1 2 3 4  1 2 3 4  1 2 3 4  1 2 3 4  1 2 3 4  1 2 3 4]; 
        effortBlock2 = effortBlock1;
%         effortBlock1 = [ repmat(1,1,nTrialByBlock/4) , repmat(2,1,nTrialByBlock/4) , repmat(3,1,nTrialByBlock/4), repmat(4,1,nTrialByBlock/4)]; 
%         effortBlock2 = [ repmat([repmat(1,1,nTrialByBlock/16) , repmat(2,1,nTrialByBlock/16) , repmat(3,1,nTrialByBlock/16), repmat(4,1,nTrialByBlock/16)],1,4)]; 
        
    % only currentReward offered
%         currentReward = 1; % it could be 2 or 3 or 4
        rewardBlock1 = [ 1 1 1 1  2 2 2 2  3 3 3 3  4 4 4 4 1 1 1 1  2 2 2 2  3 3 3 3  4 4 4 4];  
        rewardBlock2 = rewardBlock1;  
        sideBlock = [repmat(1,1,nTrialByBlock/2) , repmat(2,1,nTrialByBlock/2)];

        sideBlock = [ 1 1 1 1  1 1 1 1 1 1 1 1  1 1 1 1 2 2 2 2  2 2 2 2 2 2 2 2  2 2 2 2 ];
        
% alternative specifications

% pseudo-random permutation across options sequence
    e1 = []; r1 = [];e2 = []; r2 = []; side = [];
    for iBlock = 1:(trialNumber/nTrialByBlock)*10
        index = randperm(nTrialByBlock); 
%         index = [1:nTrialByBlock]; 
        
        e1 = [e1, effortBlock1(index) ];
        r1 = [r1, rewardBlock1(index)];
        e2 = [e2, effortBlock2(index) ];
        r2 = [r2, rewardBlock2(index)];
        side = [side, sideBlock(index)];
    end


% --- Performance Criterions --- %
omissionCriterion = 30; % maximal n° of successive omissions authorized
withdrawalCriterion = 10; % maximal n° of successive failed trials which are repeated
anticipationCriterion = 0.10*100; % (%) of force level exerted before the time of choice
reactionTimeCriterion = effortLvls(1)*100; % (%) of force level exerted after starting onsert
% effortPrecision = 0.8;
effortReleaseCriterion = 0;
calib_increment = 10;


% ------- Temporal parameters (sec.)------- %

starting_duration = 3; % time before the experiment start running
% spotduration=1;
min_redSpotDuration=0.75;
max_redSpotDuration=2.25; 
% min_cueduration=0.5;
% max_cueduration=1.5;  
anticipationduration=0.1;
responseduration=2.5;
% intertrialduration=3;
min_intertrialduration=1;
max_intertrialduration=2;
failure_delay = 4;
punishTrialDuration = 0.5;
punishDuration = 1;

redSpotDuration=rand(1,10*numel(trials))*(max_redSpotDuration-min_redSpotDuration)+min_redSpotDuration;
intertrialduration=rand(1,10*numel(trials))*(max_intertrialduration-min_intertrialduration)+min_intertrialduration;


% ------- Pre-Allocating variables to save ------- %
nt_block = 4;
rewards=nan(2,10*numel(trials));
efforts=nan(2,10*numel(trials));
rewards(:,1:nt_block) = [r1(1:nt_block) ; r2(1:nt_block)]; % initial trial's factors
efforts(:,1:nt_block) = [e1(1:nt_block) ; e2(1:nt_block)];

preparetime=nan(1,10*numel(trials));
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
% choice=zeros(1,10*numel(trials));
choice=side;
error=ones(1,10*numel(trials));
repetition=zeros(1,10*numel(trials));


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

initGUI;



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
    
% Adaptive Control of Design
    if  isnan(rewards(:,ntrial))
%         % display intructions
%         Screen(window,'Flip'); 
%         Screen('TextSize', window, 20);
%         [width,hight]=RectSize(Screen('TextBounds',window,'répéter bloc'));
%         Screen('DrawText',window,'répéter bloc',x-width/2-x/2,y+300,[100 100 100]);
%         [width,hight]=RectSize(Screen('TextBounds',window,'bloc suivant'));
%         Screen('DrawText',window,'bloc suivant',x-width/2+x/2,y+300,[100 100 100]);
%         Screen(window,'Flip');
%         
%         % monitor response
%         [~, keyCode] = KbWait;
%         nt_block = 4;
%         if  keyCode(key.left)==1
%             side = [ side, nan(1,nt_block)];
%             rewards = [ rewards, nan(2,nt_block)];
%             efforts = [ efforts, nan(2,nt_block)];
%             
%             side(:,ntrial:end) =  [ side(:,ntrial-nt_block:ntrial-1) , side(:,ntrial:end-nt_block)];
%             r1 = [ r1, nan(1,nt_block)];
%             r1(ntrial:end) = [ r1(ntrial-nt_block:ntrial-1) , r1(ntrial:end-nt_block) ];
%             r2 = [ r2, nan(1,nt_block)];
%             r2(ntrial:end) = [ r2(ntrial-nt_block:ntrial-1) , r2(ntrial:end-nt_block) ];
%             e1 = [ e1, nan(1,nt_block)];
%             e1(ntrial:end) = [ e1(ntrial-nt_block:ntrial-1) , e1(ntrial:end-nt_block) ];
%             e2 = [ e2, nan(1,nt_block)];
%             e2(ntrial:end) = [ e2(ntrial-nt_block:ntrial-1) , e2(ntrial:end-nt_block) ];
%             
%         elseif  keyCode(key.right)==1
%         else
%             task_exit = 1;
%         end

        rewards(:,ntrial:ntrial+nt_block-1) =  [r1(ntrial:ntrial+nt_block-1) ; r2(ntrial:ntrial+nt_block-1)];
        efforts(:,ntrial:ntrial+nt_block-1) =  [e1(ntrial:ntrial+nt_block-1) ; e2(ntrial:ntrial+nt_block-1)];
        choice=side;

    end
        
    
    
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
            [gui]=updateGUI(gui,s,grip);
            gripbaseline.grip{s,ntrial}(i)=grip;
            gripbaseline.time{s,ntrial}(i)=Tgrip;

            [keyisdown, secs, keycode] = KbCheck;
            if keyisdown==1 && keycode(KbName('ESCAPE'))
                task_exit = 1;
            elseif keyisdown==1 && keycode(KbName('SPACE'))
                withdrawal = withdrawalCriterion  ;
            end
        end
        WaitSecs('UntilTime',backgroundtime+(1/fs)*i); % sample at 50Hz (for MIE device)
    end
    for s=1:2
         baseline(s,ntrial) = mean(gripbaseline.grip{s,ntrial});
         gui.b(s) =  baseline(s,ntrial) ;
         [gui]=updateGUI(gui,s,grip);
    end
    
    % select side (only one if forced choice, two if free choice)
    if choice(ntrial)==0
        side2display = [1:2];
    else
        side2display =  choice(ntrial);
    end
    
    % Display preparatory red spot
    Screen('DrawTexture',window,background,[],rect_bg);
    for s = side2display
        Screen('DrawTexture',window,pic_spot_red,[],rect_spot{s});
    end
    preparetime(ntrial) = Screen(window,'Flip');
    gui.offer(end) = 1;

    
    % Monitor anticipation error
    i=0;level = [0 0];
    while GetSecs<(preparetime(ntrial) + redSpotDuration(ntrial) + anticipationduration)
        i=i+1;
        for s = 1:2
            [grip,Tgrip]=ReadGripValue('SerialMBB',Handle,s);
            [gui]=updateGUI(gui,s,grip);
            gripanticipation.grip{s,ntrial}(i)=grip - baseline(s,ntrial);
            gripanticipation.time{s,ntrial}(i)=Tgrip;
            level(s)=positivepart((grip-baseline(s,ntrial))/calibration(s,ntrial)*100);
            if level(s)>=anticipationCriterion % anticipation error
                reactiontime(ntrial)= GetSecs-preparetime(ntrial);
                error(ntrial)= 2 ;
%                 choice(ntrial)=s;
            end
        end
        if error(ntrial)== 2  % interruption precoce if effort onset too early
            break
        end
            % Display go signal
            Screen('DrawTexture',window,background,[],rect_bg);
                % Reward/effort information display during preparation
                % phase
                 for s = side2display
                    if GetSecs>( preparetime(ntrial) + redSpotDuration(ntrial) )
                         Screen('DrawTexture',window,pic_spot_green,[],rect_spot{s});
                    else
                         Screen('DrawTexture',window,pic_spot_red,[],rect_spot{s});
                    end
                 end

            % Display effort scale
            for iSide = side2display
                Display_ChoiceRE(window,pic_incentive,xEffortScale(iSide),yEffortScale,level(iSide),rewards(iSide,ntrial),effortLvls(efforts(iSide,ntrial)), effortPrecision, 4 );
                Display_RewardRE(window,pic_incentive,xReward(iSide),yEffortScale,rewards(iSide,ntrial));
            end  
            Screen(window,'Flip'); 

        WaitSecs('UntilTime',backgroundtime+(1/fs)*i); % sample at 50Hz (for MIE device)
    end
    responsetime(ntrial)=GetSecs-anticipationduration;
    gui.period(end) = 1;

%% Monitor performance
    if error(ntrial)~=2
        i=0; level = [0 0]; stop_monitoring=0; feedback = 0;pf = [0 0];
        while GetSecs<(responsetime(ntrial)+responseduration+reward_unit_duration*gain(ntrial))
            i=i+1;
            for s = 1:2
                [grip,Tgrip]=ReadGripValue('SerialMBB',Handle,[s]);
                [gui]=updateGUI(gui,s,grip);
                gripdata.grip{s,ntrial}(i)=grip - baseline(s,ntrial);
                gripdata.time{s,ntrial}(i)=Tgrip;
                level(s)=positivepart((grip-baseline(s,ntrial))/calibration(s,ntrial)*100);
                if (level(s) >= reactionTimeCriterion ) && isnan(reactiontime(ntrial))
                    reactiontime(ntrial)= GetSecs-preparetime(ntrial);
                    error(ntrial) = 4 ; % precision error (lower /upper to the range of imposed force)
                    if repetition(ntrial)~=2 % define choice for this only if not a repetion of the option
%                         choice(ntrial) = find(level == max(level)) ; 
                    end
                    gripdata.responseOnset{s,ntrial}(i)=1;
                else
                    gripdata.responseOnset{s,ntrial}(i)=0;
                end
                % Display go signal
                Screen('DrawTexture',window,background,[],rect_bg);
                for s = side2display
                    Screen('DrawTexture',window,pic_spot_green,[],rect_spot{s});
                end
                % Display effort scale
                for iSide = side2display
                   pf(iSide) = max( [ level(iSide) ,  pf(iSide) ] ) ;
%                    Display_ChoiceRE(window,pic_incentive,xEffortScale(iSide),yEffortScale, pf(iSide),rewards(iSide,ntrial),effortLvls(efforts(iSide,ntrial)), effortPrecision, error(ntrial) );
                   Display_ChoiceRE(window,pic_incentive,xEffortScale(iSide),yEffortScale, level(iSide),rewards(iSide,ntrial),effortLvls(efforts(iSide,ntrial)), effortPrecision, error(ntrial) );
                   Display_RewardRE(window,pic_incentive,xReward(iSide),yEffortScale,rewards(iSide,ntrial));
                end                    
            end
            
           
            % Monitor effort duration
            for iSide = 1:2
                if (level(iSide) >= 100*effortLvls(efforts(iSide,ntrial))) 
                    forceduration(iSide,ntrial) = forceduration(iSide,ntrial)+0.04;
                    leveltime(ntrial) = GetSecs - responsetime(ntrial);
%                     Screen('DrawTexture',window,pic_spot_blue,[],rect_spot{s});

                    % upper bound
                    if level(iSide) >= 100*(effortLvls(efforts(iSide,ntrial))+ effortPrecision)  
                        error(ntrial) = 5 ; % precision error (upper to the range of imposed force)
                        stop_monitoring = 1;
                        Screen('DrawTexture',window,pic_spot_green,[],rect_spot{s});
%                             % punishment delay
%                             Screen('DrawTexture',window,background,[],rect_bg);
%                             for iSide = side2display
%                                Screen('DrawTexture',window,pic_spot_green,[],rect_spot{iSide});
%                                Display_ChoiceRE(window,pic_incentive,xEffortScale(iSide),yEffortScale, level(iSide),rewards(iSide,ntrial),effortLvls(efforts(iSide,ntrial)), effortPrecision, error(ntrial) );
%                             end  
%                             Screen(window,'Flip'); 
%                             WaitSecs(punishDuration);
                        break;
                    end
                    
                    
                elseif stop_monitoring == 0
                         forceduration(iSide,ntrial) = 0;
                end
            end
            
            if stop_monitoring == 0 && (forceduration(1,ntrial) >= effort_unit_duration || forceduration(2,ntrial) >= effort_unit_duration ) 
                if find(level == max(level)) ~= choice(ntrial) ;
                    error(ntrial)=3; % change of mind errors
                    break;
                elseif error(ntrial) ~= 5
                    error(ntrial)=0; % do not consider anymore response as a precision error but as a correct response
                end 
                stop_monitoring = 1;
                gui.response(end) = 1;
            end;
            
            if stop_monitoring == 1  % exit recording criterion
%                 break
                % Reward Feedback:
                if choice(ntrial)~=0 && feedback==0   && error(ntrial) == 0                                     % if a choice is made
                    if forceduration(choice(ntrial),ntrial) >= effort_unit_duration ...
                            && (level(choice(ntrial)) < 100*(effortLvls(efforts(choice(ntrial),ntrial)) - effortReleaseCriterion ))...    % be under the  lower bound
%                             && (level(choice(ntrial)) < 100*(effortLvls(efforts(choice(ntrial),ntrial)) + effortPrecision  ))   % be under the upper bound
%                             && (level(choice(ntrial)) >= 100*(effortLvls(efforts(choice(ntrial),ntrial))) ) %  be over the  lower bound
%                         Display_RewardRE(window,pic_incentive,xReward(choice(ntrial)),yEffortScale,rewards(choice(ntrial),ntrial));    
%                          Screen('DrawTexture',window,pic_spot_blue,[],rect_spot{s});
                         if error(ntrial) ==0  
                            gain(ntrial) = rewardLvls(rewards(choice(ntrial),ntrial));
                            feedbacktime(ntrial)=GetSecs;
                            outp(888,00000001);
                            gui.io = 1;
                            feedback = 1;
                         end
                    end
                elseif error(ntrial) == 5    
                        outp(888,0);
                        gui.io = 0;
                        break;
                elseif feedback==1   
                    Screen('DrawTexture',window,pic_spot_blue,[],rect_spot{s});
                    if GetSecs> feedbacktime(ntrial) + reward_unit_duration*gain(ntrial)
                        outp(888,0);
                        gui.io = 0;
                        break;
                    end
                end      
                
            end
            
            
            Screen(window,'Flip'); 
            WaitSecs('UntilTime',responsetime(ntrial)+(1/fs)*i); % sample at 50Hz (for MIE device)
        end
        outp(888,0);
        gui.io = 0;
        total=total+gain(ntrial);

        % Punishment if overshoot: trial stop + delay
        if error(ntrial) == 5
            % punishment delay
            Screen('DrawTexture',window,background_lateralizationError,[],rect_bg);
            for iSide = side2display
               Screen('DrawTexture',window,pic_spot_green,[],rect_spot{iSide});
               Display_ChoiceRE(window,pic_incentive,xEffortScale(iSide),yEffortScale, level(iSide),rewards(iSide,ntrial),effortLvls(efforts(iSide,ntrial)), effortPrecision, error(ntrial) );
            end  
            Screen(window,'Flip'); 
            WaitSecs(punishDuration);
        end

        % Data to save
        if choice(ntrial)~=1 && choice(ntrial)~=2;choice(ntrial)=0;end
%         if choice(ntrial)~=side(ntrial) && choice(ntrial)~=0 ;error(ntrial)=3;end % lateralization errors
        if choice(ntrial)~=0;
            force(ntrial)=positivepart(max(gripdata.grip{choice(ntrial),ntrial}));
            perf(ntrial)= force(ntrial)/calibration(choice(ntrial),ntrial)*100;
            sumforce(ntrial)= sum(positivepart(gripdata.grip{choice(ntrial),ntrial})./responseduration); 
            sumperf(ntrial)= sumforce(ntrial)/calibration(choice(ntrial),ntrial)*100;
        end
            
    
    end
    
    % update design
    updateDesign_choiceRE;

    
    gui.error = [gui.error , error(ntrial) ];  
    gui.choice = [gui.choice , choice(ntrial) ]; 
    gui.reward = [gui.reward , gain(ntrial) ]; 
        

    
    
    
     
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

ordinalRewardLeft = rewards(1,1:ntrial);
ordinalRewardRight = rewards(2,1:ntrial);
ordinalEffortLeft = efforts(1,1:ntrial);
ordinalEffortRight = efforts(2,1:ntrial);
cardinalRewardLeft = rewardLvls(ordinalRewardLeft);
cardinalRewardRight = rewardLvls(ordinalRewardRight);
cardinalEffortLeft = effortLvls(ordinalEffortLeft);
cardinalEffortRight = effortLvls(ordinalEffortRight);

data=[(1:ntrial);rewards(1,1:ntrial);rewards(2,1:ntrial);efforts(1,1:ntrial);efforts(2,1:ntrial);
      cardinalRewardLeft ; cardinalRewardRight ; cardinalEffortLeft ;cardinalEffortRight ;
      choice(1:ntrial);error(1:ntrial);gain(1:ntrial);repetition(1:ntrial);
      force(1:ntrial);sumforce(1:ntrial);forceduration(1:ntrial);perf(1:ntrial);sumperf(1:ntrial);baseline(1:ntrial);
      reactiontime(1:ntrial);
      feedbacktime(1:ntrial);leveltime(1:ntrial);
      intertrialduration(1:ntrial); redSpotDuration(1:ntrial);]'; 
data = mat2dataset(data);
data = set(data,'VarNames', {'trialNumber','ordinalRewardLeft','ordinalRewardRight','ordinalEffortLeft','ordinalEffortRight',...
                                'cardinalRewardLeft','cardinalRewardRight','cardinalEffortLeft','cardinalEffortRight',...
                                'sideChoice','errorType','gain','isRepeatedTrial',...
                                'force','cumulativeForce','forceDuration','perf','cumulativePerf','baseline',...
                                'reactionTime','feedbackTimeMarker','perfTimeMarker','intertrialTime','offerTime'}); 
                            



end
