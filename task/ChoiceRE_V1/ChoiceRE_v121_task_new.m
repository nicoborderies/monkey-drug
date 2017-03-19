function [data,gripdata,gripbaseline] = ChoiceRE_v121_task(root,calib,reward_unit_duration,effort_unit_duration)
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
% END
% SAVE data
% PROCESS offline analysis



% clear all;
% close all;



%% 1] Initialize external parameters     
    % directory (every directory called should be in the root dir)
        taskdir = [root '\task'];
        picdir= [taskdir '\display'];
        
    % devices
        % grip (setting the interface, manual parametrization)
            % offset = 100
            % lower bound left = 60/ lower bound right = 50 
            MatOs_multichannel;
            clear Handle;
            Handle=InitializeGrip('SerialMBB',[1 2]); % 2 grips

        % computer devices
        init_computer;

        % serial port
            config_io; % installation
            outp(888,0); % set to zero
            
    % random number generator
         rand('state',sum(100*clock));

    % images
        cd(picdir);

        background_correct = Screen('MakeTexture',window,imread('bg_100.bmp'));
        background_error = Screen('MakeTexture',window,imread('bg_lateralization.bmp'));

        pic_incentive = Screen('MakeTexture',window,imread(['reward_drop2.bmp']));

        pic_spot_red=Screen('MakeTexture',window,imread('spot_1_red.bmp'));
        pic_spot_green=Screen('MakeTexture',window,imread('spot_1_green.bmp'));
        pic_spot_blue=Screen('MakeTexture',window,imread('spot_1_blue.bmp'));

        rect_spot_red=CenterRectOnPoint(Screen('Rect',pic_spot_red),xSpot,ySpot);
        rect_spot_green=CenterRectOnPoint(Screen('Rect',pic_spot_green),xSpot,ySpot);
        rect_spot_blue=CenterRectOnPoint(Screen('Rect',pic_spot_blue),xSpot,ySpot);
        rect_bg=CenterRectOnPoint(Screen('Rect',background),x,y);


%% 2] Initialize task parameters

    cd(root);
    
    % design
        [design,data] = init_design(calib);


%% 3 ] Execute Task 

    cd(taskdir);

    % control starting
        % display
        Screen('TextSize', window, 20);
        [width,hight]=RectSize(Screen('TextBounds',window,'appuyer sur une touche pour commencer'));
        Screen('DrawText',window,'appuyer sur une touche pour commencer',x-width/2,y+300,[100 100 100]);
        Screen(window,'Flip');
        
        % monitor response
        KbWait;
        Screen(window,'Flip');
        WaitSecs(starting_duration);
   
    %% trial iteration
        while task_exit==0    % exit condition

            % increment
                ntrial = ntrial +1; % trial 

            %% Prepare trial
                % adapt design


                % display background
                    if ntrial==1 
                        Screen('DrawTexture',window,background_correct,[],rect_bg);
                    else
                        switch error(ntrial-1)
                            case 0
                                Screen('DrawTexture',window,background_correct,[],rect_bg);
                            otherwise
                                Screen('DrawTexture',window,background_error,[],rect_bg);
                        end
                    end

                % mesure time
                    backgroundtime=Screen(window,'Flip');


                % inter-trial interval (baseline monitoring)
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

            % display preparatory red spot
                Screen('DrawTexture',window,background,[],rect_bg);
                Screen('DrawTexture',window,pic_spot_red,[],rect_spot_red);
                prepareTime = Screen(window,'Flip');


            %% Option display
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
                            if repetition(ntrial)~=2 % define choice for this only if not a repetion of the option
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
                            end
                        end      
                        total=total+gain(ntrial);

                    end


                    Screen(window,'Flip'); 
                    WaitSecs('UntilTime',responsetime(ntrial)+0.04*i); % sample at 50Hz (for MIE device)
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
                    trials(end+1) = trials(end)+1;
                    calib(end+1) = calib(end);
            end



            if ntrial ~= trials(end)
                    rewards(:,ntrial+1) = [r1(ntrial+1-err) ; r2(ntrial+1-err)]; 
                    efforts(:,ntrial+1) = [e1(ntrial+1-err) ; e2(ntrial+1-err)];
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
      reactiontime(1:ntrial);responsetime(1:ntrial);feedbacktime(1:ntrial);leveltime(1:ntrial)]'; 
data = mat2dataset(data);
data = set(data,'VarNames', {'trialNumber','ordinalRewardLeft','ordinalRewardRight','ordinalEffortLeft','ordinalEffortRight',...
                                'cardinalRewardLeft','cardinalRewardRight','cardinalEffortLeft','cardinalEffortRight',...
                                'sideChoice','errorType','gain','isRepeatedTrial',...
                                'force','cumulativeForce','forceDuration','perf','cumulativePerf','baseline',...
                                'reactionTime','responseTimeMarker','feedbackTimeMarker','perfTimeMarker'}); 
                            



end
