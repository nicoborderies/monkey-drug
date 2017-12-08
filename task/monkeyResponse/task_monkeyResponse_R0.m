function [] = task_monkeyResponse_R0(subName,nSession,experimenterName,varargin)
%task_monkeyResponse_R0 - 
%
% specifications: 
%       - task structure: 
%       - experimental conditions: 
%       - randomization:
%       - calibration: 
%       - start criterion:
%       - stop criterion:
%       - trial structure: 
%       - trial update criterion:
%       - minimal ntrial: 
%
%
% Syntax: 
%
% Inputs:
%    subName - subject identification number (string)
%    session - session number (double)
%    options:
%       fullscreen - logical flag to display full-screen the task,  default value=1 (logical)
%
% Outputs:
%
% Pseudocode:
%   setup OPTIONS
%   declare PARAMETERS
%   setup METADATA
%   setup CONFIGURATIONS
%   load STIMULI
%   make CONDITIONS
%   while START_CRITERION
%       display START_STIMULI
%       get RESPONSE
%   end
%   while ~STOP_CRITERION
%       display CONDITION_STIMULI
%       get RESPONSE
%       check TRIAL_VALIDITY
%       send FEEDBACK
%       if STOP_CRITERION_1
%           check STOP_CRITERION_2
%           if STOP_CRITERION_2
%               break
%           end
%       end
%       update CONDITION
%       display TRIAL_SUMMARY
%   end
%   save DATA & METADATA
%   display TASK_SUMMARY
%
%
% Example: 
%  
%
% Requirements: 
%   Subfunctions: 
%   MAT-files: 
%   MATLAB products: MATLAB, Statistics and Machine Learning Toolbox,
%                    Psychtoolbox
%
% See also:
%
% Author:  Nicolas Borderies
% email address: nico.borderies@gmail.com 
% March 2017; Last revision: March 2017

%% Configuration
% -----------------------------------------------

% Check Options
%-----------------------------------------------
% defaults options
taskName = 'monkeyResponse_R0';
optionList = {'fullscreen','rewardControllerFlag','dynamoFlag'};
subList = {'aliosha','bob','dracula','esmeralda'};
fullscreen=1;
rewardControllerFlag=1;
dynamoFlag=1;
Parameters = struct;

% update to specified options
inputList = varargin;
nInput = numel(inputList);
for iarg=1:nInput
   if ischar(inputList{iarg})
       if ismember(inputList{iarg},optionList)
           eval([ 'Parameters.options.' inputList{iarg} ' = inputList{iarg+1};' ]);
       end
   end
end
% nOptions = numel(fieldnames(Parameters.options));

% Get Parameters
%-----------------------------------------------

% Set Interface Configuration
%-----------------------------------------------
Interface = monkeyResponse_app ;

% Get Metadata
%-----------------------------------------------
if exist('subName')~=1
    subName=input('subject name?','s');
elseif ~ismember(subName,subList)
    error(' "subid" should be one of these names: aliosha, bob, dracula, esmeralda ');
end

if exist('nSession')~=1
    nSession=input('session number ?');
elseif ~isnumeric(nSession) && round(nSession)==nSession
    error(' "nSession" should be an integer number ');
end

if exist('experimenterName')~=1
    experimenterName=input('experimenter name ?','s');
elseif ~ischar(experimenterName)
    error(' "experimenterName" should be a string ');
end

clck = clock;
time = [num2str(clck(2)) '_' num2str(clck(3)) '_' num2str(clck(1)) '_' num2str(clck(4)) 'h' num2str(clck(5)) 'min' ];

Metadata = struct;
Metadata.taskName = taskName;
Metadata.subName = subName;
Metadata.nSession = nSession;
Metadata.experimenterName = experimenterName;
Metadata.time = time;
Metadata.duration = '0h0min';

% Set Directory Configuration
%-----------------------------------------------
% define script directory
Metadata.taskdir=pwd;
cd ..
uproot=pwd;

% define group & subject result directory
resultdir=[Metadata.taskdir '\resultats'];
if exist(resultdir,'dir')~=7
    mkdir resultats
end
cd(resultdir);
Metadata.subdir=[resultdir '\' Metadata.subName ];
if exist(Metadata.subdir,'dir')~=7
    mkdir(resultdir,['\' Metadata.subName]);
end
Metadata.dataFileName = strcat(Metadata.taskName,'_',Metadata.subName,'_n',num2str(Metadata.nSession),'_',Metadata.time);

% subdirectories
%%% define image directory
cd(Metadata.taskdir);
Metadata.imgdir=[Metadata.taskdir '\images'];
addpath(genpath(Metadata.imgdir));
%%% define subfunctions directory
% cd(Metadata.taskdir);
% Metadata.subfunc=[Metadata.taskdir '\subfunctions'];
% addpath(genpath(Metadata.subfunc));
%%% define serial port controller directory
cd(Metadata.taskdir);
portdir=[Metadata.taskdir '\IO32'];
addpath(genpath(portdir));
%%% define dynamo directory
cd(Metadata.taskdir);
portdir=[Metadata.taskdir '\newMatlabOscillo'];
addpath(genpath(portdir));

% Set Screen Configuration
%-----------------------------------------------
screenid = Screen('Screens');
Screen('Preference', 'SkipSyncTests', 1); %--> or see 'help SyncTrouble'

% open testing window
Parameters.display.window = Screen('OpenWindow',2,[0 0 0],[]); % testing window

% configure
Screen('TextSize', Parameters.display.window, 40);
Screen('TextFont', Parameters.display.window, 'arial');
[Parameters.display.W, Parameters.display.H]=Screen('WindowSize',Parameters.display.window);
Parameters.display.xcenter = Parameters.display.W/2;
Parameters.display.ycenter  = Parameters.display.H/2;

% Reset random number generator
%-----------------------------------------------
rand('state',sum(100*clock));

% Set device configurations
%-----------------------------------------------
Device = struct;

% keyboard
[~,~,keycode] = KbCheck;
DisableKeysForKbCheck(find(keycode==1));

KbName('UnifyKeyNames');
Device.key.left = KbName('LeftArrow');
Device.key.right = KbName('RightArrow');
Device.key.space = KbName('Space') ;
Device.key.escape = KbName('ESCAPE') ;
Device.key.return =  KbName('Return');

% Reward controller (Parallel Port Com)
if rewardControllerFlag
    config_io;
    outp(888,0); % put all pins to zero
%     Device.reward.isOpen = 0;
    Device.reward.time2volume = 0.9;
    Device.reward.volumeUnit = 2;
    Device.reward.timeUnit = Device.reward.volumeUnit/Device.reward.time2volume;
end

% Dynamometer 
if dynamoFlag
    %%% user-controlled configuration
    % MatOs_multichannel;
    % disp('check that the signal for both channels is around: 100 au. (+/-10)');
    %%% automatic configuration
    Handle = InitializeGrip('SerialMBB',[1 2]); % 2 grips
    Device.dynamo = Handle; % 2 grips
end

% Video Camera

% Motion Sensor

% Microphone Sensor

% Audio speaker
InitializePsychSound(1);
Device.audio.freq = 1000;
Device.audio.duration = 0.10;
Device.audio.intensity = 1;
Device.audio.handle = PsychPortAudio('Open', [], 1, 1, [], 2);
PsychPortAudio('Volume', Device.audio.handle, Device.audio.intensity);
%%% condition sound
soundfile = 'xylophone_sound.wav';
[s,f] = audioread(soundfile);s = s';
sound1 = PsychPortAudio('CreateBuffer', Device.audio.handle, s);
%%% reward sound
Device.audio.beep = MakeBeep(Device.audio.freq, Device.audio.duration);
sound2 = PsychPortAudio('CreateBuffer', Device.audio.handle, [Device.audio.beep; Device.audio.beep]);


% Load Stimuli
%-----------------------------------------------
Stimuli = struct;

cd(Metadata.imgdir); % enter img dir
% images

% spots
Stimuli.xspot = Parameters.display.xcenter;
Stimuli.yspot = Parameters.display.ycenter;
Stimuli.wspot = Parameters.display.W*0.025;
Stimuli.hspot = Parameters.display.W*0.025;
Stimuli.backgroundImg = ones(1600,900,3)*0.5;
Stimuli.redspotCol = [255 0 0] ;
Stimuli.greenspotCol = [0 255 0] ;
Stimuli.bluespotCol = [0 0 255] ;
Stimuli.backgroundRect = CenterRectOnPoint(Screen('Rect',Parameters.display.window),Stimuli.xspot,Stimuli.yspot);
Stimuli.redspotRect = [Stimuli.xspot-Stimuli.wspot/2, Stimuli.yspot-Stimuli.hspot/2, Stimuli.xspot+Stimuli.wspot/2, Stimuli.yspot+Stimuli.hspot/2 ] ;
Stimuli.greenspotRect = [Stimuli.xspot-Stimuli.wspot/2, Stimuli.yspot-Stimuli.hspot/2, Stimuli.xspot+Stimuli.wspot/2, Stimuli.yspot+Stimuli.hspot/2 ] ;
Stimuli.bluespotRect = [Stimuli.xspot-Stimuli.wspot/2, Stimuli.yspot-Stimuli.hspot/2, Stimuli.xspot+Stimuli.wspot/2, Stimuli.yspot+Stimuli.hspot/2 ] ;

% option positions

% force visual feedback  

cd(Metadata.taskdir); % exit

% Make experimental design
%-----------------------------------------------

% durations 
Parameters.durations.itiRange = [0.5 1];
Parameters.durations.conditionRange = [0.5 1];
Parameters.durations.responsenRange = [ 1 ];
Parameters.durations.repetitionPenalty = 0.1;
Parameters.durations.responseRewardDelay = 0.05;

% conditions
%%% sample size
Parameters.conditions.nMinTrial = 100;
Parameters.conditions.nMaxTrial = Parameters.conditions.nMinTrial*20;
Parameters.conditions.nTrialByBlock = 10;

%%% variables levels

%%% variables creation
conditionNames = {'nCondition','itiDuration','conditionDuration','responseDuration'};
designNames = [ {'nTrial','nRepetition'} , conditionNames];

Conditions = array2table(nan(Parameters.conditions.nMaxTrial,numel(conditionNames)),'VariableNames',conditionNames);
Design = array2table(nan(Parameters.conditions.nMaxTrial,numel(designNames)),'VariableNames',designNames);


%%% randomization
n = Parameters.conditions.nMaxTrial;
Design.nTrial = [1:n]';
Design.nRepetition = zeros(n,1);

%%% --- durations
Conditions.nCondition = [1:n]';
% iti
dur = Parameters.durations.itiRange(1) + range(Parameters.durations.itiRange)*rand(n,1);
Conditions.itiDuration = dur;
% condition
dur = Parameters.durations.conditionRange(1) + range(Parameters.durations.conditionRange)*rand(n,1);
Conditions.conditionDuration = dur;
% response
dur = repmat(Parameters.durations.responsenRange,n,1);
Conditions.responseDuration = dur;

% go/nogo condition
Conditions.isGo = randi([0 1],n,1);
Conditions.conditionDuration = (1 - Conditions.isGo).*Conditions.conditionDuration;
Conditions.responseDuration = (Conditions.isGo).*Conditions.responseDuration;

Design.nCondition(1) = Conditions.nCondition(1);
Design.itiDuration(1) = Conditions.itiDuration(1);
Design.conditionDuration(1) = Conditions.conditionDuration(1);
Design.responseDuration(1) = Conditions.responseDuration(1);
Design.isGo(1) = Conditions.isGo(1);


% experimental criterions
Parameters.criterions.responseThreshold = 0.1;
Parameters.criterions.maxRepetition = 20;
Parameters.criterions.maxNumberOfTrial = 400;
Parameters.criterions.dynamoCalib = 200;
Parameters.criterions.ntrial2updateCalib = 10;
Parameters.criterions.updateCalibModulation = 1;


% Data preparation
%-----------------------------------------------
behaviorNames = {'rewardOutcome','error','participation','correctResponse','peakForce','responseTime'};
Behaviors = array2table(nan(Parameters.conditions.nMaxTrial,numel(behaviorNames)),'VariableNames',behaviorNames);

n = Parameters.conditions.nMaxTrial;
Behaviors.rewardOutcome = zeros(n,1);
Behaviors.error = nominal(repmat('miss',n,1));
Behaviors.participation = zeros(n,1);
Behaviors.correctResponse = zeros(n,1);


Dynamo = struct('forceSignal',[],'time',[],'baselineSignal',[],'nTrial',[],...
                'iti',[],'condition',[],'response',[]);

% Initialization
%-----------------------------------------------

Counters.iCondition = 1;
Counters.iTrial = 0;
Counters.nParticipation = 0;
Counters.nCorrect = 0;
Counters.vReward = 0;

Events.startTask = 0;
Events.pauseTask = 0;
Events.stopTask = 0;
Events.onsetCondition = 0;
Events.onsetResponse = 0;
Events.response = 0;
Events.dynamoBaseline = [0;0];

Interface = set_properties(Interface,Metadata,Parameters,Counters,Events,...
Behaviors,Conditions,Design,Device,Dynamo);


%% Check Start Criterion
% ----------------------------------------------

while Events.startTask == 0
    % check keyboard response
    clc;
    disp('Press any key to start the task');
    WaitSecs(1);
    KbWait;
    Events.startTask = 1;
end
WaitSecs(1);
clck = clock;

%% Run task 
% ---------------------------------------------

while Events.stopTask == 0
    
    Counters.iTrial = Counters.iTrial+1;
    
    % inter-trial interval
    Screen(Parameters.display.window,'Flip');
    set_periodLamp(Interface,'iti');
    Events.onsetITI = GetSecs;
    while GetSecs < Events.onsetITI + Design.itiDuration(Counters.iTrial) % checking loop
        
        % monitor dynamo
        if dynamoFlag
            signal = nan(2,1);
            for s = 1:2
                [signal(s),t]=ReadGripValue('SerialMBB',Device.dynamo,s);
            end
            Dynamo.forceSignal = [Dynamo.forceSignal , signal];
            Dynamo.baselineSignal = [Dynamo.baselineSignal , Events.dynamoBaseline ];
            Dynamo.time = [Dynamo.time , t];
            Dynamo.nTrial = [Dynamo.nTrial , Counters.iTrial];
            Dynamo.iti = [Dynamo.iti , 1];
            Dynamo.condition = [Dynamo.condition , 0];
            Dynamo.response = [Dynamo.response , 0];
            force = signal - Events.dynamoBaseline ;
            if force(1)>=Parameters.criterions.responseThreshold*Parameters.criterions.dynamoCalib...
            || force(2)>=Parameters.criterions.responseThreshold*Parameters.criterions.dynamoCalib
                Events.response = 1;
                if Counters.iTrial>1
                    Behaviors.error(Counters.iTrial-1) = 'late';
                    Behaviors.participation(Counters.iTrial-1) = 1;
                    Behaviors.responseTime(Counters.iTrial-1) = GetSecs - Events.onsetITI ;
                end
            end
        end
        
    end
    Events.dynamoBaseline = mean(Dynamo.forceSignal(:,Dynamo.nTrial==Counters.iTrial & Dynamo.iti==1),2);

    
    % display condition stimuli
    Screen('FillRect', Parameters.display.window, Stimuli.redspotCol, Stimuli.redspotRect);
    Screen(Parameters.display.window,'Flip');
    % audio & experimenter feedback
    set_periodLamp(Interface,'condition');
    PsychPortAudio('FillBuffer', Device.audio.handle, sound1);
    PsychPortAudio('Start', Device.audio.handle, 1, 0, 0);
    % onset
    Events.response = 0;
    Events.onsetCondition = GetSecs;
    while GetSecs < Events.onsetCondition + Design.conditionDuration(Counters.iTrial) % checking loop
        
        % monitor keyboard response
        [keyisdown, secs, keycode] = KbCheck;
        if keyisdown==1 && keycode(Device.key.space)==1
            Events.response = 1;
            Behaviors.error(Counters.iTrial) = 'early';
            Behaviors.participation(Counters.iTrial) = 1;
            Behaviors.responseTime(Counters.iTrial) = GetSecs - Events.onsetCondition ;
        end
        
        % monitor dynamo
        if dynamoFlag
            signal = nan(2,1);
            for s = 1:2
                [signal(s),t]=ReadGripValue('SerialMBB',Device.dynamo,s);
            end
            Dynamo.forceSignal = [Dynamo.forceSignal , signal];
            Dynamo.baselineSignal = [Dynamo.baselineSignal , Events.dynamoBaseline ];
            Dynamo.time = [Dynamo.time , t];
            Dynamo.nTrial = [Dynamo.nTrial , Counters.iTrial];
            Dynamo.iti = [Dynamo.iti , 0];
            Dynamo.condition = [Dynamo.condition , 1];
            Dynamo.response = [Dynamo.response , 0];
            force = signal - Events.dynamoBaseline ;
            if force(1)>=Parameters.criterions.responseThreshold*Parameters.criterions.dynamoCalib...
            || force(2)>=Parameters.criterions.responseThreshold*Parameters.criterions.dynamoCalib
                Events.response = 1;
                Behaviors.error(Counters.iTrial) = 'early';
                Behaviors.participation(Counters.iTrial) = 1;
                Behaviors.responseTime(Counters.iTrial) = GetSecs - Events.onsetCondition ;
            end
        end
        
    end
    
    % get response
    if Events.response == 0 && Design.isGo(Counters.iTrial)==1
        set_periodLamp(Interface,'response');
        Events.responseFB = 0;
        Events.waitResponse = 1;
        Events.onsetResponse = GetSecs;
        while Events.waitResponse % checking loop

            % display response onset signal
            Screen('FillRect', Parameters.display.window, Stimuli.greenspotCol, Stimuli.greenspotRect);

            % monitor keyboard response
            [keyisdown, secs, keycode] = KbCheck;
            if keyisdown==1 && keycode(Device.key.space)==1
                Events.response = 1;
                Behaviors.error(Counters.iTrial) = 'correct';
                Behaviors.participation(Counters.iTrial) = 1;
                Behaviors.correctResponse(Counters.iTrial) = 1;
                Behaviors.responseTime(Counters.iTrial) = GetSecs - Events.onsetResponse ;
            elseif keyisdown==1 && keycode(Device.key.escape)==1
                Events.stopTask = 1;
            elseif keyisdown==1 && keycode(Device.key.return)==1
                Events.pauseTask = 1;
            end
            
            % monitor dynamo
            if dynamoFlag
                Events.dynamoBaseline = mean(Dynamo.forceSignal(:,Dynamo.nTrial==Counters.iTrial & Dynamo.iti==1),2);
                signal = nan(2,1);
                for s = 1:2
                    [signal(s),t]=ReadGripValue('SerialMBB',Device.dynamo,s);
                end
                Dynamo.forceSignal = [Dynamo.forceSignal , signal];
                Dynamo.baselineSignal = [Dynamo.baselineSignal , Events.dynamoBaseline ];
                Dynamo.time = [Dynamo.time , t];
                Dynamo.nTrial = [Dynamo.nTrial , Counters.iTrial];
                Dynamo.iti = [Dynamo.iti , 0];
                Dynamo.condition = [Dynamo.condition , 0];
                Dynamo.response = [Dynamo.response , 1];
                force = signal - Events.dynamoBaseline ;
                if (force(1)>=Parameters.criterions.responseThreshold*Parameters.criterions.dynamoCalib...
                || force(2)>=Parameters.criterions.responseThreshold*Parameters.criterions.dynamoCalib)...
                && Events.response == 0
                    Events.response = 1;
                    Behaviors.error(Counters.iTrial) = 'correct';
                    Behaviors.participation(Counters.iTrial) = 1;
                    Behaviors.correctResponse(Counters.iTrial) = 1;
                    Behaviors.responseTime(Counters.iTrial) = GetSecs - Events.onsetResponse ;
                end
            end

            % feedback to response: visual & audio
            if Behaviors.correctResponse(Counters.iTrial)==1
               % audio-visual feedback
               Screen('FillRect', Parameters.display.window, Stimuli.bluespotCol, Stimuli.bluespotRect);
               set_periodLamp(Interface,'reward');
               if Events.responseFB == 0
                   PsychPortAudio('FillBuffer', Device.audio.handle, sound2);
                   PsychPortAudio('Start', Device.audio.handle, 1, 0, 0);
                   Events.responseFB = 1;
               end
            end
            Screen(Parameters.display.window,'Flip');
           
            % update monitoring condition
            %%% fixed monitoring duration 
%             Events.waitResponse = ( GetSecs < Events.onsetResponse + Design.responseDuration(Counters.iTrial) ) ; 
            %%% monitoring until response
            Events.waitResponse = ( (GetSecs < Events.onsetResponse + Design.responseDuration(Counters.iTrial)) & Events.response == 0 ) ; 

        end
    end
    
    % feedback after response phase: reward
    WaitSecs(Parameters.durations.responseRewardDelay);
    if rewardControllerFlag
        if Behaviors.correctResponse(Counters.iTrial)==1
           Behaviors.rewardOutcome(Counters.iTrial) = 1;
           outp(888,00000001);
           WaitSecs(Device.reward.timeUnit*Behaviors.rewardOutcome(Counters.iTrial));
           outp(888,0);
           Counters.vReward = Counters.vReward + Device.reward.time2volume*Device.reward.timeUnit*Behaviors.rewardOutcome(Counters.iTrial) ; 
        end
    end

    % compute performance criterions
    if Behaviors.participation(Counters.iTrial)==1 || Design.isGo(Counters.iTrial)==0
        Design.nRepetition(Counters.iTrial+1) = 0;
        Counters.iCondition = Counters.iCondition+1;
        force = Dynamo.forceSignal(:,Dynamo.nTrial==Counters.iTrial)-Dynamo.baselineSignal(:,Dynamo.nTrial==Counters.iTrial);
        Behaviors.peakForce(Counters.iTrial) = max(max(force));
    else
        Design.nRepetition(Counters.iTrial+1) = Design.nRepetition(Counters.iTrial)+1;
    end
    
    % update conditions
        % trial-by-trial
        Design.nCondition(Counters.iTrial+1) = Counters.iCondition;
        penalty = Parameters.durations.repetitionPenalty*Design.nRepetition(Counters.iTrial+1);
        Design.itiDuration(Counters.iTrial+1) = Conditions.itiDuration(Counters.iCondition) + penalty;
        Design.conditionDuration(Counters.iTrial+1) = Conditions.conditionDuration(Counters.iCondition);
        Design.responseDuration(Counters.iTrial+1) = Conditions.responseDuration(Counters.iCondition);
        Design.isGo(Counters.iTrial+1) = Conditions.isGo(Counters.iCondition);
        % initial calibration update (to do)
        ncorrect = sum(Behaviors.correctResponse(1:Counters.iTrial));
        if ncorrect >= Parameters.criterions.ntrial2updateCalib
            Parameters.criterions.dynamoCalib = Parameters.criterions.updateCalibModulation*nanmax(Behaviors.peakForce(1:Counters.iTrial));
        end
        
    % display summary
    Counters.nParticipation = sum(Behaviors.participation);
    Counters.nCorrect = sum(Behaviors.correctResponse);
    clck2 = clock;
    clck2 = clck2 - clck;
    duration = [ num2str(clck2(4)) 'h' num2str(clck2(5)) 'min' ];
    Metadata.duration = duration;
    Interface = set_properties(Interface,Metadata,[],Counters,Events,...
                               Behaviors,Conditions,Design,[],Dynamo);
    Interface = display_dynamo(Interface);
    drawnow;
    [Interface,Metadata,Parameters,Counters,Events,...
                  Behaviors,Conditions,Design,Device,Dynamo] = update_properties(Interface);

   % check pause criterion
   if Events.pauseTask==1
       Screen(Parameters.display.window,'Flip');
       disp('PAUSE');
       disp('tap "dbcont" to resume ');
       exit=0;
       while exit==0
%             pause('off');
            keyboard;
            [keyisdown, secs, keycode] = KbCheck;
            if keycode(Device.key.return)==1
                exit=1;
            elseif keycode(Device.key.escape)==1
                exit=1;
                Events.stopTask = 1;
            end
       end
       Events.stopTask = Interface.Events.stopTask;
       Events.pauseTask=0;
       disp('RESUME');
   end
    
    % check stop criterions
        % maximal number of repetition
%         stopCondition = Design.nRepetition(Counters.iTrial+1) > Parameters.criterions.maxRepetition ;
        % maximal number of trials 
        stopCondition = Design.nTrial(Counters.iTrial) >= Parameters.criterions.maxNumberOfTrial ;
        % conditional interruption
        if stopCondition 
            Screen(Parameters.display.window,'Flip');
            WaitSecs(1);
            disp('check monkey willingness to consume reward');
            disp('press ''SPACE'' to continue ');
            disp('press ''ESCAPE'' to stop');
            response=0;
            while response==0
                [response, secs, keycode] = KbCheck;
                if response==1 && keycode(Device.key.escape)==1
                    Events.stopTask = 1;
                end
            end
        end        
    
end

%% Close devices
% ---------------------------------------------
Screen(Parameters.display.window,'Flip');
Screen('CloseAll');
if dynamoFlag
    CloseGripDevice('SerialMBB',Device.dynamo);
end

%% Save data & metadata
% ---------------------------------------------

cd(Metadata.subdir);

nmax = Counters.iTrial;
Conditions = Conditions(1:nmax,:);
Design = Design(1:nmax,:);
Behaviors = Behaviors(1:nmax,:);
save(Metadata.dataFileName,'Metadata','Parameters','Conditions','Design','Behaviors','Dynamo');
disp('Data saving complete.');

cd(Metadata.taskdir);


%% Summarize task data
% ---------------------------------------------

% Summary = struct;
% Summary.n_trial = ;
% Summary.n_condition = ;
% Summary.n_participation = ;
% Summary.n_correct = ;
% Summary.participation_rate = ;
% Summary.correct_rate = ;
% disp(Summary);

%% Stop all processes
%---------------------------------------------
disp('Press any key to terminate all processes.');
WaitSecs(1);
pause;
delete(Interface);
clc;
clear all;
close all;

end