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
%   while ~STOP_CRITERION_1
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
optionList = {'fullscreen'};
subList = {'aliosha','bob','dracula','esmeralda'};
fullscreen=1;
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

% Parameters.display
% Parameters.device.key



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
    experimenterName=input('session number ?');
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

% Set Directory Configuration
%-----------------------------------------------
% define script directory
Metadata.taskdir=pwd;
cd ..
uproot=pwd;

% define group & subject result directory
resultdir=[uproot '\resultats'];
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

% Set Keyboard Configuration
%-----------------------------------------------
Device = struct;
[~,~,keycode] = KbCheck;
DisableKeysForKbCheck(find(keycode==1));

KbName('UnifyKeyNames');
Device.key.left = KbName('LeftArrow');
Device.key.right = KbName('RightArrow');
Device.key.space =KbName('Space') ;
Device.key.escape =KbName('ESCAPE') ;

% Reset random number generator
%-----------------------------------------------
rand('state',sum(100*clock));

% Set device configurations
%-----------------------------------------------

% Reward controller (Parallel Port Com)
config_io;
outp(888,0); % put all pins to zero

% Dynamometer 

% Video Camera

% Motion Sensor

% Microphone Sensor


% Load Stimuli
%-----------------------------------------------
Stimuli = struct;

cd(Metadata.imgdir); % enter img dir
% images
% background=Screen('MakeTexture',window,imread('bg_100.bmp'));
% rect_bg=CenterRectOnPoint(Screen('Rect',background),x,y);

% spots
Stimuli.xspot = Parameters.display.xcenter;
Stimuli.yspot = Parameters.display.ycenter;
Stimuli.wspot = Parameters.display.W*0.05;
Stimuli.hspot = Parameters.display.W*0.05;
Stimuli.backgroundImg = ones(1600,900,3)*0.5;
Stimuli.redspotImg = cat(3,ones(80,80),zeros(80,80),zeros(80,80)) ;
Stimuli.greenspotImg = cat(3,zeros(80,80),ones(80,80),zeros(80,80)) ;
Stimuli.bluespotImg = cat(3,zeros(80,80),zeros(80,80),ones(80,80)) ;

% option positions

% force visual feedback  

cd(Metadata.taskdir); % exit

% Make experimental design
%-----------------------------------------------

% durations 
Parameters.durations.itiRange = [0.5 1.5];
Parameters.durations.conditionRange = [1 2];
Parameters.durations.responsenRange = [ 3 ];

% conditions
%%% sample size
Parameters.conditions.nMinTrial = 100;
Parameters.conditions.nMaxTrial = Parameters.conditions.nMinTrial*20;

%%% variables levels

%%% variables creation
conditionNames = {'nCondition','itiDuration','conditionDuration','responseDuration'};
designNames = [ {'nTrial','nRepetition'} , conditionNames];

Conditions = array2table(nan(Parameters.conditions.nMaxTrial,numel(conditionNames)),'VariableNames',conditionNames);
Design = array2table(nan(Parameters.conditions.nMaxTrial,numel(designNames)),'VariableNames',designNames);

%%% randomization

% experimental criterions
% Parameters.criterions.responseThreshold = ;
Parameters.criterions.maxRepetition = 20;


% Data preparation
%-----------------------------------------------
behaviorNames = {'rewardOutcome','error','participation','correctResponse','peakForce','responseTime'};
Behaviors = array2table(nan(Parameters.conditions.nMaxTrial,numel(behaviorNames)),'VariableNames',behaviorNames);

Dynamo = struct('forceSignal',[],'time',[],'baselineSignal',[],'nTrial',[],...
                'iti',[],'condition',[],'response',[]);

% Initialization
%-----------------------------------------------
Counters.iCondition = 0;
Counters.iTrial = 0;

Events.startTask = 0;
Events.stopTask = 0;
% Events.correctTrial = 0;
% Events.omission = 0;

%% Check Start Criterion
% ----------------------------------------------

while Events.startTask == 0
    % check keyboard response
    disp('Press any key to start the task\n');
    WaitSecs(1);
    KbWait;
    Events.startTask = 1;
end
WaitSecs(1);


%% Run task 
% ---------------------------------------------

while Events.stopTask == 0
    
    Counters.iTrial = Counters.iTrial+1;
    
end



end