function [design,data] = init_design(calib)
%% init_design

design =  struct;

% experimental factors
    % levels
        design.levels.reward = [ 1 , 2 , 3 , 4]; % sec. of open volume == 0.85 ml
        design.levels.effort = [ 0.10 , 0.40 , 0.70 , 1.00]; % maximal force
        design.levels.laterality = [1 , 2 , 3 , 4];
        
        design.calibration = calib;


    % trial number
        nTrialByBlock = 256 ; design.nTrialByBlock=nTrialByBlock;
        nBlock = 10 ; design.nBlock=nBlock;
        design.trialNumber = design.nTrialByBlock*design.nBlock; trialNumber = design.trialNumber;
        
    % design matrix
        design.matrix = mat2dataset([1:design.trialNumber]','VarNames','trialNumber');
        
        % controlled factors
             % orthogonal reward & effort magnitude on side 1 & side 2 
             effortBlock1 = [ repmat(1,1,nTrialByBlock/4) , repmat(2,1,nTrialByBlock/4) , repmat(3,1,nTrialByBlock/4), repmat(4,1,nTrialByBlock/4)]; 
             effortBlock2 = [ repmat([repmat(1,1,nTrialByBlock/16) , repmat(2,1,nTrialByBlock/16) , repmat(3,1,nTrialByBlock/16), repmat(4,1,nTrialByBlock/16)],1,4)]; 

             currentReward = 1; % it could be 2 or 3 or 4
             rewardBlock1 = [ repmat([repmat(1,1,nTrialByBlock/64) , repmat(2,1,nTrialByBlock/64) , repmat(3,1,nTrialByBlock/64), repmat(4,1,nTrialByBlock/64)],1,16)];
             rewardBlock2 = [ repmat([1,2,3,4],1,64)];  % similar rewards for both options

             sideBlock = [repmat(1,1,nTrialByBlock/2) , repmat(2,1,nTrialByBlock/2)];
             
             % random permuted sequence
             e1 = []; r1 = [];e2 = []; r2 = []; side = [];
            for iBlock = 1:(trialNumber/nTrialByBlock)
                index = randperm(nTrialByBlock); 

                e1 = [e1, effortBlock1(index) ];
                r1 = [r1, rewardBlock1(index)];
                e2 = [e2, effortBlock2(index) ];
                r2 = [r2, rewardBlock2(index)];
                side = [side, sideBlock(index)];
            end
            
            design.matrix.ordinalEffortLeft(:,1) = e1';
            design.matrix.ordinalEffortRight(:,1) = e2';
            design.matrix.ordinalRewardLeft(:,1) = r1';
            design.matrix.ordinalRewardRight(:,1) = r2';
            design.matrix.side(:,1) = side';
            
            design.matrix.calibration = [ repmat(calib(1),trialNumber,1) , repmat(calib(2),trialNumber,1)];
            
    % criterions
        design.criterions.maxOmission = 30; % maximal n° of successive omissions authorized
        design.criterions.maxFailure = 0; % maximal n° of successive failed trials which are repeated
        design.criterions.anticipationForce  = 20; % (%) of force level exerted before the time of choice
        design.criterions.responsetimeForce  = 20; % (%) of force level exerted after starting onsert
        design.criterions.precisionForce  = 0.20;
        design.criterions.releaseForce = 0;
        
    % timing parameters (sec)
        design.timing.starting_duration = 3; % time before the experiment start running
        design.timing.min_redSpotDuration=0.5;
        design.timing.max_redSpotDuration=2.5; 
        design.timing.anticipationduration=0.1;
        design.timing.responseduration=2.5;
        design.timing.min_intertrialduration=0.5;
        design.timing.max_intertrialduration=1.5;
        design.timing.failure_delay = 4;
        
        % timing factors
            design.matrix.redSpotDuration=rand(trialNumber,1)*(design.timing.max_redSpotDuration-design.timing.min_redSpotDuration)+design.timing.min_redSpotDuration;
            design.matrix.intertrialduration=rand(trialNumber,1)*(design.timing.max_intertrialduration-design.timing.min_intertrialduration)+design.timing.min_intertrialduration;

% initialize counter  variables
    design.count.totalGain=0; %total gain counter
    design.count.task_exit=0;  %exit flag to end task
    design.count.err=0;   %error counter
    design.count.omission = 0; % omission counter
    design.count.withdrawal = 0; % withdrawal counter
    design.count.ntrial = 0; %starting trial
    
    
% pre-allocate data matrix
    data.responsetime=nan(trialNumber,1);
    data.reactiontime=nan(trialNumber,1);
    data.feedbacktime=nan(trialNumber,1);
    data.leveltime   =nan(trialNumber,1);

    data.force=nan(trialNumber,1);
    data.sumforce=nan(trialNumber,1);
    data.forceduration=zeros(trialNumber,2);
    data.baseline=nan(trialNumber,1);
    data.perf=nan(trialNumber,1);
    data.sumperf=nan(trialNumber,1);
    data.gain=zeros(trialNumber,1);
    data.choice=zeros(trialNumber,1);
    % choice=side;
    data.error=nan(trialNumber,1);
    data.repetition=zeros(trialNumber,1);

    data.gripdata=nan(trialNumber,100);
    data.gripbaseline=nan(trialNumber,100);
    data.gripanticipation=nan(trialNumber,100);
    
    data = struct2dataset(data);
    

end