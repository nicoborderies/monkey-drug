%% updateDesign_choiceRE

% instructions conditionnal to error type
    switch error(ntrial)
        case 0 % correct
            % reset counters
                omission = 0;
                withdrawal = 0;
            repetition(ntrial+1)=0;

        case 1 % omission
            err=err+1;
            omission = omission + 1; 
            repetition(ntrial+1)= 1 + (repetition(ntrial)==2);

        case 2 % anticipation
            err=err+1;
            repetition(ntrial+1)= 1 + (repetition(ntrial)==2);
%             repetition(ntrial+1)=2; % TO COMMENT

        case 3 % laterality error
            err=err+1;
            repetition(ntrial+1)=1;

        case 4 % undershoot failure
            err=err+1;
            withdrawal = withdrawal + 1;
            repetition(ntrial+1)=2;

        case 5 % overshoot failure
            err=err+1;
            withdrawal = withdrawal + 1;
            repetition(ntrial+1)=2;
    end

% stopping criterions
    if withdrawal >= withdrawalCriterion
        repetition(ntrial+1)=0;
    end
    if ntrial == trials(end) && repetition(ntrial+1)==0;
        task_exit = 1;
    end
    
% adaptation of next trial
    if repetition(ntrial+1)~=0 % TO COMMENT
%     if repetition(ntrial+1)==2 % TO COMMENT
        trials(end+1) = trials(end)+1;
        calib(end+1) = calib(end);
        
        side = [ side, NaN];
        side(ntrial+1:end) = [ side(ntrial) , side(ntrial+1:end-1) ];
        choice = [ choice, NaN];
        choice(ntrial+1:end) = side(ntrial+1:end);
        
        r1 = [ r1, NaN];
        r1(ntrial+1:end) = [ r1(ntrial) , r1(ntrial+1:end-1) ];
        r2 = [ r2, NaN];
        r2(ntrial+1:end) = [ r2(ntrial) , r2(ntrial+1:end-1) ];
        e1 = [ e1, NaN];
        e1(ntrial+1:end) = [ e1(ntrial) , e1(ntrial+1:end-1) ];
        e2 = [ e2, NaN];
        e2(ntrial+1:end) = [ e2(ntrial) , e2(ntrial+1:end-1) ];
        rewards = [ rewards, nan(2,1)];
        rewards(:,ntrial+1:end) = [ [r1(ntrial+1);r2(ntrial+1)] , rewards(:,ntrial+1:end-1) ];
        efforts = [ efforts, nan(2,1)];
        efforts(:,ntrial+1:end) = [ [e1(ntrial+1);e2(ntrial+1)] , efforts(:,ntrial+1:end-1) ];
    end
