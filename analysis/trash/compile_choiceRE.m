function [bigData,bigResult] = compile_choiceRE(directory,subjectList,sessionList)

%% define default variables
fileName = 'choiceRE_';




%% import

% initialize
    bigData = struct;
    bigResult = struct;
    experimentData= []; 
    
% subject iteration
for iSub = 1:numel(subjectList)
    % init
    if exist('subjectData'); clear subjectData; end
        
        
    % session iteration
    for iSess = 1:numel(sessionList.(subjectList{iSub}))
        
        % load
        sessionName = sessionList.(subjectList{iSub}){iSess};
        file = dir([ directory filesep '*' fileName '*' subjectList{iSub} '*' sessionName '*.mat' ]);
        load([ directory filesep file.name]);

        % adapt format/ add session variables
            bigData.subjects{iSub}.sessions{iSess} = data ;
%             bigResult.subjects{iSub}.sessions{iSess} = result ;
            if class(data)== 'dataset'; data = dataset2table(data) ; end
            data.session = repmat(sessionName,height(data),1);
            data.session = nominal(data.session);
            
             % process force signal
            [ data ] = processForce_choiceRE( data, gripdata );
            
            
            % concatenate & save sessions
            if ~exist('subjectData')
                subjectData = data ;
            else
                [subjectData] = mergetable(subjectData,data);
            end

    end
    
    
    % concatenate & save subjects
        bigData.subjects{iSub}.interSession = subjectData ;
        experimentData = [experimentData ; subjectData];
        
   

end

% save experiment
    bigData.interSubject = experimentData ;





end