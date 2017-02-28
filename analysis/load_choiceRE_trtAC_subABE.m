%% load_choiceRE_trtAC_subABE

clc;
clear all;
close all;

%% Directory Configuration

[root,vbadir,analysisdir,datadir,resultdir] = setPathChoiceRE;
ethogramdir = [root '\monkey_pharma_choice\ethogram\data' ];

%%
% load multisession design
cd(resultdir);
design = readtable('drug_schedule.xlsx','Sheet',4);
% reformat
    % suppress irrelevant quote
    for is = 1:numel(design.session)
        design.session{is} = strrep(design.session{is},'''','');
    end
    
    trtList = {'clonidine','placebo','atomoxetine'};
    design.treatment = nominal(design.treatment);
    design.treatment = reordercats(design.treatment,trtList);
    

%% Compilation

cd(datadir);
fileName = 'choiceRE_';
start = 137;

for is = start:numel(design.session)
    
    % display
    fprintf('loading session %d / %d \n' ,is,numel(design.session));
    
    % load task file
    dataname = [ datadir filesep '*' fileName '*' design.subject{is} '*' design.session{is} '*.mat' ];
    file = dir(dataname);
    if isempty(file);
        warning(['cannot find file : ' dataname ]);
    end
    session = load([ datadir filesep file.name]);
    
    % load ethogram file
    if design.ethogram(is)==1
        dataname = [ ethogramdir filesep '*ethogram*' design.subject{is} '*' design.session{is} '*.mat' ];
        file = dir(dataname);
        if isempty(file);
            warning(['cannot find file : ' dataname ]);
        end
        behavior = load([ ethogramdir filesep file.name]);
    end
    
    
    % format
    if isa(session.data,'dataset'); session.data = dataset2table(session.data) ; end

    % preprocess
        % force
        [ session.data ] = processForce_choiceRE( session.data, session.gripdata );
        
        % ethogram
        if design.ethogram(is)==1
            [session.ethogram] = processEthogram(behavior.ethogram);
        end
    
    % add session information
    session.data.session   = nominal(repmat(design.session{is},height(session.data),1));
    session.data.treatment = repmat(design.treatment(is),height(session.data),1);
    if design.ethogram(is)==1
        session.ethogram.subject = repmat(design.subject(is),height(session.ethogram),1);
        session.ethogram.session = nominal(repmat(design.session(is),height(session.ethogram),1));
        session.ethogram.treatment = repmat(design.treatment(is),height(session.ethogram),1);
    end

    % concatenate
    if exist('data')~=1
        data = session.data ;
    else
        [data] = mergetable(data,session.data);
    end
    if design.ethogram(is)==1
        if exist('ethogram')~=1
            ethogram = session.ethogram ;
        else
            [ethogram] = mergetable(ethogram,session.ethogram);
        end
    end
end

%% Save
    cd([resultdir]);
    date = clock;
    strDate = [ num2str(date(3)) '_' num2str(date(2)) '_' num2str(date(1)) ];
    
    analysisname = 'choiceRE_trtNA_ALL_';
    save([ analysisname '_all_' strDate ],'data','ethogram','design');  
    
    


