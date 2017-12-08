%% analyze_monkeydrug_session

% Reset
clear all;
clc;


%% Parameters definitions

% folders
datadir = 'B:\nicolas.borderies\projets\monkey_pharma_choice\data\monkeydrug-2017';
codedir = 'B:\nicolas.borderies\projets\monkey_pharma_choice\code\analysis\monkeydrug-2017';
addpath(genpath(datadir));
addpath(genpath(codedir));

% session parameters
subjectName = input('subject name? (aliosha,esmeralda)','s');
sessionNumber = input('session number?');

% parameters
option=struct;
option.dispFig = 1;
option.dispTab = 1;


%% Preparation

% load data file
cd(datadir);
subdir = [datadir filesep subjectName ];
cd(subdir);
filetemplate = [ 'monkey_drug*' subjectName '*n' num2str(sessionNumber) '*.mat'];
f = dir(filetemplate);
% - catch exception
if isempty(f) ; error(['File {' filetemplate '} not found in the current folder!']);end
load(f.name);

%% Analysis


stat = struct; % store result into a struct
f = cell(0); % store figure into a cell array

% data selection
subset = (Design.phase=='test');
subset_2 = (subset & Design.nRepetition==0);

% data curation
Behaviors.correctResponse(Behaviors.participation==0) = NaN ; 
Behaviors.peakForce(Behaviors.error ~= 'correct' &  Behaviors.error ~= 'undershoot' ) = NaN ; 

% 1 - general descriptive stats
% --- computation
stat.descriptive.n_trial = numel(Behaviors.participation);
stat.descriptive.n_participation = sum(Behaviors.participation==1);
stat.descriptive.rate_participation = sum(Behaviors.participation==1)/stat.descriptive.n_trial;
stat.descriptive.n_response = sum(Behaviors.response==1);
stat.descriptive.n_correct = sum(Behaviors.correctResponse==1);
stat.descriptive.rate_correct = sum(Behaviors.correctResponse==1)/stat.descriptive.n_participation;
stat.descriptive.mean_force = nanmean(Behaviors.peakForce);
stat.descriptive.mean_responseTime = nanmean(Behaviors.responseTime);
if option.dispTab; disp(stat.descriptive);end

% 2 - distribution of response types
% --- computation
responseType = unique(Behaviors.error);
n = histcounts(Behaviors.error)';
stat.descriptive.responsedist = table(responseType,n);
disp(stat.descriptive.responsedist);
% --- display
if option.dispTab
    f{1} = figure;
    p = pie(Behaviors.error);
end

% 3 - dependance to experimental factor
% -- compute conditional statistics (all trials)
outputList = {'participation','correctResponse','cumulativeForce'};
for iVar = 1:numel(outputList)
    
    % 3.1 - difficulty manipulation
    [group,difficultyLevel] = findgroups(Design.difficulty(subset));
    varName = outputList{iVar};
    outputVar = Behaviors.(varName)(subset);
    outputStat = splitapply(@nanmean,outputVar,group);
    outputTab = table(difficultyLevel,outputStat,'VariableNames',{'difficultyLevel',varName});
    disp(outputTab);
    
    % 3.2 - incentive manipulation
    [group,incitationLevel] = findgroups(Design.incitation(subset));
    varName = outputList{iVar};
    outputVar = Behaviors.(varName)(subset);
    outputStat = splitapply(@nanmean,outputVar,group);
    outputTab = table(incitationLevel,outputStat,'VariableNames',{'incitationLevel',varName});
    disp(outputTab);

end
% % -- compute conditional statistics (only new trials)
% for iVar = 1:numel(outputList)
%     [group,difficultyLevel] = findgroups(Design.difficulty(subset_2));
%     varName = outputList{iVar};
%     outputVar = Behaviors.(varName)(subset_2);
%     outputStat = splitapply(@nanmean,outputVar,group);
%     outputTab = table(difficultyLevel,outputStat,'VariableNames',{'difficultyLevel',varName});
%     disp(outputTab);
% end
% -- display histograms
f{2} = figure;hold on
for iD = difficultyLevel'
    histogram(Behaviors.cumulativeForce(Design.difficulty==iD),10,'Normalization','probability','FaceColor',[0 0 1-1*(iD/max(difficultyLevel))]);
end
xlabel('cumulative force (au.)');
ylabel('probability density (%)');


% 4 - force dynamical profile
% --- parameters
dt = mean(diff(Dynamo.time));
window = 2;
n_window = round(window/dt);
% --- align
response_onset = find(Dynamo.response>circshift(Dynamo.response,1))';
side = Design.targetSide(Dynamo.nTrial);
indexside = ([1:numel(side)]-1)*2 + side';
force = [Dynamo.forceSignal(indexside) , nan(1,n_window)];
time = [Dynamo.time , nan(1,n_window)];
method = 'peak';
switch method
    case 'epoch'
        % --- 1: alignement on response epoch onset
        force_pulse = cell2mat(arrayfun(@(x) force(x:x+n_window),response_onset,'un',0));
        time_pulse = cell2mat(arrayfun(@(x) time(1,x:x+n_window),response_onset,'un',0));
    case 'peak'
        % --- 2: alignement on response acceleration peak
        onset_finder = @(y) find(diff(y)==nanmax(diff(y)),1,'first') - 50;
        % onset_finder = @(y) find((y)==nanmax((y)),1,'first') - 100;
        onset_pulse = cell2mat(arrayfun(@(x) onset_finder(force(x:x+n_window)),response_onset,'un',0));
        onset_pulse = response_onset + onset_pulse;
        force_pulse = cell2mat(arrayfun(@(x) force(x:x+n_window),onset_pulse,'un',0));
        time_pulse = cell2mat(arrayfun(@(x) time(1,x:x+n_window),onset_pulse,'un',0));
end
% --- correction
force_pulse = force_pulse - force_pulse(:,1); % baseline
force_pulse = force_pulse/nanmax(nanmax(force_pulse)); % normalisation
time_pulse = time_pulse - time_pulse(:,1); % onset-substraction
% --- average
subset_pulse = (Behaviors.error~='early');
subset_pulse = find(Behaviors.error(subset_pulse)=='correct' | Behaviors.error(subset_pulse)=='undershoot');
% -- conditionalize
f{3} = figure;hold on; handle = {};
for iD = difficultyLevel'
    subset_condition = find(Behaviors.correctResponse(subset_pulse)==1 & Design.difficulty(subset_pulse)==iD);
    force_response = nanmean(force_pulse(subset_condition,:),1);
    force_variability = sem(force_pulse(subset_condition,:),1);
    time_response =  nanmean(time_pulse(subset_condition,:),1);
    col_condition = [1*(1 - iD/max(difficultyLevel)) 1*(1 - iD/max(difficultyLevel)) 1];
    [h,hp] = boundedline( time_response , force_response , force_variability , 'alpha','transparency',0.5); 
    handle = [handle,h];
    set(h,'Color',col_condition,'LineWidth',2);set(hp,'FaceColor',col_condition);
end
legend(handle,{'level 1','level 2','level 3','level 4'});
ylabel('force (%fmax)');
xlabel('time (sec.)');





