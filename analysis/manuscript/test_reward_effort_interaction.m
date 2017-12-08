%% test_reward_effort_interaction


%% data definition
% ------------------

% -- 3-level factor
subject = data.subject;
% -- 2-level factor
sessionNumber = data.sessionNumber;
treatment = data.treatment;
treatment = reordercats(treatment,trtList);
% -- 1-level factor
sum_r = round(data.ordinalRewardLeft + data.ordinalRewardRight);
sum_e = round(data.ordinalEffortLeft + data.ordinalEffortRight);
dr = round(data.ordinalRewardRight - data.ordinalRewardLeft);
de = round(data.ordinalEffortRight - data.ordinalEffortLeft);
abs_dr = abs(dr); abs_dr(dr==0)=NaN;
abs_de = abs(de); abs_de(de==0)=NaN;
r_chosen = round(data.chosenOrdinalReward); r_chosen(isnan(r_chosen)) = 0;
e_chosen = round(data.chosenCardinalEffort,1);  e_chosen(isnan(e_chosen)) = 0;
nt = data.trialNumber;
nnt = data.normalizedTrialNumber;
choice_t = (data.sideChoice-1)*2 - 1 ;
choice_t = circshift(choice_t,1);
choice_t(nt==1) = NaN;
% -- measurements
participation = selectionParticipate; % participation
choice = data.sideChoice-1; % choice
choice(participation==0)=NaN;
choice(e_chosen==0)=NaN;
choiceR = double(sign(choice*2-1)==sign(dr));
choiceR(participation==0)=NaN;
choiceR(e_chosen==0)=NaN;
choiceR(sign(dr)==0)=NaN;
choiceE = double(sign(choice*2-1)==sign(de));
choiceE(participation==0)=NaN;
choiceE(e_chosen==0)=NaN;
choiceE(sign(de)==0)=NaN;
force = data.perf; % force
force(participation==0)=NaN;
force(e_chosen==0)=NaN;
g1 = findgroups(subject,sessionNumber,e_chosen); % force, corrected for required effort
force_effort = splitapply(@nanmean,force,g1);
g2 = findgroups(subject,sessionNumber);
force_sub = splitapply(@nanmean,force,g2);
force_corrected = force - force_effort(g1) + force_sub(g2);
rt = data.responseTime+eps; % response time
rt(participation==0)=NaN;
rt(e_chosen==0)=NaN;
g1 = findgroups(subject,sessionNumber,dr,de); 
choice_cond = splitapply(@nanmean,choice,g1);
choice_error = (choice -  choice_cond(g1)).^2 ;
rt_cond = splitapply(@nanmean,rt,g1);
rt_error = (rt -  rt_cond(g1)).^2 ;
g1 = findgroups(subject,sessionNumber,r_chosen,e_chosen);
force_cond = splitapply(@nanmean,force,g1);
force_error = (force -  force_cond(g1)).^2 ;
force_error_scaled = ((force -  force_cond(g1))./force_cond(g1)).^2 ;

% data preparation
BETA = [];
T = [];
SUB = [];
SESS = [];

%% statistical inference
% ----------------------


for isub = 1:nsub

    % selection
    % - subject
    subname = subjectList{isub};
    select_subject = ismember(subject, subname);  
    [~,subtrt] = ismember(unique(treatment),trtList);
    subsess = unique(session(select_subject));

    for is = subsess'

        % selection
        % - session
        select_session = (session==(is));
        % - trials
        select_trials = selectionChosen;
        % - global
%         selection = select_subject & select_session & select_trials;
        selection = select_subject & select_session;

        % choice
%             % data preparation
%             predictor = [dr,de,choice_t];
%             predictor = predictor(selection,:);
%             predictor = nanzscore(predictor);
%             y = choice(selection);
% 
%             % 1-level fit
%             formula = ' choice ~ 1 + dR + dE + dR:dE + choice_t_1';
%             varnames = {'dR';'dE';'choice_t_1';'choice'};
%             stat = fitglm(predictor,y,formula,'VarNames',varnames,...
%                           'Distribution','binomial','link','logit');
%             coef = stat.Coefficients;
        
        % participation
            % data preparation
            predictor = [sum_r,sum_e,nt];
            predictor = predictor(selection,:);
            predictor = nanzscore(predictor);
            y = participation(selection);

            % 1-level fit
            formula = ' participation ~ 1 + sum_r + sum_e + sum_r:sum_e + nt';
            varnames = {'sum_r';'sum_e';'nt';'participation'};
            stat = fitglm(predictor,y,formula,'VarNames',varnames,...
                          'Distribution','binomial','link','logit');
            coef = stat.Coefficients;
            disp(coef);
        
        
        % store
        BETA = [ BETA ; coef.Estimate' ];
        T = [ T ; unique(treatment(selection)) ];
        SUB = [ SUB ; unique(subject(selection)) ];
        SESS = [ SESS ; unique(session(selection)) ];
        
    end

end

for ibeta = 1:size(BETA,2)
    
    % data preparation
    predictor = table( (T=='clonidine') , (T=='placebo') , (T=='atomoxetine'), SUB,...
                       'VariableNames',{'clonidine','placebo','atomoxetine','subject'});
    y = table(BETA(:,ibeta),'VariableNames',{'beta'});
    % 2-3-level fit
    formula = ' beta ~ clonidine + atomoxetine + (clonidine|subject) + (atomoxetine|subject)';
    stat = fitglme([predictor,y],formula);
    coef = stat.Coefficients;
    disp(coef);

end

