%% Atomoxetine_analysis
%

% Requirements:
%   Script: do_choiceRE_NA_ALL_analysis
%   Subfunctions: 
%   Data-files: monkeydrug_NA_dataset
%   Matlab-version:
%   Matlab-toolbox: 
%
% See also: 
%
% Author: Nicolas Borderies
% email address: nico.borderies@gmail.com 
% February 2017; Last revision: 


%% CLONIDINE-PLACEBO comparisons
% -------------------------------
trt2comp = [2 3];

% between-treatments comparisons of choices (HR/HE/ correct choice)
% choice HE
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'accuracy_E','choice = HE (%)',trt2comp,'RFX');
% choice HR
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'accuracy_R','choice = HR (%)',trt2comp,'RFX');
                            
% choice HRHE
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'choice_HRHE','choice = HRHE (%)',trt2comp,'RFX');
        
% between-treatments comparisons of participation
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'participation_rate','participation(%)',trt2comp,'RFX');
                            
% between-treatments comparisons of force peaks
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'peak_force','peak force(%)',trt2comp,'RFX');
                            
% between-treatments comparisons of total reward
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'session_reward','total reward(ml)',trt2comp,'RFX');
    
% between-treatments comparisons of response time
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'responsetime','response time (s)',trt2comp,'RFX');

% between-treatments comparisons of session duration
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'session_time','session duration (ntrial)',trt2comp,'RFX');


%% CLONIDINE-PLACEBO multivariate comparisons
% -------------------------------

% ethogram behavior frequencies
varnames = {'AFFILIATIVE','AGGRESSIVE','ALERT','AVOIDANCE','FOOD CONSUMPTION','FOOD SEARCH',...
                'OBJECT MANIPULATION','STEREOTYPING','RESTING','SELF INJURING','SELF PROTECTING','SUBMISSION'};
[sessiondata,stat] =  multcomp_monkeydrug(design,sessiondata,'behaviorFreq',...
                                varnames,...
                                trt2comp,'RFX');

% ethogram aggregate metrics
varnames = {'diversity','stability','positivity','activity'};
[sessiondata,stat] =  multcomp_monkeydrug(design,sessiondata,varnames,...
                                varnames,...
                                trt2comp,'RFX');

%% CLONIDINE-PLACEBO comparisons of GLM model parameters
% -------------------------------

% logistic regression of participation + MFX comparison of treatments
formula = ' participation ~ 1 + sumR + sumE  + ntrial + part_t ';
effectName = {'beta_0','sumR','sumE','ntrial','part_t'};
[data,stat] =  do_glm_monkeydrug(design,data,'participation',...
                                formula,effectName,...
                                'logit','selectSubject',trt2comp,'RFX','MAP');
sessiondata = [sessiondata , array2table(stat.parameters,'VariableNames',stat.paramNames)];


% logistic regression of choice + MFX comparison of treatments
formula = 'choice ~ 1 + dR + dE + dR:ntrial + dE:ntrial + choice_t';
effectName = {'b_right','dR','dE','choice_t','dRxntrial','dExntrial'};
[data,stat] =  do_glm_monkeydrug(design,data,'choice',...
                                formula,effectName,...
                                'logit','selectSubject & selectionChosen',trt2comp,'RFX','MAP');
sessiondata = [sessiondata , array2table(stat.parameters,'VariableNames',stat.paramNames)];

                           
% regression of force + MFX comparison of treatments
formula = 'force ~ 1 + side + side:echoice + rchoice*echoice + ntrial';
effectName = {'k_0','rchoice','echoice','side','ntrial','rchoicexechoice ','echoicexside'};
[data,stat] =  do_glm_monkeydrug(design,data,'force',...
                                formula,effectName,...
                                'identity','selectSubject & selectionChosen',trt2comp,'RFX','MLH');
sessiondata = [sessiondata , array2table(stat.parameters,'VariableNames',stat.paramNames)];

                                                      
%% CLONIDINE-PLACEBO comparisons of bivariate associations
% -----------------------------------------------

                            
                            
%% Between-session metric correlations
% -----------------------------------------------
                             
% B0 participation tendency - Echosen weight correlations
[stat] = scattercomp_monkeydrug(design,sessiondata,'echoice','echoice','beta_0','beta_0',...
                                                      'withinConditionConfound',...
                                                      [],[],trt2comp);                        
                            
                            
                            
                            
                            