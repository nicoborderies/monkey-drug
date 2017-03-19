%% Clonidine_analysis
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


%% CLONIDINE-PLACEBO univariate comparisons
% -------------------------------
trt2comp = [1 2];

% between-treatments comparisons of choices (HR/HE/ correct choice)
% choice HE
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'accuracy_E','choice = HE (%)',[1 2],'RFX');
% choice HR
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'accuracy_R','choice = HR (%)',[1 2],'RFX');
                            
% choice HRHE
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'choice_HRHE','choice = HRHE (%)',[1 2],'RFX');
        
% between-treatments comparisons of participation
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'participation_rate','participation(%)',[1 2],'RFX');
                            
% between-treatments comparisons of force peaks
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'peak_force','peak force(%)',[1 2],'RFX');
                            
% between-treatments comparisons of total reward
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'session_reward','total reward(ml)',[1 2],'RFX');
    
% between-treatments comparisons of response time
[stat,model]  = barcomp_monkeydrug(sessiondata,...
                                'responsetime','response time (s)',[1 2],'RFX');

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
                                'logit','selectSubject',[1 2],'RFX','MAP');
sessiondata = [sessiondata , array2table(stat.parameters,'VariableNames',stat.paramNames)];


% logistic regression of choice + MFX comparison of treatments
formula = 'choice ~ 1 + dR + dE + dR:ntrial + dE:ntrial + choice_t';
effectName = {'b_right','dR','dE','choice_t','dRxntrial','dExntrial'};
[data,stat] =  do_glm_monkeydrug(design,data,'choice',...
                                formula,effectName,...
                                'logit','selectSubject & selectionChosen',[1 2],'RFX','MAP');
sessiondata = [sessiondata , array2table(stat.parameters,'VariableNames',stat.paramNames)];

                            
% manova on regression coefficients by treatment
[data,stat] =  do_glm_monkeydrug(design,data,'choice',...
                                formula,effectName,...
                                'logit','selectSubject & selectionChosen',[1 2],'MANOVA','MLH');

% logistic regression of choice + MFX comparison of treatments (with interaction term dR*dE)
formula = 'choice ~ 1 + dR*dE + dR:ntrial + dE:ntrial + choice_t';
effectName = {'b_right','dR','dE','choice_t','dRxdE','dRxntrial','dExntrial'};
[data,stat] =  do_glm_monkeydrug(design,data,'choice',...
                                formula,effectName,...
                                'logit','selectSubject & selectionChosen',[1 2],'RFX','MAP');

% logistic regression of choice + MFX comparison of treatments (with interaction term sumR*dE)
formula = 'choice ~ 1 + dR + dE + sumR:dE + dR:ntrial + dE:ntrial + choice_t';
effectName = {'b_right','dR','dE','choice_t','sumRxdE','dRxntrial','dExntrial'};
[data,stat] =  do_glm_monkeydrug(design,data,'choice',...
                                formula,effectName,...
                                'logit','selectSubject & selectionChosen',[1 2],'RFX','MAP');
                           
% regression of force + MFX comparison of treatments
formula = 'force ~ 1 + side + side:echoice + rchoice*echoice + ntrial';
effectName = {'k_0','rchoice','echoice','side','ntrial','rchoice*echoice ','echoice*side'};
[data,stat] =  do_glm_monkeydrug(design,data,'force',...
                                formula,effectName,...
                                'identity','selectSubject & selectionChosen',[1 2],'RFX','MLH');
                                                      
%% CLONIDINE-PLACEBO comparisons of bivariate associations
% -----------------------------------------------

% force-velocity relationship comparison
stat = linecomp_monkeydrug(design,data,'peak_force','peak force(%)','peak_velocity','peak force(%)','peak_force_bin','glm',...
                            10,0,[0 1],[],'selectSubject & selectionChosen',[1 2]);  
                            
% choiceHE-ntrial relationship comparison
stat = linecomp_monkeydrug(design,data,'nt','session progression (ntrial)','choiceE','choice = HE (%)','nnt','glm',...
                            10,0,[0 1],[],'selectSubject & selectionChosen',[1 2]);  

% choiceHE-sumR relationship comparison
stat = linecomp_monkeydrug(design,data,'sumR','sumR','choiceE','choice = HE (%)','sumR','glm',...
                            7,0,[0 1],[],'selectSubject & selectionChosen',[1 2]); 
                            
%% Between-session metric correlations
% -----------------------------------------------

% succes rate - choiceHE correlations
[stat] = scattercomp_monkeydrug(design,sessiondata,'success_force','correct force (%)','choice_E','choice = HE (%)',...
                                                      'betweenConditionConfound',...
                                                      [0 1],[0 1],[1 2]);  
                            
% effort sensitivity - participation bias correlations           
[stat] = scattercomp_monkeydrug(design,sessiondata,'dE','dE','beta_0','beta_0',...
                                                      'withinConditionConfound',...
                                                      [],[],[1 2]);       
                            
% B0 participation tendency - Echosen weight correlations
[stat] = scattercomp_monkeydrug(design,sessiondata,'dE','dE','k_0','k_0',...
                                                      'withinConditionConfound',...
                                                      [],[],trt2comp);  
                            
                            
                            
                            
                            
                            
                            