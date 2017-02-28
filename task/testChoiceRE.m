 %% test choiceRE
clear all;
clc;

%% dummy inputs
     subid='999';
    experimenter='999';
    treatment='999';
    calib(1)=200;
    calib(2)=200;
    reward_unit_duration = 0.5;
    
    effort_unit_duration = 0.2;
    
%% Directory initialization
    root='C:\Users\seb.grousset\Desktop\Nicolas\PharmacoSinges\MotivEffort V5';
    resultdir=[root '\data'];
    
%% test
        [data,gripdata,gripbaseline] = ChoiceRE_v124_task(root,calib,reward_unit_duration,effort_unit_duration);
