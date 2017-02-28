%% load_choiceRE_CLO_BobAliosha


%% Directory Configuration
%__________________________________________________________________________

clc;
clear all;
close all;

[root,vbadir,analysisdir,datadir,resultdir] = setPathChoiceRE;

%% Analysis
% today
today = {'12_2015'};

% compilation
   subjectList = {'bob','aliosha'};
   sessionList.bob = {'5_9_2016','5_10_2016','5_11_2016',...                                        % placebo
                      '5_30_2016','5_31_2016','6_1_2016','6_3_2016',...                             % clonidine
                      '6_7_2016','6_9_2016','6_10_2016',...                                         % placebo
                      '6_13_2016','6_14_2016','6_15_2016','6_16_2016','6_17_2016',...               % clonidine
                      '6_21_2016','6_22_2016','6_24_2016',...                                       % placebo
                      '6_28_2016','6_29_2016','6_30_2016'};                                         % clonidine
                  
   sessionList.aliosha = {'5_3_2016','5_4_2016','5_5_2016','5_6_2016',...                           % placebo
                      '5_30_2016','5_31_2016','6_1_2016','6_3_2016',...                             % clonidine
                      '6_7_2016','6_9_2016',...                                                     % placebo
                      '6_13_2016','6_14_2016','6_15_2016','6_16_2016','6_17_2016',...               % clonidine
                      '6_21_2016','6_24_2016',...                                                   % placebo
                      '6_28_2016','6_29_2016','6_30_2016'};                                         % clonidine
                  
%% Preprocessing
    [bigData,bigResult] = compile_choiceRE(datadir,subjectList,sessionList);
    data = bigData.interSubject;
    
    


