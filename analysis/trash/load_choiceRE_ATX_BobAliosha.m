%% load_choiceRE_ATX_BobAliosha


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
   sessionList.bob = {'2_22_2016','2_23_2016','2_24_2016','2_25_2016','2_26_2016',...
                      '2_29_2016','3_1_2016','3_2_2016','3_3_2016','3_4_2016',...
                      '3_29_2016','3_30_2016','3_31_2016','4_1_2016',...
                      '4_11_2016','4_13_2016','4_14_2016','4_15_2016',...
                      '4_25_2016','4_26_2016','4_27_2016','4_28_2016','4_29_2016',...
                      '5_9_2016','5_10_2016','5_11_2016'};                      
   sessionList.aliosha = {'2_22_2016','2_23_2016','2_24_2016','2_25_2016','2_26_2016',...
                      '2_29_2016','3_1_2016','3_2_2016','3_3_2016','3_4_2016',...
                      '3_29_2016','3_30_2016','3_31_2016','4_1_2016',...
                      '4_11_2016','4_13_2016','4_14_2016','4_15_2016',...
                      '4_25_2016','4_26_2016','4_27_2016','4_28_2016','4_29_2016',...
                      '5_3_2016','5_4_2016','5_5_2016'};  
                  
%% Preprocessing
    [bigData,bigResult] = compile_choiceRE(datadir,subjectList,sessionList);
    data = bigData.interSubject;
    
    

