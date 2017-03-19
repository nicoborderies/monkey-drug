%% load_choiceRE_ATX_Aliosha052016


%% Directory Configuration
%__________________________________________________________________________
% 
% clc;
% clear all;
% close all;

[root,vbadir,analysisdir,datadir,resultdir] = setPathChoiceRE;

%% Analysis
% today
today = {'12_2015'};

% compilation
   subjectList = {'aliosha'};                  
   sessionList.aliosha = {'5_3_2016','5_4_2016','5_5_2016','5_6_2016'};  
                  
%% Preprocessing
    [bigData,bigResult] = compile_choiceRE(datadir,subjectList,sessionList);
    data = bigData.interSubject;
    
    
%% Core analysis

%         processChoiceRE_2(dataTable);
%         displayChoiceRE_2(dataTable);


%         export_fig([ resultdir filesep 'choiceRE_' subjectList{iSub} '_' 'participation' '_' today{1}  ],'-png','-m3','-a2',fig.participation);
%         export_fig([ resultdir filesep 'choiceRE_' subjectList{iSub} '_' 'force' '_' today{1}   ],'-png','-m3','-a2',fig.force);
%         export_fig([ resultdir filesep 'choiceRE_' subjectList{iSub} '_' 'choice' '_' today{1}  ],'-png','-m3','-a2',fig.choice);


