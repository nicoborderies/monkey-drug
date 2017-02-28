%% load_choiceRE_Esmeralda


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
   subjectList = {'esmeralda'};
   sessionList.esmeralda = {'6_20_2016'};  
 
%% Preprocessing
    [bigData,~] = compile_choiceRE(datadir,subjectList,sessionList);
    data = bigData.interSubject;
    

