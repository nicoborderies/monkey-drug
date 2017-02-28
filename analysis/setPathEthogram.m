function [root,vbadir,analysisdir,datadir,resultdir] = setPathEthogram

% define path
    % External Hard Drive 
    hardDrive = 'H:'; % icm
%     hardDrive = 'F:'; % pc
    icmDisk = 'C:\Users\nicolas.borderies\Documents';
    disk = icmDisk;
    
    root = [ disk '\Nicolas_Sauvegarde\Projets scientifiques'];
    vbadir = [ disk '\Nicolas_Sauvegarde\MATLAB\GitHub\VBA-toolbox'];

    projectdir = [root filesep 'Projet PNH\Pharmaco-Motiv-actual'];
    datadir=[projectdir filesep 'ethogram\data'];
    analysisdir = [projectdir filesep 'analysis'];
    resultdir = [projectdir filesep 'results'];
    
% set path
    addpath(genpath(datadir));
    addpath(genpath(analysisdir));
    addpath(genpath(resultdir));
    addpath(genpath(vbadir));