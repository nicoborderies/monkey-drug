function [root,vbadir,analysisdir,datadir,resultdir] = setPathChoiceRE

% define path
    % External Hard Drive 
    hardDrive = 'H:'; % icm
%     hardDrive = 'F:'; % pc
    icmDisk = 'B:';
    disk = icmDisk;
    root = [ disk '\nicolas.borderies\Projets scientifiques'];
    vbadir = [ disk '\nicolas.borderies\MATLAB\GitHub\VBA-toolbox'];

    projectdir = [root filesep 'monkey_pharma_choice'];
    datadir=[projectdir filesep 'data'];
    analysisdir = [projectdir filesep 'analysis'];
    resultdir = [projectdir filesep 'results'];
    
% set path
%     addpath(genpath(datadir));
%     addpath(genpath(analysisdir));
%     addpath(genpath(resultdir));
%     addpath(genpath(vbadir));