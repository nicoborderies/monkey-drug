clear all;
close all;



%% Software Configurations 
% try
    
% ------- Identification ------- %
% subid=input('Identifiant du sujet? ( esmeralda, bob)   :','s');
% nsession=input('session du test? (1, 2 ou 3)   :');
% taskClock = clock;
% taskTime = [ num2str(taskClock(3)) '_' num2str(taskClock(2)) '_' num2str(taskClock(1)) '__' num2str(taskClock(4)) 'h' num2str(taskClock(5)) 'min'  ];
% experimenter=input('Nom de l''experimentateur? (ex: nicolas, sophie, caroline...)   :','s');
% treatment=input('Traitement? (ex: free, placebo, atomoxetine, l-dopa...) :','s');
% calib=input('calibration de la pince? (UA)   :');
% trialNumber=input('Nombre d''essais? (multiple de 30) :');
% resultname=strcat(['choiceRE_',subid,'_sess',num2str(nsession),'_' taskTime]);


% ------- Directory Configuration ------- %
% root='C:\Users\seb.grousset\Desktop\Nicolas\PharmacoSinges\MotivEffort V5';
% resultdir=[root '\data'];
% taskdir = [root '\task'];
% picdir= [taskdir '\display'];

MatOs;
