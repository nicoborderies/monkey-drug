%% ChoiceRE_v5.1.21 (the one that attract luck)
%

% reset
clear all;

%% Session identification
    subid=input('Identifiant du sujet? ( esmeralda, bob)   :','s');
    experimenter=input('Nom de l''experimentateur? (ex: nicolas, sophie, caroline...)   :','s');
    treatment=input('Traitement? (ex: free, placebo, atomoxetine, l-dopa...) :','s');
    calib(1)=input('calibration de la pince (canal 1 = pince gauche)? (UA)   :');
    calib(2)=input('calibration de la pince (canal 2 = pince droite)? (UA)   :');
    reward_unit_duration = input('unit� de r�compenses (sec.)? :');
    effort_unit_duration = input('unit� d''effort (sec.)? :');

%% Directory initialization
    root='C:\Users\seb.grousset\Desktop\Nicolas\PharmacoSinges\MotivEffort V5';
    resultdir=[root '\data'];


%% Behavioral Tasks
%   1) Calibration Phase
        % automatic calibration
%             try
%                 [calib_data,calib_gripdata,calib_gripbaseline] = ChoiceRE_v121_calibration(root,calib,reward_unit_duration,effort_unit_duration);
%             catch err % enable error reporting & close on going interface
%                 Screen('CloseAll');
%                 Handle=InitializeGrip('SerialMBB',[1 2]);
%                 CloseGripDevice('SerialMBB',Handle);    
%             end
%             sideCalib = {'calibrationLeft','calibrationRight'};
%             for iSide=1:2
%                 try           
% %                     newcalib(iSide) = max(calib_data.(sideCalib{iSide})(calib_data.(sideCalib{iSide})~=max(calib_data.(sideCalib{iSide}))));
%                     newcalib(iSide) = calib_data.(sideCalib{iSide})(end);
% 
%                 catch
%                     newcalib(iSide)  = calib(iSide) ;
%                     fprintf( '%s hasn''t changed\n',sideCalib{iSide});
%                 end
%             end
%             disp(newcalib);
            

%         % manual claibration
                newcalib = calib;
    


     
    %% 2) Test Phase
        try % test for syntax error
           
            % task execution & recording
            switch subid
                case 'esmeralda'
                    [data,gripdata,gripbaseline] = ChoiceRE_v121_task(root,newcalib,reward_unit_duration,effort_unit_duration);
                case 'bob'
                    [data,gripdata,gripbaseline] = ChoiceRE_v121_task(root,newcalib,reward_unit_duration,effort_unit_duration);
                case 'aliosha'
                    [data,gripdata,gripbaseline] = ChoiceRE_v121_task(root,newcalib,reward_unit_duration,effort_unit_duration);
                otherwise
                    warning('invalid subject name %s',subid);
            end
            
            % saving
                % add session information
                     subject = repmat({subid},numel(data.trialNumber),1);
                     experimenter = repmat({experimenter},numel(data.trialNumber),1);
                     treatment = repmat({treatment},numel(data.trialNumber),1);
                     calibLeft = repmat(newcalib(1),numel(data.trialNumber),1);
                     calibRight = repmat(newcalib(2),numel(data.trialNumber),1);
                     rewardUnit = repmat(reward_unit_duration,numel(data.trialNumber),1);
                     effortUnit = repmat(effort_unit_duration,numel(data.trialNumber),1);
                     sessionData = mat2dataset([subject , experimenter , treatment]);
                     sessionData = set(sessionData,'VarNames', {'subject','experimenter','treatment'}); 
                     sessionData = [ sessionData , mat2dataset([ calibLeft , calibRight , rewardUnit , effortUnit  ]) ];
                     sessionData = set(sessionData,'VarNames', {'subject','experimenter','treatment',...
                                                            'calibLeft','calibRight','rewardUnit','effortUnit'}); 
                    data = [ sessionData , data] ;


                % matfile
                     taskClock = clock;
                     taskTime = [ num2str(taskClock(2)) '_' num2str(taskClock(3)) '_' num2str(taskClock(1)) '__' num2str(taskClock(4)) 'h' num2str(taskClock(5)) 'min'  ];
                     resultname=strcat(['choiceRE_',subid,'_',taskTime]);
                     cd(resultdir);
                     save(resultname,'data','gripdata','gripbaseline');
                     save(resultname,'calib_data','-append');
            
            
            
        catch err % enable error reporting & close on going interface
            Screen('CloseAll');
            Handle=InitializeGrip('SerialMBB',[1 2]);
            CloseGripDevice('SerialMBB',Handle);    
        end
        
%% Process Data    
[result,data] = processChoiceRE(data);
save(resultname,'data','result','-append');


%% Display Result
[fig] = displayChoiceRE(result,data);

% reset
% clear all;
    