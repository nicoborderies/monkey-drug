%% load_ethogram_ATX_BobAliosha


%% Directory Configuration
%__________________________________________________________________________

clc;
clear all;
close all;
%%
[root,vbadir,analysisdir,datadir,resultdir] = setPathEthogram;

%% List sessions

% compilation
   subjectList = {'bob','aliosha'};
   sessionList.bob = {'3_29_2016','3_30_2016','3_31_2016','4_1_2016',...
                       '4_6_2016','4_7_2016','4_8_2016','4_11_2016','4_14_2016','4_15_2016',...
                       '4_25_2016','4_26_2016','4_28_2016'};                      
   sessionList.aliosha = {'3_29_2016','3_30_2016','3_31_2016','4_1_2016',...
                        '4_6_2016','4_7_2016','4_8_2016','4_11_2016','4_14_2016','4_15_2016',...
                       '4_25_2016','4_26_2016','4_28_2016'};  
                  
%% Preprocessing

    id =0;
    vNames={'subject','session','trt','behavior','time'};
    for iv = 1:numel(vNames)
        eval([vNames{iv} '= [];']);
    end
    
    % subject iteration
    for iSub = 1:numel(subjectList)
        % session iteration
        for iSess = 1:numel(sessionList.(subjectList{iSub}))

            % load
            sessionName = sessionList.(subjectList{iSub}){iSess};
            file = dir([ datadir filesep '*ethogram*' subjectList{iSub} '*' sessionName '*.mat' ]);
            load([ datadir filesep file.name]);
            
            % format
                behavCat = fieldnames(ethogram); behavCat = behavCat(1:end-1);
                behavCat = nominal(behavCat);

                indCat = zeros(numel(behavCat),numel(ethogram.time));
                cat = nan(1,numel(ethogram.time));
                for icat=1:numel(behavCat)
                   indCat(icat,find(ethogram.(char(behavCat(icat)))==1)) =  find(ethogram.(char(behavCat(icat)))==1);
                end
                for i=1:numel(ethogram.time)
                    try
                        cat(i) = find(indCat(:,i)>0);
                    catch
                        break
                    end
                end
                cat = cat(~isnan(cat));
                cat2 = [0,cat(1:end-1)];
                sw = (cat~=cat2);
                ind = find(sw==1);
                for ii = 2:numel(ind)
                    id=id+1;
                    behavior = [ behavior ; behavCat(cat(ind(ii-1))) ];
                    time =   [ time ;  ethogram.time(ind(ii))- ethogram.time(ind(ii-1)) ];
                    subject = [ subject ; nominal(subid) ];
                    trt =   [ trt ;  nominal(treatment) ];
                    session = [ session ; nominal(sessionName) ];
                end
            
        end
    end
    
    data = table(subject,session,trt,behavior,time);
    data.Properties.VariableNames = vNames;
    
