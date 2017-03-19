%% Ethogram_analysis
%   this script execute the following analysis:
%

% Requirements:
%   Script: do_choiceRE_NA_ALL_analysis
%   Subfunctions: 
%   Data-files: monkeydrug_NA_dataset
%   Matlab-version:
%   Matlab-toolbox: 
%
% See also: 
%
% Author: Nicolas Borderies
% email address: nico.borderies@gmail.com 
% March 2017; Last revision: 

%%

    % variables
    trt2comp = [2 3];
    ytext = 'number';
    vartext = 'diversity';
    xtext = vartext;
    aggregate = 0;
    nbehavior = 12;
    nbin = [ nbehavior, ntrt , maxSess ];
    if aggregate 
        nbin = [ ntrt , maxSess ]; 
        nbehavior = 1;
    end
    nbin2 = [ nbehavior, ntrt , nsub ];
    behaviorList = getlevels(ethogram.behavior);
    xlist =    {'AFFILIATIVE','AGGRESSIVE','ALERT','AVOIDANCE','FOOD CONSUMPTION','FOOD SEARCH',...
                'OBJECT MANIPULATION','STEREOTYPING','RESTING','SELF INJURING','SELF PROTECTING','SUBMISSION'};

    Y = nan(nbin2);
    Z = nan(nbin2);
    T = [];
    BEHAVIORFREQ = [];
    DIVERSITY = [];
    STABILITY = [];
    POSITIVITY = [];
    ACTIVITY = [];
    SUB = [];
    SESSION = [];
    
        
    for isub = 1:nsub
        
        % select
        sub = subjectList{isub};
        selectSubject = ismember(ethogram.subject, sub);  
        selection = selectSubject   ;

        % variables
        ysub = nan(nbin);
        session = ethogram.sessionNumber(selection);
        treatment = ethogram.treatment(selection);
        behavior = ethogram.behavior(selection);
        time = ethogram.time(selection);
        valence = double( behavior=='AFFILIATIVE' |...
                    behavior=='FOOD_CONSUMPTION' |...
                    behavior=='OBJECT_MANIPULATION' |...
                    behavior=='SELF_PROTECTION' |...
                    behavior=='SUBMISSION' ) ;
        valence(valence==0) = -1;
        valence(behavior=='RESTING') = 0;

        % stats
        [~,subtrt] = ismember(unique(treatment),trtList);
        subsess = unique(session);
        
        for is = subsess'

            % metric
            freq = histcounts(double(behavior(session==(is))),'Normalization','Probability','BinLimits',[1,nbehavior]);
            totaltime = tapply(time(session==(is)),{behavior(session==(is))},@nansum);
            diversity = numel(unique(behavior(session==(is))));
            stability = nanmean(time(session==(is)));
            positivity  = nanmean(time( session==(is) & valence==1 ));
            activity = mean(~(behavior(session==(is))=='RESTING'));
            
            ib = double(unique(behavior(session==(is))));
            trtname = unique(treatment(session==(is)));
            itrt = double(trtname);
%             Y(ib,itrt,is) = totaltime;  
            if aggregate 
                eval(['vary = ' vartext ';']);
                ysub(itrt,is) = vary;
            else
            ysub(:,itrt,is) = freq;     
            end     
                BEHAVIORFREQ = [ BEHAVIORFREQ ; freq];
                DIVERSITY = [ DIVERSITY ; diversity ];
                STABILITY = [ STABILITY ; stability ];
                POSITIVITY = [ POSITIVITY ; positivity ];
                ACTIVITY = [ ACTIVITY ; activity ];
                T = [ T ; trtname ];
                SUB = [ SUB ; nominal(sub) ];
        end
        SESSION = [SESSION ; subsess];


        % averaging
         dim=3;
         if aggregate; dim=2; end
         y = nanmean(ysub,dim);     Y(:,subtrt,isub) = y;
         z = sem(ysub,dim);         Z(:,subtrt,isub) = z;
        
%         % display
%         if isub==1; fig = figure; set(fig,'Name','ethogram'); end
%         subplot(1,nsub,isub);hold on ; clear h ;
%         for it=trt2comp
%             x = [1:nbehavior];
%             if aggregate; x = [1]; end
%             xx = x + (it-2)*0.25;
%             if aggregate
%                 [ h(it)  ] = barplot( xx ,y(it),z(it), col{it} );
%             else
%                 [ h(it)  ] = barplot( xx ,y(:,it),z(:,it), col{it} );
%             end
%             h(it).BarWidth = 0.15;
%             
%         end
        
        % legending
%         if isa(h,'matlab.graphics.Graphics') && isub==nsub
%             legend([h(1) h(2) h(3)],trtList);
%         end
%         ax = gca; 
%         ax.TickLength = [0 0];
%         if ~aggregate
%             ax.XTick = x;
%             ax.XTickLabel = cellstr(xlist);
%         end
%         ax.XTickLabelRotation = 45; 
%         ylabel(ytext); 
%         if aggregate; xlabel(xtext);  end
%         t = title(['monkey ' sub(1) ]); 
%         
    end

% stores in sessiondata
[~,sublist,sesslist] = findgroups(SUB,SESSION);
sessiondata.behaviorFreq = nan(height(sessiondata),12);
sessiondata.diversity = nan(height(sessiondata),1);
sessiondata.stability = nan(height(sessiondata),1);
sessiondata.positivity = nan(height(sessiondata),1);
sessiondata.activity = nan(height(sessiondata),1);
for i = 1:numel(SESSION)
    index = find(sessiondata.subject==sublist(i) & sessiondata.session_number==sesslist(i));
    sessiondata{index,'behaviorFreq'} = BEHAVIORFREQ(i,:);
    sessiondata{index,'diversity'} = DIVERSITY(i);
    sessiondata{index,'stability'} = STABILITY(i);
    sessiondata{index,'positivity'} = POSITIVITY(i);
    sessiondata{index,'activity'} = ACTIVITY(i);
end

