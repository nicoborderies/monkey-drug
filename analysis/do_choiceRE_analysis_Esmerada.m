%% do_choiceRE_analysis_Esmerada

%% local variables
 % data acquisition parameters
        forceSamplingFreq = 25; % Hz.
        rewardUnit2Volume = 0.9; % conversion rate between valve opening time (sec) & water volume (ml)
        drugName = 'atx';
        
        col = { [0 0 1]*0,...
                [0 0 1]*0.25,...
                [0 0 1]*0.75,...
                [0 0 1]*1};

    % data selection
        selectionRepetition = (data.isRepeatedTrial==1) ;
        
        selectionMissed = ( data.errorType==1);
        selectionPremature = ( data.errorType==2);
        selectionSwitch = ( data.errorType==3);
        selectionIncorrect = ( data.errorType==4 | data.errorType==5);

        selectionParticipate = ~selectionMissed ;
        selectionResponseTime =  ~selectionMissed & ~isnan(data.intertrialTime) ;
      
        selectionChosen = ( ~selectionMissed & ~selectionPremature & ~selectionSwitch  );
        selectionCorrect = ( selectionChosen  & ~selectionIncorrect );

    % subjectList 
        subjectList = unique(data.subject);
        sessionList = unique(data.session);
        
        
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        selection = selectSubject & selectionChosen;
        calib = [unique(data.calibLeft(selection)),unique(data.calibRight(selection))];
        data.perf(selection) = data.force(selection)./(calib(data.sideChoice(selection))');
    end
        
%% Force adaptation

fig = figure; set(fig,'Name','force_esmeralda');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
%                 selection = selectSubject & selectionChosen & ~selectionRepetition;
        selection = selectSubject & selectionChosen ;
        
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        reward = data.chosenCardinalReward(selection);
        correct = selectionCorrect(selection);
        novelty = (~selectionRepetition(selection));

        
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(force,{effort},@nanmean);
        err =   tapply(force,{effort},@sem);
        
        for i = 1
            ind = [1:numel(x)];
            plot(x,x,'--','Color',col{1},'LineWidth',1.5);
            h = plotSpread(force,'distributionIdx',effort,'distributionColors','k',...
                        'categoryIdx',correct,'categoryColors',{'r','g'});
            h = errbar( x ,mu(ind), err(ind) ,...
                     'Color',[0 0 0],'LineWidth',2,...
                     'horiz');
        end

        ax = gca;
%         ax.XTick = x;
%         ax.XLim = x;
%         ax.XLim = ax.YLim;
        ob = findobj('-property','MarkerSize');
        for io = 1:numel(ob); ob(io).MarkerSize = 10 ; end
        
%         ax.XTickLabels = {'left','right'};
        xlabel(texlabel('effort required (%fmax)'));
        ylabel('force peak (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end 
    
    
%% Force adaptation

fig = figure; set(fig,'Name','force2_esmeralda');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        % histrogram 
        subplot(1,numel(subjectList),iSub);hold on
%                 selection = selectSubject & selectionChosen & ~selectionRepetition;
        selection = selectSubject & selectionChosen ;
        
        session = data.session(selection);
        treatment = data.treatment(selection);
        force = data.perf(selection);
        effort = data.chosenCardinalEffort(selection);
        reward = data.chosenCardinalReward(selection);
        side = data.sideChoice(selection);
        correct = selectionCorrect(selection);
        
        x = tapply(effort,{effort,side},@nanmean);
        mu =  tapply(force,{effort,side},@nanmean);
        err =   tapply(force,{effort,side},@sem);
        
       [~,~,h(1)] = errorscat(x(:,1) ,mu(:,1), err(:,1),'k'); h(1).LineStyle ='--';
       [~,~,h(2)] = errorscat(x(:,2) ,mu(:,2), err(:,2),'k');

        ax = gca;
        legend([h(1) h(2)],{'left','right'});
%         ax.XTick = x;
%         ax.XLim = x;
%         ax.XLim = ax.YLim;
        ob = findobj('-property','MarkerSize');
        for io = 1:numel(ob); ob(io).MarkerSize = 10 ; end
        
%         ax.XTickLabels = {'left','right'};
        xlabel(texlabel('effort required (%fmax)'));
        ylabel('force peak (%fmax)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end   

    
%% Participation

fig = figure; set(fig,'Name','pariticipation_esmeralda');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        
        selection = selectSubject ;
        session = data.session(selection);
        treatment = data.treatment(selection);
        
        part = selectionChosen(selection);
        effort = data.cardinalEffortLeft(selection);
        reward = data.cardinalRewardLeft(selection);


        % effort
        x = tapply(effort,{effort},@nanmean);
        mu =  tapply(part,{effort},@nanmean);
        err =   tapply(part,{effort},@sem);
        
        subplot(2,numel(subjectList),iSub);hold on
        for i = 1
           [~,~,h(i)] = errorscat(x ,mu, err,col{i});
        end

        ax = gca;
        xlabel(texlabel('effort required (%fmax)'));
        ylabel('participation (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
        
        % reward
        x = tapply(reward,{reward},@nanmean);
        mu =  tapply(part,{reward},@nanmean);
        err =   tapply(part,{reward},@sem);
        
        subplot(2,numel(subjectList),numel(subjectList)+iSub);hold on
        for i = 1
           [~,~,h(i)] = errorscat(x ,mu, err,col{i});
        end

        ax = gca;
        xlabel(texlabel('reward offered '));
        ylabel('participation (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
    end  
    
%%
    fig = figure; set(fig,'Name','choice_esmeralda');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        
        selection = selectSubject & selectionChosen;
        session = data.session(selection);
        treatment = data.treatment(selection);
        
        part = selectionChosen(selection);
        de = data.deltaOrdinalEffort(selection);
        dr = data.deltaOrdinalReward(selection);
        choice = data.sideChoice(selection)-1;

        % reward
        x = tapply(dr,{dr},@nanmean);
        y =  tapply(choice,{dr},@nanmean);
        z =   tapply(choice,{dr},@sem);
        
        subplot(2,numel(subjectList),iSub);hold on
        for i = 1
           [~,~,h(i)] = errorscat(x ,y, z,col{i});
        end
        ax = gca;
        xlabel(texlabel('Delta R'));
        ylabel('right choice (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
        
          % effort
        x = tapply(de,{de},@nanmean);
        y =  tapply(choice,{de},@nanmean);
        z =   tapply(choice,{de},@sem);
        
        subplot(2,numel(subjectList),iSub+1);hold on
        for i = 1
           [~,~,h(i)] = errorscat(x ,y, z,col{i});
        end
        ax = gca;
        xlabel(texlabel('Delta E'));
        ylabel('right choice (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
        
    end  
    
%%
fig = figure; set(fig,'Name','choice2_esmeralda');
    for iSub = 1:numel(subjectList)
        sub = subjectList{iSub};
        selectSubject = ismember(data.subject, sub);
        
        
        selection = selectSubject & selectionChosen;
        session = data.session(selection);
        treatment = data.treatment(selection);
        
        part = selectionChosen(selection);
        e1 = data.ordinalEffortLeft(selection);
        e2 = data.ordinalEffortRight(selection);
        choice = data.sideChoice(selection)-1;
        x1 = tapply(e1,{e1,e2},@nanmean);
        x2 = tapply(e2,{e1,e2},@nanmean);
        y =  tapply(choice,{e1,e2},@nanmean);
        z =   tapply(choice,{e1,e2},@sem);
        
        subplot(1,numel(subjectList),iSub);hold on
        for i = 1
            x = [1:numel(y(:))];
           [~,~,h(i)] = errorscat(x ,y(:), z(:),col{i});
        end
        ax = gca;
        ax.XTick = x;
        ax.XTickLabel = {'1/1','2/1','3/1','4/1',...
                         '1/2','2/2','3/2','4/2',...
                         '1/3','2/3','3/3','4/3',...
                         '1/4','2/4','3/4','4/4'};
        xlabel(texlabel('E_left / E_right'));
        ylabel('right choice (%)'); 
        yy =ylim;
        t = title(['monkey ' sub(1) ]); 
        
    end  
    
    
    