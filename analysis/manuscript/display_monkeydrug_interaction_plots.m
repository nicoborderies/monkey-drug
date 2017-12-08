%% display_monkeydrug_interaction_plots
%


% data definition
% -- 3-level factor
subject = data.subject;
% -- 2-level factor
sessionNumber = data.sessionNumber;
treatment = data.treatment;
treatment = reordercats(treatment,trtList);
% -- 1-level factor
sum_r = round(data.ordinalRewardLeft + data.ordinalRewardRight);
sum_e = round(data.ordinalEffortLeft + data.ordinalEffortRight);
dr = round(data.ordinalRewardRight - data.ordinalRewardLeft);
de = round(data.ordinalEffortRight - data.ordinalEffortLeft);
abs_dr = abs(dr); abs_dr(dr==0)=NaN;
abs_de = abs(de); abs_de(de==0)=NaN;
r_chosen = round(data.chosenOrdinalReward); r_chosen(isnan(r_chosen)) = 0;
e_chosen = round(data.chosenCardinalEffort,1);  e_chosen(isnan(e_chosen)) = 0;
nnt = data.normalizedTrialNumber;
% -- measurements
participation = selectionParticipate; % participation
choice = data.sideChoice-1; % choice
choice(participation==0)=NaN;
choice(e_chosen==0)=NaN;
choiceR = double(sign(choice*2-1)==sign(dr));
choiceR(participation==0)=NaN;
choiceR(e_chosen==0)=NaN;
choiceR(sign(dr)==0)=NaN;
choiceE = double(sign(choice*2-1)==sign(de));
choiceE(participation==0)=NaN;
choiceE(e_chosen==0)=NaN;
choiceE(sign(de)==0)=NaN;
force = data.perf; % force
force(participation==0)=NaN;
force(e_chosen==0)=NaN;
g1 = findgroups(subject,sessionNumber,e_chosen); % force, corrected for required effort
force_effort = splitapply(@nanmean,force,g1);
g2 = findgroups(subject,sessionNumber);
force_sub = splitapply(@nanmean,force,g2);
force_corrected = force - force_effort(g1) + force_sub(g2);
rt = data.responseTime+eps; % response time
rt(participation==0)=NaN;
rt(e_chosen==0)=NaN;
g1 = findgroups(subject,sessionNumber,dr,de); 
choice_cond = splitapply(@nanmean,choice,g1);
choice_error = (choice -  choice_cond(g1)).^2 ;
rt_cond = splitapply(@nanmean,rt,g1);
rt_error = (rt -  rt_cond(g1)).^2 ;
g1 = findgroups(subject,sessionNumber,r_chosen,e_chosen);
force_cond = splitapply(@nanmean,force,g1);
force_error = (force -  force_cond(g1)).^2 ;
force_error_scaled = ((force -  force_cond(g1))./force_cond(g1)).^2 ;


% -- predictions
predicted_participation = data.predicted_participation;
predicted_choice = data.predicted_choice;
predicted_force = data.predicted_force;

% variable selection
xvar = sum_e; % sum_r , sum_e, dr, de , abs_dr, abs_de, r_chosen, e_chosen, nnt
xvar_name = '|E_{right} - E_{left}|';
x2var = sum_r; % sum_r , sum_e, dr, de , abs_dr, abs_de, r_chosen, e_chosen, nnt
x2var_name = '|R_{right} - R_{left}|';
% R_{left} + R_{right} , E_{left} + E_{right},
% |R_{right} - R_{left}|, |E_{right} - E_{left}|,
% E_{chosen}, R_{chosen}
% duration (% total trial number)
nbin = 10;
xbin = sum_e;  % sum_r , sum_e, dr, de , rchosen, echosen
x2bin = sum_r;  % sum_r , sum_e, dr, de , rchosen, echosen
yvar = participation; % participation , choiceR , choiceE , force, force_corrected, choice_error, force_error, force_error_scaled, rt_error
yvar_name = 'participation (%)'; % participation (%), choice = HR (%), choice = HE (%), force peak (%fmax), response time (s)
% choice variability (au), force variability (au), scaled force variability (au), response time variability (au)
y2var = predicted_force; % predicted_participation
trt2plot = {'placebo','atomoxetine'}; % {'clonidine','placebo','atomoxetine'};
isBin = 0;
dispPrediction = 0;
effect = 'random';
disp_x = 1;

% 1-level statistics
xbin = xvar;
x2bin = x2var;
if isBin
    xbin = quantileranks(xbin,nbin);
    x2bin = quantileranks(x2bin,nbin);
end
[g1,subject,sessionNumber,xbin,x2bin] = findgroups(subject,sessionNumber,xbin,x2bin);
treatment = splitapply(@unique,treatment,g1);
xvar = splitapply(@nanmean,xvar,g1);
x2var = splitapply(@nanmean,x2var,g1);
variableNames = {'participation','choice','force','choiceR','choiceE','yvar','y2var'...
                 'predicted_participation','predicted_choice','predicted_force'};
 for ivar = 1:numel(variableNames)
    eval([  variableNames{ivar} ' = splitapply(@nanmean,' variableNames{ivar} ',g1);' ]); 
 end
 
% 2-level statistics
g2 = findgroups(subject,treatment,xbin,x2bin);
g3 = findgroups(treatment,xbin,x2bin);
 for ivar = 1:numel(variableNames)
    eval([  variableNames{ivar} '_sub = splitapply(@nanmean,' variableNames{ivar} ',g2);' ]);
    eval([  variableNames{ivar} '_group = splitapply(@nanmean,' variableNames{ivar} ',g3);' ]);
    eval([  variableNames{ivar} '_2 = ' variableNames{ivar} ' - ' variableNames{ivar} '_sub(g2) + ' variableNames{ivar} '_group(g3) ;' ]);
 end

% display
switch effect
    case 'random'
        % - random 3 level effect
        clear g;
        alpha=0.8;
        for isub = 1:numel(unique(subject))
            subset = ( ismember(treatment,trt2plot) & nominal(subject)==subjectList{isub});
            g(1,isub) = gramm('x',xvar,'y',yvar,'size',x2var,'color',treatment,'subset',subset);
            g(1,isub).set_color_options('map',vertcat(col{ismember(trtList,trt2plot)}),'lightness',100);
            g(1,isub).set_order_options('color',trt2plot);
            if isBin
                g(1,isub).stat_summary('type','sem','geom',{'errorbar','point'},'bin_in',nbin);
            else
                g(1,isub).stat_summary('type','sem','geom',{'errorbar','point'});
            end
            g(1,isub).set_names('x',xvar_name,'y',yvar_name,'color','treatment');
            g(1,isub).set_title(subjectList{isub});
            g(1,isub).axe_property('XTick',unique(xvar));
        end
        g.draw;
        hold on;
        for isub = 1:numel(unique(subject))
            if dispPrediction; g(1,isub).update('y',y2var);  else; g(1,isub).update('y',yvar); end
            if isBin
                g(1,isub).stat_summary('type','sem','geom',{'line'},'bin_in',nbin);
            else
                g(1,isub).stat_summary('type','sem','geom',{'line'});
            end
            g(1,isub).no_legend();
        end
        g.draw;
        for isub = 1:numel(unique(subject))
            axes(g(1,isub).facet_axes_handles);
            set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                                'LineWidth',1.5,'FaceAlpha',alpha);
        end
    case 'mixed'
        % - mixed 3 level effect
        clear g;
        alpha=0.8;
        subset = ( ismember(treatment,trt2plot));
        if disp_x
            g = gramm('x',xvar,'y',yvar_2,'color',treatment,'subset',subset);
        else
            g = gramm('x',ones(size(yvar_2)),'y',yvar_2,'color',treatment,'subset',subset);
        end
        g.set_color_options('map',vertcat(col{ismember(trtList,trt2plot)}),'lightness',100);
        g.set_order_options('color',trt2plot);
        if disp_x
            if isBin
                g.stat_summary('type','sem','geom',{'errorbar','point'},'bin_in',nbin);
            else
                g.stat_summary('type','sem','geom',{'errorbar','point'});
            end
            g.set_names('x',xvar_name,'y',yvar_name,'color','treatment');
            if ~isBin; g.axe_property('XTick',unique(xvar)); end
        else
            g.stat_boxplot('width',0.6);
            g.set_names('x','','y','','color','treatment');
            g.axe_property('XTick',[]);
            g.set_title(yvar_name);

        end
        g.axe_property('Color','none');
        g.draw;
        hold on;
        if disp_x
            if dispPrediction;  g.update('y',y2var_2); else; g.update('y',yvar_2); end
            if isBin
                g.stat_summary('type','sem','geom',{'line'},'bin_in',nbin);
            else
                g.stat_summary('type','sem','geom',{'line'});
            end
            g.no_legend();
            g.draw;
        end
        axes(g.facet_axes_handles); ax= gca;
        ax.XLabel.Interpreter = 'tex';
        set_all_properties('FontName','Arial Narrow','FontWeight','bold','FontSize',16,...
                            'LineWidth',1.5,'FaceAlpha',alpha);
end
                
                
% add supplemental plot
% g.update('y',choiceE_2);
% g.stat_summary('type','sem','geom',{'errorbar','point','line'},'bin_in',nbin);
% g.no_legend();
% g.draw;
% set_all_properties('XLim',[0 1]);
% set_all_properties('YLim',[0 1]);
                
% ob = gco;
% ob2 = gco;
% legend([ob ob2],'HR','HE');