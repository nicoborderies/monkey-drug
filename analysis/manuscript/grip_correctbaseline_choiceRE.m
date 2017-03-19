function [ truebaseline, gripdataCorrect ] = grip_correctbaseline_choiceRE( gripdata )


trialduration   = 3; 
threshold1      = 2; % threshold to position the baseline.
threshold2      = 5; % threshold to position the dips.

% f = figure; hold on; set(f,'Color',[1 1 1]);
ind = cellfun(@isempty,gripdata.grip(:));
indTrial = find(ind==0)';

for i=indTrial
    
    
    % first, extract baseline
    bestfit=0;
    pointbas=round(min(gripdata.grip{i}))-50; % start 20 Newton below minimal value (to have room for the threshold window)
    pointhaut=round(max(gripdata.grip{i})*.8); % end at maximal value
    baseline(i)=pointbas;
    for B=pointbas:pointhaut
        
        dist=abs(B-gripdata.grip{i});
        fit=sum(dist<threshold1);
        
        if fit>bestfit
            baseline(i)=B;
            bestfit=fit;
        end
        
    end
    
    
    % Next, find eventual dips before and after the force peak (that indicate a negative baseline)
    index_peak  = find(gripdata.grip{i}==max(gripdata.grip{i}),1,'first');
    index_min1  = find(gripdata.grip{i}(1:index_peak)==min(gripdata.grip{i}(1:index_peak)),1,'first');
    index_min2  = find(gripdata.grip{i}(index_peak:end)==min(gripdata.grip{i}(index_peak:end)),1,'first')+index_peak-1;
    
    diff1=abs(baseline(i)-gripdata.grip{i}(index_min1));
    diff2=abs(baseline(i)-gripdata.grip{i}(index_min2));
    
    if (diff1>threshold2) & (diff2>threshold2)
        min1=gripdata.grip{i}(index_min1);
        min2=gripdata.grip{i}(index_min2);
        dip(i)=1;
    else
        min1=[];
        min2=[];
        dip(i)=0;
    end
    
    % Correct Grip dynamic (no correction occurs if no dips are detected)
    gripdataCorrect{i}=gripdata.grip{i};
    if ~isnan(index_min1) & ~isnan(index_min2)
        if gripdataCorrect{i}(index_min1)<= 2*baseline(i) ...
                && max(gripdataCorrect{i}(1:index_min1))<= 0.25*max(gripdataCorrect{i}) % disp correction requires closeness from interpolated baseline & positive ascending signal
            gripdataCorrect{i}(1:index_min1)=-gripdataCorrect{i}(1:index_min1);
            dip(i)=1;
        else
            dip(i)=0;
        end
        if  gripdataCorrect{i}(index_min2)<= 2*baseline(i)
            gripdataCorrect{i}(index_min2:end)=-gripdataCorrect{i}(index_min2:end);
        end
    end
    
    switch dip(i)
        case 0
            truebaseline(i) = (baseline(i));
        case 1
            truebaseline(i) = -(baseline(i));
    end
    gripdataCorrect{i} = gripdataCorrect{i} - truebaseline(i);


%     % Plot trial grip to check errors  
% 
%     clf('reset');
%     hold on;
%     plot(gripdata.grip{i});
%     xx=xlim;yy=ylim;
%     plot([1 length(gripdata.grip{i})],[baseline(i) baseline(i)],'g');
%     if ~isempty(min1),plot([index_min1 index_min1],[yy(1) yy(2)],'r'),end;
%     if ~isempty(min1),plot([index_min2 index_min2],[yy(1) yy(2)],'r'),end;
%     text(xx(2),yy(2),num2str(i));
%     hold off;
%     pause
%     
%     % Plot trial grip to check correction  
%     clf('reset');
%     hold on;
%     h1 = plot(gripdata.grip{i},'k');
%     h2 = plot(gripdataCorrect{i},'Color',[0.5 0.5 0.5]);
%     h3 = plot(ones(1,numel(gripdata.grip{i}))*baseline(i),'r');
%     xlabel('time');
%     ylabel('raw force)');
%     title('selection of responses with correct timing');
%     legend([h1,h2,h3],{'raw signal','corrected signal','baseline'})
%     hold off;
%     pause
%     
    
end


% % Extract median true baseline
% isdip   = round(mean(dip));
% switch isdip
%     case 0
%         truebaseline=median(baseline);
%     case 1
%         truebaseline=-median(baseline);
% end


end

