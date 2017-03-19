function [tab] = processEthogram(ethogram)

        % setup
        id =0;
        vNames={'behavior','time'};
        for iv = 1:numel(vNames)
            eval([vNames{iv} '= [];']);
        end
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
        end
        
        tab = table(behavior,time);
        tab.Properties.VariableNames = vNames;

end