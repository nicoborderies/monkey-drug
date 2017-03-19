%% doEthogram_analysis

root = 'C:\Users\bastien.blain\Desktop';

% data to select
subList =  {'bob','aliosha'};
sessionList.aliosha = {'3_29_2016','3_30_2016','3_31_2016','4_1_2016'};
sessionList.bob = {'3_29_2016','3_30_2016','3_31_2016','4_1_2016'};

for isub = 1:numel(subList)
   for isess =  1:numel(sessionList.(subList{isub}))
    
       
       % loadings
       filename = ['ethogram_' subList{isub} '_' sessionList.(subList{isub}){isess} '*.mat' ];
       load([ root '\ethogram\data' filesep filename  ]);
       
       
   end
end