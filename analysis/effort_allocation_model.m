%% effort_allocation_model - choiceRE
clc;
clear all;
close all;

%%
% range
level_1 = [1 2 3 ];
level_2 = [1 2 3 4];
level_3 = [0.2 0.4 0.6 0.8];

% setup
f = figure;
f.Name = 'force_model simulations';
i=1;
j=1;
k=0.4;

for i = level_1
% for j = level_2
% for k = level_3

    % parameters
    kr = 1*level_1(1);
    ke = -(1)*level_1(1);
    kn = 0;
    krt = 0;
    ket = 0;
    kt = 0;
    k0 = 0.5*level_1(1);

    s0 = 0.02*level_1(i);
    se = 0.01*level_1(1);

    t0 = 0.2 ;
    gamma = 1 ; 
    
    % inputs/outputs
    R = j;
    E = k;
    N = 0.1;

    F = 0;

    % formula
    Fmin = 0.1;
    VR = @(R) kr.*R ;
%     cost = @(F) ke*F ;
    cost = @(F) ke.*F.^2 ;
    success = @(F)  1 - normcdf( E , F , s0 + se.*F ) ;
%     engagement = @(F) k0.*(F/E) ;
%     engagement = @(F) k0.*(1 - safepos(-F/E+E)) ;
    engagement = @(F) k0.*((F/E) - safepos(F/E-E)) ;

    DV = @(F) success(F).*VR(R) + cost(F) + engagement(F) ; 
    minusDV = @(F) - DV(F); 

    % solution
    bound = [0.01, 2.00];
   [F1,V1] = fminbnd(minusDV, bound(1) , bound(2) );
   [F2,V2] =  fminbnd(minusDV, F1 , bound(2) );
    V = [V1,V2];    F = [F1,F2];
    [V,i] = min(V);
    F = F(i);
    
    % display
    figure(f); hold on;
%     fplot(engagement,[0 1],'k');
%     fplot(success,[0 1],'k');
%     fplot(cost,[0 1],'k');
    h = fplot(DV,[0 1],'b');
    h.LineWidth = i;
    yy = [-success(F) success(F)]*VR(R)*1.2;
    plot([F F],[yy(1) DV(F)],'r');
    
    
    % legend
    ax = gca;
%     ax.YLim = yy;
    xlabel('force peak');
    ylabel('value');
    setFigProper('FontSize',20,'LineWidth',2);
    
    
end

%%

