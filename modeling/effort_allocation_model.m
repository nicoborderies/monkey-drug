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
    ke = 1*level_1(1);
    kn = 0;
    ks = 0;
    kf = 0;
    k0 = 0.1*level_1(1);

    se = 0.01*level_1(i);

    t0 = 0.2 ;
    gamma = 1 ; 
    
    % inputs/outputs
    R = j;
    E = k;
    N = 0.1;
    Emin=0.2;

    % formula
    motorcost = @(f) f^2;
    VR = @(R) kr*R*(1+ks*N) ;
    cost = @(F) ke*motorcost(F)*(1+kf*N) ;
    success = @(F) 1 - normcdf( E , F , se ) ;
%     engagement = @(F) k0.*((F/Emin) - safepos(F/Emin-1)) ;
    engagement = @(F) k0*(F/Emin) ;

    DV = @(F) success(F).*VR(R) - cost(F) + engagement(F) ; 
    minusDV = @(F) - DV(F); 
    dValue = @(F) (k0*exp(20))/(Emin*(exp(20) + exp((20*F)/Emin))) - 2*F*ke*(N*kf + 1) + (2^(1/2)*R*kr*exp(-(E - F)^2/(2*se^2))*(N*ks + 1))/(2*se*pi^(1/2)) ; 
    dValue = @(F) dValue(F)/10;

    % solution
    bound = [0.01, 2.00];
   [F1,V1] = fminbnd(minusDV, bound(1) , bound(2) );
   [F2,V2] =  fminbnd(minusDV, F1 , bound(2) );
    V = [V1,V2];    F = [F1,F2];
    [V,i] = min(V);
    F = F(i);
    
    % display
    figure(f); hold on;
    fplot(engagement,[0 1],'k');
%     fplot(success,[0 1],'k');
%     fplot(cost,[0 1],'k');
    h = fplot(DV,[0 1],'b');
    h = fplot(dValue,[0 1],'m');

    h.LineWidth = i;
    yy = [-success(F) success(F)]*VR(R)*1.2;
    plot([F F],[yy(1) DV(F)],'r');
    plot([0 1],[0 0],'k--');
    
    
    % legend
    ax = gca;
%     ax.YLim = yy;
    xlabel('force peak');
    ylabel('value');
%     setFigProper('FontSize',20,'LineWidth',2);
    
    
end

%%

