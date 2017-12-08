%% effort_allocation_model - choiceRE
function [] = effort_allocation_model()

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
    param.kr = 1*level_1(1);
    param.ke = 0.25*level_1(1);
    param.kn = 0;
    param.ks = 0;
    param.kf = 0;
    param.k0 = 0.1*level_1(1);

    param.se = 0.01*level_1(i);
    param.alpha = 0.25*level_1(i);

    param.t0 = 0.2 ;
    param.gamma = 1 ; 
    
    % inputs/outputs
    input.R = j;
    input.E = k;
    input.N = 0.1;
    input.Emin=0.2;

%     % formula
%     motorcost = @(f) f^2;
%     VR = @(R) kr*R*(1+ks*N) ;
%     cost = @(F) ke*motorcost(F)*(1+kf*N) ;
%     success = @(F) 1 - normcdf( E , F , se ) ;
% %     engagement = @(F) k0.*((F/Emin) - safepos(F/Emin-1)) ;
%     engagement = @(F) k0*(F/Emin) ;
% 
%     DV = @(F) success(F).*VR(R) - cost(F) + engagement(F) ; 
%     minusDV = @(F) - DV(F); 
%     dValue = @(F) (k0*exp(20))/(Emin*(exp(20) + exp((20*F)/Emin))) - 2*F*ke*(N*kf + 1) + (2^(1/2)*R*kr*exp(-(E - F)^2/(2*se^2))*(N*ks + 1))/(2*se*pi^(1/2)) ; 
%     dValue = @(F) dValue(F)/10;
% 
%     % solution
%     bound = [0.01, 2.00];
%    [F1,V1] = fminbnd(minusDV, bound(1) , bound(2) );
%    [F2,V2] =  fminbnd(minusDV, F1 , bound(2) );
%     V = [V1,V2];    F = [F1,F2];
%     [V,i] = min(V);
%     F = F(i);
    
    [output] = force_model(input,param);
    F = output.F_star;
    
    % display
    figure(f); hold on;
    fplot(output.activation,[0 1],'k');
    fplot(output.success,[0 1],'k');
    fplot(output.cost,[0 1],'k');
    h = fplot(output.value,[0 1],'b');
%     h = fplot(dValue,[0 1],'m');

    h.LineWidth = i;
    yy = [-output.success(F) output.success(F)]*output.benefit(input.R)*1.2;
    plot([F F],[yy(1) output.value(F)],'r');
    plot([0 1],[0 0],'k--');
    
    
    % legend
    ax = gca;
%     ax.YLim = yy;
    xlabel('force peak');
    ylabel('value');
    
    
end

end




function [output] = force_model(input,param)
    % formula
    output.motorcost = @(f) f^2;
    output.benefit = @(R) param.kr*R ;
    output.cost = @(F) param.ke*output.motorcost(F) ;
%     output.success = @(F) 1 - normcdf( E , F ,  param.se ) ;
    output.success = @(F) 1 - gamcdf( input.E , F/(param.alpha)+1 ,  1/(param.alpha) ) ;
    output.activation  = @(F)  param.k0*F*(1+ param.kn*input.N);
    output.value  = @(F) output.success(F).*output.benefit(input.R) - output.cost(F) + output.activation(F);
    output.neg_value  = @(F) - output.value(F);

%     DV = @(F) success(F).*VR(R) - cost(F) + engagement(F) ; 
%     minusDV = @(F) - DV(F); 
%     dValue = @(F) (k0*exp(20))/(Emin*(exp(20) + exp((20*F)/Emin))) - 2*F*ke*(N*kf + 1) + (2^(1/2)*R*kr*exp(-(E - F)^2/(2*se^2))*(N*ks + 1))/(2*se*pi^(1/2)) ; 
%     dValue = @(F) dValue(F)/10;
%     
%     % solution
    bound = [0.01, 2.00];
   [F1,V1] = fminbnd(output.neg_value, bound(1) , bound(2) );
   [F2,V2] =  fminbnd(output.neg_value, F1 , bound(2) );
    V = [V1,V2];    F = [F1,F2];
    [V,i] = min(V);
    output.F_star = F(i);
    

end


