function  [ gx ] = g_choiceRE(x,P,u,inG )
% INPUT
% - x : [useless]
% - P : regression coefficients ((n+1)x1)
% - u : design matrix (nx1)
% - in : [useless]
% OUTPUT
% - gx : P(response=1|u)


%% inputs / outputs

R = [ u(1) u(2) ];
E = [ u(3) u(4) ];
N = [ u(5) ];
SIDE = [ u(6) ];
SIDE = SIDE + (SIDE==0);
OMIT = (u(6)==0);
Yt = [ u(7) u(8)];
Rt = [ u(9) ];
Et = [ u(10) ];

%%  parameters
      
nP = numel(P); param = [];
for iP =1:nP
    transform = inG.transform{iP};
    param(iP) = transform(P(iP));
end

% generic parameters
kr = param(1); % reward weight
ks = param(2); % satiety weight
ke = param(3); % effort weight
kf = param(4); % fatigue weight
k0 = param(5); % activation weight
kn = param(6); % dynamic weight

% decision parameters
beta = 1; % temperature
bm = param(7); % motor bias
bp = param(8); % perseveration bias

% force parameters
se = param(9); % motor noise

% response-time parameters
t0 = param(10); % non-decision time
theta = param(11); % rt scaling factor


%% prediction

% init
gx = zeros(4,1);

% option value
    motorcost = @(f) f^2;
    
%     V0 = - k0 - bp*Yt(1);
%     V1 =  kr*R(1)*(1+ks*N) ...
%         - ke*motorcost(E(1))*(1+kf*N) ...
%         - bp*Yt(2); 
%     V2 = (kr*R(2))*(1+ks*N) ...
%        - (ke*motorcost(E(2)))*(1+kf*N) ...
%        + bp*Yt(2) ...
%        + bm; 
   
    V0 = - ( k0*(1+kf*N) + bp*Yt(1) );
    V1 =  kr*R(1) ...
        - ke*motorcost(E(1)) ...
        - bp*Yt(2); 
    V2 = (kr*R(2)) ...
       - (ke*motorcost(E(2))) ...
       + bp*Yt(2) ...
       + bm; 
   
   
    V = [V1 , V2];
    Vchosen = V(SIDE);
    Vchosen = Vchosen*(OMIT==0);

    

% participation
    DV = (V1+V2) - V0 ; 
    PART = sig( DV/beta );
    
% choice
    DV = (V2-V1) ; 
    CHOICE = sig( DV/beta );
    
% force
    % formula
    Emin=0.2;
    
    % 1) computational form
%     VR = @(R) kr*R(SIDE)*(1+ks*N) ;
%     cost = @(F) ke*motorcost(F)*(1+kf*N) ;
%     success = @(F) 1 - normcdf( E(SIDE) , F , se ) ;
%     engagement = @(F) k0.*((F/Emin) - safepos(F/Emin-1)) ;
%     engagement = @(F) k0*F ;
%     DV = @(F) success(F).*VR(R) - cost(F) + engagement(F) ; 
    
    % 1) full numerical optimisation
%     minusDV = @(F) - DV(F); 
%     bound = [0.01, 2.00];
%    [F1,V1] = fminbnd(minusDV, bound(1) , bound(2) );
%    [F2,V2] =  fminbnd(minusDV, F1 , bound(2) );
%     V = [V1,V2];    F = [F1,F2];
%     [V,i] = min(V);
%     FORCE = F(i);

%     2) analytical form
      % force threshold basal benfit + normal distributed force
%       dValue = @(F) (k0*exp(20))/(Emin*(exp(20) + exp((20*F)/Emin))) - 2*F*ke*(N*kf + 1) + (2^(1/2)*R(SIDE)*kr*exp(-(E(SIDE) - F)^2/(2*se^2))*(N*ks + 1))/(2*se*pi^(1/2)) ; 
      % constant force basal benfit + normal distributed force
%       dValue = @(F) k0-2*F*ke*(N*kf+1)+(2^(1/2)*R(SIDE)*kr*exp(-(E(SIDE)-F)^2/(2*se^2))*(N*ks+1))/(2*se*pi^(1/2));
      % constant force basal benfit + gamma distributed force
      % - with only activational dynamic
%       dValue = @(F) k0*(N*kf+1)-2*F*ke+(2^(1/2)*R(SIDE)*kr*exp(-(log(E(SIDE))-log(F))^2/(2*se^2)))/(2*F*se*pi^(1/2)) ;
%       negValue = @(F)F^2*ke-F*k0*(N*kf+1)+R(SIDE)*kr*(erf((2^(1/2)*(log(E(SIDE))-log(F)))/(2*se))/2-1/2);
      % - with full dynamics
%       dValue = @(F) k0*(N*kn+1)-2*F*ke*(N*kf+1)+(2^(1/2)*R(SIDE)*kr*exp(-(log(E(SIDE))-log(F))^2/(2*se^2))*(N*ks+1))/(2*F*se*pi^(1/2));
%       negValue = @(F) F^2*ke*(N*kf+1)-F*k0*(N*kn+1)+R(SIDE)*kr*(N*ks+1)*(erf((2^(1/2)*(log(E(SIDE))-log(F)))/(2*se))/2-1/2);
      % with full dynamics + normal distributed force + scaled variability
%       negValue = @(F) F^2*ke*(N*kf+1)-F*k0*(N*kn+1)+R(SIDE)*kr*(N*ks+1)*(erf((2^(1/2)*(E(SIDE)-F))/(2*F*se))/2-1/2);
        % with full dynamics + normal distributed force + unscaled variability
      negValue = @(F) F^2*ke*(N*kf+1)-F*k0*(N*kn+1)+R(SIDE)*kr*(N*ks+1)*(erf((2^(1/2)*(E(SIDE)-F))/(2*se))/2-1/2);
      
%     2) analytical differentiation + numerical resolution
%       FORCE = fzero(dValue,[eps 2-eps]);
      FORCE = fminbnd(negValue,eps,2-eps);

    
% log-rt
    uncertainty = - ( (1-PART).*log(1-PART) + (PART.*CHOICE).*log(PART.*CHOICE) + (PART.*(1-CHOICE)).*log(PART.*(1-CHOICE)) );
    LOGRT = t0 + theta*uncertainty;

% integrate
    gx = [PART ; CHOICE ; FORCE ; LOGRT ];

%% derivative
%     dgdx = [];
%     dgdP = zeros(numel(P),1);

end

function y=sig(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end