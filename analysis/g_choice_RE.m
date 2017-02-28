function  [ gx ] = g_choiceRE(x,P,u,in )
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
SIDE = (SIDE>0) + 1;
Yt = [ u(7) u(8)];
Rt = [ u(9) ];
Et = [ u(10) ];

%%  parameters
      
nP = numel(P); param = [];
for iP =1:nP
    transform = inG.transform{iP};
    param(iP) = transform(P(iP));
end

kr = param(1);
kr2 = param(2);
ke = param(3);
ke2 = param(4);
kn = param(5);
krt = param(6);
ket = param(7);
kt = param(8);
k0 = param(9);

s0 = param(10);
se = param(11);

t0 = param(12);
gamma = param(13);


%% prediction

% init
gx = zeros(4,1);

% option value
    V0 = - k0 + kn*N + kt*Yt(1);
    
    V1 =  kr*R(1)*(1+krt*Rt) ...
        + ke*E(1)*(1+ket*Et) ...
        - kt*Yt(2); 
    
    V2 = (kr*R(2) + kr2*R(2))*(1+krt*Rt) ...
       + (ke*E(2) + ke2*E(2) )*(1+ket*Et) ...
       + kt*Yt(2); 
   
    V = [V1 , V2];
    Vchosen = V(SIDE);
    
    

% participation
    DV = (V1+V2) - V0 ; 
    PART = sig( DV );
    
% choice
    DV = (V2-V1) ; 
    CHOICE = sig( DV );
    
% force
    % formula
    VR = (kr*R(SIDE) + (SIDE-1)*kr2*R(2))*(1+krt*Rt) ;
    cost = @(F) (ke*F + (SIDE-1)*ke2*F)*(1+ket*Et) ;
%     cost = @(F) (ke*F^2 + (SIDE-1)*ke2*F^2)*(1+ket*Et) ;
    success = @(F) 1 - normcdf( E(SIDE) , F , s0 + se*F ) ;
    engagement = @(F) k0.*((F/E(SIDE)) - safepos(F/E(SIDE)-E(SIDE))) ;
    DV = @(F) success(F).*VR(R) + cost(F) + engagement(F) ; 
    
    % optimisation
    minusDV = @(F) - DV(F); 
    bound = [0.01, 2.00];
   [F1,V1] = fminbnd(minusDV, bound(1) , bound(2) );
   [F2,V2] =  fminbnd(minusDV, F1 , bound(2) );
    V = [V1,V2];    F = [F1,F2];
    [V,i] = min(V);
    FORCE = F(i);
    
% log-rt
    DV = (V1+V2) - V0 ; 
    LOGRT = t0 + gamma*DV;

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