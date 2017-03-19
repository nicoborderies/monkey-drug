%% cost_benefit_model_analytics
% analytical derivations of the cost-benefit model
%
% Nicolas Borderies - March 2017


clear all;
clc; close all;

%%
% declare variables
%%% inputs
parameterList = {'R','E','N','Emin'};
np = numel(parameterList);
for ip=1:np
   eval( [parameterList{ip} ' = sym(''' parameterList{ip} '''); '] ); 
end

%%% parameters
parameterList = {'kr','ke','ks','kf','k0','se'};
np = numel(parameterList);
for ip=1:np
   eval( [parameterList{ip} ' = sym(''' parameterList{ip} '''); '] ); 
end

%%% output
parameterList = {'F'};
np = numel(parameterList);
for ip=1:np
   eval( [parameterList{ip} ' = sym(''' parameterList{ip} '''); '] ); 
end

% functions
syms benefit(kr,R,ks,N) ;
benefit(kr,R,ks,N) = kr*R*(1+ks*N);

syms cost(ke,F,kf,N) ;
cost(ke,F,kf,N) = ke*(F^2)*(1+kf*N);

syms normCDF(x,mu,sigma) ;
normCDF(x,mu,sigma) =  0.5*(1+erf( (x-mu)/(sigma*sqrt(2)) )) ;

syms success(E,F,se)
success(E,F,se) = 1 - normCDF( E , F , se );

syms safepos(x)
safepos(x) = log(1+exp(20*x))/20;

syms engagement(k0,F,Emin)
engagement(k0,F,Emin) = k0*((F/Emin) - safepos(F/Emin-1));
% engagement(k0,F,Emin) = k0*(F/E) ;
% engagement(k0,F,Emin) = 0 ;


syms value(kr,ke,k0,ks,kf,se,R,E,N,Emin,F)
value(kr,ke,k0,ks,kf,se,R,E,N,Emin,F) = ...
    engagement(k0,F,Emin) + benefit(kr,R,ks,N)*success(E,F,se) - cost(ke,F,kf,N) ;


% differentiation
dvalue = diff(value,F);
d2value = diff(dvalue,F);

% find optimal value
%%% assumptions
parameterList = {'kr','ke','k0','ks','kf','se','R','E','N','Emin','F'};
np = numel(parameterList);
for ip=1:np
   x = parameterList{ip};
   eval( [ 'assume(' x '>0 & ' x '<3 & ' x '~=0); '] ); 
end

eqn = dvalue(kr,ke,k0,ks,kf,se,R,E,N,Emin,F)==0 ; 
% eqn = dvalue(1,1,1,1,1,0.1,2,0.4,1,0.2,F)==0 ; 
% eqn2 = d2value(kr,ke,k0,ks,kf,se,R,E,N,0.2,F)<0 ; 
Fstar = solve(eqn,F,'Real',1);

% display
dvalue = simplify(dvalue);
eval ( [ 'dV = @(F) ' char(dvalue(kr,ke,k0,ks,kf,se,R,E,N,Emin,F)) ]);
eval ( [ 'dV = @(F) ' char(dvalue(1,1,1,1,1,0.1,2,0.4,1,0.2,F)) ]);
f0 = fzero(dV,1);


% taylor approximation
approx_dvalue = taylor(dvalue,F,0.4,'Order',10);
% eqn = approx_dvalue(kr,ke,k0,ks,kf,se,R,E,N,Emin,F)==0 ; 
% eqn = approx_dvalue(1,1,1,1,1,0.1,2,0.4,1,0.2,F)==0 ; 

% Fstar = solve(eqn,F);
% eval ( [ 'approx_dV = @(F) ' char(approx_dvalue(kr,ke,k0,ks,kf,se,R,E,N,Emin,F)) ]);
eval ( [ 'approx_dV = @(F) ' char(approx_dvalue(1,1,1,1,1,0.1,2,0.4,1,0.2,F)) ]);
f02 = fzero(approx_dV,1);


