% function  [ gx,dgdx,dgdP ] = g_choice_ATX(x,P,u,in )
function  [ gx ] = g_choice_ATX(x,P,u,in )

% INPUT
% - x : [useless]
% - P : regression coefficients ((n+1)x1)
% - u : design matrix (nx1)
% - in : [useless]
% OUTPUT
% - gx : P(response=1|u)

beta = P;
indpos = [1 2 4 5 6];
indneg = [3];
beta(indpos) = safepos(beta(indpos));
beta(indneg) = -safepos(-beta(indneg));

r = u(1:2);
e = u(3:4);
nt = u(5);
trt = u(6);
rr = u(7);
ee = u(8);


%% prediction
    % decision value
       kr =  beta(1) + beta(7)*(trt==2)  ;
       krt = beta(5) + beta(8)*(trt==2)  ;
       kr = kr*( 1 / ( 1+ krt*rr ) );
       
       ke =   beta(3) + beta(9)*(trt==2);
       ket =  beta(6) + beta(10)*(trt==2);
       ke = ke*( 1 + ket*ee );
       
       vr = kr*(beta(2)*r(2)-r(1)) ;
       ve = ke*(beta(4)*e(2)-e(1)) ;
       
       temperature = 1 + beta(11)*(trt==2); 
       dv = (vr + ve)*(temperature);

    % sigmoidal transform
        gx = sig( dv );
    

%% derivative
%     dgdx = [];
%     dgdP = zeros(numel(P),1);
%     
%     
%     dgdP(1) = (1+beta(5)*nt)*(beta(2)*r(2)-r(1))*gx*(1-gx);
%     dgdP(2) = (beta(1)*(1+beta(5)*nt))*(r(2))*gx*(1-gx);
%     dgdP(3) = (1+beta(6)*nt)*(beta(4)*e(2)-e(1))*gx*(1-gx);
%     dgdP(4) = (beta(3)*(1+beta(6)*nt))*(e(2))*gx*(1-gx);
%     
%     dgdP(5) = beta(1)*nt*(beta(2)*r(2)-r(1))*gx*(1-gx);
%     dgdP(6) = beta(3)*nt*(beta(4)*e(2)-e(1))*gx*(1-gx);
%     
%     dgdP(7) = (trt==2)*(beta(2)*r(2)-r(1))*gx*(1-gx);
%     dgdP(8) = (trt==2)*nt*(beta(2)*r(2)-r(1))*gx*(1-gx);
%     dgdP(9) = (trt==2)*(beta(4)*e(2)-e(1))*gx*(1-gx);
%     dgdP(10) = (trt==2)*nt*(beta(4)*e(2)-e(1))*gx*(1-gx);
%     dgdP(11) = (trt==2)*(vr + ve)*gx*(1-gx);


end

function y=sig(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end