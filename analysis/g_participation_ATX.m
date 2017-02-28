function  [ gx,dgdx,dgdP ] = g_participation_ATX(x,P,u,in )
% INPUT
% - x : [useless]
% - P : regression coefficients ((n+1)x1)
% - u : design matrix (nx1)
% - in : [useless]
% OUTPUT
% - gx : P(response=1|u)

beta = P;
indpos = [1 2 5 6 7   9 10 11 12   13 16 17 18 19 20];
indneg = [3 4 8];
beta(indpos) = safepos(beta(indpos));
beta(indneg) = -safepos(-beta(indneg));

tmax=10;

y = 0;
p0= 0;
p1= 0;
d0 = zeros(1,tmax-1);

%% prediction
    % option value
        kr = beta(1:2)' + beta(9)*(u(6)==2);
        krt =  beta(17) + beta(19)*(u(6)==2);
        kr = kr*( 1 / ( 1+ krt*u(17) ) );
        
        ke = beta(3:4)' + beta(11)*(u(6)==2);
        ket =  beta(18) + beta(20)*(u(6)==2);
        ke = ke*( 1 + ket*u(18) );
        
       w =    [kr , ke];
       y = y + w*u(1:4); 
    
    % dynamic
        p0 = beta(5) + beta(13)*(u(6)==2) ; 
        p1 = beta(6) + beta(14)*(u(6)==2);
        k =  beta(7) + beta(15)*(u(6)==2) ;
        kf = beta(8) + beta(16)*(u(6)==2) ;
        
        for it=2:tmax
            d0(it-1) =  u(6+it)*exp(-k*(it-1)) ;
        end
        y = y + p0 ;
        y = y +  (p1*u(7))*(1+sum(d0));
%         y = y +  p0*(u(7)<0) +  p1*(u(7)>0) ;
        y = y + kf*u(5) ;

    % sigmoidal transform
        gx = sig( y );
    

%% derivative
    dgdx = [];
    dgdP = zeros(numel(P),1);
    
    for iu = 1:4
        dgdP(iu) = u(iu)*gx*(1-gx);
%         dgdP(iu+8) = (u(6)==2)*u(iu)*gx*(1-gx);
    end
    dgdP(9) = (u(6)==2)*(u(1)+u(2))*gx*(1-gx);
    dgdP(10) = 0;
    dgdP(11) = (u(6)==2)*(u(3)+u(4))*gx*(1-gx);
    dgdP(12) = 0;

    dgdP(5) = gx*(1-gx);
    dgdP(6) = u(7)*(1+sum(d0))*gx*(1-gx);
%     dgdP(5) = (u(7)<0)*gx*(1-gx);
%     dgdP(6) = (u(7)>0)*gx*(1-gx);
    
    dgdP(7) = -( u(7) +([1:tmax-1]*d0'))*gx*(1-gx);
    dgdP(8) = u(5)*gx*(1-gx);
    
    dgdP(13) = (u(6)==2)*gx*(1-gx);
    dgdP(14) = ((u(6)==2)*u(7))*(1+sum(d0))*gx*(1-gx);
%     dgdP(13) = (u(6)==2)*(u(7)<0)*gx*(1-gx);
%     dgdP(14) = (u(6)==2)*(u(7)>0)*gx*(1-gx);
    
    dgdP(15) = -(u(6)==2)*( u(7) +([1:tmax-1]*d0'))*gx*(1-gx);
    dgdP(16) = (u(6)==2)*u(5)*gx*(1-gx);
    
    dgdP(17) = (-u(17))*( 1 / ( 1 + beta(17)*u(17) )^2 )...
               *sum( beta(1:2)' + (beta(9:10)')*(u(6)==2) )*gx*(1-gx);
    dgdP(18) = u(18)...
               *sum( beta(3:4)' + (beta(11:12)')*(u(6)==2) )*gx*(1-gx);
    dgdP(19) = (-u(17)*(u(6)==2))*( 1 / ( 1 + (beta(17) + beta(19)*(u(6)==2))*u(17) )^2 )...
               *sum( beta(1:2)' + (beta(9:10)')*(u(6)==2) )*gx*(1-gx);
    dgdP(20) = u(18)*(u(6)==2)...
               *sum( beta(3:4)' + (beta(11:12)')*(u(6)==2) )*gx*(1-gx);
           
           

end

function y=sig(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end