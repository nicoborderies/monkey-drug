function [stat] = VBA_choice_ATX(x,y,paramNames,c,stat)


% define input & output
    if isequal(class(x),'table')
        predictor = x;
        x = predictor{:,1:end-1}';
        y = predictor{:,end}';
        Xnames = predictor.Properties.VariableNames(1:end-1);
    end

% function
    g_fname = @g_choice_ATX;

% priors
    nphi = 6 + 5;
    % prior mean
        priors.muPhi = zeros(nphi,1);     
%         priors.muPhi([1 2 3 4 5 6]) = 1;    

    % prior covariance
        priors.SigmaPhi = 1e0*eye(nphi); 
    options.priors = priors;        % include priors in options structure

% dimensions
    dim = struct('n',0,... % number of hidden states
        'n_theta',0 ,...   % number of evolution parameters
        'n_phi',nphi,...     % number of observation parameters
        'p',1,...          % output (data) dimension
        'n_t',size(y,2));   % number of time samples or trials
    options.dim             = dim;

% options
    
    if ~isempty(find(isnan(y)))
        options.isYout          =  zeros(dim.p,dim.n_t);   % data exclusion
        options.isYout(isnan(y))     =  1;  
    end
    y(isnan(y)) = 0;
    x(isnan(x)) = 0;
    options.DisplayWin      = 1;
    options.verbose         = 1;
    % options.extended=1;
    % options.sources(1).type = 1;    % binomial data
    % options.sources(1).out = 1;    
    options.binomial = 1;    % binomial data
    options.kernelSize = 0; 
    options.GnFigs = 0 ;

                

if nargin<5
    
    stat = struct;

    % Call inversion routine
        [posterior,out] = VBA_NLStateSpaceModel(y,x,[],g_fname,dim,options);

    % extract
        stat.beta = array2table([posterior.muPhi';diag(posterior.SigmaPhi)' ; ones(1,numel(posterior.muPhi)) ],...
                            'VariableNames',paramNames,...
                            'RowNames',{'mu','std','p'});

        % reduced model
        for i = 1:numel(posterior.muPhi)
            priors2.muPhi = zeros(nphi,1);         % prior mean on observation params
            priors2.SigmaPhi = 1e0*eye(nphi); % prior covariance on observation params
            priors2.SigmaPhi(i,i) = 0;
            [F2,~] = VBA_SavageDickey(posterior,options.priors,out.F,dim,priors2);
            dF = out.F - F2 ;
            beta{3,i} = 1./(1+exp(dF));
        end


        stat.logE = out.F;
        dF = out.F - out.diagnostics.LLH0 ;
        stat.p = 1./(1+exp(dF));
        stat.BCA = out.fit.acc;
        stat.yy = out.suffStat.gx;
        stat.corrBeta = cov2corr(posterior.SigmaPhi);
        stat.posterior = posterior;

end
    
     % contrast
    if nargin>3
        stat.contrast = array2table([ c , nan(size(c,1),1) ],...
                        'VariableNames',[ paramNames ,{'logE'}]);
        stat.contrast.logE(1) = stat.logE;
        
        for i = 1:size(c,1)
            priors2 = priors;     
            ind = c(i,:)==0;
            priors2.SigmaPhi(ind,ind) = 0 ; 
            [F2,posterior2] = VBA_SavageDickey(stat.posterior,options.priors,stat.logE,dim,priors2);
            stat.contrast.logE(i) = F2;
            if F2==nanmax(stat.contrast.logE)
                stat.beta = array2table([posterior2.muPhi';diag(posterior2.SigmaPhi)' ; ones(1,numel(posterior2.muPhi)) ],...
                            'VariableNames',paramNames,...
                            'RowNames',{'mu','std','p'});  
            end
        end
    end

    % figure
    indpos = [1 2 4 5 6];
    indneg = [3];
    for ind=1:numel(stat.beta{1,:})
        if ismember(ind,indpos)
            stat.beta{1,ind} = safepos(stat.beta{1,ind});
        elseif ismember(ind,indneg)
            stat.beta{1,ind} = -safepos(-stat.beta{1,ind});
        end
    end
    
     f = displayLogit( stat.beta ) ;
     if nargin>3
       try
         f = displayModelSelection( stat.contrast.logE );
       end
     end
     

end
