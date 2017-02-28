function [ output_args ] = tapply( x, list, fname, factorisation, nbin )

if nargin<4
    for i = 1:length(list)
         factorisation{i} = 'discrete';
    end
end
for i = 1:length(list)
    if strcmp(factorisation{i},'continuous') && nargin<5
        nbin = repmat(6,1,length(list));
    end
end

x=x(:);
c=['['];
for i = 1:length(list)
  l=list{i};
  l=l(:);
  switch factorisation{i}
      case 'discrete'
        dim(i)=numel(unique(l));
      case 'continuous'
        dim(i)=nbin(i);
  end
  
  
  if length(l) ~=length(x)
      error(['length(var)~=length(list{' num2str(i) '}'])
  end
  c=[c 'i' num2str(i) ' '];
end

c=[c ']'];

if length(list)==1
output_args=nan(1,dim);
else
output_args=nan(dim);
end

for i = 1:numel(output_args)
    eval([c '=ind2sub(dim,i);']);
    test=ones(1, length(x));
    for j = 1:length(list)
        eval(['ind = i' num2str(j) ';']);
        l=list{j};
        l=l(:);
        switch factorisation{j}
             case 'discrete'
                f=unique(l);
                test(l~=f(ind))=0;
            case 'continuous'
                f = quantile(l,[1/nbin(j):1/nbin(j):1]);
                if ind==1
                    test( l>f(ind))=0;    
                else
                    test( l<=f(ind-1) | l>f(ind))=0;
                end
        end
    end
    try
    output_args(i)=fname(x(test==1));
    end
end
end

