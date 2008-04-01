function [genes]=SAM_mod(data,idx,np,desired_fdg,varargin)
%   SAM_mod
%       Automation of the Significance Analysis of Microarrays
%
%   genes=SAM_mod(data,idx,np,desired_fdg)
%           data            Gene Expression data inwhich each row
%                           represents a gene and each column represents 
%                           a sample.
%           idx             A boolean array differentiating between the two
%                           disease states.
%           np              Number of permutations desired (default is 100)
%           desired_fdg     The desired number of Median Falsely Detected
%                           Genes
%
%           genes           The row indices of genes which show
%                           differential expression between the two classes
%
%       SAM_mod(...,'Waitbar',WaitBarHandle,[start stop])
%       
%       WaitBarHandle       The handle of a waitbar object
%       [start stop]        The proportion of the wait inwhich this 
%                           instance should cover.
%
%
%   Modified from code produced by Bilge:  A binary search is preformed to
%   automatically determine the largest possible gene list with a specific 
%   Median Number of Falsely Detected Genes

if  isempty(varargin)
    h=waitbar(0,'Processing SAM Data');
    bar_prefix=[];
    bar_start=0;
    bar_stop=1;
elseif strcmpi(varargin{1},'waitbar')
    h=varargin{2};
    bar_prefix=get(h,'UserData');
    bar_start=varargin{3}(1);
    bar_stop=varargin{3}(2);
end

[nf,ns]=size(data);
labels=unique(idx);
I1=find(idx==labels(1));
n1=length(I1);
I2=find(idx==labels(2));
n2=length(I2);
mu1=mean(data(:,I1),2);
mu2=mean(data(:,I2),2);
sigma1=std(data(:,I1),[],2);
sigma2=std(data(:,I2),[],2);
a=(1/n1+1/n2)/(n1+n2-2);
s=sqrt(a*(sum((data(:,I1)-mu2*ones(1,n1)).^2,2)+sum((data(:,I2)-mu1*ones(1,n2)).^2,2)));
q=quantile(s,(1:100)/100);
alphas=(0:.05:1);
cv=zeros(1,length(alphas));
waitbar(bar_start,h,[bar_prefix ' Calculating Alphas'])
for i=1:length(alphas)
    waitbar(bar_start+(bar_stop-bar_start)*(0.4*i/length(alphas)),h)
    alpha=alphas(i);
    dalpha=(mu1-mu2)./(s+quantile(s,alpha));
    v=zeros(1,99);
    for j=1:99
        I=find((s>=q(i))&(s<q(i+1)));
        v(j)=median(abs(dalpha-median(dalpha)));
    end
    cv(i)=std(v)/mean(v);
end
s0=quantile(s,min(alphas(find(cv==min(cv)))));
d=(mu2-mu1)./(s+s0);
[sd,sdo]=sort(d);

waitbar(bar_start+0.4*(bar_stop-bar_start),h,[bar_prefix 'Preforming Permutations'])
dp=zeros(nf,np);
for i=1:np
    waitbar(bar_start+(bar_stop-bar_start)*(0.4+.4*i/np),h)
    [nouse,o]=sort(rand(1,ns));
    pI1=o(1:n1);
    pI2=o((n1+1):end);
    pmu1=mean(data(:,pI1),2);
    pmu2=mean(data(:,pI2),2);
    psigma1=std(data(:,pI1),[],2);
    psigma2=std(data(:,pI2),[],2);
    ps=sqrt(a*(sum((data(:,pI1)-pmu2*ones(1,n1)).^2,2)+sum((data(:,pI2)-pmu1*ones(1,n2)).^2,2)));
    dp(:,i)=sort((pmu2-pmu1)./(ps+s0));
end
de=mean(dp,2);

fc=zeros(nf,1);
I=find(mu1~=0);
fc(I)=mu2(I)./mu1(I);

delta=unique(sort(abs(sd-de)));
%ntd=zeros(size(delta));
%fdr=zeros(size(delta));
lower=1;
upper=length(delta);
i=ceil((upper-lower)/2)+lower;
waitbar(bar_start+0.8*(bar_stop-bar_start),h,[bar_prefix 'Finding Gene List'])
while(upper~=i)
    i=ceil((upper-lower)/2)+lower;
    ntd=nnz(abs(sd-de)>delta(i));
    ld=max([sd(1) sd(find(sd-de<-delta(i), 1, 'last' ))]);
    ud=min([sd(end) sd(find(sd-de>delta(i), 1 ))]);
    fdr=mean(mean((dp<=ld)|(dp>ud)));
    if nf*fdr>desired_fdg    %more falsely detected genes then desired
        lower=i;    %so move to higher delta
    elseif nf*fdr<desired_fdg   %less falsely detected genes then desired
        upper=i;
    else
        break
    end
end

waitbar(bar_stop,h,[bar_prefix 'Completing SAM'])
ginds=find(abs(sd-de)>delta(i));
genes=sdo(ginds);
