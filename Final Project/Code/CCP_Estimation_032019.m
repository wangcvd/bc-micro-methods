clear

nRS=100; %# of obs to be selected
load prizes_bycomp.mat
COMPS=competitionid;
warning('off','all')
for z=1:57,
clearvars -except COMPS z nRS

COMP=COMPS(z,1);

%Load Data
load CCP_Estimation_Sample_032019
t=t(competitionid==COMP,1);
t_prime=t_prime(competitionid==COMP,1);
teamid=teamid(competitionid==COMP,1);
priscore_normal=priscore_normal(competitionid==COMP,1);
pubscore_normal=pubscore_normal(competitionid==COMP,1);
clear competitionid
nSub = size(t,1);

%Construct vector with max score of each team over time
teamid_unique=unique(teamid);
nTeams = size(teamid_unique,1);
team_maxscore = zeros(nSub,1);
for i=1:nTeams,
    aux = pubscore_normal(teamid==teamid_unique(i,1),1);
    team_maxscore(teamid==teamid_unique(i,1),1) = cummax(aux);
end

%Construct matrix with top K scores in the leaderboard -- using only one
%score per team (i.e., the teams' max score up to that point)
K=10;
Y=zeros(nSub,K);
for i=1:nSub,
    aux=[team_maxscore(1:i,1);-5*ones(K,1)];
    orderedY=-sort(-aux);
    Y(i,:)=orderedY(1:K,1)';
end

clear aux 

%Indicators for whether team is in top positions
top1 = team_maxscore>=Y(:,1);
top3 = team_maxscore>=Y(:,3);
top10 = team_maxscore>=Y(:,10);
top10=top10-top3;
top3=top3-top1;
top1=top1+0;

%Posterior on each players type
load(sprintf('%02d/%s_%02d.mat', COMP, 'density_estimates_EM', COMP))

gamma_matrix = zeros(nTeams,nTypes);
Gamma=zeros(size(t,1),nTypes);
for i=1:nTeams,
    aux = pubscore_normal(teamid==teamid_unique(i),1);
    aux = [prctile(aux,50), prctile(aux,75), max(aux)]';
    aux_prod = zeros(1,nTypes);
    for h=1:nTypes,
        aux_prod(1,h) = max(exp(-(aux-MU_1(h))'*(aux-MU_1(h))/(2*SIGMA_1(h)^2))/sqrt(2*pi*SIGMA_1(h)^2),eps);
    end
    for h=1:nTypes,
        gamma_matrix(i,h)=(PI_1(h)*aux_prod(:,h))/(aux_prod*PI_1);
    end
    Gamma(teamid==teamid_unique(i),:)=ones(sum(teamid==teamid_unique(i)),1)*gamma_matrix(i,:);
end
type = Gamma(:,1)<Gamma(:,2);
type = 1 + type;
%adjusting type calculation if sampling leads to too few players of a
%particular type
if min(sum(gamma_matrix(:,1)>gamma_matrix(:,2)), sum(gamma_matrix(:,2)>gamma_matrix(:,1)))<0.5*nTeams*min(PI_1)
    na=floor(nTeams*min(PI_1));
    [~,I]=min(PI_1);
    sorted = sort(gamma_matrix(:,I));
    type(Gamma(:,I)>sorted(end-na+1))=I;
    display('adjustment in type calculation was needed')
end

%Matrix of state variables
X = [ones(size(t,1),1),t, t.^2, t.^3, team_maxscore, team_maxscore.^2,team_maxscore.*t, Y, Y.*Y, Y.*[t*ones(1,K)],top1, top3, top10,top1.*t, top3.*t, top10.*t];

%Distribution of arrival times
last=t_prime>1.05;
timebetween=t_prime-t;
lambda=1/mean(timebetween(last==0));
ExpCDF = @(x) 1 - exp(-lambda*x);
ExpPDF = @(x) lambda*exp(-lambda*x);

%Optimize LogLikelihood function
NR = 15;
FHAT = 100;
betahat=zeros(size(X,2),1);
for i=1:NR,
rng(i);
beta0=normrnd(0,1,size(X,2),2);
options=optimset('MaxIter',10000,'MaxFunEval',10000,'Display','off','TolX',1e-10,'TolFun',1e-6);
[b1,fhat] = fminunc(@(x) Lik_CCP(lambda,X,t,t_prime,type,x),beta0,options);
%[i,fhat,FHAT]
if fhat<FHAT,
    betahat=b1;
    FHAT=fhat;
end
end

Obs_CCP_Estimation=nSub;
save(sprintf('%02d/%s_%02d.mat', COMP,'CCP_estimates', COMP),'betahat','FHAT','Obs_CCP_Estimation','-v7.3')

PrPlay = exp(X*betahat(:,1));
PrPlay = PrPlay./(1+PrPlay);
PrPlay1=PrPlay;

PrPlay = exp(X*betahat(:,2));
PrPlay = PrPlay./(1+PrPlay);
PrPlay2=PrPlay;

PrPlay = PrPlay1.*(type==1) + PrPlay2.*(type==2);
PrPlay(isnan(PrPlay)==1)=1;
% cdfplot(PrPlay)
[z, COMP]
mean(PrPlay)
end
warning('on','all')