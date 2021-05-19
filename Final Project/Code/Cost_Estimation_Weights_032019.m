clearvars -except COMP COMPS ZCOMP workers

%COMPetition
%COMP=2748;

%Load Data
load CCP_Estimation_Sample_032019
t=t(competitionid==COMP,1);
t_prime=t_prime(competitionid==COMP,1);
teamid=teamid(competitionid==COMP,1);
priscore_normal=priscore_normal(competitionid==COMP,1);
pubscore_normal=pubscore_normal(competitionid==COMP,1);
clear competitionid
nSub = size(t,1);

%Load Estimates
load(sprintf('%02d/%s_%02d.mat', COMP,'CCP_estimates', COMP))
load(sprintf('%02d/%s_%02d.mat', COMP,'density_estimates_EM', COMP))
load(sprintf('%02d/%s_%02d.mat', COMP,'prize', COMP))
load(sprintf('%02d/%s_%02d.mat', COMP,'random_sample', COMP))
load(sprintf('%02d/%s_%02d.mat', COMP,'pub_priv_conddensity_MLestimates', COMP))
load(sprintf('%02d/%s_%02d.mat', COMP,'weights', COMP))
p_f = @(X,b) exp(X*b)./(1+exp(X*b));

%Distribution of arrival times
last=t_prime>1.05;
timebetween=t_prime-t;
lambda=1/mean(timebetween(last==0));

%Distribution of timing of entrants times
teamid_unique=unique(teamid);
nTeams = size(teamid_unique,1);
entry_time = zeros(nTeams,1);
for i=1:nTeams,
entry_time(i,1)=min(t(teamid==teamid_unique(i,1),1));
end
mu = 1/mean(entry_time);

%Posterior on each players type
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
clear Gamma gamma_matrix aux aux_prod h i FHAT L_STAR Obs_CCP_Estimation

%Simulation of benefits for one observation: obs
ObsRS=RS.*[1:1:nSub]';
ObsRS(ObsRS<1)=[];
nRS=size(ObsRS,1);
BENEFITS=zeros(nRS,1);

obs = 10000;
max_entrants=nTeams;
DS = 100000;
NS=200; %number of simulations
benefits=zeros(NS,1);

tic
for w=1:nRS
obs=ObsRS(w,1);
%Construct matrix with state variables up to the time of submission obs
teamid_t=teamid(1:obs,1);
teamid_unique_t=unique(teamid_t);
nTeams_t = size(teamid_unique_t,1);
team_maxscore_t = zeros(nTeams_t,1);
pubscore_normal_t=pubscore_normal(1:obs,1);
t_t=t(1:obs,1);
t_prime_t=t_prime(1:obs,1);
type_t=team_maxscore_t;
y_max=type_t;
t_last=y_max;
t_draws=t_last;
for i=1:nTeams_t
    aux = pubscore_normal_t(teamid_t==teamid_unique_t(i,1),1);
    y_max(i,1) = max(aux);
    aux = t_t(teamid_t==teamid_unique_t(i,1),1);
    t_last(i,1)=aux(end);
    aux = t_prime_t(teamid_t==teamid_unique_t(i,1),1);
    t_draws(i,1)=aux(end);
    type_t(i,1)=max(type(teamid==teamid_unique_t(i,1),1));
end
%time of next play, score, identity of player, type
S_Obs=[t_draws,y_max,zeros(nTeams_t,1),type_t];

parfor z = 1:NS
%%%%%%%construct extended - state variables matrix
rng(z);
t_entry =  min(t_draws)+exprnd(1/mu,max_entrants,1); %simulate number of entrants and their entry times
t_entry=t_entry(t_entry<1,1);
n_e = size(t_entry,1);
rng(z);
u=unifrnd(0,1,n_e,1);
type_e = u>PI_1(1);
type_e = type_e + 1;
y_entry = -4*ones(n_e,1); %simulate scores of entrants

%time of next play, score, identity of player, type
S_t = [S_Obs;t_entry,y_entry,zeros(n_e,1),type_e];
[~, I]=sort(S_t(:,1));
S_t=[S_t(I,1), S_t(I,2), S_t(I,3), S_t(I,4)];
S_t(1,3)=1; %indicator for the team we're keeping track of

%draws to be used
rng(z);
T_draws_save = exprnd(1/lambda,DS,1);
Y_draws_save = [normrnd(MU_1(1),SIGMA_1(1),DS,1),normrnd(MU_1(2),SIGMA_1(2),DS,1)];
S_t_save = S_t;
P_draws_save = unifrnd(0,1,DS,1);

%Not play
S_t=S_t_save;
T_draws=T_draws_save;
Y_draws=Y_draws_save;
P_draws=P_draws_save;
i_notplay=0;

S_t(1,1)=1.1;
S_t(1,2) = max(S_t(1,2),Y_draws(1,S_t(1,4)));
P_draws(1)=[];
Y_draws(1,:)=[];
[~, I]=sort(S_t(:,1));
S_t=[S_t(I,1), S_t(I,2), S_t(I,3), S_t(I,4)];
ns_t=size(S_t(S_t(:,1)<1),1);
while ns_t>0
    tt = S_t(1,1);
    S_t(1,2) = max(S_t(1,2),Y_draws(1,S_t(1,4)));
    
    %constructing X vector in CCP
    %top scores
    K=10;
    aux=[S_t(:,2);-5*ones(K,1)];
    orderedY=-sort(-aux);
    Y=orderedY(1:K,1)';
    %indicators for being in top positions
    top1 = S_t(1,2)>=Y(:,1);
    top3 = S_t(1,2)>=Y(:,3);
    top10 = S_t(1,2)>=Y(:,10);
    top10=top10-top3;
    top3=top3-top1;
    top1=top1+0;
    X = [1,tt, tt.^2, tt.^3, S_t(1,2), S_t(1,2).^2,S_t(1,2).*tt, Y, Y.*Y, Y.*[tt*ones(1,K)],top1, top3, top10,top1.*tt, top3.*tt, top10.*tt];
    
    play = P_draws(1)<p_f(X,betahat(:,S_t(1,4)));
    if play == 1
        S_t(1,1) = tt+T_draws(1);
        T_draws(1)=[];
    else
        S_t(1,1)=1.1;
    end

    P_draws(1)=[];
    Y_draws(1,:)=[];
    [~, I]=sort(S_t(:,1));
    S_t=[S_t(I,1), S_t(I,2), S_t(I,3), S_t(I,4)];
    ns_t=size(S_t(S_t(:,1)<1),1);
    i_notplay=i_notplay+1;
end
i_notplay;
payoff_notplay = exp(BETA*S_t(S_t(:,3)==1,2))/sum(exp(BETA*S_t(:,2)));

%Play
S_t=S_t_save;
T_draws=T_draws_save;
Y_draws=Y_draws_save;
P_draws=P_draws_save;
S_t(1,1) = S_t(1,1)+T_draws(1);
S_t(1,2) = max(S_t(1,2),Y_draws(1,S_t(1,4)));
P_draws(1)=[];
T_draws(1)=[];
Y_draws(1,:)=[];
[~, I]=sort(S_t(:,1));
S_t=[S_t(I,1), S_t(I,2), S_t(I,3), S_t(I,4)];
ns_t=size(S_t(S_t(:,1)<1),1);

i_play=0;
while ns_t>0
    tt = S_t(1,1);
    S_t(1,2) = max(S_t(1,2),Y_draws(1,S_t(1,4)));
    
    %constructing X vector in CCP
    %top scores
    K=10;
    aux=[S_t(:,2);-5*ones(K,1)];
    orderedY=-sort(-aux);
    Y=orderedY(1:K,1)';
    %indicators for being in top positions
    top1 = S_t(1,2)>=Y(:,1);
    top3 = S_t(1,2)>=Y(:,3);
    top10 = S_t(1,2)>=Y(:,10);
    top10=top10-top3;
    top3=top3-top1;
    top1=top1+0;
    X = [1,tt, tt.^2, tt.^3, S_t(1,2), S_t(1,2).^2,S_t(1,2).*tt, Y, Y.*Y, Y.*[tt*ones(1,K)],top1, top3, top10,top1.*tt, top3.*tt, top10.*tt];
    
    play = P_draws(1)<p_f(X,betahat(:,S_t(1,4)));
    if play == 1
        S_t(1,1) = tt+T_draws(1);
        T_draws(1)=[];
    else
        S_t(1,1)=1.1;
    end

    P_draws(1)=[];
    Y_draws(1,:)=[];
    [~, I]=sort(S_t(:,1));
    S_t=[S_t(I,1), S_t(I,2), S_t(I,3), S_t(I,4)];
    ns_t=size(S_t(S_t(:,1)<1),1);
    i_play=i_play+1;
end
i_play;
payoff_play = exp(BETA*S_t(S_t(:,3)==1,2))/sum(exp(BETA*S_t(:,2)));
benefits(z,1)=payoff_play-payoff_notplay;
%[z,payoff_play-payoff_notplay]
end

BENEFITS(w,1)=mean(benefits);
end
toc

save(sprintf('%02d/%s_%02d.mat', COMP,'simulated_benefits', COMP),'BENEFITS','-v7.3')

%Load benefits
load(sprintf('%02d/%s_%02d.mat', COMP,'simulated_benefits', COMP))

%Optimize LogLikelihood function
t=t(RS==1,1);
t_prime=t_prime(RS==1,1);
type=type(RS==1,1);
last=t_prime>1.01;
class = (type-1)*2 + last+1;
weights = weight(class,1);


NR = 100;
FHAT = 100;
for i=1:NR,
rng(i);
alpha0=normrnd(0,1,1,1);
options=optimset('MaxIter',10000,'MaxFunEval',10000,'Display','off','TolX',1e-10,'TolFun',1e-6);
[a1,fhat] = fminunc(@(alpha) LikW_Cost(lambda,t,t_prime,max(BENEFITS,eps),type,alpha,weights),alpha0,options);
%[i,fhat,FHAT]
if fhat<FHAT
    alpha=exp(a1);
    FHAT=fhat;
end
end

options=optimset('MaxIter',250,'MaxFunEval',1,'Display','off','TolX',1e-10,'TolFun',1e-6);
[~,~,~,~,~,HESSIAN] = fminunc(@(x) LikW_Cost(lambda,t,t_prime,max(BENEFITS,eps),type,x,weights),log(alpha),options);
VAR               = diag(inv(HESSIAN*100));
SE                = sqrt(VAR/100);
%using delta method
SE=SE.*alpha;
[alpha,SE,alpha./SE]


save(sprintf('%02d/%s_%02d.mat', COMP,'cost_w', COMP),'alpha','SE','FHAT','-v7.3')

delete(gcp)

