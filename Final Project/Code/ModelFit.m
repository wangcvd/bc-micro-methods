clear

load prizes_bycomp.mat
COMPS=competitionid;
simulated=zeros(57,4);
actual=zeros(57,4);


for ww=1:57,
clearvars -except COMPS ww simulated actual

COMP=COMPS(ww,1);

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
load(sprintf('%02d/%s_%02d.mat', COMP,'entry_arrival', COMP))
p_f = @(X,b) exp(X*b)./(1+exp(X*b));

%Posterior on each players type
teamid_unique=unique(teamid);
nTeams = size(teamid_unique,1);
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

actual(ww,:)=[nSub,sum(type==1,1),nTeams,max(pubscore_normal)];

%Simulation of contests

obs = 10000;
max_entrants=nTeams;
DS = 100000;
NS=200; %number of simulations
%parpool(12)
%[number of submissions, num of submissions type 1, number of teams, max score]
OUTCOMES=zeros(NS,4);

tic
parfor z = 1:NS
%%%%%%%construct extended - state variables matrix
rng(z);
t_entry = exprnd(1/mu,max_entrants,1); %simulate number of entrants and their entry times
t_entry=t_entry(t_entry<1,1);
n_e = size(t_entry,1);
rng(z);
u=unifrnd(0,1,n_e,1);
type_e = u>PI_1(1);
type_e = type_e + 1;
y_entry = -4*ones(n_e,1); %simulate scores of entrants

%time of next play, score, identity of player, type
S_t = [t_entry,y_entry,zeros(n_e,1),type_e];
[~, I]=sort(S_t(:,1));
S_t=[S_t(I,1), S_t(I,2), S_t(I,3), S_t(I,4)];

%draws to be used
rng(z);
T_draws = exprnd(1/lambda,DS,1);
Y_draws = [normrnd(MU_1(1),SIGMA_1(1),DS,1),normrnd(MU_1(2),SIGMA_1(2),DS,1)];
P_draws = unifrnd(0,1,DS,1);

n_sub=0;
ns_t=size(S_t(S_t(:,1)<1),1);
num_teams=ns_t;
n_sub_1=0;
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
    n_sub=n_sub+1;
    n_sub_1=n_sub_1+(S_t(1,4)==1);
end

OUTCOMES(z,:)=[n_sub,n_sub_1,num_teams,max(S_t(:,2))];

end
ww
toc

simulated(ww,:)=mean(OUTCOMES,1);

end 

figure
x=log(actual(:,1));
y=log(simulated(:,1));
scatter(x,y)
hold on
refline(1,0)
hold off
title('total submissions (logs): simulated vs actual')
a=corr(x,y);
text(min(x)+0.1, max(y)-0.1, ['Coef of correlation = ' num2str(a)])


figure
x=log(actual(:,2));
y=log(simulated(:,2));
scatter(x,y)
hold on
refline(1,0)
hold off
title('total submissions type 1 (logs): simulated vs actual')
a=corr(x,y);
text(min(x)+0.1, max(y)-0.1, ['Coef of correlation = ' num2str(a)])

figure
x=log(actual(:,1)-actual(:,2));
y=log(simulated(:,1)-simulated(:,2));
scatter(x,y)
hold on
refline(1,0)
hold off
title('total submissions type 2 (logs): simulated vs actual')
a=corr(x,y);
text(min(x)+0.1, max(y)-0.1, ['Coef of correlation = ' num2str(a)])

figure
x=log(actual(:,1)./actual(:,3));
y=log(simulated(:,1)./simulated(:,3));
scatter(x,y)
hold on
refline(1,0)
hold off
title('av submissions per team (logs): simulated vs actual')
a=corr(x,y);
text(min(x)+0.1, max(y)-0.1, ['Coef of correlation = ' num2str(a)])

figure
x=log(actual(:,4));
y=log(simulated(:,4));
scatter(x,y)
hold on
refline(1,0)
hold off
title('overall max score (logs): simulated vs actual')
a=corr(x,y);
text(min(x)+0.1, max(y)-0.1, ['Coef of correlation = ' num2str(a)])

save('modelfit.mat','actual','simulated','-v7.3')

%exporting to csv for analysis in Stata
clear
close all

load prizes_bycomp
contests=competitionid;
clear prize_1 prize_2 prize_3 prize_4 prize_5
load modelfit

out=zeros(57,7);
%
EST = zeros(size(contests,1),8);
for i=1:size(contests,1),
    COMP=contests(i,1);
load(sprintf('%02d/%s_%02d.mat', COMP,'density_estimates_EM', COMP))
A = [actual(i,2);actual(i,1)-actual(i,2)];
S = [simulated(i,2);simulated(i,1)-simulated(i,2)];
[~, ind]=max(MU_1+3*SIGMA_1);
 out(i,:)=[COMP,actual(i,1),simulated(i,1),A(ind,1),S(ind,1),actual(i,4),simulated(i,4)];
end

csvwrite('modelfit_out.csv',out)

