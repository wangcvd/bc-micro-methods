clearvars -except COMP COMPS ZCOMP workers

%Discretizing state space
Y = [-0.5:0.01:4]';
nY = size(Y,1);
iY = @(y) fix(min(max(0,((floor(y*100)/100)-min(Y))/0.01) +1,nY));
T= [0:0.01:1]';
nT=size(T,1);

%Loading estimates of primitives
load(sprintf('%02d/%s_%02d.mat', COMP,'density_estimates_EM', COMP))
load(sprintf('%02d/%s_%02d.mat', COMP,'pub_priv_conddensity_MLestimates', COMP))
load(sprintf('%02d/%s_%02d.mat', COMP,'entry_arrival', COMP))
load(sprintf('%02d/%s_%02d.mat', COMP,'cost_w', COMP))
alpha=alpha*ones(2,1); %uncomment if cost dist is uniform across types
ExpCDF = @(x) 1 - exp(-lambda*x);
ExpPDF = @(x) lambda*exp(-lambda*x);
% alpha = alpha*ones(2,1);
H = @(x,a) x.^a;
EC = @(x,a) (a/(1+a))*x.^(a+1);
[~, id_high]=max(MU_1+3*SIGMA_1);

%Number of teams
load CCP_Estimation_Sample_032019
teamid=teamid(competitionid==COMP,1);
clear competitionid priscore_normal pubscore_normal t t_prime 
teamid_unique=unique(teamid);
nTeams = size(teamid_unique,1);
clear teamid teamid_unique

%PDF of scores
g=zeros(nY,2);
g(:,1)=normpdf(Y,MU_1(1),SIGMA_1(1));
g(:,2)=normpdf(Y,MU_1(2),SIGMA_1(2));
g(:,1)=g(:,1)/sum(g(:,1));
g(:,2)=g(:,2)/sum(g(:,2));
F=zeros(nY,nY,2);
for k=1:2
F(1,:,k)=g(:,k)';
for j=2:nY
F(j,:,k)=g(:,k)';
F(j,j,k)=sum(g(1:j,k));
F(j,1:j-1,k)=0;
end
end

Q0 = exp(BETA*Y)./sum(exp(BETA*Y));
Qinitial=Q0;
Q=Q0;

for it=1:3

%Solve model given Q0
V = zeros(nY,nT,2);
P = V;
tic
for k=1:2
V(:,nT,k)=Q;
for i=1:nT-1
   t = nT-i;
   A = (1-ExpCDF(T(end)-T(t)))*Q;
   B = zeros(nY,1);
   for tau=t+1:1:nT
       B= B+(ExpCDF(T(tau)-T(t))-ExpCDF(T(tau-1)-T(t)))*F(:,:,k)*V(:,tau,k);
   end 
   q=max((A+B-Q),eps);
   V(:,t,k) = H(q,alpha(k)).*(A+B-EC(q,alpha(k)))+(1-H(q,alpha(k))).*Q; 
   P(:,t,k)=H(q,alpha(k));
end
end
toc

%Contest simulation
max_entrants=nTeams;
DS = 100000;
NS=200;
QPRIME=zeros(nY,NS);
tic
OUTCOMES=zeros(NS,4);

parfor z=1:NS
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
    tt_index=floor(tt*100)+1;
    S_t(1,2) = max(S_t(1,2),Y_draws(1,S_t(1,4)));
    yy_index=iY(S_t(1,2));
    
    play = P_draws(1)<P(yy_index,tt_index,S_t(1,4));
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
QPRIME(:,z)=exp(BETA*Y)./(sum(exp(BETA*S_t(:,2)))+exp(BETA*Y));
OUTCOMES(z,:)=[n_sub,n_sub_1,num_teams,max(S_t(:,2))];

end
toc
Qprime=mean(QPRIME,2);
Q0=Q;
Q=Qprime;

end

outcomes=mean(OUTCOMES,1);
save(sprintf('%02d/%s_%02d.mat', COMP,'eq_noleader_w', COMP),'Q','-v7.3')

% figure
% plot(Y,Q0)
% hold on
% plot(Y,Qprime)

% figure
% plot(Y,Qinitial)
% hold on
% plot(Y,Qprime)

error=sqrt(sum((Q0-Qprime).^2))
load modelfit
load prizes_bycomp.mat
iCOMP=max((competitionid==COMP).*[1:1:57]');
simulated=simulated(iCOMP,:);
nS_nl = [outcomes(1,2);outcomes(1,1)-outcomes(1,2)];
nS_l = [simulated(1,2);simulated(1,1)-simulated(1,2)];

disp('submissions: With leaderboard & without & diff' )
[simulated(1,1), outcomes(1,1), simulated(1,1)-outcomes(1,1)]
disp('submissions high types: With leaderboard & without & diff' )
[nS_l(id_high), nS_nl(id_high), nS_l(id_high)-nS_nl(id_high)]
disp('max score: With leaderboard & without & diff' )
[simulated(1,4), outcomes(1,4), simulated(1,4)-outcomes(1,4)]
outcomes_noleader=outcomes;
save(sprintf('%02d/%s_%02d.mat', COMP,'outcomes_noleader_w', COMP),'error','outcomes_noleader','-v7.3')

delete(gcp)
