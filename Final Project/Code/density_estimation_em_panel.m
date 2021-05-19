clear

load prizes_bycomp.mat
COMPS=competitionid;

for z=1:57,
clearvars -except COMPS z
COMP=COMPS(z,1);
nTypes=2;
load starting_v2

%Load Data
load CCP_Estimation_Sample_032019
t=t(competitionid==COMP,1);
t_prime=t_prime(competitionid==COMP,1);
teamid=teamid(competitionid==COMP,1);
priscore_normal=priscore_normal(competitionid==COMP,1);
pubscore_normal=pubscore_normal(competitionid==COMP,1);
clear competitionid
NObs = size(t,1);

%Teams
player_uniq=unique(teamid);
player_sel=teamid;
N = size(player_uniq,1);
SCORES = zeros(N,3);
num_scores=3;
for i=1:N,
    aux = pubscore_normal(player_sel==player_uniq(i),1); %full vector of scores for player
    aux = [prctile(aux,50), prctile(aux,75), max(aux)]';
    SCORES(i,:)=aux;
end
SCORES_RESHAPED = reshape(SCORES,[N*num_scores 1]);

gamma_matrix = zeros(N,nTypes);
LIK=cell(NS,1);
ML=zeros(1,NS);
TH=zeros(3*nTypes,NS);
f_0=zeros(size(SCORES_RESHAPED,1),nTypes);

for j=1:NS,
    MU_0 = STARTING(1:nTypes,j);%[-0.5; 0.5]+0.1;
    SIGMA_0 = exp(STARTING((nTypes+1):2*nTypes,j));
    PI_0 = STARTING((2*nTypes+1):(3*nTypes-1),j)/(sum(STARTING((2*nTypes+1):(3*nTypes-1),j))+0.5);
    PI_0=[1-sum(PI_0); PI_0];
    THETA_0 = [MU_0;SIGMA_0;PI_0];
    THETA = THETA_0;
    
    error = 1;
    tol = 1e-8;
    for i=1:nTypes,
    f_0(:,i) = PI_0(i)*normpdf(SCORES_RESHAPED,MU_0(i),SIGMA_0(i));
    end
    L_0 = mean(log(sum(f_0,2))); %log likelihood evaluated at initial parameters
    L=[L_0];
    while error>tol,
        %Expectation step
        Gamma=zeros(NObs,nTypes);
        for i=1:N,
            aux = SCORES(i,:)';
            aux_prod = zeros(1,nTypes);
            for h=1:nTypes,
                aux_prod(1,h) = max(exp(-(aux-MU_0(h))'*(aux-MU_0(h))/(2*SIGMA_0(h)^2))/sqrt(2*pi*SIGMA_0(h)^2),eps); 
                        %compute prod of densities
            end
            for h=1:nTypes,
                gamma_matrix(i,h)=(PI_0(h)*aux_prod(:,h))/(aux_prod*PI_0); %compute posterior on type
            end

        end
        
        %Maximization step
        MU_1=zeros(nTypes,1);
        SIGMA_1=MU_1;
        for h=1:nTypes,
           MU_1(h,1)=sum(kron(ones(num_scores,1),gamma_matrix(:,h)).*SCORES_RESHAPED)/sum(kron(ones(num_scores,1),gamma_matrix(:,h)));
           SIGMA_1(h,1)=sqrt(sum(kron(ones(num_scores,1),gamma_matrix(:,h)).*((SCORES_RESHAPED-MU_1(h)).^2))/sum(kron(ones(num_scores,1),gamma_matrix(:,h))));
        end
        PI_1=mean(gamma_matrix,1)'; %PI_1 is updated at the team level (as in the Dreb and Trivedi paper)
        THETA_1 = [MU_1;SIGMA_1;PI_1];
        
        %Convergence
        for i=1:nTypes,
            f_0(:,i) = PI_1(i)*normpdf(SCORES_RESHAPED,MU_1(i),SIGMA_1(i));
        end
        L_1= mean(log(sum(f_0,2)));
        if L_1<L_0,
            L_1=L_0;
            SIGMA_1=SIGMA_0;
            PI_1=PI_0;
            MU_1=MU_0;
        end
        error=(L_1-L_0)^2;
        MU_0 = MU_1;
        SIGMA_0 = SIGMA_1;
        PI_0 = PI_1;
        L_0=L_1;
        L=[L;L_1];
    end
    LIK{j,1}=L;
    TH(:,j)=THETA_1;
    ML(:,j)=L_1;
end
[~, I]=max(ML);
THETA_STAR = TH(:,I);
MU_1 = THETA_STAR(1:nTypes,1);
SIGMA_1 = THETA_STAR((nTypes+1):2*nTypes,1);
PI_1 = THETA_STAR((2*nTypes+1):3*nTypes,1);
iter=size(LIK{I,1},1);
L_STAR=max(LIK{I,1});
save(sprintf('%02d/%s_%02d.mat', COMP,'density_estimates_EM', COMP),'nTypes','MU_1','SIGMA_1','PI_1','L_STAR','-v7.3')

[z, COMP]
PI_1'

end
