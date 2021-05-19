clear

nRS=100; %# of obs to be selected
load prizes_bycomp.mat
COMPS=competitionid;

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

%Posterior on each players type
load(sprintf('%02d/%s_%02d.mat', COMP,'density_estimates_EM', COMP))
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

p=unifrnd(0,1,nSub,1);
last= t_prime>1.05;

class = (type-1)*2 + last+1;
A=grpstats(ones(size(class)),class,'sum');
nS_Class = zeros(4,1);
nS_Class(2) = min(10,A(2));
nS_Class(4) = min(10,A(4));
nS_Class(1)=50-nS_Class(2);
nS_Class(3)=50-nS_Class(4);

%rate of last submission
rate(1)=mean(1./grpstats(ones(size(class(type==1,1))),teamid(type==1,1),'sum'));
rate(2)=mean(1./grpstats(ones(size(class(type==2,1))),teamid(type==2,1),'sum'));

RS=zeros(size(class));
for w=1:1:4
    if nS_Class(w,1)==0
        continue
    end
sorted = sort(p(class==w,1));
RS(class==w,1)=p(class==w,1)>=sorted(end-nS_Class(w)+1);
end

% trs=type(RS==1,1);
% minT=min(sum([trs==1, trs==2]));
% h=h+1;
[z, nS_Class',sum(RS)]

Q=[PI_1(1)*(1-rate(1)),PI_1(1)*rate(1),PI_1(2)*(1-rate(2)),PI_1(2)*rate(2)]';
H = nS_Class./sum(nS_Class);
weight = Q./H;
save(sprintf('%02d/%s_%02d.mat', COMP,'random_sample', COMP),'RS','-v7.3')
save(sprintf('%02d/%s_%02d.mat', COMP,'weights', COMP),'weight','-v7.3')

end


