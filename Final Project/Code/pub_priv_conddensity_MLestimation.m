clear

load prizes_bycomp.mat
COMPS=competitionid;

for z=1:57,
clearvars -except COMPS z 

COMP=COMPS(z,1);

%Load Data
load CCP_Estimation_Sample_032019
t=t(competitionid==COMP,1);
t_prime=t_prime(competitionid==COMP,1);
teamid=teamid(competitionid==COMP,1);
pri=priscore_normal(competitionid==COMP,1);
pub=pubscore_normal(competitionid==COMP,1);
clear competitionid
N = size(t,1);

%Assume score_{priv} = \alpha + score_{pub} \beta + \varepsilon,
%where \varepsilon is distributed according to a Gumbel distribution
%(see https://en.wikipedia.org/wiki/Gumbel_distribution)
%We normalize the dispersion parameter to 1

%Density: f(\varepsilon) = exp(-(\varepsilon + exp(-\varepsilon))) 
%Log likelihood of an obs: l_i = -(\varepsilon + exp(-\varepsilon))
theta0=[0;1];
f = @(theta) sum(((pri-theta(1)-pub*theta(2))+exp(-(pri-theta(1)-pub*theta(2)))))/N;

NS=20;
Initial=normrnd(0,1,2,NS);
options=optimset('MaxIter',2500,'MaxFunEval',10000,'Display','off','TolX',1e-10,'TolFun',1e-6);
TH = zeros(2,NS);
FH = zeros(1,NS);

for i=1:NS,
    [thetahat, fhat]=fminsearch(@(theta) f(theta),Initial(:,i),options);
    TH(:,i)=thetahat;
    FH(:,i)=fhat;
end
[~, I]=min(FH);
thetahat = TH(:,I);
fhat = FH(:,I);

BETA = thetahat(2);

save(sprintf('%02d/%s_%02d.mat', COMP,'pub_priv_conddensity_MLestimates', COMP),'BETA','-v7.3')

end
