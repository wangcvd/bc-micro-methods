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
priscore_normal=priscore_normal(competitionid==COMP,1);
pubscore_normal=pubscore_normal(competitionid==COMP,1);
clear competitionid
nSub = size(t,1);

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

save(sprintf('%02d/%s_%02d.mat', COMP,'entry_arrival', COMP),'mu','lambda','-v7.3')

[z, COMP]
end