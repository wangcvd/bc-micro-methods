clear
close all

load prizes_bycomp.mat
PRIZE= prize_1+prize_2+prize_3+prize_4+prize_5;
COMPS=competitionid;
SUBS=[];
SUBSH=[];
MAXS=[];
out=[];
for z=1:1
    COMP=COMPS(z,1);
    
    %if isfile(sprintf('%02d/%s_%02d.mat', COMP,'outcomes_noleader_w', COMP))==1
        load(sprintf('%02d/%s_%02d.mat', COMP,'outcomes_noleader_w', COMP))
        load(sprintf('%02d/%s_%02d.mat', COMP,'cost_w', COMP))
        expcost=alpha./(1+alpha);
        outcomes=outcomes_noleader;
        load modelfit
        simulated=simulated(z,:);
        load(sprintf('%02d/%s_%02d.mat', COMP,'density_estimates_EM', COMP))
        [~, id_high]=max(MU_1+3*SIGMA_1);
        nS_nl = [outcomes(1,2);outcomes(1,1)-outcomes(1,2)];
        nS_l = [simulated(1,2);simulated(1,1)-simulated(1,2)];
        SUBS=[SUBS;z,COMP,simulated(1,1), outcomes(1,1), simulated(1,1)-outcomes(1,1)];
        SUBSH=[SUBSH;z,COMP,nS_l(id_high), nS_nl(id_high), nS_l(id_high)-nS_nl(id_high)];
        MAXS=[MAXS;z,COMP,simulated(1,4), outcomes(1,4), simulated(1,4)-outcomes(1,4)];
        out=[out;[z,simulated(1,1), outcomes(1,1),nS_l(id_high), nS_nl(id_high),simulated(1,4), outcomes(1,4),expcost,PRIZE(z,1),COMP]];
    %end

end

csvwrite('outputstata_strategic.csv',out)