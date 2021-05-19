clear
set more off


use 2445



*FIGURE 1A 

hist timetodeadline_norm, fraction bin(10) scheme(sj) graphregion(color(white)) bgcolor(white) ytitle(, size(large)) xtitle(, size(large)) xlabel(,labsize(large)) ylabel(,labsize(large))
graph export Figures/SubmissionsTime.pdf, replace


*FIGURE 1B

preserve
collapse (sum) k_tc, by(competitionid firstsubmissioncompteam_normcat)
bysort competitionid : egen total_cumul_teams=total(k_tc)
bysort competitionid  (firstsubmissioncompteam_normcat): gen cumul_teams_aux=k_tc if _n==1
bysort competitionid (firstsubmissioncompteam_normcat): replace cumul_teams_aux = (cumul_teams[_n-1] +  k_tc) if _n>1
gen cumul_teams=cumul_teams_aux/total_cumul_teams
gen firstsubmissioncompteam_normcat2=firstsubmissioncompteam_normcat/10

lpoly cumul_teams firstsubmissioncompteam_normcat2, noscatter ci title(" ") graphregion(color(white)) bgcolor(white) legend(off) scheme(sj) ytitle(, size(large)) xtitle(, size(large)) xlabel(,labsize(large)) ylabel(,labsize(large))
graph export Figures/fig_cumteams.pdf, replace
restore


*FIGURE 2 

lpoly timebetweensubmissions timetodeadline_norm, noscatter title(" ") ci graphregion(color(white)) bgcolor(white) legend(off) scheme(sj) mcolor(gs12) ytitle(, size(large)) xtitle(, size(large)) xlabel(,labsize(large)) ylabel(,labsize(large))
graph export Figures/fig_timebtwsubmissions.pdf, replace


*TABLE 3

quietly: reghdfe pubscore_normal , absorb( teamid) vce(robust)
quietly: estimates store p1
quietly:reghdfe pubscore_normal  sub_number, absorb( teamid) vce(robust)
quietly: estimates store p2
quietly: reghdfe pubscore_normal  if inc75==1, absorb( teamid) vce(robust)
quietly: estimates store p3
quietly:reghdfe pubscore_normal sub_number  if inc75==1, absorb( teamid) vce(robust)
quietly: estimates store p4

esttab p1  p2 p3 p4, r2 se 


*TABLE 4: 

*Max pub score over time at team--competition level
gen maxpubscoret_team=.
bysort competitionid teamid (datesubmitted_stand): replace maxpubscoret_team=pubscore_normal if _n==1
bysort competitionid teamid (datesubmitted_stand): replace maxpubscoret_team=max(pubscore_normal,maxpubscoret_team[_n-1]) if _n>1

*Deviations from max public score
gen dev_leader=maxpubscoret-maxpubscoret_team
bysort competitionid: egen sdaux=sd(dev_leader)
gen dev_leader_st=dev_leader/sdaux
bysort competitionid teamid (datesubmitted_stand): gen lastsubmission_team = 1 if _n==_N

replace lastsubmission_team=0 if lastsubmission_team==.

*Last submission on deviation from leader
gen inter_dev_75=dev_leader_st*inc75
quietly: reg lastsubmission_team dev_leader_st, vce(robust)
estimates store p1
quietly: reg lastsubmission_team dev_leader_st inter_dev_75 inc75, vce(robust)
test dev_leader_st+inter_dev_75=0
estadd scalar asdf=r(p)
estimates store p2
esttab p1 p2, r2 se



*TABLE 5

**Drastic submissions on participation
gen timetodeadline_norm_part= round(timetodeadline_norm*1000)/1000
bysort competitionid timetodeadline_norm_part: gen Nsubms_part=_N

gen drastic01=0
bys competitionid: replace drastic01=1 if changemaxpub>=0.01 &  timetodeadline_norm>=0.2
gen time_dummy01 = drastic01*timetodeadline_norm_part

gen  dummy_after01=0
gen  dummy_RD01=.
gen time_around_disruptive=.
scalar timewindow = 0.05
levelsof competitionid, local(levelsc)
foreach comp of local levelsc{
levelsof time_dummy01 if competitionid==`comp', local(levels) 
foreach t of local levels{
qui replace dummy_after01=1 if timetodeadline_norm_part>=(`t'+ 0.001) & timetodeadline_norm_part<=(`t' + `=timewindow') & `t'>0 & competitionid==`comp'
qui replace dummy_RD01=round(1000*(`t'-timetodeadline_norm_part))/1000 if timetodeadline_norm_part>=(`t'- `=timewindow') & timetodeadline_norm_part<=(`t' + `=timewindow') & `t'>0 & competitionid==`comp'
 replace time_around_disruptive = timetodeadline_norm_part-`t'  if timetodeadline_norm_part>=(`t' - `=timewindow') & timetodeadline_norm_part<=(`t' + `=timewindow') & `t'>0 & competitionid==`comp'
}
}
gen LogN=log(Nsubms_part)
gen timesq=timetodeadline_norm*timetodeadline_norm
bysort competitionid timetodeadline_norm_part: gen k_part=1 if _n==1


**** Decomposition by types top 100
gen top100=0
replace top100=1 if ranking<=100
bysort competitionid timetodeadline_norm_part top100: gen Nsubs_top100=_N
gen LogN_100=log(Nsubs_top100)
gen inter100=dummy_after01*top100
bysort competitionid timetodeadline_norm_part top100: gen k_part100=1 if _n==1
**** Decomposition by types top 50
gen top50=0
replace top50=1 if ranking<=50
bysort competitionid timetodeadline_norm_part top50: gen Nsubs_top50=_N
gen LogN_50=log(Nsubs_top50)
gen inter50=dummy_after01*top50
bysort competitionid timetodeadline_norm_part top50: gen k_part50=1 if _n==1
**** Decomposition by types top 10
gen top10=0
replace top10=1 if ranking<=10
bysort competitionid timetodeadline_norm_part top10: gen Nsubs_top10=_N
gen LogN_10=log(Nsubs_top10)
gen inter10=dummy_after01*top10
bysort competitionid timetodeadline_norm_part top10: gen k_part10=1 if _n==1
**** Decomposition by types top 5
gen top5=0
replace top5=1 if ranking<=5
bysort competitionid timetodeadline_norm_part top5: gen Nsubs_top5=_N
gen LogN_5=log(Nsubs_top5)
gen inter5=dummy_after01*top5
bysort competitionid timetodeadline_norm_part top5: gen k_part5=1 if _n==1
**** Decomposition by types top 3
gen top3=0
replace top3=1 if ranking<=3
bysort competitionid timetodeadline_norm_part top3: gen Nsubs_top3=_N
gen LogN_3=log(Nsubs_top3)
gen inter3=dummy_after01*top3
bysort competitionid timetodeadline_norm_part top3: gen k_part3=1 if _n==1
**** Decomposition by inc75
bysort competitionid timetodeadline_norm_part inc75: gen Nsubs_inc75=_N
gen LogN_inc75=log(Nsubs_inc75)
gen inter_inc75=dummy_after01*inc75
bysort competitionid timetodeadline_norm_part inc75: gen k_inc75=1 if _n==1

quietly: reg LogN dummy_after01 timetodeadline_norm_part timesq if dummy_RD01!=. & k_part==1, r 
estimates store p1
quietly: reg LogN_inc75 dummy_after01 inter_inc75 inc75 timetodeadline_norm_part timesq if dummy_RD01!=. & k_inc75==1, r 
test dummy_after01+inter_inc75=0
estadd scalar asdf=r(p)
estimates store p2
quietly: reg LogN_50 dummy_after01 inter50 top50 timetodeadline_norm_part timesq if dummy_RD01!=. & k_part50==1, r 
test dummy_after01+inter50=0
estadd scalar asdf=r(p)
estimates store p3
quietly: reg LogN_10 dummy_after01 inter10 top10 timetodeadline_norm_part timesq if dummy_RD01!=. & k_part10==1, r 
test dummy_after01+inter10=0
estadd scalar asdf=r(p)
estimates store p4
esttab p1 p2 p3 p4, r2 se 


