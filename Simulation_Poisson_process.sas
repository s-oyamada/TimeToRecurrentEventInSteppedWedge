/*

******************************************************************************************************************************************************************************************
  Program name : Simulation_Poisson_process.sas
  Author : Shunsuke Oyamada
  SAS version : 9.4
  Description : Simulation of time-to-recurrent-event generated based on (Mixed-) Poisson process
                in a stepped wedge cluster randomized trial using open cohort design
  Reference : Oyamada S, Chiu SW, Yamaguchi T, 2022.
              Comparison of statistical models for estimating intervention effects based on time-to-recurrent-event in stepped wedge cluster randomized trial using open cohort design.
              Under submission to journal.
  Notes : This simulation code contains only the analyses with stratification by clusters.
******************************************************************************************************************************************************************************************

Copyright (c) 2022 Shunsuke Oyamada

*/


%macro sim_poisson
(cluster_size=,      /*sample size in each cluster*/
 cluster_num=,       /*number of clusters (= number of steps)*/
 switch_distance=,   /*distance between switches (= length between steps)*/
 intv_effect=,       /*true intervention effect (regression coefficient)*/
 scale=,             /*scale parameter of exponential distribution*/
 sim_num=,           /*number of simulations*/
 sigma_c=,           /*variance on the random effect representing the variation between clusters*/
 sigma_s=,           /*variance on the random effect representing the variation between subjects*/
 entry=,             /*timing of trial entry*/
 follow=,            /*follow-up period*/
 entry_range=,       /*acceptable range of entry to the trial, stepend or studyend*/
 seed_e=,            /*seed of pseudo-random number to be used to determine the point of trial entry*/
 seed_r1=,           /*seed of the pseudo-random number to be used for the 1st time-to-recurrent-event*/
 seed_r2=,           /*seed of the pseudo-random number to be used for the 2nd time-to-recurrent-event*/
 seed_r3=,           /*seed of the pseudo-random number to be used for the 3rd time-to-recurrent-event*/
 seed_c=,            /*seed of pseudo-random number to be used to determine the cluster effect*/
 seed_s=             /*seed of pseudo-random number to be used to determine the inter-individual variability*/
);


/*
**************************
  Data generation process
**************************
*/
data ttre_ (drop= _:);
  steplength = &switch_distance;
  retain _prn1 &seed_e _prn2 &seed_r1 _prn3 &seed_r2 _prn4 &seed_r3 _prn5 &seed_c _prn6 &seed_s;
  do sim = 1 to &sim_num;
    do cluster = 1 to &cluster_num;
    call rannor(_prn5, _cerror);
    cerror = &sigma_c*_cerror;
      do id = 1 to &cluster_size;
        call ranuni(_prn1,_e);
        call ranuni(_prn2,_u1);
        call ranuni(_prn3,_u2);
        call ranuni(_prn4,_u3);
        call rannor(_prn6, _serror);
        serror = &sigma_s*_serror;

        *generate SWCRT scheme;
        studybegin = '01feb2017'd;
        stepend    = '27jan2018'd;
        studyend   = stepend + (&switch_distance * &follow);

        entrydate  = studybegin+int((&entry_range-studybegin)*(_e/&entry));
        switchdate = intnx(steplength,studybegin,cluster,'same');
        switchtime = max(switchdate-entrydate, 0);

        *generate censoring time;
        call streaminit(1234);
        c        = rand("weibull", 1.7191, 272.1588/*Inversion of 0.003674*/);
        dropdate = entrydate+int(c);

        *generate true survival time (exponential-distributed baseline hazard function);
        _logu1  = -log(_u1);
        _logu2  = -log(_u2);
        _logu3  = -log(_u3);
        _linpre = cerror + serror;
        _scalee = &scale*exp(_linpre);
        _nom1   = _logu1-_scalee*(switchtime)+&scale*exp(_linpre+&intv_effect)*(switchtime);
        _nom2   = _logu2-_scalee*(switchtime)+&scale*exp(_linpre+&intv_effect)*(switchtime);
        _nom3   = _logu3-_scalee*(switchtime)+&scale*exp(_linpre+&intv_effect)*(switchtime);
        _denom  = &scale*exp(_linpre+&intv_effect);
        if _logu1 < _scalee*(switchtime) then t1 = (_logu1/_scalee); else t1 = (_nom1/_denom);
        if _logu2 < _scalee*(switchtime) then t2 = (_logu2/_scalee); else t2 = (_nom2/_denom);
        if _logu3 < _scalee*(switchtime) then t3 = (_logu3/_scalee); else t3 = (_nom3/_denom);

        output;
        format studybegin stepend studyend switchdate dropdate entrydate date9.;
      end;
    end;
  end;
run;

data ttre_2_;
  set ttre_(in=in1 keep=studybegin stepend studyend switchdate dropdate switchtime entrydate cluster t1 sim id rename=(t1=t))
      ttre_(in=in2 keep=studybegin stepend studyend switchdate dropdate switchtime entrydate cluster t2 sim id rename=(t2=t))
      ttre_(in=in3 keep=studybegin stepend studyend switchdate dropdate switchtime entrydate cluster t3 sim id rename=(t3=t))
  ;
run;

proc sort data=ttre_2_; by sim cluster id t; run;

data ttre_3_;
  set ttre_2_;
  by sim cluster id;
  if first.id then visit=0;
  visit+1;
  recdate = entrydate+int(t);
run;

proc sort data=ttre_3_; by sim cluster id studybegin stepend studyend switchdate dropdate switchtime entrydate visit; run;

proc transpose data=ttre_3_ out=ttre_4_ prefix=recdate;
  var recdate;
  by sim cluster id studybegin stepend studyend switchdate dropdate switchtime entrydate;
  id visit;
run;

data ttre_5_;
  set ttre_4_(drop=_NAME_);
  *calculate censoring indicator;
  *1st record;
  if recdate1 > min(dropdate, studyend) then do;
    Status1    = 0;
    obenddate1 = min(dropdate, studyend);
  end;
  else do;
    Status1    = 1;
    obenddate1 = recdate1;
  end;
  *2nd record;
  if recdate2 > min(dropdate, studyend) then do;
    Status2    = 0;
    obenddate2 = min(dropdate, studyend);
  end;
  else do;
    Status2    = 1;
    obenddate2 = recdate2;
  end;
  *3rd record;
  if recdate3 > min(dropdate, studyend) then do;
    Status3    = 0;
    obenddate3 = min(dropdate, studyend);
  end;
  else do;
    Status3    = 1;
    obenddate3 = recdate3;
  end;

  *calculate observed survival time;
  survtime1 = obenddate1-entrydate;
  survtime2 = (obenddate2-obenddate1)+survtime1;
  survtime3 = (obenddate3-obenddate2)+survtime2;

  *counting process style;
  TStart1 = 0;
  TStop1 = survtime1;
  TStart2 = survtime1;
  TStop2 = survtime2;
  TStart3 = survtime2;
  TStop3 = survtime3;

  output;
  format studybegin stepend studyend switchdate dropdate recdate1 obenddate1 recdate2 obenddate2 recdate3 obenddate3 entrydate date9.;
run;

data ttre;
  set ttre_5_(in=in1 keep=studybegin stepend studyend switchdate dropdate switchtime entrydate cluster recdate1 obenddate1 survtime1 TStart1 TStop1 Status1 sim id
                     rename=(recdate1=recdate obenddate1=obenddate survtime1=survtime TStart1=TStart TStop1=TStop Status1=Status))
      ttre_5_(in=in2 keep=studybegin stepend studyend switchdate dropdate switchtime entrydate cluster recdate2 obenddate2 survtime2 TStart2 TStop2 Status2 sim id
                     rename=(recdate2=recdate obenddate2=obenddate survtime2=survtime TStart2=TStart TStop2=TStop Status2=Status))
      ttre_5_(in=in3 keep=studybegin stepend studyend switchdate dropdate switchtime entrydate cluster recdate3 obenddate3 survtime3 TStart3 TStop3 Status3 sim id
                     rename=(recdate3=recdate obenddate3=obenddate survtime3=survtime TStart3=TStart TStop3=TStop Status3=Status))
  ;
  if in1 then visit=1;
  if in2 then visit=2;
  if in3 then visit=3;
run;

proc sort data=ttre; by sim cluster id visit; run;

data ttre;
  set ttre;
  if obenddate < studybegin then delete;
run;


/*
**********************************************************
  AG(Andersen-Gill) model with stratification by clusters
**********************************************************
*/
ods output parameterestimates=ag_strt;
ods listing close;
proc phreg data = ttre covs;
  class cluster;
  model (TStart, TStop)*Status(0)=z_t/rl ties=exact;
  if TStop < switchtime then z_t=0; else z_t=1;
  by sim;
  strata cluster;
run;

data ag_strt_est;
  set ag_strt;
  where parameter='z_t';
  variable = parameter;
  model    = 'ag';
  true     = &intv_effect;
  bias     = estimate-true;
  mse      = bias**2;
  cover    = (hrlowercl<= exp(true)<= hruppercl);
  drop parameter;
run;


/*
******************************************************
  Data processing for the implementation of PWP model
******************************************************
*/
proc sort data=ttre; by sim cluster id visit; run;
data ttre_pwp(keep=TStart TStop Status id cluster visit switchtime Gaptime sim obenddate switchdate entrydate);
  retain LastStatus;
  set ttre;
  by sim cluster id;
  if first.id then LastStatus=1;
  if (Status=0 and LastStatus=0) then delete;
  LastStatus=Status;
  Gaptime=TStop-TStart;
run;


/*
***********************************************************************************
  PWP(Prentice-Williams-Peterson) Total-Time model with stratification by clusters
***********************************************************************************
*/
ods output parameterestimates=pwp_tt_strt;
ods listing close;
proc phreg data = ttre_pwp;
  model (TStart, TStop)*Status(0)=z_t/rl ties=exact;
  if TStop < switchtime then z_t=0; else z_t=1;
  strata Visit cluster;
  by sim;
run;

data pwp_tt_strt_est;
  set pwp_tt_strt;
  where parameter='z_t';
  variable = parameter;
  model    = 'pwp_tt';
  true     = &intv_effect;
  bias     = estimate-true;
  mse      = bias**2;
  cover    = (hrlowercl<= exp(true)<= hruppercl);
  drop parameter;
run;


/*
*********************************************************************************
  PWP(Prentice-Williams-Peterson) Gap-Time model with stratification by clusters
*********************************************************************************
*/
ods output parameterestimates=pwp_gt_strt;
ods listing close;
proc phreg data = ttre_pwp;
  model Gaptime*Status(0)=z_t/rl ties=exact;
  if Gaptime < switchtime then z_t=0; else z_t=1;
  strata Visit cluster;
  by sim;
run;

data pwp_gt_strt_est;
  set pwp_gt_strt;
  where parameter='z_t';
  variable = parameter;
  model    = 'pwp_gt';
  true     = &intv_effect;
  bias     = estimate-true;
  mse      = bias**2;
  cover    = (hrlowercl<= exp(true)<= hruppercl);
  drop parameter;
run;


/*
************************************************
  Calculate control-to-intervention ratio ratio
************************************************
*/
data ttre_pwp; set ttre_pwp;
  controltime = max(min(obenddate,switchdate)-entrydate, 0);
  trttime = max(obenddate-max(entrydate,switchdate), 0);
run;

proc summary data=ttre_pwp;
  var controltime trttime; by sim;
  output out=ttre_cti sum=;
run;

data ttre_cti; set ttre_cti;
  ctiratio=controltime/trttime;
run;


/*
******************************************************************
  Calculate average observed number of recurrences per individual
******************************************************************
*/
data ttre_anr_; set ttre_pwp;
  where Status=1;
run;

proc freq data=ttre_anr_;
  tables Status / noprint out=ttre_anr;
run;

data ttre_anr; set ttre_anr;
  anr=count/(&cluster_size*&cluster_num*&sim_num);
run;


/*
**********************************
  Calculate censoring proportions
**********************************
*/
proc sort data=ttre_pwp; by sim cluster id Status; run;

proc sort data=ttre_pwp out=ttre_cens nodupkey; by sim cluster id; run;


/*
*********
  Output
*********
*/
ods listing;
title1 "Poisson, n=&cluster_size, m=&cluster_num, steplength=&switch_distance, sigma_c=&sigma_c, sigma_s=&sigma_s" ;
title2 "entry=&entry, follow=&follow, entry_range=&entry_range, intervention_effect=&intv_effect, scale=&scale, sim_num=&sim_num" ;


title3 'AG(Andersen-Gill) model with stratification by clusters';
proc means data = ag_strt_est n mean;
  var estimate bias mse cover;
  output out=ag_strt_summary n= mean= / autoname; run;
run;


title3 'PWP(Prentice-Williams-Peterson) Total-Time model with stratification by clusters';
proc means data = pwp_tt_strt_est n mean;
  var estimate bias mse cover;
  output out=pwp_tt_strt_summary n= mean= / autoname; run;
run;


title3 'PWP(Prentice-Williams-Peterson) Gap-Time model with stratification by clusters';
proc means data = pwp_gt_strt_est n mean;
  var estimate bias mse cover;
  output out=pwp_gt_strt_summary n= mean= / autoname; run;
run;


title3 "Control-to-intervention ratio";
proc means data = ttre_cti n mean;
    var ctiratio;
run;


title3 "Average observed number of recurrences per individual";
proc means data=ttre_anr mean;
var anr;
run;


title3 "Censoring proportions";
proc freq data=ttre_cens;
tables Status / missing;
run;


title;
quit;

%mend sim_poisson;
