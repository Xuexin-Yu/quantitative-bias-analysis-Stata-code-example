*****************************************************************
Here, we provide example Stata code to calculate the simulation estimates and intervals for the three-way interaction term between loneliness duration, years of follow-up, and gender, with sensitivity values following a trapezoidal distribution between 0.25 and 0.30 among men and 0.35-0.40 among women, and specificity values were constantly hold at 0.92-1.00 for both genders with a trapezoidal uniform distribution.
*****************************************************************

******************************************
****Step 1*********************************
****bias parameter specification***********
*****Trapezoidal distribution, N=10000****
*******************************************

gen rownumber=_n
	   
****1.1 sensitivity value for men**********	   
gen MIN=0.25
gen MAX=0.30 
gen DL=MIN+0.015
gen DU=MAX-0.015
***here, we specified the trapezoidal distribution for the selected bias range set (min dl du max)***

gen X=(runiform()*(MAX+DU-DL-MIN) + (MIN+DL))/2 if rownumber<10001
gen se_m=. 
replace  se_m = X if DL <X & X<DU  & rownumber<10001
replace  se_m = MIN+sqrt((DL-MIN)*(2*X-MIN-DL)) if X<DL & rownumber<10001
replace se_m= MAX-sqrt((MAX+DU-2*X)*(MAX-DU)) if X>DU & rownumber<10001
drop MIN MAX DL DU X 
****1.2 sensitivity value for women********
gen MIN=0.35
gen MAX=0.40
gen DL=MIN+0.015
gen DU=MAX-0.015
gen X=(runiform()*(MAX+DU-DL-MIN) + (MIN+DL))/2 if rownumber<10001
gen se_f=. 
replace  se_f = X if DL <X & X<DU  & rownumber<10001
replace  se_f = MIN+sqrt((DL-MIN)*(2*X-MIN-DL)) if X<DL & rownumber<10001
replace se_f= MAX-sqrt((MAX+DU-2*X)*(MAX-DU)) if X>DU & rownumber<10001
drop MIN MAX DL DU X 
****1.3 specificity value for men**********
gen MIN=0.92 
gen MAX=1.00
gen DL=MIN+0.02
gen DU=MAX-0.02
gen X=(runiform()*(MAX+DU-DL-MIN) + (MIN+DL))/2 if rownumber<10001
gen sp_m=. 
replace  sp_m = X if DL <X & X<DU  & rownumber<10001
replace  sp_m = MIN+sqrt((DL-MIN)*(2*X-MIN-DL)) if X<DL & rownumber<10001
replace sp_m= MAX-sqrt((MAX+DU-2*X)*(MAX-DU)) if X>DU & rownumber<10001
drop MIN MAX DL DU X 

****1.4 specificity value for women**********
gen MIN=0.92 
gen MAX=1.00
gen DL=MIN+0.02
gen DU=MAX-0.02
gen X=(runiform()*(MAX+DU-DL-MIN) + (MIN+DL))/2 if rownumber<10001
gen sp_f=. 
replace  sp_f = X if DL <X & X<DU  & rownumber<10001
replace  sp_f = MIN+sqrt((DL-MIN)*(2*X-MIN-DL)) if X<DL & rownumber<10001
replace sp_f= MAX-sqrt((MAX+DU-2*X)*(MAX-DU)) if X>DU & rownumber<10001
drop MIN MAX DL DU X 


*******************************************
****Step 2*********************************
****PPV & NPV generation*******************
*******************************************
******************************************

*****Based observed data in the 2x2 table, we calculated the NPV and PPV at each time point from 1996-2004 for men and women separately***************
****ed, ued, eud, and ueud reprents cells in the distribution of self-reported loneliness (exposed [e] vs. unexposed [ue]) and outcome (disease [d] vs. un-disease [ud]) in the 2x2 table at each time pint from 1996-2004*********
program define npvppv
****Men, 1996********************************
gen ed = 161
gen ued= 1532 
gen eud= 159 
gen ueud= 1534 

gen ed_t=(ed-(1-sp_m)*(ed + ued))/(se_m-(1-sp_m)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_m)*(eud+ueud))/(se_m-(1-sp_m))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_m*ed_t)	
generate t_neg_d=(sp_m*ued_t)
generate f_neg_d=((1-se_m)*ed_t)	
generate f_pos_d=((1-sp_m)*ued_t)

generate t_pos_ud=(se_m*eud_t) 	
generate t_neg_ud=(sp_m*ueud_t)
generate f_neg_ud=((1-se_m)*eud_t)	
generate f_pos_ud=((1-sp_m)*ueud_t)

generate PPV_md1=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_md1=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_mud1=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_mud1=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
****Men, 1998********************************
gen ed = 249
gen ued= 1442 
gen eud= 179
gen ueud= 1516

gen ed_t=(ed-(1-sp_m)*(ed + ued))/(se_m-(1-sp_m)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_m)*(eud+ueud))/(se_m-(1-sp_m))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_m*ed_t)	
generate t_neg_d=(sp_m*ued_t)
generate f_neg_d=((1-se_m)*ed_t)	
generate f_pos_d=((1-sp_m)*ued_t)

generate t_pos_ud=(se_m*eud_t) 	
generate t_neg_ud=(sp_m*ueud_t)
generate f_neg_ud=((1-se_m)*eud_t)	
generate f_pos_ud=((1-sp_m)*ueud_t)

generate PPV_md2=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_md2=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_mud2=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_mud2=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
****Men, 2000********************************
gen ed = 226
gen ued= 1467 
gen eud= 159
gen ueud= 1534

gen ed_t=(ed-(1-sp_m)*(ed + ued))/(se_m-(1-sp_m)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_m)*(eud+ueud))/(se_m-(1-sp_m))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_m*ed_t)	
generate t_neg_d=(sp_m*ued_t)
generate f_neg_d=((1-se_m)*ed_t)	
generate f_pos_d=((1-sp_m)*ued_t)

generate t_pos_ud=(se_m*eud_t) 	
generate t_neg_ud=(sp_m*ueud_t)
generate f_neg_ud=((1-se_m)*eud_t)	
generate f_pos_ud=((1-sp_m)*ueud_t)

generate PPV_md3=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_md3=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_mud3=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_mud3=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
****Men, 2002********************************
gen ed = 264
gen ued= 1428
gen eud= 164
gen ueud= 1530

gen ed_t=(ed-(1-sp_m)*(ed + ued))/(se_m-(1-sp_m)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_m)*(eud+ueud))/(se_m-(1-sp_m))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_m*ed_t)	
generate t_neg_d=(sp_m*ued_t)
generate f_neg_d=((1-se_m)*ed_t)	
generate f_pos_d=((1-sp_m)*ued_t)

generate t_pos_ud=(se_m*eud_t) 	
generate t_neg_ud=(sp_m*ueud_t)
generate f_neg_ud=((1-se_m)*eud_t)	
generate f_pos_ud=((1-sp_m)*ueud_t)

generate PPV_md4=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_md4=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_mud4=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_mud4=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
****Men, 2004********************************
gen ed = 294
gen ued= 1396 
gen eud= 136
gen ueud= 1560

gen ed_t=(ed-(1-sp_m)*(ed + ued))/(se_m-(1-sp_m)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_m)*(eud+ueud))/(se_m-(1-sp_m))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_m*ed_t)	
generate t_neg_d=(sp_m*ued_t)
generate f_neg_d=((1-se_m)*ed_t)	
generate f_pos_d=((1-sp_m)*ued_t)

generate t_pos_ud=(se_m*eud_t) 	
generate t_neg_ud=(sp_m*ueud_t)
generate f_neg_ud=((1-se_m)*eud_t)	
generate f_pos_ud=((1-sp_m)*ueud_t)

generate PPV_md5=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_md5=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_mud5=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_mud5=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
******Women, 1996*********************
gen ed= 489
gen ued= 2332 
gen eud= 373
gen ueud= 2452

gen ed_t=(ed-(1-sp_f)*(ed + ued))/(se_f-(1-sp_f)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_f)*(eud+ueud))/(se_f-(1-sp_f))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_f*ed_t)	
generate t_neg_d=(sp_f*ued_t)
generate f_neg_d=((1-se_f)*ed_t)	
generate f_pos_d=((1-sp_f)*ued_t)

generate t_pos_ud=(se_f*eud_t) 	
generate t_neg_ud=(sp_f*ueud_t)
generate f_neg_ud=((1-se_f)*eud_t)	
generate f_pos_ud=((1-sp_f)*ueud_t)

generate PPV_fd1=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_fd1=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_fud1=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_fud1=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
******Women, 1998*********************
gen ed = 642
gen ued= 2186
gen eud= 401
gen ueud= 2417

gen ed_t=(ed-(1-sp_f)*(ed + ued))/(se_f-(1-sp_f)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_f)*(eud+ueud))/(se_f-(1-sp_f))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_f*ed_t)	
generate t_neg_d=(sp_f*ued_t)
generate f_neg_d=((1-se_f)*ed_t)	
generate f_pos_d=((1-sp_f)*ued_t)

generate t_pos_ud=(se_f*eud_t) 	
generate t_neg_ud=(sp_f*ueud_t)
generate f_neg_ud=((1-se_f)*eud_t)	
generate f_pos_ud=((1-sp_f)*ueud_t)

generate PPV_fd2=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_fd2=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_fud2=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_fud2=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
******Women, 2000*********************
gen ed = 672
gen ued= 2151
gen eud= 429
gen ueud= 2394

gen ed_t=(ed-(1-sp_f)*(ed + ued))/(se_f-(1-sp_f)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_f)*(eud+ueud))/(se_f-(1-sp_f))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_f*ed_t)	
generate t_neg_d=(sp_f*ued_t)
generate f_neg_d=((1-se_f)*ed_t)	
generate f_pos_d=((1-sp_f)*ued_t)

generate t_pos_ud=(se_f*eud_t) 	
generate t_neg_ud=(sp_f*ueud_t)
generate f_neg_ud=((1-se_f)*eud_t)	
generate f_pos_ud=((1-sp_f)*ueud_t)

generate PPV_fd3=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_fd3=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_fud3=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_fud3=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
******Women, 2002*********************
gen ed = 704
gen ued= 2117
gen eud= 379
gen ueud= 2446

gen ed_t=(ed-(1-sp_f)*(ed + ued))/(se_f-(1-sp_f)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_f)*(eud+ueud))/(se_f-(1-sp_f))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_f*ed_t)	
generate t_neg_d=(sp_f*ued_t)
generate f_neg_d=((1-se_f)*ed_t)	
generate f_pos_d=((1-sp_f)*ued_t)

generate t_pos_ud=(se_f*eud_t) 	
generate t_neg_ud=(sp_f*ueud_t)
generate f_neg_ud=((1-se_f)*eud_t)	
generate f_pos_ud=((1-sp_f)*ueud_t)

generate PPV_fd4=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_fd4=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_fud4=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_fud4=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*
******Women, 2004*********************
gen ed = 745
gen ued= 2063 
gen eud= 400
gen ueud= 2438

gen ed_t=(ed-(1-sp_f)*(ed + ued))/(se_f-(1-sp_f)) 
////cell A*
generate ued_t=(ed +ued)-ed_t 
///***cell B
generate eud_t=(eud-(1-sp_f)*(eud+ueud))/(se_f-(1-sp_f))
 ///Cell C//
generate ueud_t=(eud+ueud)-eud_t 
////cell D/////

generate t_pos_d=(se_f*ed_t)	
generate t_neg_d=(sp_f*ued_t)
generate f_neg_d=((1-se_f)*ed_t)	
generate f_pos_d=((1-sp_f)*ued_t)

generate t_pos_ud=(se_f*eud_t) 	
generate t_neg_ud=(sp_f*ueud_t)
generate f_neg_ud=((1-se_f)*eud_t)	
generate f_pos_ud=((1-sp_f)*ueud_t)

generate PPV_fd5=t_pos_d/(t_pos_d+f_pos_d)
generate NPV_fd5=t_neg_d/(t_neg_d+f_neg_d)
generate PPV_fud5=t_pos_ud/(t_pos_ud+f_pos_ud)
generate NPV_fud5=t_neg_ud/(t_neg_ud+f_neg_ud) 

drop e* t* f* u*

end

npvppv

*******************************************
****Step 3*********************************
****Record-level correction****************
****Repeat modeling analyses, N=10000******
*******************************************

program define correction
sort PPV_md1

gen threeway=. 
//////threeway interaction between loneliness (continous), years of follow-up, and gender****\
gen threeway1=.
gen threeway2=.
gen threeway3=.

//////threeway interaction between loneliness (categorical), years of follow-up, and gender****\

forvalues i=1/10000{

gen rflonenew=.
	

***among men*****
*****binomail trail by outcome status, i.e., 50% memory scores*****
replace rflonenew = (rbinomial(1,PPV_md1[`i'])) if rflone==1 & ragender==1 & year==1996 & disease==1
replace rflonenew = (rbinomial(1,PPV_md2[`i'])) if rflone==1 & ragender==1 & year==1998 & disease==1
replace rflonenew = (rbinomial(1,PPV_md3[`i'])) if rflone==1 & ragender==1 & year==2000 & disease==1
replace rflonenew = (rbinomial(1,PPV_md4[`i'])) if rflone==1 & ragender==1 & year==2002 & disease==1
replace rflonenew = (rbinomial(1,PPV_md5[`i'])) if rflone==1 & ragender==1 & year==2004 & disease==1

replace rflonenew = (rbinomial(1,PPV_mud1[`i'])) if rflone==1 & ragender==1 & year==1996 & disease==0
replace rflonenew = (rbinomial(1,PPV_mud2[`i'])) if rflone==1 & ragender==1 & year==1998 & disease==0
replace rflonenew = (rbinomial(1,PPV_mud3[`i'])) if rflone==1 & ragender==1 & year==2000 & disease==0
replace rflonenew = (rbinomial(1,PPV_mud4[`i'])) if rflone==1 & ragender==1 & year==2002 & disease==0
replace rflonenew = (rbinomial(1,PPV_mud5[`i'])) if rflone==1 & ragender==1 & year==2004 & disease==0


replace rflonenew = (rbinomial(1,1-NPV_md1[`i'])) if rflone==0 & ragender==1 & year==1996 & disease==1
replace rflonenew = (rbinomial(1,1-NPV_md2[`i'])) if rflone==0 & ragender==1 & year==1998 & disease==1
replace rflonenew = (rbinomial(1,1-NPV_md3[`i'])) if rflone==0 & ragender==1 & year==2000 & disease==1
replace rflonenew = (rbinomial(1,1-NPV_md4[`i'])) if rflone==0 & ragender==1 & year==2002 & disease==1
replace rflonenew = (rbinomial(1,1-NPV_md5[`i'])) if rflone==0 & ragender==1 & year==2004 & disease==1

replace rflonenew = (rbinomial(1,1-NPV_mud1[`i'])) if rflone==0 & ragender==1 & year==1996 & disease==0
replace rflonenew = (rbinomial(1,1-NPV_mud2[`i'])) if rflone==0 & ragender==1 & year==1998 & disease==0
replace rflonenew = (rbinomial(1,1-NPV_mud3[`i'])) if rflone==0 & ragender==1 & year==2000 & disease==0
replace rflonenew = (rbinomial(1,1-NPV_mud4[`i'])) if rflone==0 & ragender==1 & year==2002 & disease==0
replace rflonenew = (rbinomial(1,1-NPV_mud5[`i'])) if rflone==0 & ragender==1 & year==2004 & disease==0


*****among women*****
replace rflonenew = (rbinomial(1,PPV_fd1[`i'])) if rflone==1 & ragender==2 & year==1996 & disease==1
replace rflonenew = (rbinomial(1,PPV_fd2[`i'])) if rflone==1 & ragender==2 & year==1998 & disease==1
replace rflonenew = (rbinomial(1,PPV_fd3[`i'])) if rflone==1 & ragender==2 & year==2000 & disease==1
replace rflonenew = (rbinomial(1,PPV_fd4[`i'])) if rflone==1 & ragender==2 & year==2002 & disease==1
replace rflonenew = (rbinomial(1,PPV_fd5[`i'])) if rflone==1 & ragender==2 & year==2004 & disease==1

replace rflonenew = (rbinomial(1,PPV_fud1[`i'])) if rflone==1 & ragender==2 & year==1996 & disease==0
replace rflonenew = (rbinomial(1,PPV_fud2[`i'])) if rflone==1 & ragender==2 & year==1998 & disease==0
replace rflonenew = (rbinomial(1,PPV_fud3[`i'])) if rflone==1 & ragender==2 & year==2000 & disease==0
replace rflonenew = (rbinomial(1,PPV_fud4[`i'])) if rflone==1 & ragender==2 & year==2002 & disease==0
replace rflonenew = (rbinomial(1,PPV_fud5[`i'])) if rflone==1 & ragender==2 & year==2004 & disease==0


replace rflonenew = (rbinomial(1,1-NPV_fd1[`i'])) if rflone==0 & ragender==2 & year==1996 & disease==1
replace rflonenew = (rbinomial(1,1-NPV_fd2[`i'])) if rflone==0 & ragender==2 & year==1998 & disease==1
replace rflonenew = (rbinomial(1,1-NPV_fd3[`i'])) if rflone==0 & ragender==2 & year==2000 & disease==1
replace rflonenew = (rbinomial(1,1-NPV_fd4[`i'])) if rflone==0 & ragender==2 & year==2002 & disease==1
replace rflonenew = (rbinomial(1,1-NPV_fd5[`i'])) if rflone==0 & ragender==2 & year==2004 & disease==1

replace rflonenew = (rbinomial(1,1-NPV_fud1[`i'])) if rflone==0 & ragender==2 & year==1996 & disease==0
replace rflonenew = (rbinomial(1,1-NPV_fud2[`i'])) if rflone==0 & ragender==2 & year==1998 & disease==0
replace rflonenew = (rbinomial(1,1-NPV_fud3[`i'])) if rflone==0 & ragender==2 & year==2000 & disease==0
replace rflonenew = (rbinomial(1,1-NPV_fud4[`i'])) if rflone==0 & ragender==2 & year==2002 & disease==0
replace rflonenew = (rbinomial(1,1-NPV_fud5[`i'])) if rflone==0 & ragender==2 & year==2004 & disease==0


replace rflonenew=0 if year>2004
replace rflonenew=0 if year<1996
egen lone37new=sum(rflonenew), by(hhidpn)
replace lone37new=3 if lone37new>3 & lone37new<6
replace lone37new=. if lone37c==.

****coding loneliness as continous variable**************


****pooled analyses-continous*********
mixed  zmemimp c.year_center##c.lone37new##i.ragender c.year_center#c.year_center  bage_center  i.race  i.raeduc i.bmarfi  i.brwork i.bwealthfi  bisolationfi  bcesd bradla  if year>2002 & bage>49 & lone37c<. || hhidpn: c.year_center , var cov(un) reml
local threeway=el(r(table),1,11)
display `threeway'
replace threeway=`threeway' in `i'


******coding loneliness as categorical variable****************
mixed  zmemimp c.year_center##i.lone37new##i.ragender c.year_center#c.year_center  bage_center  i.race  i.raeduc i.bmarfi  i.brwork i.bwealthfi  bisolationfi  bcesd bradla  if year>2002 & bage>49 & lone37c<.  || hhidpn: c.year_center , var cov(un) reml
local threeway1=el(r(table),1,25)
display `threeway1'
local threeway2=el(r(table),1,27)
display `threeway2'
local threeway3=el(r(table),1,29)
display `threeway3'
replace threeway1=`threeway1' in `i'
replace threeway2=`threeway2' in `i'
replace threeway3=`threeway3' in `i'



drop lone37new rflonenew

}

*Using standard error estimates from conventional analysis to calculate simulation interval incorporating random error
replace threeway= threeway-(.0024931 * rnormal(0,1))

replace threeway1=threeway1-(.0063986 * rnormal(0,1))
replace threeway2=threeway2-(.0089531   * rnormal(0,1))
replace threeway3=threeway3-(.0088428 * rnormal(0,1))

end   

correction





