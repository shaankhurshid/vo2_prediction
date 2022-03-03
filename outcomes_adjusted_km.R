# Depends
library(data.table)
library(plyr)
library(stringr)
library(prodlim)
library(Cairo)
library(cmprsk)
library(glmnet)
library(ggplot2)

### LINKING
# Load models for inference
outcomes_set <- fread(file='01-24-2022_c3po_pclr_infer_MGB-c3po_last_stfu_ecg.tsv')
outcomes_mgh <- outcomes_set[s3_path=='ecg_mgh_hd5s']
outcomes_bwh <- outcomes_set[s3_path=='ecg_bwh_hd5s']

# Load phenotype files
wide <- fread(file='hf-wide-2021-05-08.tsv')
af <- fread(file='wide_with_af_051021.csv')
mi <- fread(file='wide_with_mi_stroke_051021.csv')

# Load linker
linker <- fread(file='mgh_bwh_mrn_empi_linker.txt')
bwh_linker <- linker[Hospital=='BWH']
mgh_linker <- linker[Hospital=='MGH']

# Link ECGs to wide file
## BWH
setkey(bwh_linker,MRN); setkey(outcomes_bwh,sample_id)
outcomes_bwh[bwh_linker,linker_id := i.LINKER_ID]

## MGH
setkey(mgh_linker,MRN); setkey(outcomes_mgh,sample_id)
outcomes_mgh[mgh_linker,linker_id := i.LINKER_ID]

outcomes_set <- rbind(outcomes_bwh,outcomes_mgh) # 520868 - 401007 = 119861

### Variable ascertainment and cleanup
# Cast continuous variables
numerics <- c('start_fu','last_encounter','Dem.Date_Of_Death_Age.no_filter','mi_age')
for (j in numerics){set(outcomes_set,j=j,value=as.numeric(str_extract(outcomes_set[[j]],'\\d+')))}

# Add AF and MI outcomes
setkey(mi,id); setkey(af,id); setkey(outcomes_set,linker_id)
outcomes_set[mi,':='(mi_age = i.mi_age, global_censor_age = i.global_censor_age)]
outcomes_set[af,af_age := i.af_age]

#### COVARIATES

# BMI
# Recover a few extra heights
wide[str_detect(start_fu_Height,"\\d+\\'\\d+"),
     start_fu_Height := as.numeric(str_extract(start_fu_Height,"^\\d+"))*5+
       as.numeric(str_extract(start_fu_Height,"\\d+$"))]
wide[str_detect(start_fu_Height,"\\d+'"),
     start_fu_Height := as.numeric(str_extract(start_fu_Height,"^\\d+"))*5]

# Cast continuous variables
numerics <- c('start_fu_Height','start_fu_Weight')
for (j in numerics){set(wide,j=j,value=as.numeric(wide[[j]]))}

# Harmonize units and QC
## Wt
wide[start_fu_Weight_units %in% c('pound','Pounds','Lbs/Oz',''),wt_kg := start_fu_Weight/2.2] # Pounds assumption if not kg
wide[start_fu_Weight_units %in% c('Kilograms','kg'),wt_kg := start_fu_Weight]
wide[c(wt_kg < 20 | wt_kg > 450), wt_kg := NA] # ultimately 127988 NAs
wide[,id := as.numeric(id)]

## Ht
wide[c(start_fu_Height_units %in% c('Inches','inch')),ht_cm := as.numeric(start_fu_Height)*2.54]
wide[c((start_fu_Height_units=='') | is.na(start_fu_Height_units)),ht_cm := as.numeric(str_detect(start_fu_Height,'\\d+'))*2.54]
wide[start_fu_Height_units %in% c('Centimeters','cm'),ht_cm := as.numeric(start_fu_Height)]
wide[start_fu_Height_units=='Feet',ht_cm := ((as.numeric(start_fu_Height)/12)*2.54)]
wide[c(ht_cm < 91 | ht_cm > 305), ht_cm := NA] # ultimately 76409 NAs remain
wide[,id := as.numeric(id)]

## Remove outliers in BMI prior to ht/wt recovery
wide[,':='(start_fu_BMI_clean = ifelse(start_fu_BMI < 15 | start_fu_BMI > 120,NA,start_fu_BMI))]

# Recover ht/wt from BMI where possible
wide[c(is.na(wt_kg) & !is.na(start_fu_BMI_clean) & !is.na(ht_cm)), wt_kg := start_fu_BMI_clean*(ht_cm/100)**2]
wide[c(is.na(ht_cm) & !is.na(start_fu_BMI_clean) & !is.na(wt_kg)), ht_cm := sqrt(wt_kg/start_fu_BMI_clean)*100]

# Recover BMI from ht/wt where possible
wide[c(!is.na(ht_cm) & is.na(start_fu_BMI_clean) & !is.na(wt_kg)), start_fu_BMI_clean := wt_kg/(ht_cm/100)**2]

# Reapply QC
wide[c(ht_cm < 91 | ht_cm > 305), ht_cm := NA] 
wide[c(wt_kg < 20 | wt_kg > 450), wt_kg := NA] 
wide[c(start_fu_BMI_clean < 15 | start_fu_BMI_clean > 120), start_fu_BMI_clean := NA] 

# Convert units for PCLR model
wide[,':='(ht_in = ht_cm*0.393701,wt_lb = wt_kg/2.2)]

# Cast ages
numerics <- c('start_fu','last_encounter','Dem.Date_Of_Death_Age.no_filter',
              'inpatient_hf_age')
for (j in numerics){set(wide,j=j,value=as.numeric(str_extract(wide[[j]],'\\d+')))}

# Remove missing start_fu
wide <- wide[!is.na(start_fu)] #520868

# Remove missing sex
wide <- wide[Dem.Gender.no_filter != ''] #520868 - 16 = 520852

# Merge variables
setkey(wide,id); setkey(outcomes_set,linker_id)
outcomes_set[wide,':='(age = i.start_fu/365.25, sex = ifelse(i.Dem.Gender.no_filter == 'Female',1,0),
                       ht = i.ht_in, wt = i.wt_lb, bmi = i.start_fu_BMI_clean,
                       start_fu = i.start_fu, last_encounter = i.last_encounter,
                       hf_age = i.inpatient_hf_age, death_age = i.Dem.Date_Of_Death_Age.no_filter)]

# Filter to complete data
outcomes_set <- outcomes_set[!is.na(bmi) & !is.na(ht) & !is.na(wt)] # 119861 - 34184 = 85677

# Remove in CPET dataset
not_holdout <- fread(file='not_holdout_021822.csv')
setkey(not_holdout,MRN)
not_holdout[mgh_linker,linker_id := i.LINKER_ID]
outcomes_set <- outcomes_set[!(linker_id %in% not_holdout$linker_id)] # 85677 - 252 = 85425

# Remove missing last encounter
outcomes_set <- outcomes_set[!is.na(last_encounter)] #85425 - 707 = 84718

#### OUTCOMES
## MI
outcomes_set[,':='(mi = ifelse(!is.na(mi_age),1,0))]
outcomes_set[,':='(prev_mi = ifelse(c(mi==1 & (mi_age <= start_fu)),1,0))]
outcomes_set[,':='(incd_mi = ifelse(mi==1 & (mi_age > start_fu) & !is.na(global_censor_age) & (mi_age <= global_censor_age),1,0))]
outcomes_set[,':='(time_to_mi = ifelse(mi == 1,pmin(mi_age - start_fu,global_censor_age - start_fu)/365.25,
                                       pmin((last_encounter - start_fu),(global_censor_age - start_fu),
                                            (Dem.Date_Of_Death_Age.no_filter-start_fu),na.rm=T)/365.25))]
## AF
outcomes_set[,':='(af = ifelse(!is.na(af_age),1,0))]
outcomes_set[,':='(prev_af = ifelse(c(af==1 & (af_age <= start_fu)),1,0))]
outcomes_set[,':='(incd_af = ifelse(af==1 & (af_age > start_fu) & !is.na(global_censor_age) & (af_age <= global_censor_age),1,0))]
outcomes_set[,':='(time_to_af = ifelse(af == 1,pmin(af_age - start_fu,global_censor_age - start_fu)/365.25,
                                       pmin((last_encounter - start_fu),(global_censor_age - start_fu),
                                            (Dem.Date_Of_Death_Age.no_filter-start_fu),na.rm=T)/365.25))]
## HF
outcomes_set[,':='(hf = ifelse(!is.na(hf_age),1,0))]
outcomes_set[,':='(prev_hf = ifelse(c(hf==1 & (hf_age <= start_fu)),1,0))]
outcomes_set[,':='(incd_hf = ifelse(hf==1 & (hf_age > start_fu) & !is.na(global_censor_age) & (hf_age <= global_censor_age),1,0))]
outcomes_set[,':='(time_to_hf = ifelse(hf == 1,pmin(hf_age - start_fu,global_censor_age - start_fu)/365.25,
                                       pmin((last_encounter - start_fu),(global_censor_age - start_fu),
                                            (Dem.Date_Of_Death_Age.no_filter-start_fu),na.rm=T)/365.25))]
## Death
outcomes_set[,':='(death = ifelse(!is.na(death_age),1,0))]
outcomes_set[,':='(prev_death = ifelse(c(death==1 & (death_age <= start_fu)),1,0))]
outcomes_set[,':='(incd_death = ifelse(death==1 & (death_age > start_fu),1,0))]
outcomes_set[,':='(incd_death_before_censor = ifelse(death==1 & (death_age > start_fu) & !is.na(global_censor_age) & (death_age <= global_censor_age),1,0))]
outcomes_set[,':='(time_to_death = ifelse(death == 1,(death_age - start_fu)/365.25,
                                          pmin((last_encounter - start_fu),(global_censor_age - start_fu),
                                               (death_age-start_fu),na.rm=T)/365.25))]

#### INFERENCE
load('lm_pclr_basic.RData')
load('cv.lm_pclr_basic.RData')

# Generate X matrices
# Define variable space
## Impute TT bike for all
outcomes_set[,':='(tt_bike = 1, tt_row = 0, tt_tread = 0)]

## Prepare for scaling
basic_covars <- c('age','sex','bmi','tt_bike','tt_tread','tt_row')
pclr_only <- names(outcomes_set)[str_detect(names(outcomes_set),'c3po_pclr')]  
pclr_basic <- c(basic_covars,pclr_only)

# Standardize continuous variables for shrinkage
## Don't scale the all zero columns (generates NAs)
pclr_zero <- outcomes_set[,lapply(.SD,mean),.SDcols=pclr_only]
pclr_zero <- names(pclr_zero)[pclr_zero == 0]

## Now scale
for (j in names(outcomes_set)[c((names(outcomes_set) %in% pclr_basic[!(pclr_basic %in% c('sex','tt_bike','tt_tread','tt_row'))])
                                & (!names(outcomes_set) %in% pclr_zero))]){
  set(outcomes_set,j=j,value=scale(outcomes_set[[j]]))
}

# Generate reduced X matrix
x_pclr_basic = data.matrix(outcomes_set[,.SD,.SDcols=pclr_basic])

# Now infer
outcomes_set[,pclr_basic_fitted := predict(lm_pclr_basic,cv.lm_pclr_basic$lambda.min,newx=x_pclr_basic,type='link')]

# Sets
af_set <- outcomes_set[prev_af==0 & time_to_af > 0]
mi_set <- outcomes_set[prev_mi==0 & time_to_mi > 0]
hf_set <- outcomes_set[prev_hf==0 & time_to_hf > 0]
death_set <- outcomes_set[prev_death==0 & time_to_death > 0]

# Standardize exposure variables
af_set[,':='(vo2_std = (pclr_basic_fitted - mean(pclr_basic_fitted))/sd(pclr_basic_fitted),
             vo2_less14 = ifelse(pclr_basic_fitted < 14,1,0),
             vo2_less10pct = ifelse(pclr_basic_fitted < quantile(pclr_basic_fitted,0.10),1,0))]
mi_set[,':='(vo2_std = (pclr_basic_fitted - mean(pclr_basic_fitted))/sd(pclr_basic_fitted),
             vo2_less14 = ifelse(pclr_basic_fitted < 14,1,0),
             vo2_less10pct = ifelse(pclr_basic_fitted < quantile(pclr_basic_fitted,0.10),1,0))]
hf_set[,':='(vo2_std = (pclr_basic_fitted - mean(pclr_basic_fitted))/sd(pclr_basic_fitted),
             vo2_less14 = ifelse(pclr_basic_fitted < 14,1,0),
             vo2_less10pct = ifelse(pclr_basic_fitted < quantile(pclr_basic_fitted,0.10),1,0))]
death_set[,':='(vo2_std = (pclr_basic_fitted - mean(pclr_basic_fitted))/sd(pclr_basic_fitted),
                vo2_less14 = ifelse(pclr_basic_fitted < 14,1,0),
                vo2_less10pct = ifelse(pclr_basic_fitted < quantile(pclr_basic_fitted,0.10),1,0))]

# Plot density of true VO2
ggplot() + geom_density(data=death_set,aes(x=death_set$pclr_basic_fitted),fill="#8c96c6",alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,70,10),expand=c(0,0),limits=c(0,70)) +
  scale_y_continuous(breaks=seq(0,0.05,0.01),expand=c(0,0),limits=c(0,0.05)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted Peak VO2 (mL/kg/min)',y='Density') 
ggsave('mgb_predicted_vo2_density.pdf',
       height=2,width=2.5,units='in',scale=4)

# Models
af_set2 <- survSplit(Surv(time_to_af,incd_af) ~ ., data=af_set, cut=10,episode='tgroup')
mod_af <- coxph(Surv(time_to_af,incd_af) ~ vo2_std + age:strata(tgroup) + sex + bmi,data=af_set2)
mod_af2 <- coxph(Surv(time_to_af,incd_af) ~ vo2_less14 + age:strata(tgroup) + sex + bmi,data=af_set2)
mod_af3 <- coxph(Surv(time_to_af,incd_af) ~ vo2_less10pct + age:strata(tgroup) + sex + bmi,data=af_set2)

mod_mi <- coxph(Surv(time_to_mi,incd_mi) ~ vo2_std + age + sex + bmi,data=mi_set)
mod_mi2 <- coxph(Surv(time_to_mi,incd_mi) ~ vo2_less14 + age + sex + bmi,data=mi_set)
mod_mi3 <- coxph(Surv(time_to_mi,incd_mi) ~ vo2_less10pct + age + sex + bmi,data=mi_set)

hf_set2 <- survSplit(Surv(time_to_hf,incd_hf) ~ ., data=hf_set, cut=c(1,2,3,4),episode='tgroup')
mod_hf <- coxph(Surv(time_to_hf,incd_hf) ~ vo2_std + age:strata(tgroup) + sex + bmi,data=hf_set2)
mod_hf2 <- coxph(Surv(time_to_hf,incd_hf) ~ vo2_less14 + age:strata(tgroup) + sex + bmi,data=hf_set2)
mod_hf3 <- coxph(Surv(time_to_hf,incd_hf) ~ vo2_less10pct + age:strata(tgroup) + sex + bmi,data=hf_set2)

death_set2 <- survSplit(Surv(time_to_death,incd_death) ~ ., data=death_set,cut=seq(2,18,2),episode='tgroup')
mod_death <- coxph(Surv(time_to_death,incd_death) ~ vo2_std + age:strata(tgroup) + sex + bmi:strata(tgroup),data=death_set2)
mod_death2 <- coxph(Surv(time_to_death,incd_death) ~ vo2_less14 + age:strata(tgroup) + sex + bmi:strata(tgroup),data=death_set2)
mod_death3 <- coxph(Surv(time_to_death,incd_death) ~ vo2_less10pct + age:strata(tgroup) + sex + bmi:strata(tgroup),data=death_set2)

# KM curves
## AF
af_set[,vo2_graphical := factor(ifelse(vo2_less10pct==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_af <- prodlim(Hist(time_to_af,incd_af)~vo2_graphical,data=af_set)

# Plot
CairoPDF(file='km_af_less10.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_af,"cuminc",ylim=c(0,0.40),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.40,0.10),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c('#ef3b2c','#02818a'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.40*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Below 10th percentile","At or above 10th percentile"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.8,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=seq(0,10,2), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.2,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.8) # descriptor for N at risk
dev.off()

## MI
mi_set[,vo2_graphical := factor(ifelse(vo2_less10pct==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_mi <- prodlim(Hist(time_to_mi,incd_mi)~vo2_graphical,data=mi_set)

# Plot
CairoPDF(file='km_mi_less10.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_mi,"cuminc",ylim=c(0,0.20),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.20,0.05),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c('#ef3b2c','#02818a'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.20*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Below 10th percentile","At or above 10th percentile"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.8,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=seq(0,10,2), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.1,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.8) # descriptor for N at risk
dev.off()

## HF
hf_set[,vo2_graphical := factor(ifelse(vo2_less10pct==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_hf <- prodlim(Hist(time_to_hf,incd_hf)~vo2_graphical,data=hf_set)

# Plot
CairoPDF(file='km_hf_less10.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_hf,"cuminc",ylim=c(0,0.40),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.40,0.10),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c('#ef3b2c','#02818a'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.40*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Below 10th percentile","At or above 10th percentile"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.8,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=seq(0,10,2), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.2,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.8) # descriptor for N at risk
dev.off()

## Death
death_set[,vo2_graphical := factor(ifelse(vo2_less10pct==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_death <- prodlim(Hist(time_to_death,incd_death)~vo2_graphical,data=death_set)

# Plot
CairoPDF(file='km_death_less10.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_death,"cuminc",ylim=c(0,0.40),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.40,0.10),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c('#ef3b2c','#02818a'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.40*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Below 10th percentile","At or above 10th percentile"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.8,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=seq(0,10,2), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.2,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.8) # descriptor for N at risk
dev.off()

# KM curves
## AF
af_set[,vo2_graphical := factor(ifelse(vo2_less14==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_af <- prodlim(Hist(time_to_af,incd_af)~vo2_graphical,data=af_set)

# Plot
CairoPDF(file='km_af_less14.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_af,"cuminc",ylim=c(0,0.40),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.40,0.10),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c('#ef3b2c','#02818a'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.40*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.8,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=seq(0,10,2), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.2,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.8) # descriptor for N at risk
dev.off()

## MI
mi_set[,vo2_graphical := factor(ifelse(vo2_less14==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_mi <- prodlim(Hist(time_to_mi,incd_mi)~vo2_graphical,data=mi_set)

# Plot
CairoPDF(file='km_mi_less14.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_mi,"cuminc",ylim=c(0,0.20),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.20,0.05),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c('#ef3b2c','#02818a'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.20*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.8,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=seq(0,10,2), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.1,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.8) # descriptor for N at risk
dev.off()

## HF
hf_set[,vo2_graphical := factor(ifelse(vo2_less14==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_hf <- prodlim(Hist(time_to_hf,incd_hf)~vo2_graphical,data=hf_set)

# Plot
CairoPDF(file='km_hf_less14.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_hf,"cuminc",ylim=c(0,0.40),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.40,0.10),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c('#ef3b2c','#02818a'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.40*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.8,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=seq(0,10,2), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.2,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.8) # descriptor for N at risk
dev.off()

## Death
death_set[,vo2_graphical := factor(ifelse(vo2_less14==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_death <- prodlim(Hist(time_to_death,incd_death)~vo2_graphical,data=death_set)

# Plot
CairoPDF(file='km_death_less14.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_death,"cuminc",ylim=c(0,0.40),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.40,0.10),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c('#ef3b2c','#02818a'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.40*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.8,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=seq(0,10,2), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.2,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.8) # descriptor for N at risk
dev.off()

### Death adjust
## AF
# Generate status variable (1 = AF, 2 = died, 0 = censored)
af_set[,died_without_af := ifelse(incd_death_before_censor==1 & incd_af==0,1,0)]
af_set[,failure := ifelse(incd_af==1,1,
                          ifelse(died_without_af==1,2,0))]

vo2_af <- boot_cmp(time='time_to_af',status='failure',data=af_set,runs=200,exposure='vo2_std')
x_matrix <- model.matrix(~age + sex + bmi + vo2_std,data=af_set)[,-1]
hr_vo2_af <- crr(ftime=af_set$time_to_af,fstatus=af_set$failure,
                 failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_af <- c(hr_vo2_af,hr_vo2_af+1.96*sd(vo2_af),hr_vo2_af-1.96*sd(vo2_af))
ci_vo2_af <- exp(-ci_vo2_af)

vo2_less14_af <- boot_cmp(time='time_to_af',status='failure',data=af_set,runs=200,exposure='vo2_less14')
x_matrix <- model.matrix(~age + sex + bmi + vo2_less14,data=af_set)[,-1]
hr_vo2_less14_af <- crr(ftime=af_set$time_to_af,fstatus=af_set$failure,
                        failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_less14_af <- c(hr_vo2_less14_af,hr_vo2_less14_af-1.96*sd(vo2_less14_af),hr_vo2_less14_af+1.96*sd(vo2_less14_af))
ci_vo2_less14_af <- exp(ci_vo2_less14_af)

vo2_less10pct_af <- boot_cmp(time='time_to_af',status='failure',data=af_set,runs=200,exposure='vo2_less10pct')
x_matrix <- model.matrix(~age + sex + bmi + vo2_less10pct,data=af_set)[,-1]
hr_vo2_less10pct_af <- crr(ftime=af_set$time_to_af,fstatus=af_set$failure,
                           failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_less10pct_af <- c(hr_vo2_less10pct_af,hr_vo2_less10pct_af-1.96*sd(vo2_less10pct_af),hr_vo2_less10pct_af+1.96*sd(vo2_less10pct_af))
ci_vo2_less10pct_af <- exp(ci_vo2_less10pct_af)

## MI
# Generate status variable (1 = AF, 2 = died, 0 = censored)
mi_set[,died_without_mi := ifelse(incd_death_before_censor==1 & incd_mi==0,1,0)]
mi_set[,failure := ifelse(incd_mi==1,1,
                          ifelse(died_without_mi==1,2,0))]

vo2_mi <- boot_cmp(time='time_to_mi',status='failure',data=mi_set,runs=200,exposure='vo2_std')
x_matrix <- model.matrix(~age + sex + bmi + vo2_std,data=mi_set)[,-1]
hr_vo2_mi <- crr(ftime=mi_set$time_to_mi,fstatus=mi_set$failure,
                 failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_mi <- c(hr_vo2_mi,hr_vo2_mi+1.96*sd(vo2_mi),hr_vo2_mi-1.96*sd(vo2_mi))
ci_vo2_mi <- exp(-ci_vo2_mi)

vo2_less14_mi <- boot_cmp(time='time_to_mi',status='failure',data=mi_set,runs=200,exposure='vo2_less14')
x_matrix <- model.matrix(~age + sex + bmi + vo2_less14,data=mi_set)[,-1]
hr_vo2_less14_mi <- crr(ftime=mi_set$time_to_mi,fstatus=mi_set$failure,
                        failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_less14_mi <- c(hr_vo2_less14_mi,hr_vo2_less14_mi-1.96*sd(vo2_less14_mi),hr_vo2_less14_mi+1.96*sd(vo2_less14_mi))
ci_vo2_less14_mi <- exp(ci_vo2_less14_mi)

vo2_less10pct_mi <- boot_cmp(time='time_to_mi',status='failure',data=mi_set,runs=200,exposure='vo2_less10pct')
x_matrix <- model.matrix(~age + sex + bmi + vo2_less10pct,data=mi_set)[,-1]
hr_vo2_less10pct_mi <- crr(ftime=mi_set$time_to_mi,fstatus=mi_set$failure,
                           failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_less10pct_mi <- c(hr_vo2_less10pct_mi,hr_vo2_less10pct_mi-1.96*sd(vo2_less10pct_mi),hr_vo2_less10pct_mi+1.96*sd(vo2_less10pct_mi))
ci_vo2_less10pct_mi <- exp(ci_vo2_less10pct_mi)

## AF
# Generate status variable (1 = AF, 2 = died, 0 = censored)
hf_set[,died_without_hf := ifelse(incd_death_before_censor==1 & incd_hf==0,1,0)]
hf_set[,failure := ifelse(incd_hf==1,1,
                          ifelse(died_without_hf==1,2,0))]

vo2_hf <- boot_cmp(time='time_to_hf',status='failure',data=hf_set,runs=200,exposure='vo2_std')
x_matrix <- model.matrix(~age + sex + bmi + vo2_std,data=hf_set)[,-1]
hr_vo2_hf <- crr(ftime=hf_set$time_to_hf,fstatus=hf_set$failure,
                 failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_hf <- c(hr_vo2_hf,hr_vo2_hf+1.96*sd(vo2_hf),hr_vo2_hf-1.96*sd(vo2_hf))
ci_vo2_hf <- exp(-ci_vo2_hf)

vo2_less14_hf <- boot_cmp(time='time_to_hf',status='failure',data=hf_set,runs=200,exposure='vo2_less14')
x_matrix <- model.matrix(~age + sex + bmi + vo2_less14,data=hf_set)[,-1]
hr_vo2_less14_hf <- crr(ftime=hf_set$time_to_hf,fstatus=hf_set$failure,
                        failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_less14_hf <- c(hr_vo2_less14_hf,hr_vo2_less14_hf-1.96*sd(vo2_less14_hf),hr_vo2_less14_hf+1.96*sd(vo2_less14_hf))
ci_vo2_less14_hf <- exp(ci_vo2_less14_hf)

vo2_less10pct_hr <- boot_cmp(time='time_to_hf',status='failure',data=hf_set,runs=200,exposure='vo2_less10pct')
x_matrix <- model.matrix(~age + sex + bmi + vo2_less10pct,data=hf_set)[,-1]
hr_vo2_less10pct_hf <- crr(ftime=hf_set$time_to_hf,fstatus=hf_set$failure,
                           failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
ci_vo2_less10pct_hf <- c(hr_vo2_less10pct_hf,hr_vo2_less10pct_hf-1.96*sd(vo2_less10pct_hf),hr_vo2_less10pct_hf+1.96*sd(vo2_less10pct_hf))
ci_vo2_less10pct_hf <- exp(ci_vo2_less10pct_hf)

### ADJUSTED KMs
## AF
# Relevel
af_set[,vo2_graphical := factor(ifelse(vo2_less10pct==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_af <- coxph(Surv(time_to_af,incd_af) ~ strata(vo2_graphical) + age + sex + bmi,data=af_set)

weights <- data.frame(vo2_graphical = levels(af_set$vo2_graphical),
                      age = rep(mean(af_set[sex==1]$age),2),
                      sex = rep(1,2),
                      bmi = rep(mean(af_set[sex==1]$bmi),2))

# Plot
fit <- survfit(mod_af,newdata=weights)

CairoPDF(file='km_af_less10_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 10th percentile","At or above 10th percentile"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# MALE
weights2 <- data.frame(vo2_graphical = levels(af_set$vo2_graphical),
                       age = rep(mean(af_set[sex==0]$age),2),
                       sex = rep(0,2),
                       bmi = rep(mean(af_set[sex==0]$bmi),2))

# Plot
fit <- survfit(mod_af,newdata=weights2)

CairoPDF(file='km_af_less10_adjusted_m.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 10th percentile","At or above 10th percentile"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

## MI
# Relevel
mi_set[,vo2_graphical := factor(ifelse(vo2_less10pct==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_mi <- coxph(Surv(time_to_mi,incd_mi) ~ strata(vo2_graphical) + age + sex + bmi,data=mi_set)

weights2 <- data.frame(vo2_graphical = levels(mi_set$vo2_graphical),
                       age = rep(mean(af_set[sex==1]$age),2),
                       sex = rep(1,2),
                       bmi = rep(mean(mi_set[sex==1]$bmi),2))

# Plot
fit <- survfit(mod_mi,newdata=weights)

CairoPDF(file='km_mi_less10_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 10th percentile","At or above 10th percentile"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# MALE
weights2 <- data.frame(vo2_graphical = levels(mi_set$vo2_graphical),
                       age = rep(mean(mi_set[sex==0]$age),2),
                       sex = rep(0,2),
                       bmi = rep(mean(mi_set[sex==0]$bmi),2))

# Plot
fit <- survfit(mod_mi,newdata=weights2)

CairoPDF(file='km_mi_less10_adjusted_m.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 10th percentile","At or above 10th percentile"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

## HF
# Relevel
hf_set[,vo2_graphical := factor(ifelse(vo2_less10pct==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_hf <- coxph(Surv(time_to_hf,incd_hf) ~ strata(vo2_graphical) + age + sex + bmi,data=hf_set)

weights2 <- data.frame(vo2_graphical = levels(hf_set$vo2_graphical),
                       age = rep(mean(af_set[sex==1]$age),2),
                       sex = rep(1,2),
                       bmi = rep(mean(hf_set[sex==1]$bmi),2))

# Plot
fit <- survfit(mod_hf,newdata=weights)

CairoPDF(file='km_hf_less10_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 10th percentile","At or above 10th percentile"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# MALE
weights2 <- data.frame(vo2_graphical = levels(hf_set$vo2_graphical),
                       age = rep(mean(hf_set[sex==0]$age),2),
                       sex = rep(0,2),
                       bmi = rep(mean(hf_set[sex==0]$bmi),2))

# Plot
fit <- survfit(mod_hf,newdata=weights2)

CairoPDF(file='km_hf_less10_adjusted_m.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 10th percentile","At or above 10th percentile"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

## Death
# Relevel
death_set[,vo2_graphical := factor(ifelse(vo2_less10pct==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_death <- coxph(Surv(time_to_death,incd_death) ~ strata(vo2_graphical) + age + sex + bmi,data=death_set)

weights2 <- data.frame(vo2_graphical = levels(death_set$vo2_graphical),
                       age = rep(mean(af_set[sex==1]$age),2),
                       sex = rep(1,2),
                       bmi = rep(mean(death_set[sex==1]$bmi),2))

# Plot
fit <- survfit(mod_death,newdata=weights)

CairoPDF(file='km_death_less10_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 10th percentile","At or above 10th percentile"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# MALE
weights2 <- data.frame(vo2_graphical = levels(death_set$vo2_graphical),
                       age = rep(mean(death_set[sex==0]$age),2),
                       sex = rep(0,2),
                       bmi = rep(mean(death_set[sex==0]$bmi),2))

# Plot
fit <- survfit(mod_death,newdata=weights2)

CairoPDF(file='km_death_less10_adjusted_m.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 10th percentile","At or above 10th percentile"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

### ADJUSTED KMs
## AF
# Relevel
af_set[,vo2_graphical := factor(ifelse(vo2_less14==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_af <- coxph(Surv(time_to_af,incd_af) ~ strata(vo2_graphical) + age + sex + bmi,data=af_set)

weights <- data.frame(vo2_graphical = levels(af_set$vo2_graphical),
                      age = rep(mean(af_set[sex==1]$age),2),
                      sex = rep(1,2),
                      bmi = rep(mean(af_set[sex==1]$bmi),2))

# Plot
fit <- survfit(mod_af,newdata=weights)

CairoPDF(file='km_af_less14_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# MALE
weights2 <- data.frame(vo2_graphical = levels(af_set$vo2_graphical),
                       age = rep(mean(af_set[sex==0]$age),2),
                       sex = rep(0,2),
                       bmi = rep(mean(af_set[sex==0]$bmi),2))

# Plot
fit <- survfit(mod_af,newdata=weights2)

CairoPDF(file='km_af_less14_adjusted_m.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

## MI
# Relevel
mi_set[,vo2_graphical := factor(ifelse(vo2_less14==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_mi <- coxph(Surv(time_to_mi,incd_mi) ~ strata(vo2_graphical) + age + sex + bmi,data=mi_set)

weights2 <- data.frame(vo2_graphical = levels(mi_set$vo2_graphical),
                       age = rep(mean(af_set[sex==1]$age),2),
                       sex = rep(1,2),
                       bmi = rep(mean(mi_set[sex==1]$bmi),2))

# Plot
fit <- survfit(mod_mi,newdata=weights)

CairoPDF(file='km_mi_less14_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# MALE
weights2 <- data.frame(vo2_graphical = levels(mi_set$vo2_graphical),
                       age = rep(mean(mi_set[sex==0]$age),2),
                       sex = rep(0,2),
                       bmi = rep(mean(mi_set[sex==0]$bmi),2))

# Plot
fit <- survfit(mod_mi,newdata=weights2)

CairoPDF(file='km_mi_less14_adjusted_m.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

## HF
# Relevel
hf_set[,vo2_graphical := factor(ifelse(vo2_less14==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_hf <- coxph(Surv(time_to_hf,incd_hf) ~ strata(vo2_graphical) + age + sex + bmi,data=hf_set)

weights2 <- data.frame(vo2_graphical = levels(hf_set$vo2_graphical),
                       age = rep(mean(af_set[sex==1]$age),2),
                       sex = rep(1,2),
                       bmi = rep(mean(hf_set[sex==1]$bmi),2))

# Plot
fit <- survfit(mod_hf,newdata=weights)

CairoPDF(file='km_hf_less14_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# MALE
weights2 <- data.frame(vo2_graphical = levels(hf_set$vo2_graphical),
                       age = rep(mean(hf_set[sex==0]$age),2),
                       sex = rep(0,2),
                       bmi = rep(mean(hf_set[sex==0]$bmi),2))

# Plot
fit <- survfit(mod_hf,newdata=weights2)

CairoPDF(file='km_hf_less14_adjusted_m.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

## Death
# Relevel
death_set[,vo2_graphical := factor(ifelse(vo2_less14==1,'Below','At/Above'),levels=c('Below','At/Above'))]
mod_death <- coxph(Surv(time_to_death,incd_death) ~ strata(vo2_graphical) + age + sex + bmi,data=death_set)

weights2 <- data.frame(vo2_graphical = levels(death_set$vo2_graphical),
                       age = rep(mean(af_set[sex==1]$age),2),
                       sex = rep(1,2),
                       bmi = rep(mean(death_set[sex==1]$bmi),2))

# Plot
fit <- survfit(mod_death,newdata=weights)

CairoPDF(file='km_death_less14_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# MALE
weights2 <- data.frame(vo2_graphical = levels(death_set$vo2_graphical),
                       age = rep(mean(death_set[sex==0]$age),2),
                       sex = rep(0,2),
                       bmi = rep(mean(death_set[sex==0]$bmi),2))

# Plot
fit <- survfit(mod_death,newdata=weights2)

CairoPDF(file='km_death_less14_adjusted_m.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.25),xlim=c(0,10),
     col=c('#ef3b2c','#02818a'),
     lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=seq(0,10,2),cex.axis=1.5,labels=as.character(seq(0,10,2)))
axis(2,at=seq(0,0.25,0.05),las=2,cex.axis=1.5,labels=c("0 %","5 %","10 %","15 %","20 %","25 %"))

mtext("Cumulative risk (%)",side=2,line=6,at=0.125,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.25*1.05,legend=c("Below 14 mL/kg/min","At or above 14 mL/kg/min"),
       col=c('#ef3b2c','#02818a'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

