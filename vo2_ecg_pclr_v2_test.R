# Script to assess value of ECG with VO2

# Depends
library(data.table)
library(stringr)
library(rms)
library(ggplot2)
library(glmnet)
library(epiR)
library(pROC)
library(nricens)
source('vo2_functions.R')

######## STEP 1: Pre-processing
# Load ECG / PCLR data
cpet_pclr <- fread(file='01-03-2022_c3po_pclr_infer_cpet.csv')

# Friendly names
## CPET ECG
setnames(cpet_pclr,
         c('VO2 (ml/kg/min)','VO2 Max (l/min)','VO2 % Predicted in Decimal Form (Formula)','Sex    1=male 2=Female',
           'Wt            lbs','Ht                 in','BMI (Formula)','Test         1=Bike 2=Treadmill 3=rower','Age'),
         c('vo2','vo2_max','vo2_pct','sex','wt','ht','bmi','test_type','age'))
cpet_pclr[,vo2_pct := vo2_pct * 100]

## PCLR
names(cpet_pclr)[str_detect(names(cpet_pclr),'^\\d+')] <- paste0('pclr_',names(cpet_pclr)[str_detect(names(cpet_pclr),'^\\d+')])

# Filter out missingness
## Figure out what is missing
na_sum <- function(x){sum(is.na(x))}
missing_ct <- apply(cpet_pclr,FUN=na_sum,MARGIN=2)
names(cpet_pclr)[which(missing_ct>0)] # ht, wt, bmi, vo2, PR interval, paxis
## Unify heights and weights
cpet_pclr[,":="(ht = ifelse(!is.na(ht),ht,heightin),
                wt = ifelse(!is.na(wt),wt,weightlbs))]
## Recover heights and weight from BMI when missing
cpet_pclr[,":="(ht = ifelse(!is.na(ht),ht,sqrt((weightlbs/2.2)/bmi)*39.3701),
                wt = ifelse(!is.na(wt),wt,(((heightin*0.0254)**2)*bmi)*2.2))]
## Use height and weight to calculate BMI
cpet_pclr[,":="(bmi = ifelse(!is.na(bmi),bmi,(wt/2.2)/(ht*0.0254)**2))]
## Remove missing
### Ht/wt/bmi
cpet_pclr <- cpet_pclr[!is.na(ht) & !is.na(wt) & !is.na(bmi)] # N = 3
### ECG params 
cpet_pclr <- cpet_pclr[c(!is.na(ventricularrate_md) & !is.na(qrsduration_md) # N = 11
                         & !is.na(printerval_md) & !is.na(qtinterval_md))]
### Target
cpet_pclr <- cpet_pclr[c(!is.na(vo2))] # N = 12
## Clean up variables
cpet_pclr <- cpet_pclr[,!c('error','patientage','CPP CPET Subject #','V1','s3_path')]
## Time to ECG
cpet_pclr[,":="(ecg_date = as.Date(substr(datetime,1,10),format='%Y-%m-%d'),
                cpet_date = as.Date(Date,format='%m/%d/%y'))]
cpet_pclr[,':='(cpet_to_ecg = (as.numeric(ecg_date) - as.numeric(cpet_date)))]
## Create dummy variables for test type
cpet_pclr[,':='(tt_bike = ifelse(test_type==1,1,0),
                tt_tread = ifelse(test_type==2,1,0),
                tt_row = ifelse(test_type==3,1,0))]
## Recode sex
cpet_pclr[,sex := ifelse(sex==1,0,1)]

# Isolate original holdout set
holdout <- fread(file='holdout.csv')
not_holdout <- cpet_pclr[!(MRN %in% holdout$MRN)]
holdout <- cpet_pclr[(MRN %in% holdout$MRN)]

# Split into train/test
set.seed(1)
not_holdout[,sample_var := runif(nrow(not_holdout))]
train <- not_holdout[sample_var <= 0.75]
test <- not_holdout[sample_var > 0.75]

# Save out not holdout
#write.csv(not_holdout,file='/data/arrhythmia/skhurshid/vo2/not_holdout_021822.csv',row.names=F)

# Generate X matrices
# Define variable space
outcomes_list <- c('vo2','vo2_max','vo2_pct')
basic_covars <- c('age','sex','bmi','tt_bike','tt_tread','tt_row')
ecg_basic <- c('age','sex','bmi',
               'tt_bike','tt_tread','tt_row',
               'ventricularrate_md','qrsduration_md','printerval_md','qtinterval_md')
pclr_only <- names(cpet_pclr)[str_detect(names(cpet_pclr),'c3po_pclr')]  
pclr_basic <- c(basic_covars,pclr_only)
all <- c(ecg_basic,pclr_only)

# Add additional comparators
## Unitize for OG equations
holdout[,':='(ht_cm = ht*2.54, wt_kg = wt/2.2)]
## Calculate predicted weight_cm for Wasserman equations
holdout[,predicted_wt_kg := ifelse(sex==0,ht_cm*0.79 - 60.7,ht_cm*0.65-42.8)]
## Calculate cycle factor for Wasserman equations
holdout[,cycle_factor := ifelse(sex==0,50.72-0.372*age,22.78-0.17*age)]
## Calculate Wasserman age
holdout[,wasserman_age := ifelse(age < 30,30,age)]
## Wasserman           
holdout[,wasserman := ifelse(sex==0,
                             ifelse(wt_kg >= predicted_wt_kg,
                                    ((0.0337*ht_cm-0.000165*wasserman_age*ht_cm-1.963+0.006*(wt_kg-predicted_wt_kg))/wt_kg)*1000,
                                    ((0.0337*ht_cm-0.000165*wasserman_age*ht_cm-1.963+0.014*(wt_kg-predicted_wt_kg))/wt_kg)*1000),
                             ((0.001*ht_cm*(14.783-0.11*wasserman_age)+0.006*(wt_kg-predicted_wt_kg))/wt_kg)*1000)]
holdout[,wasserman := ifelse(tt_bike==1,wasserman,wasserman*1.11)]

## Jones           
holdout[,jones := ifelse(tt_bike==1,
                         ifelse(sex==0,(-3.76+0.034*ht_cm-0.028*age+0.022*wt_kg)*1000/wt_kg,
                                (-2.26+0.025*ht_cm-0.018*age+0.010*wt_kg)*1000/wt_kg),
                         ifelse(sex==0,(4.2-0.032*age)*1000/wt_kg,(2.6-0.014*age)*1000/wt_kg))]
## Silva           
holdout[,silva := ifelse(tt_bike==1,
                         ifelse(sex==0,45.2-0.35*age-10.9*1-0.15*(wt_kg/2.2)+0.68*(ht_cm*0.393701)-0.46*2,
                                45.2-0.35*age-10.9*2-0.15*(wt_kg/2.2)+0.68*(ht_cm*0.393701)-0.46*2),
                         ifelse(sex==0,45.2-0.35*age-10.9*1-0.15*(wt_kg/2.2)+0.68*(ht_cm*0.393701)-0.46*1,
                                45.2-0.35*age-10.9*2-0.15*(wt_kg/2.2)+0.68*(ht_cm*0.393701)-0.46*1))]

# Standardize continuous variables for shrinkage
## Don't scale the all zero columns (generates NAs)
pclr_zero <- cpet_pclr[,lapply(.SD,mean),.SDcols=pclr_only]
pclr_zero <- names(pclr_zero)[pclr_zero == 0]

## Calculate and save scaling factors from training set
## relevant columns
scaled <- names(cpet_pclr)[c((names(cpet_pclr) %in% all[!(all %in% c('sex','tt_bike','tt_tread','tt_row',
                                                                     'wasserman','jones','silva'))])
                             & (!names(cpet_pclr) %in% pclr_zero))]
scale_factors <- data.table(var = scaled, 
                            means = lapply(train[,.SD,.SDcols=scaled],mean), 
                            vars = lapply(train[,.SD,.SDcols=scaled],sd))

## Now scale training set
for (j in names(cpet_pclr)[c((names(cpet_pclr) %in% all[!(all %in% c('sex','tt_bike','tt_tread','tt_row',
                                                                     'wasserman','jones','silva'))])
                             & (!names(cpet_pclr) %in% pclr_zero))]){
  set(train,j=j,value=scale(train[[j]]))
}

## Scale holdout and test sets
for (col in scaled){
  test[,col] <- ((test[,..col] - scale_factors[var==col]$means)/(scale_factors[var==col]$vars))
  holdout[,col] <- ((holdout[,..col] - scale_factors[var==col]$means)/(scale_factors[var==col]$vars))
}

# Generate reduced X matrices
x_basic = data.matrix(train[,.SD,.SDcols=basic_covars])
x_ecg_basic = data.matrix(train[,.SD,.SDcols=ecg_basic])
x_pclr = data.matrix(train[,.SD,.SDcols=pclr_only])
x_pclr_basic = data.matrix(train[,.SD,.SDcols=pclr_basic])

# Generate reduced X matrices
x_basic_holdout = data.matrix(holdout[,.SD,.SDcols=basic_covars])
x_ecg_basic_holdout = data.matrix(holdout[,.SD,.SDcols=ecg_basic])
x_pclr_holdout = data.matrix(holdout[,.SD,.SDcols=pclr_only])
x_pclr_basic_holdout = data.matrix(holdout[,.SD,.SDcols=pclr_basic])

######## STEP 2: Fit models
## Model 1: Basic
lm_basic <- glmnet(x=x_basic,y=train$vo2,family='gaussian',alpha=0.5)
set.seed(1)
cv.lm_basic <- cv.glmnet(x=x_basic,y=train$vo2,family='gaussian',alpha=0.5,
                         lambda=c(1,0.8,0.6,0.4,0.2,0.1,0.01,0.001))
coef_basic <- data.frame(coefs=coef(lm_basic,s=cv.lm_basic$lambda.min)[,1])

## Model 2: Basic + ECG
lm_ecg_basic <- glmnet(x=x_ecg_basic,y=train$vo2,family='gaussian',alpha=0.5)
set.seed(1)
cv.lm_ecg_basic <- cv.glmnet(x=x_ecg_basic,y=train$vo2,family='gaussian',alpha=0.5,
                             lambda=c(1,0.8,0.6,0.4,0.2,0.1,0.01,0.001))
coef_ecg_basic <- data.frame(coefs=coef(lm_ecg_basic,s=cv.lm_ecg_basic$lambda.min)[,1])

## Model 3: PCLR only
lm_pclr <- glmnet(x=x_pclr,y=train$vo2,family='gaussian',alpha=0)
set.seed(1)
cv.lm_pclr <- cv.glmnet(x=x_pclr,y=train$vo2,family='gaussian',alpha=0)
coef_pclr <- data.frame(coefs=coef(lm_pclr,s=cv.lm_pclr$lambda.min)[,1])

## Model 4: Basic + PCLR
lm_pclr_basic <- glmnet(x=x_pclr_basic,y=train$vo2,family='gaussian',alpha=1)
set.seed(1)
cv.lm_pclr_basic <- cv.glmnet(x=x_pclr_basic,y=train$vo2,family='gaussian',alpha=1)
coef_pclr_basic <- data.frame(coefs=coef(lm_pclr_basic,s=cv.lm_pclr_basic$lambda.min)[,1])

## Obtain fitted values in holdout set
holdout[,basic_fitted := predict(lm_basic,cv.lm_basic$lambda.min,newx=x_basic_holdout,type='link')]
holdout[,ecg_basic_fitted := predict(lm_ecg_basic,cv.lm_ecg_basic$lambda.min,newx=x_ecg_basic_holdout,type='link')]
holdout[,pclr_fitted := predict(lm_pclr,cv.lm_pclr$lambda.min,newx=x_pclr_holdout,type='link')]
holdout[,pclr_basic_fitted := predict(lm_pclr_basic,cv.lm_pclr_basic$lambda.min,newx=x_pclr_basic_holdout,type='link')]

# Save out holdout
#write.csv(holdout,file='/data/arrhythmia/skhurshid/vo2/vo2_holdout.csv',row.names=F)

vo2_mgh <- compare_linear_vo2(data=holdout,truth='vo2',
                              inference=c('basic_fitted','ecg_basic_fitted',
                                          'pclr_fitted','pclr_basic_fitted'),
                              out_path='/plots/',plot=TRUE,write_out=TRUE,
                              indicator='vo2_holdout')

pclr_v_ecg <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + pclr_basic_fitted,
                    data=holdout)

# Save out models
# save(lm_pclr_basic,file='/data/arrhythmia/skhurshid/vo2/lm_pclr_basic.RData')
# save(cv.lm_pclr_basic,file='/data/arrhythmia/skhurshid/vo2/cv.lm_pclr_basic.RData')

# MAE bootstrap
boot_mae <- ci_mae(y='vo2',x1='basic_fitted',
                   data=holdout,runs=1000)
raw_error <- mean(with(holdout,abs(vo2-basic_fitted)))
ci <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='pclr_basic_fitted',
                   data=holdout,runs=1000)
raw_error <- mean(with(holdout,abs(vo2-pclr_basic_fitted)))
ci <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='ecg_basic_fitted',
                   data=holdout,runs=1000)
raw_error <- mean(with(holdout,abs(vo2-ecg_basic_fitted)))
ci <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

# MAE bootstrap
boot_mae <- compare_mae(y='vo2',x1='ecg_basic_fitted',x2='pclr_basic_fitted',
                        data=holdout,runs=1000)
raw_error <- with(holdout,mean(abs(vo2-ecg_basic_fitted)-mean(abs(vo2-pclr_basic_fitted))))
z <- raw_error/sd(boot_mae[,3])
p <- 2*(1-pnorm(abs(z)))

# Analyze vs. old school models
vo2_os <- compare_linear_vo2(data=holdout,truth='vo2',
                             inference=c('jones','wasserman',
                                         'silva','pclr_basic_fitted'),
                             out_path='/plots/',plot=TRUE,write_out=TRUE,
                             indicator='vo2_holdout')

pclr_v_jones <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + jones, data=holdout)
pclr_v_wasserman <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + wasserman, data=holdout)
pclr_v_silva <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + silva, data=holdout)

boot_mae <- ci_mae(y='vo2',x1='jones',
                   data=holdout,runs=1000)
raw_error <- mean(with(holdout,abs(vo2-jones)))
ci_jones <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='wasserman',
                   data=holdout,runs=1000)
raw_error <- mean(with(holdout,abs(vo2-wasserman)))
ci_wasserman <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='silva',
                   data=holdout,runs=1000)
raw_error <- mean(with(holdout,abs(vo2-silva)))
ci_silva <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

## Comparisons
boot_mae <- compare_mae(y='vo2',x1='jones',x2='pclr_basic_fitted',
                        data=cpet_bwh,runs=1000)
raw_error <- with(cpet_bwh,mean(abs(vo2-jones)-mean(abs(vo2-pclr_basic_fitted))))
z <- raw_error/sd(boot_mae[,3])
jones_p <- 2*(1-pnorm(abs(z)))

boot_mae <- compare_mae(y='vo2',x1='wasserman',x2='pclr_basic_fitted',
                        data=cpet_bwh,runs=1000)
raw_error <- with(cpet_bwh,mean(abs(vo2-wasserman)-mean(abs(vo2-pclr_basic_fitted))))
z <- raw_error/sd(boot_mae[,3])
wasserman_p <- 2*(1-pnorm(abs(z)))

boot_mae <- compare_mae(y='vo2',x1='silva',x2='pclr_basic_fitted',
                        data=cpet_bwh,runs=1000)
raw_error <- with(cpet_bwh,mean(abs(vo2-silva)-mean(abs(vo2-pclr_basic_fitted))))
z <- raw_error/sd(boot_mae[,3])
silva_p <- 2*(1-pnorm(abs(z)))

# Targeted Plots
### BASIC
pdf(file='basic_vo2_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$basic_fitted,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,80),ylim=c(0,80),
     pch=19,col='#0000008C')
axis(1,at=seq(0,80,10),cex.axis=1.6)
axis(2,at=seq(0,80,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,81,81,lty=5)
dev.off()

### ECG BASIC
pdf(file='ecg_basic_vo2_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$ecg_basic_fitted,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,80),ylim=c(0,80),
     pch=19,col='#d7301f8C')
axis(1,at=seq(0,80,10),cex.axis=1.6)
axis(2,at=seq(0,80,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,81,81,lty=5)
dev.off()

### PCLR
pdf(file='pclr_basic_vo2_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$pclr_basic_fitted,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,80),ylim=c(0,80),
     pch=19,col='#2c7fb88C')
axis(1,at=seq(0,80,10),cex.axis=1.6)
axis(2,at=seq(0,80,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,81,81,lty=5)
dev.off()

# Targeted Bland-Altman Plots
### BASIC
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','basic_fitted')]
holdout[,diff_pair := vo2 - basic_fitted]
setkey(holdout,mean_pair)

pdf('ba_vo2_basic.pdf',height=4,width=4,pointsize=5)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#0000008C"

# Plot
plot(x=holdout$mean_pair,y=holdout$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,70),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,71,upper,lty=5,lwd=1,col='black')
segments(0,lower,71,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,71,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-32,at=seq(0,70,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-30,30,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=55,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=55,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### ECG BASIC
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','ecg_basic_fitted')]
holdout[,diff_pair := vo2 - ecg_basic_fitted]
setkey(holdout,mean_pair)

pdf('ba_vo2_ecg_basic.pdf',height=4,width=4,pointsize=5)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#d7301f8C"

# Plot
plot(x=holdout$mean_pair,y=holdout$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,70),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,71,upper,lty=5,lwd=1,col='black')
segments(0,lower,71,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,71,mean_diff,lty=1,lwd=1,col='#2c7fb8')

# Axes
axis(1,cex.axis=1.6,pos=-32,at=seq(0,70,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-30,30,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=55,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=55,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### PCLR BASIC
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2_max','pclr_basic_fitted')]
holdout[,diff_pair := vo2 - pclr_basic_fitted]
setkey(holdout,mean_pair)

pdf('ba_vo2_pclr_basic.pdf',height=4,width=4,pointsize=5)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#2c7fb88C"

# Plot
plot(x=holdout$mean_pair,y=holdout$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,40),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,41,upper,lty=5,lwd=1,col='black')
segments(0,lower,41,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,71,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-32,at=seq(0,40,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-30,30,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=31,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=31,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

# Targeted Density Plots
### ECG BASIC
x <- list(v1=holdout$ecg_basic_fitted,v2=holdout$vo2)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,7,1),expand=c(0,0),limits=c(0,7)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1),expand=c(0,0),limits=c(0,0.5)) +
  scale_fill_manual(values=c("#2c7fb8","#31a354"),name='',labels=c('ECG Basic VO2','True VO2')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Peak VO2 (mL/kg/min)',y='Density') 
ggsave('ecg_basic_density_vo2_max_compare.pdf',
       height=2,width=2.5,units='in',scale=4)

### PCLR BASIC
x <- list(v1=holdout$pclr_basic_fitted,v2=holdout$vo2)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,7,1),expand=c(0,0),limits=c(0,7)) +
  scale_y_continuous(breaks=seq(0,0.5,0.1),expand=c(0,0),limits=c(0,0.5)) +
  scale_fill_manual(values=c("#cb181d","#31a354"),name='',labels=c('PCLR Basic VO2','True VO2')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Peak VO2 (mL/kg/min)',y='Density') 
ggsave('pclr_basic_density_vo2_max_compare.pdf',
       height=2,width=2.5,units='in',scale=4)

# Sens/spec cutoff
## Create cats
holdout[, ':=' (pred_ecg_less14 = ifelse(ecg_basic_fitted < 14,1,0),
                true_less14 = ifelse(vo2 < 14,1,0),
                pred_pclr_less14 = ifelse(pclr_basic_fitted < 14,1,0),
                pred_basic_less14 = ifelse(basic_fitted < 14,1,0),
                true_less25pct = ifelse(vo2 < quantile(holdout$vo2,0.25),1,0),
                pred_pclr_less25pct = ifelse(pclr_basic_fitted < quantile(holdout$vo2,0.25),1,0),
                pred_basic_less25pct = ifelse(basic_fitted < quantile(holdout$vo2,0.25),1,0),
                pred_ecg_less25pct = ifelse(ecg_basic_fitted < quantile(holdout$vo2,0.25),1,0))]

## Create quartiles
holdout[, ':=' (true_quartiles = quantilize(vo2,4),
                pred_quartile = quantilize(pclr_basic_fitted,4))]

## Run table
## 14
basic_table <- table2x2(data=holdout,disease='true_less14',
                        test='pred_basic_less14',key='MRN')
basic <- epi.tests(basic_table)

ecg_table <- table2x2(data=holdout,disease='true_less14',
                      test='pred_ecg_less14',key='MRN')
ecg <- epi.tests(ecg_table)

pclr_table <- table2x2(data=holdout,disease='true_less14',
                       test='pred_pclr_less14',key='MRN')
pclr <- epi.tests(pclr_table)

## 25th percentile
basic_table_25 <- table2x2(data=holdout,disease='true_less25pct',
                           test='pred_basic_less25pct',key='MRN')
basic_25 <- epi.tests(basic_table_25)

ecg_table_25 <- table2x2(data=holdout,disease='true_less25pct',
                         test='pred_ecg_less25pct',key='MRN')
ecg_25 <- epi.tests(ecg_table_25)

pclr_table_25 <- table2x2(data=holdout,disease='true_less25pct',
                          test='pred_pclr_less25pct',key='MRN')
pclr_25 <- epi.tests(pclr_table_25)

# Make ROCs
basic_roc <- roc(response=holdout$true_less14,predictor=holdout$pred_basic_less14,ci=TRUE)
ecg_roc <- roc(response=holdout$true_less14,predictor=holdout$pred_ecg_less14,ci=TRUE)
pclr_roc <- roc(response=holdout$true_less14,predictor=holdout$pred_pclr_less14,ci=TRUE)

# Calculate SEs
pclr_basic_se <- boot_ci(status='true_less14',response1='pred_pclr_less14',response2='pred_basic_less14',data=holdout,runs=1000)
pclr_ecg_se <- boot_ci(status='true_less14',response1='pred_pclr_less14',response2='pred_ecg_less14',data=holdout,runs=1000)

## Perform z-tests
z_pclr_basic <- ((pclr_roc$auc - basic_roc$auc) / sd(pclr_basic_se))
p_pclr_basic <- 2*(1-pnorm(abs(z_pclr_basic)))

z_pclr_ecg <- ((pclr_roc$auc - ecg_roc$auc) / sd(pclr_ecg_se))
p_pclr_ecg <- 2*(1-pnorm(abs(z_pclr_ecg)))

## NRI
pclr_basic_reclass <- nribin(event=holdout$true_less14,
                             p.std=holdout$pred_basic_less14,p.new=holdout$pred_pclr_less14,
                             cut=0.5,niter=100)

pclr_ecg_reclass <- nribin(event=holdout$true_less14,
                           p.std=holdout$pred_ecg_less14,p.new=holdout$pred_pclr_less14,
                           cut=0.5,niter=100)

pclr_ecg_reclass <- nribin(event=holdout$true_less14,
                           p.std=holdout$pred_ecg_less14,p.new=holdout$pred_pclr_less14,
                           cut=0.5,niter=100)

pclr_basic_reclass25 <- nribin(event=holdout$true_less25pct,
                               p.std=holdout$pred_basic_less25pct,p.new=holdout$pred_pclr_less25pct,
                               cut=0.5,niter=100)

pclr_ecg_reclass25 <- nribin(event=holdout$true_less25pct,
                             p.std=holdout$pred_ecg_less25pct,p.new=holdout$pred_pclr_less25pct,
                             cut=0.5,niter=100)

## Error bars
# Define vars
holdout[,':='(vo2_decile = quantilize(vo2,10),
              abs_error = pclr_basic_fitted - vo2,
              rel_error = ((pclr_basic_fitted - vo2)/vo2)*100)]

# Plot absolute
y_vals <- holdout[,mean(abs_error),by='vo2_decile']
setkey(y_vals,vo2_decile)

pdf(file='abs_error_holdout.pdf',height=5,width=8,pointsize=5)
par(mar=c(10,7,1,2),oma=c(1,1,1,1))

coords <- barplot(height=y_vals$V1,col='#2c7fb88C',border=NA,
                  bty='n',xaxt='n',yaxt='n',xlim=c(0,12),
                  ylim=c(-10,10),xlab='',ylab='')

axis(1,at=coords,cex.axis=2.4,labels=rep('',length(coords)))
axis(2,cex.axis=2.4,at=seq(-10,10,5),las=2)

mtext("Error (mL/kg/min)",2,line=4.5,cex=2.2)
mtext("Decile of Peak VO2",1,line=8.5,cex=2.2)

labels <- c('<13.5','13.5-18.9','19.0-22.6','22.7-28.0','28.1-33.4',
            '33.5-38.5','38.6-41.7','41.8-46.0','46.1-50.6','>50.6')

text(x = coords,
     y = par("usr")[3] - 0.6,
     labels = labels,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.8)

dev.off()

# Plot relative
y_vals <- holdout[,mean(rel_error),by='vo2_decile']
setkey(y_vals,vo2_decile)

pdf(file='rel_error_holdout.pdf',height=5,width=8,pointsize=5)
par(mar=c(10,7,1,2),oma=c(1,1,1,1))

coords <- barplot(height=y_vals$V1,col='#d7301f8C',border=NA,
                  bty='n',xaxt='n',yaxt='n',xlim=c(0,12),
                  ylim=c(-50,50),xlab='',ylab='')

axis(1,at=coords,cex.axis=2.4,labels=rep('',length(coords)))
axis(2,cex.axis=2.4,at=seq(-50,50,25),las=2)

mtext("Error (%)",2,line=4.5,cex=2.2)
mtext("Decile of Peak VO2",1,line=8.5,cex=2.2)

labels <- c('<13.5','13.5-18.9','19.0-22.6','22.7-28.0','28.1-33.4',
            '33.5-38.5','38.6-41.7','41.8-46.0','46.1-50.6','>50.6')

text(x = coords,
     y = par("usr")[3] - 2.4,
     labels = labels,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.8)

dev.off()

# Bradley-Blackwood
## BASIC
holdout[,diff_pair := vo2 - basic_fitted]
holdout[,avg_pair := (vo2+basic_fitted)/2]

mod <- lm(holdout$diff_pair ~ holdout$avg_pair)
sy2 <- sum((holdout$diff_pair)**2)
ssr <- sum(mod$residuals**2)
df <- nrow(holdout)-2
msr <- ssr/df
bb <- (sy2-ssr)/(2*msr)
p_basic <- 1-pf(bb,df1=2,df2=df)

## ECG
holdout[,diff_pair := vo2 - ecg_basic_fitted]
holdout[,avg_pair := (vo2+ecg_basic_fitted)/2]

mod <- lm(holdout$diff_pair ~ holdout$avg_pair)
sy2 <- sum((holdout$diff_pair)**2)
ssr <- sum(mod$residuals**2)
df <- nrow(holdout)-2
msr <- ssr/df
bb <- (sy2-ssr)/(2*msr)
p_ecg <- 1-pf(bb,df1=2,df2=df)

## PCLR
holdout[,diff_pair := vo2 - pclr_basic_fitted]
holdout[,avg_pair := (vo2+pclr_basic_fitted)/2]

mod <- lm(holdout$diff_pair ~ holdout$avg_pair)
sy2 <- sum((holdout$diff_pair)**2)
ssr <- sum(mod$residuals**2)
df <- nrow(holdout)-2
msr <- ssr/df
bb <- (sy2-ssr)/(2*msr)
p_pclr <- 1-pf(bb,df1=2,df2=df)

# Targeted Plots for old school comparisons
## Original
### JONES
pdf(file='jones_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$jones,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,80),ylim=c(0,80),
     pch=19,col='#0000008C')
axis(1,at=seq(0,80,10),cex.axis=1.6)
axis(2,at=seq(0,80,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,81,81,lty=5)
dev.off()

### JONES
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','jones')]
holdout[,diff_pair := vo2 - jones]
setkey(holdout,mean_pair)
nrow(holdout[mean_pair <= 0])

png('ba_vo2_jones_holdout.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#0000008C"

# Plot
plot(x=holdout[mean_pair > 0]$mean_pair,y=holdout[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,60),ylim=c(-40,40),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,61,upper,lty=5,lwd=1,col='black')
segments(0,lower,61,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,61,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-42,at=seq(0,60,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-40,40,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=45,y=-36,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=45,y=-38,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### WASSERMAN
pdf(file='wasserman_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$wasserman,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,80),ylim=c(0,80),
     pch=19,col='#d7301f8C')
axis(1,at=seq(0,80,10),cex.axis=1.6)
axis(2,at=seq(0,80,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,81,81,lty=5)
dev.off()

### WASSERMAN
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','wasserman')]
holdout[,diff_pair := vo2 - wasserman]
setkey(holdout,mean_pair)
nrow(holdout[mean_pair <= 0])

png('ba_vo2_wasserman_holdout.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#d7301f8C"

# Plot
plot(x=holdout[mean_pair > 0]$mean_pair,y=holdout[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,60),ylim=c(-40,40),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,61,upper,lty=5,lwd=1,col='black')
segments(0,lower,61,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,61,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-42,at=seq(0,60,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-40,40,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=45,y=-36,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=45,y=-38,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### FRIEND
pdf(file='friend_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$silva,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,80),ylim=c(0,80),
     pch=19,col='#8dd3c78C')
axis(1,at=seq(0,80,10),cex.axis=1.6)
axis(2,at=seq(0,80,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,81,81,lty=5)
dev.off()

### FRIEND
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','silva')]
holdout[,diff_pair := vo2 - silva]
setkey(holdout,mean_pair)

png('ba_vo2_friend_holdout.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#8dd3c78C"

# Plot
plot(x=holdout[mean_pair > 0]$mean_pair,y=holdout[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,80),ylim=c(-60,20),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,81,upper,lty=5,lwd=1,col='black')
segments(0,lower,81,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,81,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-62,at=seq(0,80,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-60,20,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=60,y=-56,labels=paste0("Limits of Agreement: ",round(lower,2),' to ',round(upper,2)),cex=1.4)
text(x=60,y=-58,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

## RECALIBRATED
### JONES
pdf(file='jones_recal_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$jones_recal,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,80),
     pch=19,col='#0000008C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### JONES
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','jones_recal')]
holdout[,diff_pair := vo2 - jones_recal]
setkey(holdout,mean_pair)
nrow(holdout[mean_pair <= 0])

png('ba_vo2_jones_recal_holdout.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#0000008C"

# Plot
plot(x=holdout[mean_pair > 0]$mean_pair,y=holdout[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,51,upper,lty=5,lwd=1,col='black')
segments(0,lower,51,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,51,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-32,at=seq(0,50,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-30,30,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=35,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=35,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### WASSERMAN
pdf(file='wasserman_recal_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$wasserman_recal,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,80),
     pch=19,col='#d7301f8C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### WASSERMAN
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','wasserman_recal')]
holdout[,diff_pair := vo2 - wasserman_recal]
setkey(holdout,mean_pair)
nrow(holdout[mean_pair <= 0])

png('ba_vo2_wasserman_recal_holdout.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#d7301f8C"

# Plot
plot(x=holdout[mean_pair > 0]$mean_pair,y=holdout[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,51,upper,lty=5,lwd=1,col='black')
segments(0,lower,51,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,51,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-32,at=seq(0,50,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-30,30,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=15,y=26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=15,y=28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### FRIEND
pdf(file='friend_recal_holdout.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=holdout$silva_recal,y=holdout$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,80),
     pch=19,col='#8dd3c78C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,71,71,lty=5)
dev.off()

### FRIEND
holdout[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','silva_recal')]
holdout[,diff_pair := vo2 - silva_recal]
setkey(holdout,mean_pair)
nrow(holdout[mean_pair <= 0])

png('ba_vo2_friend_recal_holdout.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#8dd3c78C"

# Plot
plot(x=holdout[mean_pair > 0]$mean_pair,y=holdout[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(holdout$diff_pair)
upper <- mean(holdout$diff_pair)+1.96*sd(holdout$diff_pair)
lower <- mean(holdout$diff_pair)-1.96*sd(holdout$diff_pair)
segments(0,upper,51,upper,lty=5,lwd=1,col='black')
segments(0,lower,51,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,51,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-32,at=seq(0,50,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-30,30,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(holdout[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(holdout)
text(x=35,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=35,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

############# Secondary/subset analyses
### LBBB
no_lbbb <- holdout[lbbb==FALSE]
cor_no_lbbb <- cor.test(no_lbbb$pclr_basic_fitted,no_lbbb$vo2)
boot_mae <- ci_mae(y='vo2',x1='pclr_basic_fitted',
                   data=no_lbbb,runs=1000)
raw_error <- mean(with(no_lbbb,abs(vo2-pclr_basic_fitted)))
ci_no_lbbb <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

no_lbbb[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','pclr_basic_fitted')]
no_lbbb[,diff_pair := vo2 - pclr_basic_fitted]
upper <- mean(no_lbbb$diff_pair)+1.96*sd(no_lbbb$diff_pair)
lower <- mean(no_lbbb$diff_pair)-1.96*sd(no_lbbb$diff_pair)
agree_no_lbbb <- c(lower,upper)

### RBBB
no_rbbb <- holdout[rbbb==FALSE]
cor_no_rbbb <- cor.test(no_rbbb$pclr_basic_fitted,no_rbbb$vo2)
boot_mae <- ci_mae(y='vo2',x1='pclr_basic_fitted',
                   data=no_rbbb,runs=1000)
raw_error <- mean(with(no_rbbb,abs(vo2-pclr_basic_fitted)))
ci <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

no_rbbb[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','pclr_basic_fitted')]
no_rbbb[,diff_pair := vo2 - pclr_basic_fitted]
upper <- mean(no_rbbb$diff_pair)+1.96*sd(no_rbbb$diff_pair)
lower <- mean(no_rbbb$diff_pair)-1.96*sd(no_rbbb$diff_pair)
agree_no_rbbb <- c(lower,upper)

### LVH
no_lvh <- holdout[lvh==FALSE]
cor_no_lvh <- cor.test(no_lvh$pclr_basic_fitted,no_lvh$vo2)
boot_mae <- ci_mae(y='vo2',x1='pclr_basic_fitted',
                   data=no_lvh,runs=1000)
raw_error <- mean(with(no_lvh,abs(vo2-pclr_basic_fitted)))
ci_no_lvh <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

no_lvh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','pclr_basic_fitted')]
no_lvh[,diff_pair := vo2 - pclr_basic_fitted]
upper <- mean(no_lvh$diff_pair)+1.96*sd(no_lvh$diff_pair)
lower <- mean(no_lvh$diff_pair)-1.96*sd(no_lvh$diff_pair)
agree_no_lvh <- c(lower,upper)

### AF
no_af <- holdout[af==FALSE]
cor_no_af <- cor.test(no_af$pclr_basic_fitted,no_af$vo2)
boot_mae <- ci_mae(y='vo2',x1='pclr_basic_fitted',
                   data=no_af,runs=1000)
raw_error <- mean(with(no_af,abs(vo2-pclr_basic_fitted)))
ci_no_af <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

no_af[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','pclr_basic_fitted')]
no_af[,diff_pair := vo2 - pclr_basic_fitted]
upper <- mean(no_af$diff_pair)+1.96*sd(no_af$diff_pair)
lower <- mean(no_af$diff_pair)-1.96*sd(no_af$diff_pair)
agree_no_af <- c(lower,upper)

### Low EF
## Load low EFs
mgh_lvef <- fread(file='lvef_mgh_tc.csv')
no_lowef <- holdout[!(MRN %in% mgh_lvef[lvef<50]$MRN)]
cor_no_lowef <- cor.test(no_lowef$pclr_basic_fitted,no_lowef$vo2)
boot_mae <- ci_mae(y='vo2',x1='pclr_basic_fitted',
                   data=no_lowef,runs=1000)
raw_error <- mean(with(no_lowef,abs(vo2-pclr_basic_fitted)))
ci_no_lowef <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

no_lowef[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','pclr_basic_fitted')]
no_lowef[,diff_pair := vo2 - pclr_basic_fitted]
upper <- mean(no_lowef$diff_pair)+1.96*sd(no_lowef$diff_pair)
lower <- mean(no_lowef$diff_pair)-1.96*sd(no_lowef$diff_pair)
agree_no_lowef <- c(lower,upper)