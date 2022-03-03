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

####### PARAMS
set_alpha <- 1
set_alpha2 <- 0.5

######## STEP 1: Pre-processing
# Load ECG / PCLR data
cpet_pclr <- fread(file='/data/arrhythmia/skhurshid/vo2/01-03-2022_c3po_pclr_infer_cpet.csv')

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

# Plot density of time to ECG
ggplot() + geom_histogram(data=cpet_pclr,aes(x=cpet_pclr$cpet_to_ecg),
                          binwidth=(1 + 3.322*log(nrow(cpet_pclr))),fill="#f03b20",alpha=0.55) +
  scale_x_continuous(breaks=seq(-150,100,50),expand=c(0.01,0),limits=c(-150,100)) +
  scale_y_continuous(breaks=seq(0,1400,200),expand=c(0,0),limits=c(0,1400)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position='None',
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.6,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x=expression(paste("Days from CPET to ECG")),y='Count')
ggsave(filename='/data/arrhythmia/skhurshid/vo2/days_cpet_to_ecg.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# Plot density of true VO2
ggplot() + geom_density(data=cpet_pclr,aes(x=cpet_pclr$vo2),fill="#2c7fb8",alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,85,10),expand=c(0,0),limits=c(0,85)) +
  scale_y_continuous(breaks=seq(0,0.03,0.01),expand=c(0,0),limits=c(0,0.03)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Peak VO2 (mL/kg/min)',y='Density') 
ggsave('/data/arrhythmia/skhurshid/vo2/all_true_vo2_density.pdf',
       height=2,width=2.5,units='in',scale=4)

# Save out analysis set
#write.csv(cpet_pclr,file='/data/arrhythmia/skhurshid/vo2/cpet_pclr_021822.csv',row.names=F)

# Isolate original holdout set
holdout <- fread(file='/data/arrhythmia/skhurshid/vo2/holdout.csv')
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
  
# Standardize continuous variables for shrinkage
## Don't scale the all zero columns (generates NAs)
pclr_zero <- cpet_pclr[,lapply(.SD,mean),.SDcols=pclr_only]
pclr_zero <- names(pclr_zero)[pclr_zero == 0]

## Now scale
for (j in names(cpet_pclr)[c((names(cpet_pclr) %in% all[!(all %in% c('sex','tt_bike','tt_tread','tt_row'))])
                               & (!names(cpet_pclr) %in% pclr_zero))]){
  set(train,j=j,value=scale(train[[j]]))
  set(test,j=j,value=scale(test[[j]]))
  set(holdout,j=j,value=scale(holdout[[j]]))
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
lm_basic <- glmnet(x=x_basic,y=train$vo2,family='gaussian',alpha=set_alpha2)
set.seed(1)
cv.lm_basic <- cv.glmnet(x=x_basic,y=train$vo2,family='gaussian',alpha=set_alpha2,
                         lambda=c(1,0.8,0.6,0.4,0.2,0.1,0.01,0.001))
plot(cv.lm_basic)
coef(lm_basic,s=cv.lm_basic$lambda.min)

## Model 2: Basic + ECG
lm_ecg_basic <- glmnet(x=x_ecg_basic,y=train$vo2,family='gaussian',alpha=set_alpha2)
set.seed(1)
cv.lm_ecg_basic <- cv.glmnet(x=x_ecg_basic,y=train$vo2,family='gaussian',alpha=set_alpha2,
                             lambda=c(1,0.8,0.6,0.4,0.2,0.1,0.01,0.001))
plot(cv.lm_ecg_basic)
coef(lm_ecg_basic,s=cv.lm_ecg_basic$lambda.min)

## Model 3: PCLR only
lm_pclr <- glmnet(x=x_pclr,y=train$vo2,family='gaussian',alpha=set_alpha)
set.seed(1)
cv.lm_pclr <- cv.glmnet(x=x_pclr,y=train$vo2,family='gaussian',alpha=set_alpha)
plot(cv.lm_pclr)
coef(lm_pclr,s=cv.lm_pclr$lambda.min)

## Model 4: Basic + PCLR
lm_pclr_basic <- glmnet(x=x_pclr_basic,y=train$vo2,family='gaussian',alpha=set_alpha)
set.seed(1)
cv.lm_pclr_basic <- cv.glmnet(x=x_pclr_basic,y=train$vo2,family='gaussian',alpha=set_alpha)
plot(cv.lm_pclr_basic)
coef(lm_pclr_basic,s=cv.lm_pclr_basic$lambda.min)

## Obtain fitted values in holdout set
holdout[,basic_fitted := predict(lm_basic,cv.lm_basic$lambda.min,newx=x_basic_holdout,type='link')]
holdout[,ecg_basic_fitted := predict(lm_ecg_basic,cv.lm_ecg_basic$lambda.min,newx=x_ecg_basic_holdout,type='link')]
holdout[,pclr_fitted := predict(lm_pclr,cv.lm_pclr$lambda.min,newx=x_pclr_holdout,type='link')]
holdout[,pclr_basic_fitted := predict(lm_pclr_basic,cv.lm_pclr_basic$lambda.min,newx=x_pclr_basic_holdout,type='link')]

vo2 <- compare_linear_vo2(data=holdout,truth='vo2',
                      inference=c('basic_fitted','ecg_basic_fitted',
                                  'pclr_fitted','pclr_basic_fitted'),
                      out_path='/data/arrhythmia/skhurshid/vo2/',plot=TRUE,write_out=TRUE,
                      indicator='vo2_alpha1_holdout')

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
raw_error <- mean(with(holdout,abs(ecg_basic_fitted-pclr_basic_fitted)))
z <- raw_error/sd(boot_mae[,3])
p <- 2*(1-pnorm(abs(z)))

# Targeted Plots
### BASIC
pdf(file='/data/arrhythmia/skhurshid/vo2/basic_vo2_holdout.pdf',height=4,width=4,pointsize=5)
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
pdf(file='/data/arrhythmia/skhurshid/vo2/ecg_basic_vo2_holdout.pdf',height=4,width=4,pointsize=5)
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
pdf(file='/data/arrhythmia/skhurshid/vo2/pclr_basic_vo2_holdout.pdf',height=4,width=4,pointsize=5)
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

png('/data/arrhythmia/skhurshid/vo2/ba_vo2_basic.png',pointsize=6,res=300,
    height=1200,width=1200)
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
mtext(side=1,"Mean of peak VO2 (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in peak VO2 (mL/kg/min)",cex=1.8,line=2.5)

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

png('/data/arrhythmia/skhurshid/vo2/ba_vo2_ecg_basic.png',pointsize=6,res=300,
    height=1200,width=1200)
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
mtext(side=1,"Mean of peak VO2 (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in peak VO2 (mL/kg/min)",cex=1.8,line=2.5)

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

png('/data/arrhythmia/skhurshid/vo2/ba_vo2_pclr_basic.png',pointsize=6,res=300,
    height=1200,width=1200)
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
mtext(side=1,"Mean of peak VO2 (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in peak VO2 (mL/kg/min)",cex=1.8,line=2.5)

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
ggsave('/data/arrhythmia/skhurshid/vo2/ecg_basic_density_vo2_max_compare.pdf',
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
ggsave('/data/arrhythmia/skhurshid/vo2/pclr_basic_density_vo2_max_compare.pdf',
       height=2,width=2.5,units='in',scale=4)

# Sens/spec cutoff
## Create cats
holdout[, ':=' (pred_ecg_less14 = ifelse(ecg_basic_fitted < 14,1,0),
                true_less14 = ifelse(vo2 < 14,1,0),
                pred_pclr_less14 = ifelse(pclr_basic_fitted < 14,1,0),
                pred_basic_less14 = ifelse(basic_fitted < 14,1,0))]

## Run table
basic_table <- table2x2(data=holdout,disease='true_less14',
                     test='pred_basic_less14',key='MRN')
basic <- epi.tests(basic_table)

ecg_table <- table2x2(data=holdout,disease='true_less14',
                        test='pred_ecg_less14',key='MRN')
ecg <- epi.tests(ecg_table)

pclr_table <- table2x2(data=holdout,disease='true_less14',
                        test='pred_pclr_less14',key='MRN')
pclr <- epi.tests(pclr_table)

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