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

## relevant columns
scaled <- names(cpet_pclr)[c((names(cpet_pclr) %in% all[!(all %in% c('sex','tt_bike','tt_tread','tt_row'))])
                             & (!names(cpet_pclr) %in% pclr_zero))]
scale_factors <- data.table(var = scaled, 
                            means = lapply(train[,.SD,.SDcols=scaled],mean), 
                            vars = lapply(train[,.SD,.SDcols=scaled],sd))

# Write out scale factors
fwrite(scale_factors,file='scale_factors.csv',row.names=F)

## Now scale training set
for (j in names(cpet_pclr)[c((names(cpet_pclr) %in% all[!(all %in% c('sex','tt_bike','tt_tread','tt_row'))])
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
cv.lm_pclr_basic <- cv.glmnet(x=x_pclr_basic,y=train$vo2,family='gaussian',alpha=1)
coef_pclr_basic <- data.frame(coefs=coef(lm_pclr_basic,s=cv.lm_pclr_basic$lambda.min)[,1])

### Write out models
## Lambda values
lambdas <- data.table(model=c('basic','ecg_basic','pclr','pclr_basic'),
                      lambda=c(cv.lm_basic$lambda.min,cv.lm_ecg_basic$lambda.min,
                               cv.lm_pclr$lambda.min,cv.lm_pclr_basic$lambda.min))
write.csv(lambdas,file='/Volumes/arrhythmia/skhurshid/vo2/lambdas.csv',row.names=F)

## Models
save(lm_basic,file='basic.RData')
save(lm_ecg_basic,file='ecg_basic.RData')
save(lm_pclr,file='pclr.RData')
save(lm_pclr_basic,file='pclr_basic.RData')

######## STEP 2: BWH VALIDATION
# Load ECG / PCLR data
cpet_bwh <- fread(file='092322_bwh_ecg_inference.tsv')

# Remove one study with erroneous VO2 value
cpet_bwh <- cpet_bwh[peak_vo2 > 1]

# Pick the CPET that has the closest time to an ECG
cpet_bwh <- cpet_bwh[,.SD[which.min(abs_diff)],by='linker_id']

cpet_bwh[,c('Report_Date_Time_Age','Report_Date','datetime','abs_diff')]
# Friendly names
## CPET ECG
setnames(cpet_bwh,
         c('peak_vo2','Gender','Report_Date_Time_Age'),
         c('vo2','sex','age'))

# Load covariates
cpet_bwh_covar <- fread(file='cpet_bwh_covars.csv')

# Recodes and joins
setkey(cpet_bwh_covar,linker_id); setkey(cpet_bwh,linker_id)
cpet_bwh[cpet_bwh_covar,':='(bmi = i.bmi, ht=i.ht, wt=i.wt)]
cpet_bwh[,sex := ifelse(sex=="Male",0,1)]
cpet_bwh[,age := as.numeric(str_extract(age,'\\d+'))/365.25]

# Load modality/submax
modality <- fread(file='c3po_cpet_vo2_modality.csv')
modality[,age := as.numeric(str_extract(Report_Date_Time_Age,'\\d+'))/365.25]
cpet_bwh[,age2 := as.numeric(round(age,2))]
modality[,age2 := as.numeric(round(age,2))]
setkey(modality,linker_id,age2); setkey(cpet_bwh,linker_id,age2)
modality <- modality[cpet_bwh,nomatch=0]

# Add modality info
cpet_bwh[modality,':='(modality = i.Modality,
                       submax = i.submax)]

# Remove submax
cpet_bwh <- cpet_bwh[submax==0]

# Modality
cpet_bwh[,':='(tt_bike = ifelse(modality == 'Cycle',1,0),
               tt_row = 0,
               tt_tread = ifelse(modality == 'Treadmill',1,0))]

## Remove missing
### Ht/wt/bmi
cpet_bwh <- cpet_bwh[!is.na(bmi)]
### ECG params 
cpet_bwh <- cpet_bwh[c(!is.na(ventricularrate_md) & !is.na(qrsduration_md)
                       & !is.na(printerval_md) & !is.na(qtinterval_md))]

# Plot density of time to ECG
ggplot() + geom_histogram(data=cpet_bwh,aes(x=cpet_bwh$diff),
                          binwidth=(1 + 3.322*log(nrow(cpet_bwh))),fill="#2c7fb8",alpha=0.55) +
  scale_x_continuous(breaks=seq(-150,100,50),expand=c(0.01,0),limits=c(-150,100)) +
  scale_y_continuous(breaks=seq(0,600,100),expand=c(0,0),limits=c(0,600)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position='None',
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.6,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x=expression(paste("Days from CPET to ECG")),y='Count')
ggsave(filename='days_cpet_to_ecg_bwh.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# Write out
write.csv(cpet_bwh,file='cpet_bwh_final_062622.csv',row.names=F)

# Add additional comparators
## Calculate predicted weight for Wasserman equations
cpet_bwh[,predicted_wt := ifelse(sex==0,ht*0.79 - 60.7,ht*0.65-42.8)]
## Calculate cycle factor for Wasserman equations
cpet_bwh[,cycle_factor := ifelse(sex==0,50.72-0.372*age,22.78-0.17*age)]
## Calculate Wasserman age
cpet_bwh[,wasserman_age := ifelse(age < 30,30,age)]
## Wasserman           
cpet_bwh[,wasserman := ifelse(sex==0,
                              ifelse(wt >= predicted_wt,
                                     ((0.0337*ht-0.000165*wasserman_age*ht-1.963+0.006*(wt-predicted_wt))/wt)*1000,
                                     ((0.0337*ht-0.000165*wasserman_age*ht-1.963+0.014*(wt-predicted_wt))/wt)*1000),
                              ((0.001*ht*(14.783-0.11*wasserman_age)+0.006*(wt-predicted_wt))/wt)*1000)]
cpet_bwh[,wasserman := ifelse(modality=='Cycle',wasserman,wasserman*1.11)]

## Jones           
cpet_bwh[,jones := ifelse(modality=='Cycle',
                          ifelse(sex==0,(-3.76+0.034*ht-0.028*age+0.022*wt)*1000/wt,
                                 (-2.26+0.025*ht-0.018*age+0.010*wt)*1000/wt),
                          ifelse(sex==0,(4.2-0.032*age)*1000/wt,(2.6-0.014*age)*1000/wt))]
## Silva           
cpet_bwh[,silva := ifelse(modality=='Cycle',
                          ifelse(sex==0,45.2-0.35*age-10.9*1-0.15*(wt*2.2)+0.68*(ht*0.393701)-0.46*2,
                                 45.2-0.35*age-10.9*2-0.15*(wt*2.2)+0.68*(ht*0.393701)-0.46*2),
                          ifelse(sex==0,45.2-0.35*age-10.9*1-0.15*(wt/2.2)+0.68*(ht*0.393701)-0.46*1,
                                 45.2-0.35*age-10.9*2-0.15*(wt*2.2)+0.68*(ht*0.393701)-0.46*1))]

# Standardize continuous variables for shrinkage
## Don't scale the all zero columns (generates NAs)
pclr_zero <- cpet_bwh[,lapply(.SD,mean),.SDcols=pclr_only]
pclr_zero <- names(pclr_zero)[pclr_zero == 0]

## Scale holdout and test sets
for (col in scaled){
  cpet_bwh[,col] <- ((cpet_bwh[,..col] - scale_factors[var==col]$means)/(scale_factors[var==col]$vars))
}

# Generate reduced X matrices
x_basic_bwh = data.matrix(cpet_bwh[,.SD,.SDcols=basic_covars])
x_ecg_basic_bwh = data.matrix(cpet_bwh[,.SD,.SDcols=ecg_basic])
x_pclr_bwh = data.matrix(cpet_bwh[,.SD,.SDcols=pclr_only])
x_pclr_basic_bwh = data.matrix(cpet_bwh[,.SD,.SDcols=pclr_basic])

## Obtain fitted values in holdout set
cpet_bwh[,basic_fitted := predict(lm_basic,cv.lm_basic$lambda.min,newx=x_basic_bwh,type='link')]
cpet_bwh[,ecg_basic_fitted := predict(lm_ecg_basic,cv.lm_ecg_basic$lambda.min,newx=x_ecg_basic_bwh,type='link')]
cpet_bwh[,pclr_fitted := predict(lm_pclr,cv.lm_pclr$lambda.min,newx=x_pclr_bwh,type='link')]
cpet_bwh[,pclr_basic_fitted := predict(lm_pclr_basic,cv.lm_pclr_basic$lambda.min,newx=x_pclr_basic_bwh,type='link')]

# Mean shifted
cpet_bwh[,basic_fitted_recal := basic_fitted - (mean(basic_fitted) - mean(vo2))]
cpet_bwh[,ecg_basic_fitted_recal := ecg_basic_fitted - (mean(ecg_basic_fitted) - mean(vo2))]
cpet_bwh[,pclr_fitted_recal := pclr_fitted - (mean(pclr_fitted) - mean(vo2))]
cpet_bwh[,pclr_basic_fitted_recal := pclr_basic_fitted - (mean(pclr_basic_fitted) - mean(vo2))]

cpet_bwh[,jones_recal := jones - (mean(jones) - mean(vo2))]
cpet_bwh[,wasserman_recal := wasserman - (mean(wasserman) - mean(vo2))]
cpet_bwh[,silva_recal := silva - (mean(silva) - mean(vo2))]

# Save out dataset
write.csv(cpet_bwh,file='bwh_validation_062623.csv',row.names=F)

# Analyze
vo2 <- compare_linear_vo2_bwh(data=cpet_bwh,truth='vo2',
                              inference=c('basic_fitted','ecg_basic_fitted',
                                          'pclr_fitted','pclr_basic_fitted'),
                              out_path='/plots/',plot=TRUE,write_out=TRUE,
                              indicator='vo2_bwh')

vo2_recal <- compare_linear_vo2_bwh(data=cpet_bwh,truth='vo2',
                                    inference=c('basic_fitted_recal','ecg_basic_fitted_recal',
                                                'pclr_fitted_recal','pclr_basic_fitted_recal'),
                                    out_path='/plots/',plot=TRUE,write_out=TRUE,
                                    indicator='vo2_bwh_recal')

pclr_v_basic <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + basic_fitted,
                      data=cpet_bwh)

pclr_v_ecg <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + pclr_basic_fitted,
                    data=cpet_bwh)

pclr_v_ecg_recal <- cocor(formula = ~ vo2 + ecg_basic_fitted_recal | vo2 + pclr_basic_fitted_recal,
                          data=cpet_bwh)

# Analyze vs. old school models
vo2_os <- compare_linear_vo2_bwh(data=cpet_bwh,truth='vo2',
                                 inference=c('jones','wasserman',
                                             'silva','pclr_basic_fitted'),
                                 out_path='/plots/',plot=TRUE,write_out=TRUE,
                                 indicator='vo2_bwh')

pclr_v_jones <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + jones, data=cpet_bwh)
pclr_v_wasserman <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + wasserman, data=cpet_bwh)
pclr_v_silva <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + silva, data=cpet_bwh)

boot_mae <- ci_mae(y='vo2',x1='jones',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-jones)))
ci_jones <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='wasserman',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-wasserman)))
ci_wasserman <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='silva',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-silva)))
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

## RECAL
vo2_recal <- compare_linear_vo2_bwh(data=cpet_bwh,truth='vo2',
                                    inference=c('basic_fitted_recal','ecg_basic_fitted_recal',
                                                'pclr_fitted_recal','pclr_basic_fitted_recal'),
                                    out_path='/plots/',plot=TRUE,write_out=TRUE,
                                    indicator='vo2_bwh_recal')

pclr_v_basic <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + basic_fitted,
                      data=cpet_bwh)

pclr_v_ecg <- cocor(formula = ~ vo2 + ecg_basic_fitted | vo2 + pclr_basic_fitted,
                    data=cpet_bwh)

pclr_v_ecg_recal <- cocor(formula = ~ vo2 + ecg_basic_fitted_recal | vo2 + pclr_basic_fitted_recal,
                          data=cpet_bwh)

# Analyze vs. old school models
vo2_os_recal <- compare_linear_vo2_bwh(data=cpet_bwh,truth='vo2',
                                       inference=c('jones_recal','wasserman_recal',
                                                   'silva_recal','pclr_basic_fitted_recal'),
                                       out_path='/plots/',plot=TRUE,write_out=TRUE,
                                       indicator='vo2_bwh_recal')

boot_mae <- ci_mae(y='vo2',x1='jones_recal',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-jones_recal)))
ci_jones_recal <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='wasserman_recal',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-wasserman_recal)))
ci_wasserman_recal <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='silva_recal',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-silva_recal)))
ci_silva_recal <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

## Comparisons
boot_mae <- compare_mae(y='vo2',x1='jones_recal',x2='pclr_basic_fitted_recal',
                        data=cpet_bwh,runs=1000)
raw_error <- with(cpet_bwh,mean(abs(vo2-jones_recal)-mean(abs(vo2-pclr_basic_fitted_recal))))
z <- raw_error/sd(boot_mae[,3])
jones_p_recal <- 2*(1-pnorm(abs(z)))

boot_mae <- compare_mae(y='vo2',x1='wasserman_recal',x2='pclr_basic_fitted_recal',
                        data=cpet_bwh,runs=1000)
raw_error <- with(cpet_bwh,mean(abs(vo2-wasserman_recal)-mean(abs(vo2-pclr_basic_fitted_recal))))
z <- raw_error/sd(boot_mae[,3])
wasserman_p_recal <- 2*(1-pnorm(abs(z)))

boot_mae <- compare_mae(y='vo2',x1='silva_recal',x2='pclr_basic_fitted_recal',
                        data=cpet_bwh,runs=1000)
raw_error <- with(cpet_bwh,mean(abs(vo2-silva_recal)-mean(abs(vo2-pclr_basic_fitted_recal))))
z <- raw_error/sd(boot_mae[,3])
silva_p_recal <- 2*(1-pnorm(abs(z)))

# Save out models
save(lm_pclr_basic,file='lm_pclr_basic.RData')
save(cv.lm_pclr_basic,file='cv.lm_pclr_basic.RData')

# MAE bootstrap
boot_mae <- ci_mae(y='vo2',x1='basic_fitted',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-basic_fitted)))
ci_basic <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='ecg_basic_fitted',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-ecg_basic_fitted)))
ci_ecg <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='pclr_basic_fitted',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-pclr_basic_fitted)))
ci_pclr <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

# MAE bootstrap recal
boot_mae <- ci_mae(y='vo2',x1='basic_fitted_recal',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-basic_fitted_recal)))
ci_basic_recal <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='ecg_basic_fitted_recal',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-ecg_basic_fitted_recal)))
ci_ecg_recal <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

boot_mae <- ci_mae(y='vo2',x1='pclr_basic_fitted_recal',
                   data=cpet_bwh,runs=1000)
raw_error <- mean(with(cpet_bwh,abs(vo2-pclr_basic_fitted_recal)))
ci_pclr_recal <- c(raw_error,raw_error - 1.96*sd(boot_mae),raw_error + 1.96*sd(boot_mae))

# MAE bootstrap
boot_mae <- compare_mae(y='vo2',x1='ecg_basic_fitted',x2='pclr_basic_fitted',
                        data=cpet_bwh,runs=1000)
raw_error <- with(cpet_bwh,mean(abs(vo2-ecg_basic_fitted)-mean(abs(vo2-pclr_basic_fitted))))
z <- raw_error/sd(boot_mae[,3])
p <- 2*(1-pnorm(abs(z)))

# MAE bootstrap refitted
boot_mae <- compare_mae(y='vo2',x1='ecg_basic_fitted_recal',x2='pclr_basic_fitted_recal',
                        data=cpet_bwh,runs=1000)
raw_error <- with(cpet_bwh,mean(abs(vo2-ecg_basic_fitted_recal)-mean(abs(vo2-pclr_basic_fitted_recal))))
z <- raw_error/sd(boot_mae[,3])
p <- 2*(1-pnorm(abs(z)))

# Plot density of true VO2
ggplot() + geom_density(data=cpet_bwh,aes(x=cpet_bwh$vo2),fill="#bdbdbd") + 
  scale_x_continuous(breaks=seq(0,60,10),expand=c(0,0),limits=c(0,60)) +
  scale_y_continuous(breaks=seq(0,0.08,0.02),expand=c(0,0),limits=c(0,0.08)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Peak VO2 (mL/kg/min)',y='Density') 
ggsave('all_true_vo2_density_bwh.pdf',
       height=2,width=2.5,units='in',scale=4)

# Plot density of true VO2 in men
ggplot() + geom_density(data=cpet_bwh[sex==0],aes(x=cpet_bwh[sex==0]$vo2),fill="#2c7fb8",alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,60,10),expand=c(0,0),limits=c(0,60)) +
  scale_y_continuous(breaks=seq(0,0.08,0.02),expand=c(0,0),limits=c(0,0.08)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Peak VO2 (mL/kg/min)',y='Density') 
ggsave('m_true_vo2_density_bwh.pdf',
       height=2,width=2.5,units='in',scale=4)

# Plot density of true VO2 in women
ggplot() + geom_density(data=cpet_bwh[sex==1],aes(x=cpet_bwh[sex==1]$vo2),fill="#e31a1c",alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,60,10),expand=c(0,0),limits=c(0,60)) +
  scale_y_continuous(breaks=seq(0,0.08,0.02),expand=c(0,0),limits=c(0,0.08)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Peak VO2 (mL/kg/min)',y='Density') 
ggsave('f_true_vo2_density_bwh.pdf',
       height=2,width=2.5,units='in',scale=4)

# Targeted Plots
### BASIC
pdf(file='basic_vo2_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh[basic_fitted >= 0]$basic_fitted,y=cpet_bwh[basic_fitted >= 0]$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#0000008C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### BASIC REFITTED
pdf(file='basic_vo2_bwh_refitted.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh[basic_fitted_recal >= 0]$basic_fitted_recal,y=cpet_bwh[basic_fitted_recal >= 0]$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#0000008C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### ECG BASIC
pdf(file='ecg_basic_vo2_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh[ecg_basic_fitted >= 0]$ecg_basic_fitted,y=cpet_bwh[ecg_basic_fitted >= 0]$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#d7301f8C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### ECG BASIC REFITTED
pdf(file='ecg_basic_vo2_bwh_refitted.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh[ecg_basic_fitted_recal >= 0]$ecg_basic_fitted_recal,y=cpet_bwh[ecg_basic_fitted_recal >= 0]$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#d7301f8C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### PCLR
pdf(file='pclr_basic_vo2_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh[pclr_basic_fitted >= 0]$pclr_basic_fitted,y=cpet_bwh[pclr_basic_fitted >= 0]$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#2c7fb88C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### PCLR
pdf(file='pclr_basic_vo2_bwh_refitted.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh[pclr_basic_fitted_recal >= 0]$pclr_basic_fitted_recal,y=cpet_bwh[pclr_basic_fitted_recal >= 0]$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#2c7fb88C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### PCLR
pdf(file='pclr_basic_vo2_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh[pclr_basic_fitted >= 0]$pclr_basic_fitted,y=cpet_bwh[pclr_basic_fitted >= 0]$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#2c7fb88C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### PCLR
pdf(file='pclr_basic_vo2_bwh_refitted.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh[pclr_basic_fitted_recal >= 0]$pclr_basic_fitted_recal,y=cpet_bwh[pclr_basic_fitted_recal >= 0]$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#2c7fb88C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

######### Targeted Bland-Altman Plots
### BASIC
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','basic_fitted_recal')]
cpet_bwh[,diff_pair := vo2 - basic_fitted_recal]
setkey(cpet_bwh,mean_pair)

png('ba_vo2_basic_bwh_recal.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#0000008C"
nrow(cpet_bwh[mean_pair <= 0])

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=37,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=37,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### ECG BASIC
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','ecg_basic_fitted_recal')]
cpet_bwh[,diff_pair := vo2 - ecg_basic_fitted_recal]
setkey(cpet_bwh,mean_pair)

png('ba_vo2_ecg_basic_bwh_recal.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#d7301f8C"
nrow(cpet_bwh[mean_pair <= 0])

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
segments(0,upper,51,upper,lty=5,lwd=1,col='black')
segments(0,lower,51,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,51,mean_diff,lty=1,lwd=1,col='#2c7fb8')

# Axes
axis(1,cex.axis=1.6,pos=-32,at=seq(0,50,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-30,30,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=37,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=37,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### PCLR BASIC
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','pclr_basic_fitted_recal')]
cpet_bwh[,diff_pair := vo2 - pclr_basic_fitted_recal]
setkey(cpet_bwh,mean_pair)
nrow(cpet_bwh[mean_pair <= 0])

png('ba_vo2_pclr_basic_bwh_recal.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#2c7fb88C"

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=37,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=37,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

#####
### BASIC
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','basic_fitted')]
cpet_bwh[,diff_pair := vo2 - basic_fitted]
setkey(cpet_bwh,mean_pair)

pdf('ba_vo2_basic_bwh.pdf',height=4,width=4,pointsize=5)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#0000008C"
nrow(cpet_bwh[mean_pair <= 0])

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=37,y=26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=37,y=24,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### ECG BASIC
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','ecg_basic_fitted')]
cpet_bwh[,diff_pair := vo2 - ecg_basic_fitted]
setkey(cpet_bwh,mean_pair)

pdf('ba_vo2_ecg_basic_bwh.pdf',height=4,width=4,pointsize=5)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#d7301f8C"
nrow(cpet_bwh[mean_pair <= 0])

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
segments(0,upper,51,upper,lty=5,lwd=1,col='black')
segments(0,lower,51,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,51,mean_diff,lty=1,lwd=1,col='#2c7fb8')

# Axes
axis(1,cex.axis=1.6,pos=-32,at=seq(0,50,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-30,30,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=37,y=26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=37,y=24,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### PCLR BASIC
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','pclr_basic_fitted')]
cpet_bwh[,diff_pair := vo2 - pclr_basic_fitted]
setkey(cpet_bwh,mean_pair)
nrow(cpet_bwh[mean_pair <= 0])

pdf('ba_vo2_pclr_basic_bwh.pdf',height=4,width=4,pointsize=5)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#2c7fb88C"

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=37,y=26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=37,y=24,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

# Targeted Density Plots
### ECG BASIC
x <- list(v1=cpet_bwh$ecg_basic_fitted,v2=holdout$vo2)
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
cpet_bwh[, ':=' (pred_ecg_less14 = ifelse(ecg_basic_fitted_recal < 14,1,0),
                 true_less14 = ifelse(vo2 < 14,1,0),
                 pred_pclr_less14 = ifelse(pclr_basic_fitted_recal < 14,1,0),
                 pred_basic_less14 = ifelse(basic_fitted_recal < 14,1,0),
                 true_less25pct = ifelse(vo2 < 20.9,1,0),
                 pred_pclr_less25pct = ifelse(pclr_basic_fitted_recal < 20.9,1,0),
                 pred_basic_less25pct = ifelse(basic_fitted_recal < 20.9,1,0),
                 pred_ecg_less25pct = ifelse(ecg_basic_fitted_recal < 20.9,1,0))]

cpet_bwh[, ':=' (pred_ecg_less14_norecal = ifelse(ecg_basic_fitted < 14,1,0),
                 true_less14_norecal = ifelse(vo2 < 14,1,0),
                 pred_pclr_less14_norecal = ifelse(pclr_basic_fitted < 14,1,0),
                 pred_basic_less14_norecal = ifelse(basic_fitted < 14,1,0),
                 true_less25pct_norecal = ifelse(vo2 < 20.9,1,0),
                 pred_pclr_less25pct_norecal = ifelse(pclr_basic_fitted < 20.9,1,0),
                 pred_basic_less25pct_norecal = ifelse(basic_fitted < 20.9,1,0),
                 pred_ecg_less25pct_norecal = ifelse(ecg_basic_fitted < 20.9,1,0))]

## Create quartiles
cpet_bwh[, ':=' (true_quartiles = quantilize(vo2,4),
                 pred_quartile = quantilize(pclr_basic_fitted,4))]

## Run table
## 14
basic_table <- table2x2(data=cpet_bwh,disease='true_less14',
                        test='pred_basic_less14',key='linker_id')
basic <- epi.tests(basic_table)

ecg_table <- table2x2(data=cpet_bwh,disease='true_less14',
                      test='pred_ecg_less14',key='linker_id')
ecg <- epi.tests(ecg_table)

pclr_table <- table2x2(data=cpet_bwh,disease='true_less14',
                       test='pred_pclr_less14',key='linker_id')
pclr <- epi.tests(pclr_table)

## 25th percentile
basic_table_25 <- table2x2(data=cpet_bwh,disease='true_less25pct',
                           test='pred_basic_less25pct',key='linker_id')
basic_25 <- epi.tests(basic_table_25)

ecg_table_25 <- table2x2(data=cpet_bwh,disease='true_less25pct',
                         test='pred_ecg_less25pct',key='linker_id')
ecg_25 <- epi.tests(ecg_table_25)

pclr_table_25 <- table2x2(data=cpet_bwh,disease='true_less25pct',
                          test='pred_pclr_less25pct',key='linker_id')
pclr_25 <- epi.tests(pclr_table_25)

# Make ROCs
basic_roc <- roc(response=cpet_bwh$true_less14,predictor=cpet_bwh$pred_basic_less14,ci=TRUE)
ecg_roc <- roc(response=cpet_bwh$true_less14,predictor=cpet_bwh$pred_ecg_less14,ci=TRUE)
pclr_roc <- roc(response=cpet_bwh$true_less14,predictor=cpet_bwh$pred_pclr_less14,ci=TRUE)

basic_roc25 <- roc(response=cpet_bwh$true_less25pct,predictor=cpet_bwh$pred_basic_less25pct,ci=TRUE)
ecg_roc25 <- roc(response=cpet_bwh$true_less25pct,predictor=cpet_bwh$pred_ecg_less25pct,ci=TRUE)
pclr_roc25 <- roc(response=cpet_bwh$true_less25pct,predictor=cpet_bwh$pred_pclr_less25pct,ci=TRUE)

# Calculate SEs
pclr_basic_se <- boot_ci(status='true_less14',response1='pred_pclr_less14',response2='pred_basic_less14',data=cpet_bwh,runs=1000)
pclr_ecg_se <- boot_ci(status='true_less14',response1='pred_pclr_less14',response2='pred_ecg_less14',data=cpet_bwh,runs=1000)

## Perform z-tests
z_pclr_basic <- ((pclr_roc$auc - basic_roc$auc) / sd(pclr_basic_se))
p_pclr_basic <- 2*(1-pnorm(abs(z_pclr_basic)))

z_pclr_ecg <- ((pclr_roc$auc - ecg_roc$auc) / sd(pclr_ecg_se))
p_pclr_ecg <- 2*(1-pnorm(abs(z_pclr_ecg)))

# Calculate SEs
pclr_basic_se25 <- boot_ci(status='true_less25pct',response1='pred_pclr_less25pct',response2='pred_basic_less25pct',data=cpet_bwh,runs=1000)
pclr_ecg_se25 <- boot_ci(status='true_less25pct',response1='pred_pclr_less25pct',response2='pred_ecg_less25pct',data=cpet_bwh,runs=1000)

## Perform z-tests
z_pclr_basic25 <- ((pclr_roc25$auc - basic_roc25$auc) / sd(pclr_basic_se25))
p_pclr_basic25 <- 2*(1-pnorm(abs(z_pclr_basic25)))

z_pclr_ecg25 <- ((pclr_roc25$auc - ecg_roc25$auc) / sd(pclr_ecg_se25))
p_pclr_ecg25 <- 2*(1-pnorm(abs(z_pclr_ecg25)))

## NRI
pclr_basic_reclass <- nribin(event=cpet_bwh$true_less14,
                             p.std=cpet_bwh$pred_basic_less14,p.new=cpet_bwh$pred_pclr_less14,
                             cut=0.5,niter=100)

pclr_ecg_reclass <- nribin(event=cpet_bwh$true_less14,
                           p.std=cpet_bwh$pred_ecg_less14,p.new=cpet_bwh$pred_pclr_less14,
                           cut=0.5,niter=100)

pclr_basic_reclass25 <- nribin(event=cpet_bwh$true_less25pct,
                               p.std=cpet_bwh$pred_basic_less25pct,p.new=cpet_bwh$pred_pclr_less25pct,
                               cut=0.5,niter=100)

pclr_ecg_reclass25 <- nribin(event=cpet_bwh$true_less25pct,
                             p.std=cpet_bwh$pred_ecg_less25pct,p.new=cpet_bwh$pred_pclr_less25pct,
                             cut=0.5,niter=100)

## Error bars
# Define vars
cpet_bwh[,':='(vo2_decile = quantilize(vo2,10),
               abs_error = pclr_basic_fitted - vo2,
               rel_error = ((pclr_basic_fitted - vo2)/vo2)*100)]

# Plot absolute
y_vals <- cpet_bwh[,mean(abs_error),by='vo2_decile']
setkey(y_vals,vo2_decile)

pdf(file='abs_error_bwh.pdf',height=5,width=8,pointsize=5)
par(mar=c(10,7,1,2),oma=c(1,1,1,1))

coords <- barplot(height=y_vals$V1,col='#2c7fb88C',border=NA,
                  bty='n',xaxt='n',yaxt='n',xlim=c(0,12),
                  ylim=c(-10,10),xlab='',ylab='')

axis(1,at=coords,cex.axis=2.4,labels=rep('',length(coords)))
axis(2,cex.axis=2.4,at=seq(-10,10,5),las=2)

mtext("Error (mL/kg/min)",2,line=4.5,cex=2.2)
mtext("Decile of Peak VO2",1,line=8.5,cex=2.2)

labels <- c('<10.0','10.0-11.9','12.0-13.2','13.3-14.8','14.9-16.4',
            '16.5-18.0','18.1-19.8','19.9-22.3','22.4-26.6','>26.6')

text(x = coords,
     y = par("usr")[3] - 0.5,
     labels = labels,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.8)

dev.off()

# Plot relative
y_vals <- cpet_bwh[,mean(rel_error),by='vo2_decile']
setkey(y_vals,vo2_decile)

pdf(file='rel_error_bwh.pdf',height=5,width=8,pointsize=5)
par(mar=c(10,7,1,2),oma=c(1,1,1,1))

coords <- barplot(height=y_vals$V1,col='#d7301f8C',border=NA,
                  bty='n',xaxt='n',yaxt='n',xlim=c(0,12),
                  ylim=c(-75,100),xlab='',ylab='')

axis(1,at=coords,cex.axis=2.4,labels=rep('',length(coords)))
axis(2,cex.axis=2.4,at=seq(-75,100,25),las=2)

mtext("Error (%)",2,line=4.5,cex=2.2)
mtext("Decile of Peak VO2",1,line=8.5,cex=2.2)

labels <- c('<10.0','10.0-11.9','12.0-13.2','13.3-14.8','14.9-16.4',
            '16.5-18.0','18.1-19.8','19.9-22.3','22.4-26.6','>26.6')

text(x = coords,
     y = par("usr")[3] - 3.8,
     labels = labels,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.8)

dev.off()

## Error bars
# Define vars
cpet_bwh[,':='(abs_error_recal = pclr_basic_fitted_recal - vo2,
               rel_error_recal = ((pclr_basic_fitted_recal - vo2)/vo2)*100)]

# Plot absolute
y_vals <- cpet_bwh[,mean(abs_error_recal),by='vo2_decile']
setkey(y_vals,vo2_decile)

pdf(file='abs_error_bwh_recal.pdf',height=5,width=8,pointsize=5)
par(mar=c(10,7,1,2),oma=c(1,1,1,1))

coords <- barplot(height=y_vals$V1,col='#2c7fb88C',border=NA,
                  bty='n',xaxt='n',yaxt='n',xlim=c(0,12),
                  ylim=c(-10,10),xlab='',ylab='')

axis(1,at=coords,cex.axis=2.4,labels=rep('',length(coords)))
axis(2,cex.axis=2.4,at=seq(-10,10,5),las=2)

mtext("Error (mL/kg/min)",2,line=4.5,cex=2.2)
mtext("Decile of Peak VO2",1,line=8.5,cex=2.2)

labels <- c('<10.0','10.0-11.9','12.0-13.2','13.3-14.8','14.9-16.4',
            '16.5-18.0','18.1-19.8','19.9-22.3','22.4-26.6','>26.6')

text(x = coords,
     y = par("usr")[3] - 0.7,
     labels = labels,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.8)

dev.off()

# Plot relative
y_vals <- cpet_bwh[,mean(rel_error_recal),by='vo2_decile']
setkey(y_vals,vo2_decile)

pdf(file='rel_error_bwh_recal.pdf',height=5,width=8,pointsize=5)
par(mar=c(10,7,1,2),oma=c(1,1,1,1))

coords <- barplot(height=y_vals$V1,col='#d7301f8C',border=NA,
                  bty='n',xaxt='n',yaxt='n',xlim=c(0,12),
                  ylim=c(-50,50),xlab='',ylab='')

axis(1,at=coords,cex.axis=2.4,labels=rep('',length(coords)))
axis(2,cex.axis=2.4,at=seq(-50,50,25),las=2)

mtext("Error (%)",2,line=4.5,cex=2.2)
mtext("Decile of Peak VO2",1,line=8.5,cex=2.2)

labels <- c('<10.0','10.0-11.9','12.0-13.2','13.3-14.8','14.9-16.4',
            '16.5-18.0','18.1-19.8','19.9-22.3','22.4-26.6','>26.6')

text(x = coords,
     y = par("usr")[3] - 2.4,
     labels = labels,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.8)

dev.off()

## BASIC
cpet_bwh[,diff_pair := vo2 - basic_fitted]
cpet_bwh[,avg_pair := (vo2+basic_fitted)/2]

mod <- lm(cpet_bwh$diff_pair ~ cpet_bwh$avg_pair)
sy2 <- sum((cpet_bwh$diff_pair)**2)
ssr <- sum(mod$residuals**2)
df <- nrow(cpet_bwh)-2
msr <- ssr/df
bb <- (sy2-ssr)/(2*msr)
p_basic_bwh <- 1-pf(bb,df1=2,df2=df)

## ECG
cpet_bwh[,diff_pair := vo2 - ecg_basic_fitted]
cpet_bwh[,avg_pair := (vo2+ecg_basic_fitted)/2]

mod <- lm(cpet_bwh$diff_pair ~ cpet_bwh$avg_pair)
sy2 <- sum((cpet_bwh$diff_pair)**2)
ssr <- sum(mod$residuals**2)
df <- nrow(cpet_bwh)-2
msr <- ssr/df
bb <- (sy2-ssr)/(2*msr)
p_ecg_bwh <- 1-pf(bb,df1=2,df2=df)

## PCLR
cpet_bwh[,diff_pair := vo2 - pclr_basic_fitted]
cpet_bwh[,avg_pair := (vo2+pclr_basic_fitted)/2]

mod <- lm(cpet_bwh$diff_pair ~ cpet_bwh$avg_pair)
sy2 <- sum((cpet_bwh$diff_pair)**2)
ssr <- sum(mod$residuals**2)
df <- nrow(cpet_bwh)-2
msr <- ssr/df
bb <- (sy2-ssr)/(2*msr)
p_pclr_bwh <- 1-pf(bb,df1=2,df2=df)

# Targeted Plots for old school comparisons
## Original
### JONES
pdf(file='jones_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh$jones,y=cpet_bwh$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#0000008C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### JONES
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','jones')]
cpet_bwh[,diff_pair := vo2 - jones]
setkey(cpet_bwh,mean_pair)
nrow(cpet_bwh[mean_pair <= 0])

png('ba_vo2_jones_bwh.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#0000008C"

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=15,y=26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=15,y=28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### WASSERMAN
pdf(file='wasserman_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh$wasserman,y=cpet_bwh$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#d7301f8C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### WASSERMAN
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','wasserman')]
cpet_bwh[,diff_pair := vo2 - wasserman]
setkey(cpet_bwh,mean_pair)
nrow(cpet_bwh[mean_pair <= 0])

png('ba_vo2_wasserman_bwh.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#d7301f8C"

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=15,y=26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=15,y=28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### FRIEND
pdf(file='friend_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh$silva,y=cpet_bwh$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,60),ylim=c(0,60),
     pch=19,col='#8dd3c78C')
axis(1,at=seq(0,60,10),cex.axis=1.6)
axis(2,at=seq(0,60,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,61,61,lty=5)
dev.off()

### FRIEND
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','silva')]
cpet_bwh[,diff_pair := vo2 - silva]
setkey(cpet_bwh,mean_pair)
nrow(cpet_bwh[mean_pair <= 0])

png('ba_vo2_friend_bwh.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#8dd3c78C"

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-50,20),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
segments(0,upper,51,upper,lty=5,lwd=1,col='black')
segments(0,lower,51,lower,lty=5,lwd=1,col='black')
segments(0,mean_diff,51,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-52,at=seq(0,50,10))
axis(2,cex.axis=1.6,las=2,pos=0,at=seq(-50,20,10))
mtext(side=1,"Mean of paired VO2 values (mL/kg/min)",cex=1.8,line=2.8)
mtext(side=2,"Difference in paired VO2 values (mL/kg/min)",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=35,y=-46,labels=paste0("Limits of Agreement: ",round(lower,2),' to ',round(upper,2)),cex=1.4)
text(x=35,y=-48,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

## RECALIBRATED
### JONES
pdf(file='jones_recal_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh$jones_recal,y=cpet_bwh$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#0000008C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### JONES
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','jones_recal')]
cpet_bwh[,diff_pair := vo2 - jones_recal]
setkey(cpet_bwh,mean_pair)
nrow(cpet_bwh[mean_pair <= 0])

png('ba_vo2_jones_recal_bwh.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#0000008C"

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=35,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=35,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### WASSERMAN
pdf(file='wasserman_recal_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh$wasserman_recal,y=cpet_bwh$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#d7301f8C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,51,51,lty=5)
dev.off()

### WASSERMAN
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','wasserman_recal')]
cpet_bwh[,diff_pair := vo2 - wasserman_recal]
setkey(cpet_bwh,mean_pair)
nrow(cpet_bwh[mean_pair <= 0])

png('ba_vo2_wasserman_recal_bwh.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#d7301f8C"

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=15,y=26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=15,y=28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

### FRIEND
pdf(file='friend_recal_bwh.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=cpet_bwh$silva_recal,y=cpet_bwh$vo2,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,50),ylim=c(0,50),
     pch=19,col='#8dd3c78C')
axis(1,at=seq(0,50,10),cex.axis=1.6)
axis(2,at=seq(0,50,10),cex.axis=1.6,las=2)
mtext('Predicted peak VO2 (mL/kg/min)',1,line=3.2,cex=1.8)
mtext("True peak VO2 (mL/kg/min)",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,71,71,lty=5)
dev.off()

### FRIEND
cpet_bwh[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('vo2','silva_recal')]
cpet_bwh[,diff_pair := vo2 - silva_recal]
setkey(cpet_bwh,mean_pair)
nrow(cpet_bwh[mean_pair <= 0])

png('ba_vo2_friend_recal_bwh.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#8dd3c78C"

# Plot
plot(x=cpet_bwh[mean_pair > 0]$mean_pair,y=cpet_bwh[mean_pair > 0]$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(0,50),ylim=c(-30,30),cex=0.8)

# CI Lines
mean_diff <- mean(cpet_bwh$diff_pair)
upper <- mean(cpet_bwh$diff_pair)+1.96*sd(cpet_bwh$diff_pair)
lower <- mean(cpet_bwh$diff_pair)-1.96*sd(cpet_bwh$diff_pair)
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
agree <- nrow(cpet_bwh[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(cpet_bwh)
text(x=35,y=-26,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=35,y=-28,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

############# Secondary/subset analyses
### LBBB
no_lbbb <- cpet_bwh[lbbb==FALSE]
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
no_rbbb <- cpet_bwh[rbbb==FALSE]
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
no_lvh <- cpet_bwh[lvh==FALSE]
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
no_af <- cpet_bwh[af==FALSE]
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
bwh_lvef <- fread(file='bwh_lvef.csv')
no_lowef <- cpet_bwh[!(linker_id %in% bwh_lvef[Car.Report_Text.lvef<50]$id)] # N= 1076 - 234 = 842
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
