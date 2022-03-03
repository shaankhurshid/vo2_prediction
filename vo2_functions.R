# Bootstrap to compare MAE
compare_mae <- function(y,x1,x2,data,runs,size=nrow(data)){
  out <- list()
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    mae1 <- mean(abs(sample[,get(y)] - sample[,get(x1)]))
    mae2 <- mean(abs(sample[,get(y)] - sample[,get(x2)]))
    diff_mae <- mae1-mae2
    out[[i]] <- c(mae1,mae2,diff_mae)
    if(i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(do.call(rbind,out))
}

# Bootstrap for MAE
ci_mae <- function(y,x1,data,runs,size=nrow(data)){
  out <- list()
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    out[[i]] <- mean(abs(sample[,get(y)] - sample[,get(x1)]))
    if(i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(do.call(rbind,out))
}

# Compare VO2 max
compare_linear_vo2_max <- function(data,truth,inference,indicator,write_out=FALSE,
                           out_path=NULL,plot=FALSE){
  out <- list()
  for (i in 1:length(inference)){
    y <- data[,get(truth)]
    x <- data[,get(inference[i])]
    cor <- cor.test(x,y)
    mae <- mean(abs(x-y))
    out[[i]] <- c(inference[i],as.numeric(cor$estimate),as.numeric(cor$conf.int[1]),
                  as.numeric(cor$conf.int[2]),as.numeric(mae))
    
    if (plot==TRUE){
      pdf(paste0(out_path,inference[i],'_',indicator,'.pdf'),height=4,width=4,pointsize=5)
      par(mar=c(5,5,2,2),oma=c(1,1,1,1))
      plot(x=x,y=y,xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,8),ylim=c(0,8),
           pch=19,col='#2c7fb88C')
      axis(1,at=0:8,cex.axis=1.4)
      axis(2,at=0:8,cex.axis=1.4,las=2)
      mtext(paste0(inference[i],' predicted VO2'),1,line=4,cex=1.5)
      mtext('True VO2 max (l/min)',2,line=3.2,cex=1.5)
      text(6,2,labels=paste0('r=',round(cor$estimate,2),' (',round(cor$conf.int[1],2),
                             ' - ',round(cor$conf.int[2],2),')'),cex=1.2)
      text(6,1.7,labels=paste0('MAE=',round(mae,2)),cex=1.2)
      segments(-0.1,-0.1,8.1,8.1,lty=5)
      dev.off()
    }
  }
  output <- data.table(do.call(rbind,out))
  names(output) <- c('model','cor','lower','upper','mae')
  if (write_out==TRUE){
    write.csv(output,file=paste0(out_path,'output_',indicator,'.csv'),row.names=F)
  }
  return(output)
}

# Compare VO2
compare_linear_vo2 <- function(data,truth,inference,indicator,write_out=FALSE,
                           out_path=NULL,plot=FALSE){
  out <- list()
  for (i in 1:length(inference)){
    y <- data[,get(truth)]
    x <- data[,get(inference[i])]
    cor <- cor.test(x,y)
    mae <- mean(abs(x-y))
    out[[i]] <- c(inference[i],as.numeric(cor$estimate),as.numeric(cor$conf.int[1]),
                  as.numeric(cor$conf.int[2]),as.numeric(mae))
    
    if (plot==TRUE){
      pdf(paste0(out_path,inference[i],'_',indicator,'.pdf'),height=4,width=4,pointsize=5)
      par(mar=c(5,5,2,2),oma=c(1,1,1,1))
      plot(x=x,y=y,xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,90),ylim=c(0,90),
           pch=19,col='#2c7fb88C')
      axis(1,at=seq(0,90,10),cex.axis=1.4)
      axis(2,at=seq(0,90,10),cex.axis=1.4,las=2)
      mtext(paste0(inference[i],' predicted VO2'),1,line=4,cex=1.5)
      mtext('True VO2 (ml/kg/min)',2,line=3.2,cex=1.5)
      text(75,10,labels=paste0('r=',round(cor$estimate,2),' (',round(cor$conf.int[1],2),
                               ' - ',round(cor$conf.int[2],2),')'),cex=1.2)
      text(75,5,labels=paste0('MAE=',round(mae,2)),cex=1.2)
      segments(-0.1,-0.1,91,91,lty=5)
      dev.off()
    }
  }
  output <- data.table(do.call(rbind,out))
  names(output) <- c('model','cor','lower','upper','mae')
  if (write_out==TRUE){
    write.csv(output,file=paste0(out_path,'output_',indicator,'.csv'),row.names=F)
  }
  return(output)
}

# Compare VO2 pct
compare_linear_vo2_pct <- function(data,truth,inference,indicator,write_out=FALSE,
                           out_path=NULL,plot=FALSE){
  out <- list()
  for (i in 1:length(inference)){
    y <- data[,get(truth)]
    x <- data[,get(inference[i])]
    cor <- cor.test(x,y)
    mae <- mean(abs(x-y))
    out[[i]] <- c(inference[i],as.numeric(cor$estimate),as.numeric(cor$conf.int[1]),
                  as.numeric(cor$conf.int[2]),as.numeric(mae))
    
    if (plot==TRUE){
      pdf(paste0(out_path,inference[i],'_',indicator,'.pdf'),height=4,width=4,pointsize=5)
      par(mar=c(5,5,2,2),oma=c(1,1,1,1))
      plot(x=x,y=y,xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(0,225),ylim=c(0,225),
           pch=19,col='#2c7fb88C')
      axis(1,at=seq(0,225,25),cex.axis=1.4)
      axis(2,at=seq(0,225,25),cex.axis=1.4,las=2)
      mtext(paste0(inference[i],' predicted VO2 %'),1,line=4,cex=1.5)
      mtext('True VO2 %',2,line=3.2,cex=1.5)
      text(180,35,labels=paste0('r=',round(cor$estimate,2),' (',round(cor$conf.int[1],2),
                                ' - ',round(cor$conf.int[2],2),')'),cex=1.2)
      text(180,25,labels=paste0('MAE=',round(mae,2)),cex=1.2)
      segments(-0.1,-0.1,226,226,lty=5)
      dev.off()
    }
  }
  output <- data.table(do.call(rbind,out))
  names(output) <- c('model','cor','lower','upper','mae')
  if (write_out==TRUE){
    write.csv(output,file=paste0(out_path,'output_',indicator,'.csv'),row.names=F)
  }
  return(output)
}

# Function to make a 2x2 table
table2x2 <- function(data,disease,test,key){
  true_pos <- data[data[,get(disease)]==1]
  true_neg <- data[data[,get(disease)]==0]
  test_pos <- data[data[,get(test)]==1]
  test_neg <- data[data[,get(test)]==0]
  a <- nrow(true_pos[true_pos[,get(key)] %in% test_pos[,get(key)]])
  c <- nrow(true_pos[true_pos[,get(key)] %in% test_neg[,get(key)]])
  b <- nrow(true_neg[true_neg[,get(key)] %in% test_pos[,get(key)]])
  d <- nrow(true_neg[true_neg[,get(key)] %in% test_neg[,get(key)]])
  table <- matrix(c(a,c,b,d),ncol=2,nrow=2)
}

# Bootstrap function to obtain SE of CI
boot_ci <- function(status,response1,response2,data,runs,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    out[i] <- concordance(glm(sample[,get(status)] ~ sample[,get(response1)],data=sample))$concordance - 
                    concordance(glm(sample[,get(status)] ~ sample[,get(response2)],data=sample))$concordance
    if(i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(out)
}

# Bootstrap cmprsk
boot_cmp <- function(time,status,response,data,runs,exposure,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    x_matrix <- model.matrix(~age + sex + bmi + get(exposure),data=sample)[,-1]
    out[i] <- crr(ftime=sample[,get(time)],fstatus=sample[,get(status)],
                  failcode=1,cencode=0,cov1=x_matrix,variance=FALSE)$coef[4]
    if(i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(out)
}