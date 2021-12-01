#!/usr/bin/env Rscript
library(biomod2)

# ====== Control panel ======


datapath <- "../Input/"
trainpath <- "./Output_train/"
outpath <- "./Output_forward/"

gname <- 'Bh04v01'
gh <- 4
gv <- 1
N<-6000
SCALE <- 120
pred_var <- c('dem','twi','slp','asp','DD5','PAS','AHM')

# ============= Read input ================
df_fw <- read.csv(paste(datapath,'LC_TOPO_',gname,'.csv',sep=""))
df_fw[c('twi','slp','asp')] <- df_fw[c('twi','slp','asp')]/100 # scale
df_clm <- read.csv(paste(datapath,'CLIM_',gname,'.csv',sep=""))
df_fw['R'] <- df_fw['r']%/%SCALE
df_fw['C'] <- df_fw['c']%/%SCALE
df_fw <- merge(df_fw, df_clm, by=c("R","C"),all=TRUE)
df_fw <- df_fw[complete.cases(df_fw),!(names(df_fw) %in% c('R','C'))]
print(head(df_fw)); print(nrow(df_fw))

df_fw['r'] <- df_fw['r']+(gv-0)*N # from tile to domain, minimum gv = 0
df_fw['c'] <- df_fw['c']+(gh-1)*N # from tile to domain, minimum gh = 1

myRespXY_fw <- df_fw[,c("r","c")]
myExpl_fw <- df_fw[,pred_var]
rm(df_fw)

# =============== Forward run ==================
setwd(trainpath)
tagid <- 0 # id of end member of training
load(paste('Train',tagid,'.Rdata',sep=''))
    
    
#=========== BIOMOD prediction =============
myBiomodProjectionFuture <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                              new.env = myExpl_fw,
                                              proj.name = 'future',
                                              selected.models = 'all',
                                              binary.meth = c('TSS','ROC'),
                                              xy.new.env = myRespXY_fw,
                                              compress = FALSE,
                                              build.clamping.mask = TRUE)
    
myProjDFFuture <- as.data.frame(get_predictions(myBiomodProjectionFuture),xy=TRUE)
    
#============ Write output ==================
write.csv(myProjDFFuture,paste('../',outpath,'SUITB_',gname,'_',tagid,'.csv',sep=''))
write.csv(myRespXY_fw,paste('../',outpath,'SUITB_',gname,'_xy.csv',sep=''))




