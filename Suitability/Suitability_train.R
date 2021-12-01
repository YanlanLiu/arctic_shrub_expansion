library(biomod2)
arrayid<-0

inpath <- './Input_train/'
outpath <- './Output_train/'
pred_var <- c('dem','twi','slp','asp','DD5','PAS','AHM')
df <- read.csv(paste(inpath,'BIOMOD_train_demo.csv',sep=''))
myRespName <- paste('Train',arrayid,sep='')

print(myRespName)
print(head(df))
myResp <- as.numeric((df[,'lc0']>4) & (df[,'lc0']<8))
myRespXY <- df[,c("r","c")]
myExpl <- df[,pred_var]
rm(df)

# Print_Default_ModelingOptions()
#=========== Moddels ===============
setwd(outpath)
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
myBiomodData

myBiomodOption <- BIOMOD_ModelingOptions(RF = list(do.classif = TRUE,ntree = 10,mtry = 'default',nodesize = 5,maxnodes = NULL))

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                     models =  c('RF'),
                                     models.options = myBiomodOption,
                                     DataSplit=100,
                                     NbRunEval=1,
                                     models.eval.meth = c('TSS','ROC'),
                                     SaveObj = TRUE,
                                     rescal.all.models = FALSE,
                                     do.full.models = TRUE)

save(myBiomodModelOut,file=paste(myRespName,'.Rdata',sep=''))

eval <- get_evaluations(myBiomodModelOut,as.data.frame=T)
print(eval$Testing.data)
