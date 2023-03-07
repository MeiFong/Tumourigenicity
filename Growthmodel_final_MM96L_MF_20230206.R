library(readxl)
library(saemix)
library(gridExtra) 
library(tidyverse)
library(magrittr)
#link: https://www.r-bloggers.com/2019/09/fitting-complex-mixed-models-with-nlme-example-4/
library(nlme)
#devtools::install_github("onofriAndreaPG/aomisc")
library(aomisc)
library(emmeans)

library(nlstools)

rm(list=ls())
setwd("L:/Lab_Stats_Unit/Satomi/QIMR/MeiFong/tumour growth model")
#setwd("L:/Joint_Projects/Analysis_MF/20221220")

#MM649 <- read_excel("Mouse data_20230110_so.xlsx", sheet = "MM649_2 (2)")
MM96L <- read_excel("Mouse data_20230110_so.xlsx", sheet = "MM96L_2")

# names(MM649)
# table(MM649$Groups)
names(MM96L)
table(MM96L$Groups)

# Data prep to long
# df_long_649 <-
#   MM649 %>%
#   mutate(mouseID = paste(Groups, Mouse, sep ="_")) %>%
#   mutate(Groups = str_replace(Groups, "Dox", "D")) %>%
#   dplyr::mutate(ID = row_number()%>% as.factor()) %>%
#   pivot_longer(-c(Groups, Mouse, COD, mouseID, ID), names_to = "time",
#                names_prefix ="X", names_transform = list(time = as.integer),
#                values_to = "value") %>%
#   mutate(Treat = if_else(grepl("D",Groups ), "D", "C"),
#          Group = str_remove(Groups, "D")) %>%
#   mutate_if(is.character, as.factor) %>%
#   mutate(week = time/7) %>%
#   filter(!is.na(value))   # remove NA timepoints
# 
# table(df_long_649$Groups, df_long_649$Group)
# table(df_long_649$Groups, df_long_649$Treat)
# 
# table(df_long_649$Groups, df_long_649$Group)
# table(df_long_649$Groups, df_long_649$Treat)
# 
df_long_96 <-
  MM96L %>%
  mutate(mouseID = paste(Groups, Mouse, sep ="_")) %>%
  mutate(Groups = str_replace(Groups, "Dox", "D")) %>%
  dplyr::mutate(ID = row_number()%>% as.factor()) %>%
  pivot_longer(-c(Groups, Mouse, COD, mouseID, ID), names_to = "time",
               names_prefix ="X", names_transform = list(time = as.integer),
               values_to = "value") %>%
  mutate(Treat = if_else(grepl("D",Groups ), "D", "C"),
         Group = str_remove(Groups, "D")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(week = time/7) %>%
  filter(!is.na(value)) # remove NA timepoints

table(df_long_96$Groups, df_long_96$Group)
table(df_long_96$Groups, df_long_96$Treat)

# Visualization 
# df_long_649 %>%
#   ggplot(aes(y=value ,x = time , group= ID, color = Treat )) +
#   geom_line() + 
#   facet_grid(Group ~ Treat) +
#   geom_vline(xintercept = 13, linetype = "dashed")

df_long_96 %>%
  ggplot(aes(y=value ,x = time , group= ID, color = Treat )) +
  geom_line() +
  facet_grid(Group ~ Treat) +
  geom_vline(xintercept = 13, linetype = "dashed")


# indiviudal plot
# df_long_649 %>%
#   ggplot(aes(y=value ,x = time , group= ID, color = Treat )) +
#   geom_line() + 
#   facet_wrap(~ID) +
#   geom_vline(xintercept = 13, linetype = "dashed")



##########################################################
# subset data:
# data preparation: this will remove one sample which wasn't measured
# df_long_649_t13 <- df_long_649 %>% filter(time >=13) %>% mutate(time = time -13) 
df_long_96_t13 <- df_long_96 %>% filter(time >=13) %>% mutate(time = time -13) %>% mutate(week = time/7)


##################################################################
# Exponential growth
#########################################
# 
# mnaiv_exp_649 <- drm(value ~ time , fct =DRC.expoGrowth() , 
#                      data = df_long_649_t13   ,  
#                      curveid = Groups,
#                      pmodels = c( ~ Groups,    # initial volume at day 13
#                                   ~ Groups), # k : growth rate:relative increase per day
# ) # Transform-Both-Side (log) to account for heteroscedasticity
# summary(mnaiv_exp_649)
# coef(mnaiv_exp_649)
# plot(mnaiv_exp_649,log="")
# AIC(mnaiv_exp_649)

########################################################
# Prepare functions
###################################
# exponential function 
f.exp<-function(t,a,k){
  a*exp(k*t)
}

# # Gompertz funciton
# f.gomp3 <- function(t, b,  d, e){
#   ypred <- d*exp(-exp(-b*(t-e)))  # formula of 3P Gompertz fucntion 
#   return(ypred)
# }
# 
# # logistic function 
# f.logit3 <- function(t, b, d, e){
#   ypred <- d/(1+exp(-b*(t-e)))  # formula of 3P logistic fucntion
#   return(ypred)
# }

##################################################
# Get rough starting values 
##################################
#linearlize the formula by log
# Exponential: value ~ a*exp(k * time) => log(value) ~ log(a) + k*time
lin.exp <- lm(log(value+1) ~ time, data = df_long_96_t13)
# back-transform parameters: exp(intercept) = a, beta = k
exp.start= list(a= exp(coef(lin.exp))[1], k =coef(lin.exp)[2])

# # Gomperz : 
# xydata <- sortedXyData("time","value",df_long_96_t13) # provides the average y for each unique value of x 
# x = xydata["x"]
# y = xydata["y"]
# # Grid search - Gomperz
# #choose grid of a and b values
# bseq = seq(0,1, 0.02)
# dseq = seq(100,500,50)
# eseq = seq(10,50,2)
# 
# n.b = length(bseq)
# n.d = length(dseq)
# n.e = length(eseq)
# SSout = matrix(0,n.d*n.b*n.e,4) #matrix to save output
# cnt = 1
# for (k in 1:n.b){
#   for (j in 1:n.d){
#     for(i in 1:n.e) {
#       
#       ypred = f.gomp3(x, bseq[k],dseq[j],eseq[i]) #evaluate model w/ these parms
#       ss = sum((y-ypred)^2)  #this is our SSE objective function
#       #save values of a, b, and SSE
#       SSout[cnt,1]=bseq[k]
#       SSout[cnt,2]=dseq[j]
#       SSout[cnt,3]=eseq[i]
#       SSout[cnt,4]=ss
#       cnt = cnt +1
#     }
#   }
# }
# SSout.exp = SSout
# #find minimum SSE and associated a,b values
# mn_indx = which.min(SSout[,4])
# gomp3.start <- list(b=SSout[mn_indx,1], 
#                     d=SSout[mn_indx,2], 
#                     e=SSout[mn_indx,3] )
# 
# # Grid search - Logistic 3P
# #choose grid of a and b values
# bseq = seq(0,1, 0.02)
# dseq = seq(100,300, 20)
# eseq = seq(20,50,2)
# n.b = length(bseq)
# n.d = length(dseq)
# n.e = length(eseq)
# SSout2 = matrix(0,n.d*n.b*n.e,4) #matrix to save output
# cnt = 1
# for (k in 1:n.b){
#   for (j in 1:n.d){
#     for(i in 1:n.e) {
#       
#       ypred = f.logit3(x, bseq[k],dseq[j],eseq[i]) #evaluate model w/ these parms
#       ss = sum((y-ypred)^2)  #this is our SSE objective function
#       #save values of a, b, and SSE
#       SSout2[cnt,1]=bseq[k]
#       SSout2[cnt,2]=dseq[j]
#       SSout2[cnt,3]=eseq[i]
#       SSout2[cnt,4]=ss
#       cnt = cnt +1
#     }
#   }
# }
# SSout.logit3 = SSout
# 
# #find minimum SSE and associated parameter values
# mn_indx = which.min(SSout.exp[,4])
# gomp3.start <- list(b=SSout.exp[mn_indx,1], 
#                     d=SSout.exp[mn_indx,2], 
#                     e=SSout.exp[mn_indx,3] )
# gomp3.start
# 
# logit3.start <- list(b=SSout2[mn_indx,1], 
#                      d=SSout2[mn_indx,2], 
#                      e=SSout2[mn_indx,3] )
# logit3.start

#######################################################
# Fit naive models for starting values
###########################################
fit.nls.exp <-
  nlsList(value ~ a*exp(k*time)| Groups,
          data =df_long_96_t13 ,
          start = exp.start )


# fit.nls.gomp3 <-
#   nlsList(value ~ d*exp(-exp(-b*(time-e)))| Groups,
#           data =df_long_96_t13 ,
#           start = gomp3.start,
#           control = list(maxiter = 100000,  minFactor = 1/100000000000))
# 
# 
# fit.nls.logit3 <-
#   nlsList(value ~ d/(1+exp(-b*(time-e)))| Groups,
#           data =df_long_96_t13 ,
#           start = gomp3.start,
#           control = list(maxiter = 100000,  minFactor = 1/100000000000))

###########################################################
# Non-linear Mixed effect model
########################################################
# Exponential
# fixed effect of k
model_expo <- 
  nlme(value ~f.exp(time, a, k), data = df_long_96_t13, 
       #random = a+k ~ 1|ID,
       random = list(ID = pdDiag(a+k   ~ 1)),  # pdDiag : uncorrelated REs
       fixed = list(a ~ 1,    
                    k ~ Group*Treat),  
       #weights= varIdent(form ~1|Groups),
       weights = varPower(),
       #start = fixed.effects(model_expo), #initial estimates from naive model
       start=c(coef(fit.nls.exp)[1,1], rep(0,0),
               coef(fit.nls.exp)[1,2], rep(0,7)),
      control = nlmeControl(msMaxIter = 100000, maxIter = 100000),  #
       verbose = TRUE
  )




fixed.effects(model_expo)
nlme::intervals(model_expo)
AIC(model_expo)
VarCorr(model_expo)

# Gompertz

# model_gomp3 <- 
#   nlme(value ~ f.gomp3(time, b,d,e), data = df_long_96_t13, 
#        #  random = b +d ~ 1|ID,
#        random = list(ID = pdDiag(d~ 1)),
#        fixed = list(b ~ Groups,    
#                     d ~ Groups,
#                     e ~ 1),  
#        weights = varPower(),
#        # start = coef(mnaiv_exp_649), #initial estimates from naive model
#        start=c(coef(fit.nls.gomp3)[1,1], rep(0,7),
#                coef(fit.nls.gomp3)[1,2], rep(0,7),
#                coef(fit.nls.gomp3)[1,3], rep(0,0)),
#        control = list (maxIter = 1000, MaxEval = 1000),  #
#        verbose = TRUE
#   )
# 
# fixed.effects(model_gomp3)
# AIC(model_gomp3)
# 
# 
# # Logistic 3P
# 
# model_logistic3 <- 
#   nlme(value ~ f.logit3(time, b,d,e), data = df_long_96_t13, 
#        #random = a+k ~ 1|ID,
#        random = list(ID = pdDiag(d ~ 1)),
#        fixed = list(b ~ Groups,    
#                     d ~ Groups,
#                     e ~ 1),  
#        weights = varPower(),
#        # start = coef(mnaiv_exp_649), #initial estimates from naive model
#        start=c(coef(fit.nls.logit3)[1,1], rep(0,7),
#                coef(fit.nls.logit3)[1,2], rep(0,7),
#                coef(fit.nls.logit3)[1,3], rep(0,0)),
#        control = list (msMaxIter = 1000, msMaxEval = 1000),  #
#        verbose = TRUE
#   )
# 
# fixed.effects(model_logistic3)
# AIC(model_logistic3)


# model comparison 
anova(model_expo, model_gomp3, model_logistic3)


# marginal mean
# Pair-wise comparison
# margin.a <- emmeans(model_expo,  ~ Groups, param = "a")
# contrast.a <- contrast(margin.a, method= "revpairwise", adjust = "none") %>% as.data.frame()

margin.k <- emmeans(model_expo,  ~ Group*Treat, data = df_long_96_t13, param = "k")
#Pair-wise comparison (original code)
contrast.k <- contrast(margin.k,  method = "revpairwise", by = "Group", infer = T, adjust = "none") %>% as.data.frame()
#Control vs Treated by group
write.csv(contrast.k, file = "MM96L_tumour_vsControl.csv", row.names =F)

#contrast of contrast (D -c) between groups
contrast.k_contrast <- contrast(margin.k, method= "revpairwise", interaction = T, infer = T, adjust = "none")
write.csv(contrast.k_contrast, file = "MM96L_tumour_vsGroup.csv", row.names =F)

# constant of proportionality (k): a*exp(kt)
df_margin.k = 
  margin.k %>% as.data.frame() %>% dplyr::select(-df) %>%
  dplyr::rename(k = emmean) 

# Growth rate r: a*(1-r)^t
df_GR.r <- margin.k %>% as.data.frame() %>%
  mutate_at(c("emmean", "lower.CL", "upper.CL"), function(x) exp(x)-1) %>% # relative growth rate per day
  dplyr::select(-SE, -df) %>%
  dplyr::rename(GR.r = emmean)  

# Doubling time: log(2)/k
df_Doubling.time <-
  margin.k %>% as.data.frame() %>%
  mutate_at(c("emmean", "lower.CL", "upper.CL"), function(x) log(2)/x) %>% # relative growth rate per day
  dplyr::select(-SE, -df) %>%
  dplyr::rename(Doubling.time = emmean, lower.CL= upper.CL, upper.CL= lower.CL) 


# exponentiated contrast estimate  = relative ratio of GR between 2 groups

emmip(margin.k, ~Groups, CI=T)


# individual parameter estimate
df_k.fixed <- margin.k %>% as.data.frame() %>% dplyr::rename(k.fixed = emmean) %>% dplyr::select(Groups, k.fixed)
# df_a.fixed <- margin.a %>% as.data.frame() %>% dplyr::rename(a.fixed = emmean) %>% dplyr::select(Groups, a.fixed)
df_random <- random.effects(model_expo, augFrame = T, omitGroupingFactor=F) %>% 
  dplyr::rename(k.random = "k.(Intercept)", 
                #a.random= "a.(Intercept)"
  )


df_ind.pram <-
  merge(df_k.fixed, df_random, by ="Groups") %>%
  mutate(k.ind = k.fixed + k.random,
         # a.ind = a.fixed + a.random
  ) 



#Plot of Growth rate (per day)
df_GR.r %>%
  mutate(Treat = if_else(grepl("D",Groups ), "D", "C"),
         Group = str_remove(Groups, "D")) %>%
  ggplot(aes(y= GR.r, x = Treat)) +
  geom_point() +
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL, x= Treat, color =Treat)) +
  geom_jitter(data=df_ind.pram, aes(y=exp(k.ind)-1, x = Treat, color = Treat), alpha = 0.5, width =0.3) +
  # scale_y_continuous(limits=c(1,1.1)) +
  facet_wrap(~Group)

df_GR.r %>%
  mutate(Treat = if_else(grepl("D",Groups ), "D", "C"),
         Group = str_remove(Groups, "D")) %>%
  ggplot(aes(y= GR.r, x = Group)) +
  geom_point() +
  geom_errorbar(aes(ymax = upper.CL, ymin = lower.CL, x= Group, color =Group)) +
  geom_jitter(data=df_ind.pram, aes(y=exp(k.ind)-1, x = Group, color = Group), alpha = 0.5, width =0.3) +
  #scale_y_continuous(limits=c(1,1.1)) +
  facet_wrap(~Treat)


df_long_96_t13 %>%
  merge(., df_ind.pram %>% dplyr::select(ID, k.ind), by ="ID") %>%
  filter(Group== "shMITF") %>%
  ggplot(aes(y=value ,x = time , group= ID, color = Treat ,
             #alpha = round(k.ind, 3)
            
  )) +
  geom_line() + 
  facet_grid(Group ~ Treat) 




# fitted line by ID
df_pred.expo <- predict(model_expo, level =c(0,1)) %>% dplyr::rename(ID.1 = ID)
df_pred.gomp <- predict(model_gomp3, level =c(0,1)) %>% dplyr::rename(ID.1 = ID)


df_fitted.ID <- cbind(df_long_96_t13, df_pred.expo)

# fitted line by Groups
df_groups <- df_long_96_t13 %>% dplyr::select(Groups, time) %>% distinct() 

df_fitted.groups <-  cbind(df_long_96_t13 %>% dplyr::select("Groups", "Group", "Treat", "time") %>% distinct(),
                           predicted.groups=predict(model_expo, newdata =df_groups, level =c(0)) ) 
  

# Fitted line plot
lplot_pred.expo <-
  df_fitted.ID %>%
  ggplot(aes(y= predict.ID, x = time, group = ID, color = as.factor(Mouse))) +
  geom_line(size=0.5) +
  geom_point(aes(y=value , x=time, group = ID, color =as.factor(Mouse)),
             alpha = 0.5, size =1) +
  geom_line(data = df_fitted.groups, aes(y=predicted.groups, x= time, group = Groups),
            size =1, color ="red") +
  facet_grid(Group ~ Treat) +
  scale_y_continuous(limits=c(0,1000)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y= "Volume", x = "Days since Treatment") 

lplot_pred.expo

# Doubling time
log(2)/  as.data.frame(margin.k)$emmean


###########################################
# Save results for MM96L
path = "L:/Lab_Stats_Unit/Satomi/QIMR/MeiFong/tumour growth model/Exponential model_Growth.model/"
write.csv(df_margin.k, file = paste0(path, "parameter.k_Expomodel_96L", ".csv" ), row.names = F)
write.csv(df_GR.r, file = paste0(path, "GR.perday_Expomodel_96L", ".csv" ), row.names = F)
write.csv(df_Doubling.time, file = paste0(path, "Doubling.time_Expomodel_96L", ".csv" ), row.names = F)

write.csv(contrast.k, file = paste0(path, "Contrast_GR_Expomodel_96L", ".csv" ), row.names = F)

#dataset for plots
write.csv(df_fitted.ID, file = paste0(path, "df.input.predicted.id_Expomodel_96L", ".csv" ), row.names = F)
write.csv(df_fitted.groups, file = paste0(path, "df.predicted.groups_Expomodel_96L", ".csv" ), row.names = F)


# write.csv(margin.a %>% as.data.frame(), file = paste0(path, "Initial.vol_Expomodel", ".csv" ), row.names = F)
# write.csv(contrast.a, file = paste0(path, "Contrast_Initial.vol_Expomodel", ".csv" ), row.names = F)

ggsave(paste0(path, "lineplot_fitted_Expomodel_96L", ".png"),
       plot =lplot_pred.expo , device = "jpeg", 
       width = 6, height = 8, units = "in", dpi = 300)


#############################################
# for plot only
####################
rm(list=ls())
setwd("L:/Joint_Projects/Analysis_MF/Exponential model_Tumor.Growth/MM96L")

# input dataset and fitted data from exponential model ("predict.ID")
df_fitted.ID.96 <- read.csv("df.input.predicted.id_Expomodel_96L.csv")
# predicted line for each group
df_fitted.groups.96 <- read.csv("df.predicted.groups_Expomodel_96L.csv")

# note:
# time : calculated as days - 13 (days since day 13)
# value : tumour volume
# predict.ID: predicted value for each tumour ID


# Fitted line plot
lplot_pred.expo.96 <-
  df_fitted.ID.96 %>%
  ggplot(aes(y= predict.ID, x = time, group = ID, color = as.factor(Mouse))) +
  geom_line(size=0.5) +
  geom_point(aes(y=value , x=time, group = ID, color =as.factor(Mouse)),
             alpha = 0.5, size =1) +
  geom_line(data = df_fitted.groups.96, aes(y=predicted.groups, x= time, group = Groups),
            size =1, color ="red") +
  facet_grid(Group ~ Treat) +
  scale_y_continuous(limits=c(0,1000)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y= "Volume", x = "Days since Treatment") 

lplot_pred.expo.96
