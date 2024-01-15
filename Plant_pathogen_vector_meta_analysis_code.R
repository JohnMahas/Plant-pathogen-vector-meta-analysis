#Code used for the meta-analysis in R-syntax
#load data
datum_meta=read.csv(file.choose())
#check data
head(datum_meta)
#install and load packages
install.packages("ggalt")
install.packages("remotes")
remotes::install_github("wviechtb/metafor")
library(remotes)
library(metafor)
library(ggplot2)
#Calculate Fisher's Z using correlation coefficients
results_ri=escalc(measure="ZCOR",ri=r,ni=Sample_size,data=datum_meta)
write.csv(results_ri,"E:/AU_Classes/FISH 7350/Class_Project/Analyses/results_ri.csv")

#Calculate Fisher's Z using t,f,z,x2
results_tfzx2=escalc(measure="ZCOR",ti=t_z_f_x2,ni=Sample_size,data=datum_meta)
write.csv(results_tfzx2,"E:/AU_Classes/FISH 7350/Class_Project/Analyses/results_tfzx2_.csv")

#Calculate Fisher's Z using p values (p > 0.05 = 1)
results_pi_1=escalc(measure="ZCOR",pi=P_value_adj_1,ni=Sample_size,data=datum_meta)
write.csv(results_pi_1,"E:/AU_Classes/FISH 7350/Class_Project/Analyses/results_pi_1.csv")

#test fishers z calculated from p values (> 0.05 = 1) add random study
results_P1_random_study=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,data=datum_meta,method="REML",random=~1|Paper_Name_Short,slab=Paper_Name_Short)
summary(results_P1_random_study)
forest(results_P1_random_study)
funnel(results_P1_random_study,xlab="Effect size (Fisher's Z)")

#Fail-safe results
fsn(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,data=datum_meta)
#rank test
ranktest(results_P1_random_study,exact=FALSE)

forest(datum_meta$yi_from_ti_ri_pi_1, # These are effect sizes from each row in database
       datum_meta$vi_from_ti_ri_pi_1, # These are variances from each row in database
       annotate = FALSE,
       slab = results_P1_random_study$slab, 
       xlab = "ln(Response Ratio)", # Label for x-axis
       cex = .8, # Text side for study labels
       pch = 15, # shape of bars in forest plot
       cex.lab = 1,) # Size of x-axis label

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=scale
results_P1_random_study_ModScale=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                        mods=~Spatial_scale_m_mod,
                                        data=datum_meta,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModScale,digits=6)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=specialization
results_P1_random_study_ModSpec=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                       mods=~Specialist_generalist_both_mod-1,
                                       data=datum_meta,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModSpec)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=landscape
results_P1_random_study_ModLand=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                       mods=~Landscape_category_mod-1,
                                       data=datum_meta,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModLand)
#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=measurement_mod

results_P1_random_study_ModMeas=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                       mods=~measurement_mod-1,
                                       data=datum_meta,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModMeas)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=Vector_or_pathogen

results_P1_random_study_ModVectPath=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                           mods=~vector_or_pathogen_mod-1,
                                           data=datum_meta,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModVectPath)

#load moderator z results
datum_fishzMod=read.csv(file.choose())
head(datum_fishzMod)
fishZmod=datum_fishzMod

#transform z to r
fishZmod_r=transf.ztor(fishZmod)
write.csv(fishZmod_r,"E:/AU_Classes/FISH 7350/Class_Project/Analyses/fishZmod_r_ed.csv")

#Load pathogen alone data
datum_meta_pathogen=read.csv(file.choose())
head(datum_meta_pathogen)

#pathogen - test fishers z calculated from p values (> 0.05 = 1) add random study
results_P1_random_study_path=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,data=datum_meta_pathogen,method="REML",random=~1|Paper_Name_Short,slab=Paper_Name_Short)
summary(results_P1_random_study_path)
forest(results_P1_random_study_path)
funnel(results_P1_random_study_path,xlab="Effect size (Fisher's Z)")
summary(results_P1_random_study_p)
#Fail-safe results
fsn(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,data=datum_meta_pathogen)
#rank test
ranktest(results_P1_random_study_path,exact=FALSE)

#univariable moderator models for pathogen

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=scale
results_P1_random_study_ModScale_path=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                             mods=~Spatial_scale_m_mod,
                                             data=datum_meta_pathogen,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModScale_path,digits=6)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=specialization
results_P1_random_study_ModSpec_path=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                            mods=~Specialist_generalist_both_mod-1,
                                            data=datum_meta_pathogen,method="REML",random=~1|Paper_Name_Short)
summary(results_P1_random_study_ModSpec_path)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=landscape
results_P1_random_study_ModLand_path=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                            mods=~Landscape_category_mod-1,
                                            data=datum_meta_pathogen,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModLand_path)
#subset by transmission type
#only 4 data points for semipersistent so this was not analyzed

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=scale, subset persistent
results_P1_random_study_ModScale_path_pers=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                                  mods=~Spatial_scale_m_mod,
                                                  data=datum_meta_pathogen,subset=c(Mode_of_pathogen_transmission=="persistent"),method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModScale_path_pers,digits=6)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=scale, subset persistent
results_P1_random_study_ModScale_path_non=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                                 mods=~Spatial_scale_m_mod,
                                                 data=datum_meta_pathogen,subset=c(Mode_of_pathogen_transmission=="non-persistent"),method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModScale_path_non,digits=6)

#load moderator z results for pathogen alone
datum_fishzMod_path=read.csv(file.choose())
head(datum_fishzMod_path)
fishZmodPath=datum_fishzMod_path

#transform z to r
fishZmodPath_r=transf.ztor(fishZmodPath)
write.csv(fishZmodPath_r,"E:/AU_Classes/FISH 7350/Class_Project/Analyses/fishZmod_r_path.csv")

#Load vector alone data
datum_meta_vector=read.csv(file.choose())
head(datum_meta_vector)

#vector - test fishers z calculated from p values (> 0.05 = 1) add random study
results_P1_random_study_vect=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,data=datum_meta_vector,method="REML",random=~1|Paper_Name_Short,slab=Paper_Name_Short)
summary(results_P1_random_study_vect)
forest(results_P1_random_study_vect)
funnel(results_P1_random_study_vect,xlab="Effect size (Fisher's Z)")
summary(results_P1_random_study_p_vect)
#Fail-safe results
fsn(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,data=datum_meta_vector)
#rank test
ranktest(results_P1_random_study_vect,exact=FALSE)


#univariable moderator models for vector

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=scale
results_P1_random_study_ModScale_vect=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                             mods=~Spatial_scale_m_mod,
                                             data=datum_meta_vector,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModScale_vect,digits=6)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=specialization
results_P1_random_study_ModSpec_vect=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                            mods=~Specialist_generalist_both_mod-1,
                                            data=datum_meta_vector,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModSpec_vect)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=landscape
results_P1_random_study_ModLand_vect=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                            mods=~Landscape_category_mod-1,
                                            data=datum_meta_vector,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModLand_vect)

#test fishers z calculated from p values (> 0.05 = 1) add random study and moderator=measurement_mod

results_P1_random_study_ModMeas_vect=rma.mv(yi_from_ti_ri_pi_1,vi_from_ti_ri_pi_1,
                                            mods=~measurement_mod-1,
                                            data=datum_meta_vector,method="REML",random=~1|Paper_Name_Short)

summary(results_P1_random_study_ModMeas_vect)

#load moderator z results for vector alone
datum_fishzMod_vect=read.csv(file.choose())
head(datum_fishzMod_vect)
fishZmodvect=datum_fishzMod_vect

#transform z to r
fishZmodvect_r=transf.ztor(fishZmodvect)
write.csv(fishZmodvect_r,"E:/AU_Classes/FISH 7350/Class_Project/Analyses/fishZmod_r_vect.csv")
