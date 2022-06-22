# Jun 14, 2022, Stage: incorporating CHS and stroke outcomes.
library(survival)
library(dplyr)
library(metafor)
library(ggplot2)


cleaned_folder = "/Users/luli/Dropbox/Research/IPD/processed_CHS_stroke/"

filenames = c("aric_final", "bari_final", "Belfrail_final", "busselton_final", "CHS_data_final",
             "diabetes_final","habc_final","inchianti_final", "leiden_final", 
             "pisa_final", "vivit_final","wickham_final")

clean_data <- function(data){
  # remove rows with NAs in age, gender, and tsh
  data_rm_na = data[!is.na(data$age) & !is.na(data$gender) & !is.na(data$tsh) & !is.na(data$inc_stroke) & !is.na(data$timeinc_stroke), ]
  # impute missing outcome variables
  data_tmp = data_rm_na[,c("age", "gender", "tsh", "inc_stroke", "timeinc_stroke")]
  completedData = data_tmp
  # categorize data
  # create the following age strata:  18-49, 50-64, 65-79, and ≥ 80 years. There should be no one aged <18 in the dataset.
  if(length(which(completedData$age < 18))>0){
    completedData = completedData[-which(completedData$age < 18),]
  }
  completedData$cont_age = completedData$age 
  for (ind in 1:dim(completedData)[1]){
    tmp = completedData$age[ind]
    if(tmp > 17 & tmp < 50){
      completedData$age[ind] = 1
    } else if (tmp >= 50 & tmp < 65){
      completedData$age[ind] = 2
    } else if (tmp >= 65 & tmp < 80){
      completedData$age[ind] = 3
    } else if (tmp >= 80){
      completedData$age[ind] = 4
    }
  }
  # Exclude individuals with TSH <0.45 and >=20. The TSH values use the same units, so no conversion is necessary.
  data_tmp_3 = completedData[-c(which(completedData$tsh < 0.45), which(completedData$tsh >= 20)),]
  # Create the following TSH strata: 0.45-4.49 (Reference), 4.50-6.99, 7.00-9.99, and 10.00-19.99
  data_tmp_3$tsh_cat = ifelse(data_tmp_3$tsh >= 0.45 & data_tmp_3$tsh<=4.49, 1, 0)
  data_tmp_3$tsh_cat = data_tmp_3$tsh_cat + ifelse(data_tmp_3$tsh >= 4.50 & data_tmp_3$tsh<=6.99, 2, 0)
  data_tmp_3$tsh_cat = data_tmp_3$tsh_cat + ifelse(data_tmp_3$tsh >= 7.00 & data_tmp_3$tsh<=9.99, 3, 0)
  data_tmp_3$tsh_cat = data_tmp_3$tsh_cat + ifelse(data_tmp_3$tsh >= 10.00 & data_tmp_3$tsh<=19.99, 4, 0)
  data_tmp_3$tsh_cat = as.factor(data_tmp_3$tsh_cat)
  data_tmp_3$age = as.factor(data_tmp_3$age)
  data_tmp_3$gender = as.factor(data_tmp_3$gender)
  data_filtered_tsh = data_tmp_3
  data_filtered_tsh$tsh_cat = as.factor(data_filtered_tsh$tsh_cat)
  data_filtered_tsh$age = as.factor(data_filtered_tsh$age)
  data_filtered_tsh$gender = as.factor(data_filtered_tsh$gender)
  # Remove observations with fewer than 10 for tsh
  data_filtered_tsh=data_filtered_tsh %>% 
    group_by(age, tsh_cat) %>% 
    filter(n() >= 10)
  # Remove observations with fewer than 10 for age
  data_filtered_tsh=data_filtered_tsh %>% 
    group_by(age, gender) %>% 
    filter(n() >= 10)
  data_filtered_tsh$age <- droplevels(data_filtered_tsh$age)
  data_filtered_tsh$gender <- droplevels(data_filtered_tsh$gender)
  data_filtered_tsh$tsh_cat <- droplevels(data_filtered_tsh$tsh_cat)
  data_full_clean = data_filtered_tsh[,c("tsh_cat","age", "cont_age", "gender",
                                         "inc_stroke","timeinc_stroke")]
  return (data_full_clean)
}

# Stroke outcome analysis
# Loop through 4 age groups
for (j in 1:4){
  firstTime = TRUE
  used_indices = c()
  for (i in 1:length(filenames)){
    data = read.csv(paste(cleaned_folder, filenames[i], ".csv", sep=""))
    data_full_clean = clean_data(data)
    # filter out patients in a given age group
    data_full_clean_age = data_full_clean[which(data_full_clean$age == j), ]
    # if there are patients in this group, proceed to analysis
    if(dim(data_full_clean_age)[1] != 0){
      used_indices = c(used_indices, i)
      values_count <- sapply(lapply(data_full_clean_age[,c("gender","tsh_cat","cont_age")], unique), length)
      fit_stroke <- coxph(Surv( timeinc_stroke, inc_stroke) ~ ., data = data_full_clean_age[, c(c("gender","tsh_cat","cont_age")[values_count > 1], "inc_stroke", "timeinc_stroke")])
      if(firstTime){
        betas.local.all <- data.frame(t(fit_1$coefficients))
        var.betas.local.all <- data.frame(t(diag(vcov(fit_stroke))))
        firstTime = FALSE
      }
      else{
        betas.local.all = dplyr::bind_rows(betas.local.all, data.frame(t(fit_stroke$coefficients)))
        var.betas.local.all = dplyr::bind_rows(var.betas.local.all, data.frame(t(diag(vcov(fit_stroke)))))
      }
      var.betas.local.all[var.betas.local.all ==0] = NA
      SE.betas.local.all <- sqrt(var.betas.local.all) 
    }
  }
  studies = lapply(filenames[used_indices], function(i) {(strsplit(i,"_")[[1]][[1]])})
  
  beta_bar_RE <- rep(0, 5)
  SE_beta_bar_RE <- rep(0, 5)
  beta_bar_FE <- rep(0, 5)
  SE_beta_bar_FE <- rep(0, 5)
  vars <- vector(mode="list", length=9)
  names(vars) <- c("age2", "age3", "age4", "tsh_cat2", "tsh_cat3", "tsh_cat4", "gender1", "age1","cont_age")
  vars[[1]] <- "Age group 2 (50-64)"; vars[[2]] <- "Age group 3 (65-79)"; vars[[3]] <- "Age group 4 (>= 80)"; vars[[4]]<- "Tsh group 2 (4.50-6.99)"
  vars[[5]] <- "Tsh group 3 (7.00-9.99)"; vars[[6]] <- "Tsh group 4 (10.00-19.99)"; vars[[7]] <- "Male"; vars[[8]] <- "Age group 1 (18 - 49)";vars[[9]] <- "Age (continuous)"
  for (i in 1:5) { # 5 means five variables: gender, tsh * 3, and continuous age
    # fixed effects
    fe_res = rma(yi = betas.local.all[,i], vi = var.betas.local.all[,i], method = "FE")
    # random effects
    re_res = rma(yi = betas.local.all[,i], vi = var.betas.local.all[,i], method = "DL")
    beta_bar_FE[i] <- fe_res$beta
    SE_beta_bar_FE[i] <- fe_res$se
    beta_bar_RE[i] <- re_res$beta
    SE_beta_bar_RE[i] <- re_res$se
    # Jessie, feel free to change this output path
    # forestplot33
    png(paste0('/Users/luli/Dropbox/5_with_Lu/3_IPD_meta-analysis/3_arm_based/forest_plots_stroke/age',j,"-",colnames(betas.local.all)[i] , '.png'),
        width     = 3.25,
        height    = 2.25,
        units     = "in",
        res       = 1200,
        pointsize = 4)
    forest.rma(fe_res, 
               annotate=TRUE,
               atransf=exp,
               mlab = bquote("FE Model (" ~ I^2 == .(round(fe_res$I2,2))~"%)"),
               header=c(paste0(vars[paste0("age", j)], "\n", vars[colnames(betas.local.all)[i]][[1]]), "Hazard Ratios (95% CI)"), 
               slab = lapply(filenames[used_indices], function(i) {
                 (strsplit(i,"_")[[1]][1])
               }
               )
    )
    addpoly(re_res, 
            atransf=exp,
            row= -1.5, 
            mlab= as.expression(bquote("RE Model ("~I^2==.(round(re_res$I2,2))*"%)"))) #bquote("RE Model (" ~ I^2 == .(round(fe_res$I2,2))~"%)"))
    dev.off()
  }
}