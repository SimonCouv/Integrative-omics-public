get_saphir_pheno_associations <- function(eigengenes, y){
  
  # phenotype association models
  me_wilcox_cvd <- function(df) wilcox.test(formula=value~cvd, data=df, conf.int=TRUE)
  me_wilcox_sex <- function(df) wilcox.test(formula=value~sex, data=df, conf.int=TRUE)
  me_spear_age <- function(df) cor.test(formula = ~value+age, data=df, method="spearman")
  me_coxph_cvdw <- function(df) coxph(formula=Surv(cvdw_days,cvd)~value, data=df)
  me_lrm_cvd_sa_adj <- function(df){
    fit <- lrm(formula=cvd~value+sex+age, data = df, x=TRUE, y=TRUE)
    list(p.value = anova(fit)['value', 'P'],
         coef = fit$coefficients['value'],
         gof = resid(fit, "gof")['P'])
  }
  me_lrm_cvd <- function(df){
    fit <- lrm(formula=cvd~value, data = df, x=TRUE, y=TRUE)
    list(p.value = anova(fit)['value', 'P'],
         coef = fit$coefficients['value'],
         gof = resid(fit, "gof")['P'])
  }
  # # note on interpretation: wilcox.test estimate is positive when group 0 is higher than group 1
  # # i.e. the estimate ~ group0 - group1
  # dfa <- data.frame(x=c(rnorm(100,4,1), rnorm(150,2,1)), group=rep(c(1,0), c(100,150)))
  # dfa <- data.frame(x=c(rnorm(100,4,1), rnorm(150,2,1)), group=rep(c('1','0'), c(100,150)))
  # dfa <- data.frame(x=c(rnorm(100,4,1), rnorm(150,2,1)), group=factor(rep(c('1','0'), c(100,150))))
  # dfa <- data.frame(x=c(rnorm(100,4,1), rnorm(150,2,1)), group=factor(rep(c(1,0), c(100,150))))
  # dfa <- data.frame(x= c(rnorm(150,2,1), rnorm(100,4,1)), group=factor(rep(c(0,1), c(150,100))))
  # dfa <- data.frame(x= c(rnorm(150,2,1), rnorm(100,4,1)), group=rep(c(0,1), c(150,100)))
  # wilcox.test(formula=x~group, data = dfa, conf.int=TRUE)$estimate
  # xs <- rnorm(500,0,10)
  # dfb <- data.frame(x=xs, y=xs+rnorm(500,10,1))
  # cor.test(formula=~x+y, data=dfb, method="spearman")
  
  # browser()
  # statistical tests
  suppressWarnings(
    assoc_res <- eigengenes%>%
      rownames_to_column("sample") %>%
      setNames(sub(names(.),pattern = "#", replacement = "_")) %>%
      bind_cols(y,.) %>%
      # bind_cols(cvd[,-1]) %>%
      gather(ME, value, -one_of("id", "cvd", "cvdw_days", 
                                "age", "sex", "sample")) %>%
      # spread(sample, value) %>%
      dplyr:: group_by(ME) %>%
      nest() %>%
      dplyr::mutate(
        wilcox_cvd_res=map(data, me_wilcox_cvd),
        wilcox_sex_res=map(data, me_wilcox_sex),
        spear_age_res=map(data, me_spear_age),
        coxph_cvdw_model=map(data, me_coxph_cvdw),
        lrm_cvd_sa_adj_res = map(data, me_lrm_cvd_sa_adj),
        lrm_cvd_res = map(data, me_lrm_cvd))
  )
  
  
  
  # browser()
  pval_m <- assoc_res %>%
    dplyr::mutate(
      spear_age_p=map_dbl(spear_age_res, "p.value"),
      wilcox_sex_p=map_dbl(wilcox_sex_res, "p.value"),
      lrm_cvd_sa_adj_p = map_dbl(lrm_cvd_sa_adj_res, "p.value"),
      lrm_cvd_p = map_dbl(lrm_cvd_res, "p.value"),
      wilcox_cvd_p=map_dbl(wilcox_cvd_res, "p.value")) %>%
    dplyr::select(ME, spear_age_p:wilcox_cvd_p) %>%
    column_to_rownames("ME") %>%
    as.matrix()
  
  # BH-adjust p-values per phenotype
  pval_m <- apply(pval_m, 2, p.adjust, method="BH")
  
  # #overconservative: BH-adjust over all phenotype associations simultaneously
  # pval_m <- matrix(p.adjust(pval_m, method="BH"),
  #                  nrow = nrow(pval_m),
  #                  ncol= ncol(pval_m),
  #                  dimnames = list(rownames(pval_m),
  #                                  colnames(pval_m)))
  
  #for wilcox estimate: invert direction so positive when group 1 is higher
  est_m <- assoc_res %>%
    dplyr::mutate(
      spear_age_est=map_dbl(spear_age_res, "estimate"),
      wilcox_sex_est=map_dbl(wilcox_sex_res, "estimate")*-1,
      lrm_cvd_sa_adj_coef = map_dbl(lrm_cvd_sa_adj_res, "coef"),
      lrm_cvd_coef = map_dbl(lrm_cvd_res, "coef"),
      wilcox_cvd_est=map_dbl(wilcox_cvd_res, "estimate")*-1) %>%
    dplyr::select(ME, spear_age_est:wilcox_cvd_est) %>%
    column_to_rownames("ME") %>%
    as.matrix()  
  
  # browser()
  p_star_m <- get_p_stars(pval_m)
  
  coxph_m <- assoc_res %>%
    dplyr::mutate(cvdw_tidy=map(coxph_cvdw_model,tidy, exponentiate=FALSE)) %>%
    unnest(cvdw_tidy) %>%
    dplyr::select(ME, estimate:conf.high)%>%
    column_to_rownames("ME") %>%
    as.matrix()
  

  return(list(pval_m=pval_m, est_m=est_m, p_star_m=p_star_m, coxph_m=coxph_m))
}