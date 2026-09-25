# =============================================================================
# OBC 2.0 — Inferential Alignment Major Revision v4
# Reads the derived speaker-within-trial table from the original pipeline.
# Main revisions:
#   1. estimands defined outside model coding;
#   2. Gamma-Poisson empirical-Bayes rate calibration;
#   3. analytic variance benchmark for log-rate contrasts;
#   4. fixed-total-token breadth/depth allocation;
#   5. source-matched Proceeding breadth with fixed source-weighted target;
#   6. Monte Carlo standard errors (ADEMP-style reporting);
#   7. raw/add-0.5/EB sensitivity.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

USER_HOME <- Sys.getenv("USERPROFILE")
if (!nzchar(USER_HOME)) USER_HOME <- path.expand("~")
BASE_DIR <- file.path(USER_HOME, "Downloads", "OBC_IJCL_InferentialAlignment")
INFILE <- file.path(BASE_DIR, "03_speaker_trial_analysis_units.csv")
OUT_DIR <- file.path(USER_HOME, "Downloads", "OBC_IJCL_MAJOR_REVISION_v4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INFILE)) {
  INFILE <- file.choose()
}

set.seed(20260925)
dt <- fread(INFILE)
dt <- dt[role %chin% c("witness", "judge") & period2 %chin% c("early", "late") & tokens > 0]
dt[, role := factor(role, levels=c("witness","judge"))]
dt[, period2 := factor(period2, levels=c("early","late"))]
FEATURES <- c("modal_aux_n", "first_person_n")
CELLS <- data.table(role=c("witness","witness","judge","judge"), period2=c("early","late","early","late"))

# ---------- helpers ----------
fit_gamma_poisson <- function(y, exposure) {
  raw <- y / exposure
  m <- max(mean(raw), 1e-8)
  v <- max(var(raw), 1e-8)
  a0 <- max(m^2 / v, 0.1)
  b0 <- max(m / v, 1)
  nll <- function(z) {
    a <- exp(z[1]); b <- exp(z[2])
    -sum(lgamma(y+a)-lgamma(a)-lgamma(y+1) +
           a*(log(b)-log(b+exposure)) +
           y*(log(exposure)-log(b+exposure)))
  }
  op <- optim(log(c(a0,b0)), nll, method="L-BFGS-B",
              lower=c(-10,-10), upper=c(10,20))
  c(alpha=exp(op$par[1]), beta=exp(op$par[2]), convergence=op$convergence)
}

cell_counts <- function(N, delta_nominal) {
  nd <- as.integer(round((1+delta_nominal)*N/4))
  no <- as.integer(N/2 - nd)
  data.table(role=CELLS$role, period2=CELLS$period2,
             n=c(nd,no,no,nd),
             delta_realized=(nd-no)/(N/2))
}

contrast_vectors <- list(
  role_early      = c(-1, 0, 1, 0),
  period_witness  = c(-1, 1, 0, 0),
  interaction     = c( 1,-1,-1, 1),
  role_avg        = c(-.5,-.5,.5,.5),
  period_avg      = c(-.5,.5,-.5,.5)
)

mc_summary <- function(est, truth) {
  est <- est[is.finite(est)]
  err <- est-truth
  rmse <- sqrt(mean(err^2))
  data.table(valid=length(est), bias=mean(err), emp_sd=sd(est), rmse=rmse,
             bias_mcse=sd(est)/sqrt(length(est)),
             rmse_mcse=sd(err^2)/(sqrt(length(est))*2*rmse))
}

# ---------- EB calibration ----------
ebpars <- list()
for (feat in FEATURES) {
  dt[, paste0(feat,"_raw_rate") := get(feat)/tokens]
  dt[, paste0(feat,"_add_rate") := (get(feat)+.5)/(tokens+1)]
  dt[, paste0(feat,"_eb_rate") := NA_real_]
  for (r in c("witness","judge")) for (p in c("early","late")) {
    ii <- which(dt$role==r & dt$period2==p)
    par <- fit_gamma_poisson(dt[[feat]][ii], dt$tokens[ii])
    dt[ii, (paste0(feat,"_eb_rate")) := (get(feat)+par["alpha"])/(tokens+par["beta"])]
    ebpars[[length(ebpars)+1L]] <- data.table(feature=feat,role=r,period=p,
      alpha=par["alpha"], beta=par["beta"], prior_mean=par["alpha"]/par["beta"])
  }
}
fwrite(rbindlist(ebpars), file.path(OUT_DIR,"01_EB_parameters.csv"))

# rate diagnostics / truths
rate_diag <- rbindlist(lapply(FEATURES, function(feat) rbindlist(lapply(c("raw","add","eb"), function(s) {
  col <- paste0(feat,"_",s,"_rate")
  dt[, .(n_units=.N, zero_prop=mean(get(feat)==0), mean_tokens=mean(tokens), median_tokens=median(tokens),
         mean_rate=mean(get(col)), median_rate=median(get(col)), q95_rate=quantile(get(col),.95), var_rate=var(get(col))),
     by=.(role,period2)][, `:=`(feature=feat,scheme=s)]
}))))
fwrite(rate_diag, file.path(OUT_DIR,"02_rate_diagnostics.csv"))

get_cell_stats <- function(feature, scheme="eb") {
  col <- paste0(feature,"_",scheme,"_rate")
  z <- dt[, .(mu=mean(get(col)), sigma2=var(get(col))), by=.(role,period2)]
  z[, ord := match(paste(role,period2), paste(CELLS$role,CELLS$period2))]
  setorder(z,ord); z[]
}

get_truths <- function(feature, scheme="eb") {
  z <- get_cell_stats(feature,scheme); lm <- log(z$mu)
  rbindlist(lapply(names(contrast_vectors), function(nm) data.table(feature=feature,scheme=scheme,
    estimand=nm,truth=sum(contrast_vectors[[nm]]*lm))))
}
truths <- rbindlist(lapply(FEATURES, function(f) rbindlist(lapply(c("raw","add","eb"), function(s) get_truths(f,s)))))
fwrite(truths, file.path(OUT_DIR,"03_truth_sensitivity.csv"))

# ---------- analytic variance ----------
analytic_var <- function(feature, B, N, delta_nominal, estimand, scheme="eb") {
  cs <- get_cell_stats(feature,scheme)
  cc <- cell_counts(N,delta_nominal)
  t <- B/N
  vg <- (cs$sigma2 + cs$mu/t)/(cc$n*cs$mu^2)
  cvec <- contrast_vectors[[estimand]]
  sum((cvec^2)*vg)
}

# ---------- primary Monte Carlo ----------
sim_one <- function(feature,B,N,delta_nominal,scheme="eb") {
  cc <- cell_counts(N,delta_nominal); t <- B/N
  out <- vector("list",4)
  for (j in 1:4) {
    pool <- dt[role==CELLS$role[j] & period2==CELLS$period2[j]][[paste0(feature,"_",scheme,"_rate")]]
    rs <- sample(pool, cc$n[j], replace=TRUE)
    out[[j]] <- sum(rpois(cc$n[j], rs*t))/(cc$n[j]*t)
  }
  if (any(unlist(out)<=0)) return(rep(NA_real_,length(contrast_vectors)))
  lm <- log(unlist(out))
  vapply(contrast_vectors, function(cv) sum(cv*lm), numeric(1))
}

run_primary <- function(reps=5000L, scheme="eb") {
  res <- list(); k <- 0L
  for (feature in FEATURES) for (B in c(100000L,250000L,500000L)) for (d in c(0,.25,.50,.75,.90)) {
    cc <- cell_counts(200L,d); M <- matrix(NA_real_,nrow=reps,ncol=length(contrast_vectors)); colnames(M)<-names(contrast_vectors)
    for (rr in seq_len(reps)) M[rr,] <- sim_one(feature,B,200L,d,scheme)
    tr <- get_truths(feature,scheme)
    for (e in names(contrast_vectors)) {
      k<-k+1L; s<-mc_summary(M[,e],tr[estimand==e]$truth)
      av<-analytic_var(feature,B,200L,d,e,scheme)
      res[[k]]<-cbind(data.table(feature=feature,token_budget=B,delta_nominal=d,delta_realized=cc$delta_realized[1],estimand=e,
                                requested_reps=reps,analytic_var=av,analytic_sd=sqrt(av)),s)
    }
  }
  ans<-rbindlist(res)
  ans[, analytic_ratio:=analytic_sd/analytic_sd[delta_nominal==0],by=.(feature,token_budget,estimand)]
  ans[, rmse_ratio:=rmse/rmse[delta_nominal==0],by=.(feature,token_budget,estimand)]
  ans
}
primary_res <- run_primary(5000L,"eb")
fwrite(primary_res,file.path(OUT_DIR,"04_primary_EB_5000.csv"))

# ---------- smoothing sensitivity at 250k ----------
smooth_res <- list(); kk<-0L
for (feature in FEATURES) for (scheme in c("raw","add","eb")) for (d in c(0,.90)) {
  M<-matrix(NA_real_,nrow=2000,ncol=length(contrast_vectors)); colnames(M)<-names(contrast_vectors)
  for (rr in 1:2000) M[rr,]<-sim_one(feature,250000L,200L,d,scheme)
  tr<-get_truths(feature,scheme)
  kk<-kk+1L; s<-mc_summary(M[,"interaction"],tr[estimand=="interaction"]$truth)
  smooth_res[[kk]]<-cbind(data.table(feature=feature,scheme=scheme,delta=d,truth=tr[estimand=="interaction"]$truth,requested_reps=2000),s)
}
fwrite(rbindlist(smooth_res),file.path(OUT_DIR,"05_smoothing_sensitivity.csv"))

# ---------- fixed-budget allocation ----------
alloc_res<-list(); kk<-0L
for (feature in FEATURES) for (N in c(100L,200L,500L,1000L)) for (d in c(0,.50,.90)) {
  M<-rep(NA_real_,5000); tr<-get_truths(feature,"eb")[estimand=="interaction"]$truth
  for (rr in 1:5000) M[rr]<-sim_one(feature,250000L,N,d,"eb")["interaction"]
  kk<-kk+1L; s<-mc_summary(M,tr)
  alloc_res[[kk]]<-cbind(data.table(feature=feature,delta=d,N=N,tokens_per_unit=250000/N,requested_reps=5000,truth=tr),s)
}
fwrite(rbindlist(alloc_res),file.path(OUT_DIR,"06_fixed_budget_allocation.csv"))

# ---------- source-matched breadth ----------
# Common Proceedings per period, shared between roles.
common <- list()
for (p in c("early","late")) {
  w <- unique(dt[role=="witness" & period2==p]$proceeding)
  j <- unique(dt[role=="judge" & period2==p]$proceeding)
  common[[p]] <- intersect(w,j)
}
source_means <- list()
for (feature in FEATURES) for (r in c("witness","judge")) for (p in c("early","late")) {
  sm <- dt[role==r & period2==p & proceeding %in% common[[p]], .(source_rate=mean(get(paste0(feature,"_eb_rate")))), by=proceeding]
  source_means[[paste(feature,r,p,sep="__")]] <- sm
}
source_truth <- function(feature) {
  mu <- sapply(1:4,function(j) mean(source_means[[paste(feature,CELLS$role[j],CELLS$period2[j],sep="__")]]$source_rate))
  sum(contrast_vectors$interaction*log(mu))
}
source_sim <- function(feature,K,total_cell_exposure=62500) {
  est <- numeric(4)
  for (p in c("early","late")) {
    chosen <- sample(common[[p]],K,replace=FALSE)
    exp_per <- total_cell_exposure/K
    for (r in c("witness","judge")) {
      sm <- source_means[[paste(feature,r,p,sep="__")]]
      rates <- sm[match(chosen,proceeding)]$source_rate
      ratehat <- sum(rpois(K,rates*exp_per))/(K*exp_per)
      j <- which(CELLS$role==r & CELLS$period2==p)
      est[j] <- ratehat
    }
  }
  if(any(est<=0)) return(NA_real_)
  sum(contrast_vectors$interaction*log(est))
}
src_res<-list();kk<-0L
for(feature in FEATURES) for(K in c(3L,5L,10L,20L,40L)) {
  tr<-source_truth(feature); M<-rep(NA_real_,10000)
  for(rr in 1:10000) M[rr]<-source_sim(feature,K)
  kk<-kk+1L;s<-mc_summary(M,tr)
  src_res[[kk]]<-cbind(data.table(feature=feature,proceedings_per_period=K,truth=tr,requested_reps=10000),s)
}
fwrite(rbindlist(src_res),file.path(OUT_DIR,"07_source_matched_breadth.csv"))

# ---------- figures ----------
# Figure 1: Monte Carlo vs analytic at 250k for interaction
pdat<-primary_res[token_budget==250000 & estimand=="interaction"]
# equal variance benchmark using realized delta
pdat[, equalvar_ratio:=1/sqrt(1-delta_realized^2)]
p1<-ggplot(pdat,aes(delta_realized,rmse_ratio))+geom_line()+geom_point()+
  geom_line(aes(y=analytic_ratio),linetype=2)+geom_line(aes(y=equalvar_ratio),linetype=3)+
  facet_wrap(~feature)+labs(x="Realized delta",y="RMSE / RMSE(delta=0)",title="Monte Carlo and analytic joint-support penalties")+theme_minimal()
ggsave(file.path(OUT_DIR,"Figure_1_Benchmark.png"),p1,width=8.5,height=4.5,dpi=400)

p2dat<-fread(file.path(OUT_DIR,"06_fixed_budget_allocation.csv"))
p2<-ggplot(p2dat,aes(N,rmse,group=factor(delta),linetype=factor(delta)))+geom_line()+geom_point()+facet_wrap(~feature,scales="free_y")+
  labs(x="Speaker-within-trial units (fixed 250k tokens)",y="Interaction RMSE",linetype="delta",title="Fixed-budget breadth-depth allocation")+theme_minimal()
ggsave(file.path(OUT_DIR,"Figure_2_Allocation.png"),p2,width=8.5,height=4.5,dpi=400)

p3dat<-fread(file.path(OUT_DIR,"07_source_matched_breadth.csv"))
p3<-ggplot(p3dat,aes(proceedings_per_period,rmse))+geom_line()+geom_point()+facet_wrap(~feature,scales="free_y")+
  labs(x="Matched Proceedings per period",y="RMSE",title="Source breadth under fixed source-weighted target")+theme_minimal()
ggsave(file.path(OUT_DIR,"Figure_3_SourceBreadth.png"),p3,width=8.5,height=4.5,dpi=400)

capture.output(sessionInfo(),file=file.path(OUT_DIR,"99_sessionInfo.txt"))
cat("\nCompleted. Results: ",normalizePath(OUT_DIR,winslash="/",mustWork=FALSE),"\n",sep="")
