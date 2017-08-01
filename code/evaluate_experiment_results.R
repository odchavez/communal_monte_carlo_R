#evaluate the quality of the exeriment results
source("code/communal_monte_carlo_functions.R")



test_dat = read.csv("data/K=20/d=2_n=10000_file_num_96.csv")

#process emb. par. particles
#EXP_1_ps = load("particles/K=20d=2_n=40000_pn=1000_sn=4_gs=1_experiment_number=1.RData")
#emb_par_log_lik = experiment_log_lik(test_dat, particles, 1)
# 
#EXP_2_ps = load("particles/K=20d=2_n=40000_pn=1000_sn=4_gs=2_experiment_number=2.RData")
#EXP_2_log_lik = experiment_log_lik(test_dat, particles, 2) 

#1st round of experiments with 10 particles
load("particles/K=20d=2_n=4000_pn=10_sn=4_gs=1_experiment_number=3.RData")
emb_par_log_lik = experiment_log_lik(test_dat, particles, 3)

load("particles/K=20d=2_n=4000_pn=10_sn=4_gs=1_experiment_number=4.RData")
EXP_4_log_lik = experiment_log_lik(test_dat, particles, 4)
load("particles/K=20d=2_n=4000_pn=10_sn=4_gs=2_experiment_number=5.RData")
EXP_5_log_lik = experiment_log_lik(test_dat, particles, 5)
load("particles/K=20d=2_n=3996_pn=10_sn=4_gs=3_experiment_number=6.RData")
EXP_6_log_lik = experiment_log_lik(test_dat, particles, 6)
load("particles/K=20d=2_n=4000_pn=10_sn=4_gs=4_experiment_number=7.RData")
EXP_7_log_lik = experiment_log_lik(test_dat, particles, 7)
load("particles/K=20d=2_n=4000_pn=10_sn=4_gs=5_experiment_number=8.RData")
EXP_8_log_lik = experiment_log_lik(test_dat, particles, 8)
load("particles/K=20d=2_n=4008_pn=10_sn=4_gs=6_experiment_number=9.RData")
EXP_9_log_lik = experiment_log_lik(test_dat, particles, 9)
load("particles/K=20d=2_n=4004_pn=10_sn=4_gs=7_experiment_number=10.RData")
EXP_10_log_lik = experiment_log_lik(test_dat, particles, 10)
load("particles/K=20d=2_n=4000_pn=10_sn=4_gs=8_experiment_number=11.RData")
EXP_11_log_lik = experiment_log_lik(test_dat, particles, 11)
load("particles/K=20d=2_n=3996_pn=10_sn=4_gs=9_experiment_number=12.RData")
EXP_12_log_lik = experiment_log_lik(test_dat, particles, 12)
load("particles/K=20d=2_n=4000_pn=10_sn=4_gs=10_experiment_number=13.RData")
EXP_13_log_lik = experiment_log_lik(test_dat, particles, 13)

exp_mat_10 = cbind(emb_par_log_lik,
                   EXP_4_log_lik,
                   EXP_5_log_lik,
                   EXP_6_log_lik,
                   EXP_7_log_lik,
                   EXP_8_log_lik,
                   EXP_9_log_lik,
                   EXP_10_log_lik,
                   EXP_11_log_lik,
                   EXP_12_log_lik,
                   EXP_13_log_lik)

means_10 = colMeans(exp_mat_10)
sd_10    = apply(exp_mat_10, 2, sd)/sqrt(nrow(exp_mat_10))

#2nd round of experiments with 50 particles
load("particles/K=20d=2_n=4000_pn=50_sn=4_gs=1_experiment_number=14.RData")
emb_par_log_lik_50 = experiment_log_lik(test_dat, particles, 14)

load("particles/K=20d=2_n=4000_pn=50_sn=4_gs=1_experiment_number=15.RData")
EXP_15_log_lik = experiment_log_lik(test_dat, particles, 15)
load("particles/K=20d=2_n=4000_pn=50_sn=4_gs=2_experiment_number=16.RData")
EXP_16_log_lik = experiment_log_lik(test_dat, particles, 16)
load("particles/K=20d=2_n=3996_pn=50_sn=4_gs=3_experiment_number=17.RData")
EXP_17_log_lik = experiment_log_lik(test_dat, particles, 17)
load("particles/K=20d=2_n=4000_pn=50_sn=4_gs=4_experiment_number=18.RData")
EXP_18_log_lik = experiment_log_lik(test_dat, particles, 18)
load("particles/K=20d=2_n=4000_pn=50_sn=4_gs=5_experiment_number=19.RData")
EXP_19_log_lik = experiment_log_lik(test_dat, particles, 19)
load("particles/K=20d=2_n=4008_pn=50_sn=4_gs=6_experiment_number=20.RData")
EXP_20_log_lik = experiment_log_lik(test_dat, particles, 20)
load("particles/K=20d=2_n=4004_pn=50_sn=4_gs=7_experiment_number=21.RData")
EXP_21_log_lik = experiment_log_lik(test_dat, particles, 21)
load("particles/K=20d=2_n=4000_pn=50_sn=4_gs=8_experiment_number=22.RData")
EXP_22_log_lik = experiment_log_lik(test_dat, particles, 22)
load("particles/K=20d=2_n=3996_pn=50_sn=4_gs=9_experiment_number=23.RData")
EXP_23_log_lik = experiment_log_lik(test_dat, particles, 23)
load("particles/K=20d=2_n=4000_pn=50_sn=4_gs=10_experiment_number=24.RData")
EXP_24_log_lik = experiment_log_lik(test_dat, particles, 24)

exp_mat_50 = cbind(emb_par_log_lik_50,
                   EXP_15_log_lik,
                   EXP_16_log_lik,
                   EXP_17_log_lik,
                   EXP_18_log_lik,
                   EXP_19_log_lik,
                   EXP_20_log_lik,
                   EXP_21_log_lik,
                   EXP_22_log_lik,
                   EXP_23_log_lik,
                   EXP_24_log_lik)
means_50 = colMeans(exp_mat_50)
sd_50    = apply(exp_mat_50, 2, sd)/sqrt(nrow(exp_mat_50))

par(mfrow = c(1,2))
plot(c(1,1:10), means_10, main ="average Loglik vs number of GS (10)", xlab = "global steps", ylab = "avg Loglik")
points(c(1,1:10), means_10 + 3*sd_10)
points(c(1,1:10), means_10 - 3*sd_10)
plot(c(1,1:10), means_50, main ="average Loglik vs number of GS (50)", xlab = "global steps", ylab = "avg Loglik")
points(c(1,1:10), means_50 + 3*sd_50)
points(c(1,1:10), means_50 - 3*sd_50)

plot(apply(exp_mat_10, 2, mean), main ="Loglik vs number of GS (10)", las = 2, ylim = c(-10,0))
upper = apply(exp_mat_10, 2, quantile, probs = c(0.05,0.95))
points(upper[,1])
boxplot(exp_mat_50, main ="Loglik vs number of GS (50)", las = 2, ylim = c(-10,0))
