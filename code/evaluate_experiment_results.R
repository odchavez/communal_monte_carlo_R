#evaluate the quality of the exeriment results
source("code/communal_monte_carlo_functions.R")



test_dat = read.csv("data/K=20/d=2_n=10000_file_num_96.csv")

#process emb. par. particles
EXP_1_ps = load("particles/K=20d=2_n=400_pn=10_sn=4_gs=1_experiment_number=1.RData")
emb_par_log_lik = experiment_log_lik(test_dat, particles, 1)
 
EXP_2_ps = load("particles/K=20d=2_n=400_pn=10_sn=4_gs=2_experiment_number=2.RData")
EXP_2_log_lik = experiment_log_lik(test_dat, particles, 2) 



boxplot(emb_par_log_lik, EXP_2_log_lik)

