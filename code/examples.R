
setwd("/Users/o/Google Drive/school/Williamson Research/Implementations/python code/test_runs")

#MVN mixture Gibbs Sampler
data_size = "100000"
K         = "30"
shard_num = 1
data_file = paste('data/MVN_train_data_K=',toString(K),'_data_size=',toString(data_size),'.csv', sep = "")
dat = read.csv(data_file, header = FALSE,  nrows = 1, skip = 0)
dat


stop = FALSE
f = file(data_file, "r")
while(!stop) {
  next_line = readLines(f, n = 1)
  data_line = as.numeric(strsplit(next_line, ",")[[1]])
  ## Insert some if statement logic here
  if(length(next_line) == 0) {
    stop = TRUE
    close(f)
  }
}