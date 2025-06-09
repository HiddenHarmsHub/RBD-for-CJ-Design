source("rdb.R")

test_matrix <- read.csv("../data/C_Barnet.csv", header = FALSE)

print(test_matrix)

result <- rbd(test_matrix)



