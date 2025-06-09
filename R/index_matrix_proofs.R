test_matrix <- read.csv("data/test_grid.csv", header = FALSE)

i <- 1

print(test_matrix[, 1 : i - 1])
print(test_matrix[1 : i - 1, ])
print(test_matrix[i, ])
print(test_matrix[, i])

print(test_matrix[i, 1 : i])
print(test_matrix[, 1 : i])
print(test_matrix[1 : i, i])
print(test_matrix[i, 1 : i])

print(test_matrix[1 : i, ])
print(test_matrix[1 : i, 1 : i])
