# data1 <- read.csv("build/cell_data_20.csv")  
# plot(data1$Position, data1$Physics,col="red")
# data2 <- read.csv("build/cell_data_40.csv")  
# lines(data2$Position, data2$Physics,col="blue")
# data <- read.csv("build/cell_data_100.csv")  
# lines(data$Position, data$Physics,col="black")
# data <- read.csv("build/cell_data_15.csv")  
# lines(data$Position, data$Physics,col="red")
# https://mirrors.tuna.tsinghua.edu.cn/CRAN/


# install.packages("httpgd", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# install.packages("httpgd", dependencies=TRUE, repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

# Create a vector of file numbers

numbers<-15

data1 <- read.csv("build/cell_data_800.csv")  
plot(data1$Position, data1$Physics,col="red")
# Loop through each file number
for (number in seq(4, numbers, by = 5)) {
  # Construct the filename
  filename <- paste("build/cell_data_", number, ".csv", sep="")

  # Read data from the CSV file
  data <- read.csv(filename)

  # Plot the data using lines, with different colors for each file
  if (number == 100) {
    lines(data$Position, data$Physics, col = "black", lty = 1, type = "b", pch = 16, main = "Position vs Physics")
  } else {
  lines(data$Position, data$Physics, col = rainbow(numbers)[number], lty = 1, type = "b", pch = 16)
  }
}

# Add legend
legend("topright", legend = numbers, col = c("black", "red"), lty = 1, pch = 16, title = "File Number")
