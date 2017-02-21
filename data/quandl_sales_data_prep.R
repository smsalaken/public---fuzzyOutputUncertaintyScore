## setting working dir
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# data source: https://www.quandl.com/data/WSJ/NG_HH-Natural-Gas-Henry-Hub
WSJ.NG_HH <- read.csv("WSJ-NG_HH.csv")

#View(WSJ.NG_HH)
a <- as.character(WSJ.NG_HH$Date)
a.year <-  sapply(strsplit(a, "-") , "[", 1)
a.month <-  sapply(strsplit(a, "-") , "[", 2)
a.day <-  sapply(strsplit(a, "-") , "[", 3)

backward_horizon <- 3 # how many past sales value should be used as input?
# past sales value will be treated as inputs along with date information
data<- NULL
for (i in 5:(length(a)-backward_horizon)){
  # last column contains output
  data <- rbind(data, as.numeric(matrix(cbind(a.year[i], a.month[i], a.day[i], t(WSJ.NG_HH$Value[(i-backward_horizon):i])), nrow = 1, ncol = backward_horizon+3+1)))

}
#print(data)
colnames(data) <- c("year","month", "day", "past_sale_1", "past_sale_2", "past_sale_3", "current_sales")
write.csv(data, file = "quandl_sales_price.csv",row.names=FALSE)