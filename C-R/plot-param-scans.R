library(ggplot2)
library(gridExtra)

# gradient fill blue-red with white = 1
colfn = scale_fill_gradientn(colours = c("black", "blue", "white", "red", "black"), values = c(0, 0.5/2.5, 1/2.5, 2/2.5, 1), limits = c(0,2.5))

######################### static environment, no damage
# read file and factor-ise encoding strategy
costs.df = read.csv("costsscans-0.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative[costs.relative$b == 0 & costs.relative$k == 0,]
cost.plot.1 = ggplot(subset, aes(x=num,y=nuc,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.1

######################### static environment, oDNA damage
costs.df = read.csv("costsscans-1.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative[costs.relative$b == 0 & costs.relative$k == 0,]
cost.plot.1a = ggplot(subset, aes(x=num,y=nuc,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.1a

png("fig-1.png", height = 800, width = 800)
cost.plot.1a
dev.off()

png("fig-s1.png", height = 800, width = 800)
cost.plot.1
dev.off()

######################### dynamic environment, no damage

costs.df = read.csv("costsscans-2.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

### facet diff, no cytoplasmic degradation
subset = costs.relative[costs.relative$mu == 0 & costs.relative$nuc == 0 & costs.relative$num == 0.1,]
cost.plot.1 = ggplot(subset[subset$k == 1,], aes(x=a,y=b,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.2 = ggplot(subset[subset$a == 1,], aes(x=b,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.3 = ggplot(subset[subset$b == 1,], aes(x=a,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
grid.arrange(cost.plot.1, cost.plot.2, cost.plot.3, nrow=3)
### ^ this plot is interesting 

png("fig-s2.png", width=800, height = 600)
cost.plot.2
dev.off()

######################### dynamic environment, oDNA damage (central plot)

costs.df = read.csv("costsscans-3.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

### facet diff, no cytoplasmic degradation
subset = costs.relative[costs.relative$mu == 0 & costs.relative$nuc == 0 & costs.relative$num == 0.1,]
cost.plot.1 = ggplot(subset[subset$k == 1,], aes(x=a,y=b,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.2 = ggplot(subset[subset$a == 1,], aes(x=b,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.3 = ggplot(subset[subset$b == 1,], aes(x=a,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
grid.arrange(cost.plot.1, cost.plot.2, cost.plot.3, nrow=3)
### ^ this plot is interesting 

png("fig-2.png", width=800, height = 600)
cost.plot.2
dev.off()

png("fig-s3.png", width=800, height = 600)
grid.arrange(cost.plot.1, cost.plot.2, cost.plot.3, nrow=3)
dev.off()

# facet diff, nonzero cytoplasmic degradation
subset = costs.relative[costs.relative$mu == 0 & costs.relative$nuc == 0.1 & costs.relative$num == 0.1,]
cost.plot.1 = ggplot(subset[subset$k == 1,], aes(x=a,y=b,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.2 = ggplot(subset[subset$a == 1,], aes(x=b,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.3 = ggplot(subset[subset$b == 1,], aes(x=a,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
grid.arrange(cost.plot.1, cost.plot.2, cost.plot.3, nrow=3)

png("fig-s4.png", width=800, height = 600)
cost.plot.2
dev.off()

######################### dynamic environment, oDNA damage, different cost function

costs.df = read.csv("costsscans-4.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

### facet diff, no cytoplasmic degradation
subset = costs.relative[costs.relative$mu == 0 & costs.relative$nuc == 0 & costs.relative$num == 0.1,]
cost.plot.1 = ggplot(subset[subset$k == 1,], aes(x=a,y=b,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.2 = ggplot(subset[subset$a == 1,], aes(x=b,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.3 = ggplot(subset[subset$b == 1,], aes(x=a,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
grid.arrange(cost.plot.1, cost.plot.2, cost.plot.3, nrow=3)
### ^ this plot is interesting 

png("fig-s5.png", width=800, height = 600)
cost.plot.2
dev.off()

######################### dynamic environment, oDNA damage, different lambdas

costs.df = read.csv("costsscans-5.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative[costs.relative$mu == 0 & costs.relative$nuc == 0 & costs.relative$num == 0.1,]
cost.plot.1 = ggplot(subset[subset$k == 1,], aes(x=a,y=b,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.2 = ggplot(subset[subset$a == 1,], aes(x=b,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.3 = ggplot(subset[subset$b == 1,], aes(x=a,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
grid.arrange(cost.plot.1, cost.plot.2, cost.plot.3, nrow=3)

costs.df = read.csv("costsscans-6.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative[costs.relative$mu == 0 & costs.relative$nuc == 0 & costs.relative$num == 0.1,]
cost.plot.1a = ggplot(subset[subset$k == 1,], aes(x=a,y=b,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.2a = ggplot(subset[subset$a == 1,], aes(x=b,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.3a = ggplot(subset[subset$b == 1,], aes(x=a,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
grid.arrange(cost.plot.1, cost.plot.2, cost.plot.3, nrow=3)

png("fig-s6.png", width=800, height = 600)
grid.arrange(cost.plot.2, cost.plot.2a, nrow=2)
dev.off()

######################### noisy environments, oDNA damage
## white noise
# read file and factor-ise encoding strategy
costs.df = read.csv("costsscans-7.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative
cost.plot.1 = ggplot(subset, aes(x=num,y=nuc,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.1

## red noise
# read file and factor-ise encoding strategy
costs.df = read.csv("costsscans-8.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative
cost.plot.1a = ggplot(subset, aes(x=num,y=nuc,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.1a

png("fig-4.png", width=800, height=800)
grid.arrange(cost.plot.1, cost.plot.1a, nrow=2)
dev.off()

######################### dynamic environment, oDNA damage, different different cost function (quadratic)

costs.df = read.csv("costsscans-9.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

### facet diff, no cytoplasmic degradation
subset = costs.relative[costs.relative$mu == 0 & costs.relative$nuc == 0 & costs.relative$num == 0.1,]
cost.plot.1 = ggplot(subset[subset$k == 1,], aes(x=a,y=b,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.2 = ggplot(subset[subset$a == 1,], aes(x=b,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
cost.plot.3 = ggplot(subset[subset$b == 1,], aes(x=a,y=k,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)
grid.arrange(cost.plot.1, cost.plot.2, cost.plot.3, nrow=3)
### ^ this plot is interesting 

png("fig-s7.png", width=800, height = 600)
cost.plot.2
dev.off()

######################### noisy environments, different red noise freqs, oDNA damage
## red noise 1
# read file and factor-ise encoding strategy
costs.df = read.csv("costsscans-8.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative
cost.plot.1a = ggplot(subset, aes(x=num,y=nuc,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)

## red noise 2
# read file and factor-ise encoding strategy
costs.df = read.csv("costsscans-10.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative
cost.plot.2a = ggplot(subset, aes(x=num,y=nuc,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)

## red noise 3
# read file and factor-ise encoding strategy
costs.df = read.csv("costsscans-11.csv")
costs.df$encoding = "MT"
costs.df$encoding[costs.df$alpha==1] = "NU"

# new data frame for relative costs
costs.0 = costs.df[costs.df$encoding=="MT",]
costs.1 = costs.df[costs.df$encoding=="NU",]
costs.relative = costs.0
costs.relative$relative = costs.0$cost/costs.1$cost

subset = costs.relative
cost.plot.3a = ggplot(subset, aes(x=num,y=nuc,fill=relative)) + geom_tile() + colfn + facet_wrap(~diff)

png("fig-sx.png", width=800, height=1200)
grid.arrange(cost.plot.1a, cost.plot.2a, cost.plot.3a, nrow=3)
dev.off()
