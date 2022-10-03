library(ggplot2)
library(gridExtra)

df = read.csv("traces-gen.csv")

# static environments
subset = df[df$expt %in% c(0, 1, 2),]
g1 = ggplot(subset) + geom_path(aes(x=t, y=xmstar)) + geom_line(aes(x=t, y=xm, color=factor(alpha)))  + facet_wrap(~expt)
g2 = ggplot(subset) + geom_line(aes(x=t, y=cost, color=factor(alpha))) + facet_wrap(~expt)
g3 = ggplot(subset) + geom_path(aes(x=xmstar, y=xm, color=factor(alpha))) + facet_wrap(~expt)

png("fig-1a.png", width=600, height=800)
grid.arrange(g1, g2, g3, nrow=3)
dev.off()
		 
# static environments
subset = df[df$expt %in% c(3, 4, 5, 6),]
g1 = ggplot(subset) + geom_path(aes(x=t, y=xmstar)) + geom_line(aes(x=t, y=xm, color=factor(alpha)))  + facet_wrap(~expt)
g2 = ggplot(subset) + geom_line(aes(x=t, y=cost, color=factor(alpha))) + facet_wrap(~expt)
g3 = ggplot(subset) + geom_path(aes(x=xmstar, y=xm, color=factor(alpha))) + facet_wrap(~expt)

png("fig-2a.png", width=600, height=800)
grid.arrange(g1, g2, g3, nrow=3)
dev.off()

# noisy environments
subset = df[df$expt %in% c(7, 8, 9, 10),]
g1 = ggplot(subset) + geom_path(aes(x=t, y=xmstar)) + geom_line(aes(x=t, y=xm, color=factor(alpha)))  + facet_wrap(~expt)
g2 = ggplot(subset) + geom_line(aes(x=t, y=cost, color=factor(alpha))) + facet_wrap(~expt)
g3 = ggplot(subset) + geom_path(aes(x=xmstar, y=xm, color=factor(alpha))) + facet_wrap(~expt)

png("fig-3a.png", width=600, height=800)
grid.arrange(g1, g2, g3, nrow=3)
dev.off()
