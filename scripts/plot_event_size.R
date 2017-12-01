x <- read.table("events_number_size")
y <- data.frame(branch=x$V1,
birth=x$V2, extension=x$V5, death=-x$V8, contraction=-x$V11)
require(reshape2)
library(ggplot2)

dat <- melt(y)
dat$branch <- factor(dat$branch, levels=x$V1)


pdf("events_size.pdf", width = 4, height = 3.2, pointsize = 8)
dat1 <- subset(dat, value >= 0)
dat2 <- subset(dat, value < 0)
ggplot() +
geom_bar(data = dat1, aes(x=branch, y=value, fill=variable),stat = "identity") +
geom_bar(data = dat2, aes(x=branch, y=value, fill=variable),stat = "identity") +
scale_fill_brewer(type = "seq", palette = 1)+
xlab("Branch") +
ylab("Region size (Mbp)")+
theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
theme(axis.line = element_line(colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_rect(colour = "black"),
panel.background = element_blank())
dev.off()


