library(ggplot2)
branches <- c("Catarrhini", "Homininae",
              "Hominini", "Human", "Chimp", "Gorilla", "Rhesus",
              "Murinae", "Mouse", "Rat", "Dog")
events <- c("birth", "extension", "death", "contraction")
eventnames <- c("B", "E", "D", "C")

# violinplot
x <- data.frame(dist=0, event="", branch="", stringsAsFactors=F)
gainevents <- c(1,2)
lossevents <- c(3,4)

lineage <- c(1:4, 8:9,11)
muttypes <- list(gainevents, lossevents)
mutnames <- c("gain", "loss")
colors <- list(c("#CA2027", "#F2A481"), c("#1B75BB", "#92C4DD"))

for (k in 1:2) {
muttype <- muttypes[[k]]
for (i in lineage){
for (j in muttype) {
filename <- paste("7sp.", branches[i], ".HYPO.", events[j], ".dist2TSS", sep="")
y <- read.table(filename)$V13
xx <- data.frame(dist=log(y+1)/log(10), event=rep(eventnames[j], length(y)),
branch=rep(branches[i], length(y)), stringsAsFactors=F)
if (i==lineage[1] && j==muttype[1]) {x <- xx}
else {x <- rbind(x, xx)}
}
}

x$event <- factor(x$event, levels = eventnames[muttype])
x$branch <- factor(x$branch, levels = branches[lineage])

# Basic violin plot
dp <- ggplot(x, aes(x=branch, y=dist, fill=event)) +
geom_violin(trim=FALSE) + ylim(c(0, 7)) +
labs(title="", x="branch", y = "Distance to TSS")+
scale_fill_manual(values = colors[[k]]) + theme_classic()+
theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
theme(axis.line = element_line(colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_rect(colour = "black"),
panel.background = element_blank())

pdf(paste("Events_dist2TSS_", mutnames[k], ".pdf", sep=""), width=8, height=2)
print(dp)
dev.off()
}


