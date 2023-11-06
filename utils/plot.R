library(ggplot2)
library(reshape2)
library(dplyr)

EFclasses_performances.combine_1_to_1 <- read.delim("~/Projects/GitHub/E2P2/test/performances/EFclasses_performances.combine_1_to_1.tsv", row.names=1)
combine_1_to_1.no_na <- EFclasses_performances.combine_1_to_1[complete.cases(EFclasses_performances.combine_1_to_1),]

EFclasses_performances.original <- read.delim("~/Projects/GitHub/E2P2/test/performances/EFclasses_performances.original.tsv", row.names=1)


original.no_na <- EFclasses_performances.original[!is.na(EFclasses_performances.original$EFClass),]
original.no_na[c("Fscore")][is.na(original.no_na[c("Fscore")])] <- 0
original.no_na <- EFclasses_performances.original[complete.cases(EFclasses_performances.original), ]


EFclasses_to_prots_max_to_prots_len <- read.delim("~/Projects/GitHub/E2P2/test/performances/EFclasses_to_prots_max_to_prots_len.tsv")


p <- ggplot(EFclasses_to_prots_max_to_prots_len, aes(x = Max_Num_of_prots_of_EFs, y = cumsum(Num_of_prots)))

class_length=table(EFclasses_performances.combine_1_to_1$NumOfSeq.RPSD.)
class_length_df=as.data.frame(class_length)
p <- ggplot(class_length_df, aes(x = Var1, y = cumsum(Freq), group = 1))
p + geom_line() + geom_point()

original.no_na$range <- cut(original.no_na$NumOfSeq.RPSD., breaks = seq(0,500,50), include.lowest = T)


ggplot(original.no_na, aes(x = NumOfSeq.RPSD., y = Precision, group = NumOfSeq.RPSD.)) + geom_boxplot() + coord_cartesian(xlim = c(0, 50)) 
ggplot(original.no_na, aes(x = NumOfSeq.RPSD., y = Recall, group = NumOfSeq.RPSD.)) + geom_boxplot() + coord_cartesian(xlim = c(0, 50)) 
ggplot(original.no_na, aes(x = NumOfSeq.RPSD., y = Fscore, group = NumOfSeq.RPSD.)) + geom_boxplot() + coord_cartesian(xlim = c(0, 50)) 


####################################################################
EFclasses_performances.ec_nums <- read.delim("~/Projects/GitHub/E2P2/test/performances/EFclasses_performances.ec_nums.tsv", row.names=1)
EFclasses_performances.rxn_id <- read.delim("~/Projects/GitHub/E2P2/test/performances/EFclasses_performances.rxn_id.tsv", row.names=1)


p <- ggplot(EFclasses_performances.ec_nums, aes(x = NumOfSeq.Test., group = Group, color = Group))
p + geom_density()
p + geom_density() + scale_x_continuous(trans = "log10")

p + geom_line(aes(y = Fscore, color = "Fscore")) + geom_line(aes(y = Precision, color = "Precision")) + geom_line(aes(y = Recall, color = "Recall"))

p + geom_point() + coord_cartesian(xlim = c(0, 10)) + geom_smooth(method = "loess")


library(dplyr)

ggplot(ef.seq, aes(x="", y=Num..of.Seq.in.Test, fill=Group)) + geom_bar(stat="identity", width=1, color="white") + coord_polar("y", start=0)