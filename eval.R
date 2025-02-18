## script for evaluating, output from the DL model
library(tibble)
path <- "//Users//matga374//Desktop//log_fold1-2_x200.csv"


data <- read.csv(path,
                 header = TRUE, sep = ",")
data2 <- as_tibble(data)
ggplot(data2, aes(x=epoch)) + 
  geom_line(aes(y = train_acc, colour="training"), size=1) + 
  geom_line(aes(y = val_acc, colour="validation"), size=1) +
  ylab("accuracy") +
  scale_color_manual(name = "", values = c("training" = "steelblue", "validation" = "darkred"))

ggplot(data2, aes(x=epoch)) + 
  geom_line(aes(y = train_loss, colour="training"), size=1) + 
  geom_line(aes(y = val_loss, colour="validation"), size=1) +
  ylab("loss") +
  scale_color_manual(name = "", values = c("training" = "steelblue", "validation" = "darkred"))