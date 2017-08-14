setwd("c:/00_Lab_work/06_LFS_project/Site-specific_mutants/June1_data_analysis")

# load packages
library(ggplot2)
library(dplyr)


# read the data and check its structure
all_data <- read.csv("all_data_june1.txt")
str(all_data)

nhej_data <- subset(all_data, Category %in% c("del - no HDR", "ins - no HDR", "wild-type"))
nhej_data

ggplot(data=nhej_data, aes(x=Category, y=Percentage, fill=Orientation)) +
  geom_bar(stat = "summary", fun.y = "mean", position=position_dodge()) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8)) +
  facet_wrap(~ Knockin) + 
  theme(strip.text.x = element_text(size=16, color="red",
                                    face="bold")) + 
  scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels=c("del - no HDR" = "Deletions", "ins - no HDR" = "Insertions",
                            "wild-type" = "Wild-type")) +
  scale_y_continuous("Percentage", breaks = c(0, 10, 20, 30, 40, 50, 60, 70)) + 
  theme(
    axis.title.x = element_text(color="#993333", size= 18, face="bold", margin=margin(t = 20, unit = "pt")),
    axis.title.y = element_text(color="#993333", size= 18, face="bold", margin=margin(r = 20, unit = "pt")),
    axis.text.x = element_text(angle = 45, vjust=0.5, colour="grey20", size= 15, face="plain"),
    axis.text.y = element_text(colour="grey20",size= 15, face="plain"),
    legend.title = element_text( size = 15, face = "bold"),
    legend.text = element_text( size = 15, face = "plain"),
    legend.key.size = unit(1, "cm")
  )

ggsave("Figure4A.png", dpi = 600)

# prepare data for knock-in frequency analysis

KI_data <- subset(all_data, Category %in% c("correct HDR", "del - HDR", "ins - HDR", "unmapped - HDR"))

ggplot(data=KI_data, aes(x=Category, y=Percentage, fill=Orientation)) +
  geom_bar(stat = "summary", fun.y = "mean", position=position_dodge()) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8)) +
  geom_signif(comparisons = list(c("anti", "sense")), 
              map_signif_level=TRUE) +
  facet_wrap(~ Knockin) + 
  theme(strip.text.x = element_text(size=16, color="red",
                                    face="bold")) + 
  scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels=c("correct HDR" = "correct HDR KI", "del - HDR" = "deletions + HDR KI",
                            "ins - HDR" = "insertions + HDR KI", "unmapped - HDR" = "unmapped + HDR KI")) +
  scale_y_sqrt("Percentage, square root", breaks = c(0,0.02, 0.1, 0.2, 0.3, 0.5, 1, 1.5, 2, 2.5, 3, 3.4)) + 
  theme(
    axis.title.x = element_text(color="#993333", size= 18, face="bold", margin=margin(t = 20, unit = "pt")),
    axis.title.y = element_text(color="#993333", size= 18, face="bold", margin=margin(r = 20, unit = "pt")),
    axis.text.x = element_text(angle = 70, vjust=0.5, colour="grey20", size= 15, face="plain"),
    axis.text.y = element_text(colour="grey20",size= 15, face="plain"),
    legend.title = element_text( size = 15, face = "bold"),
    legend.text = element_text( size = 15, face = "plain"),
    legend.key.size = unit(1, "cm")
  )

ggsave("Figure4B.png", dpi = 600)
