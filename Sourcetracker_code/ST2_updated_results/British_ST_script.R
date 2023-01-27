# summarize sourcetracker results 

#### set working directory ####
setwd("/Users/sterl/OneDrive/Desktop/microArch/British/ST_Results/")

#### load data with modern dental calculus and plaque ####
data <- read.csv("ST_British_results-calculus-plaque.csv")

#### load packages ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)
library(extrafont)
library(scales)
library(reshape2)

#### reformat data ####
data.df<-as.data.frame.matrix(data)
data.long<-melt(data.df, id.vars = c("SampleID"), value.name = "Proportion")
names(data.long)[2]<-paste("Environment")

# uncomment to create a tiff image 
#tiff("British_Sourcetracker_results.tiff", units = "in", width = 6, height = 4, res = 300)

#### stacked-bar plot ####
ggplot() + geom_bar(aes(y = data.long$Proportion, x = data.long$SampleID, fill = data.long$Environment), stat = "identity") +
  labs(x = "Samples", y = "Percentage") +
  theme_bw() +
  scale_fill_manual(name = "Environment", values = c("firebrick", "lightskyblue", "thistle", "darkseagreen", "seashell4", "peachpuff"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)
        )
# dev.off()

##### Environment averages #####

mean(data$Plaque)
# [1] 52.4588

mean(data$ModernCalculus)
# [1] 3.766341

mean(data$Skin)
# [1] 1.709457

mean(data$EBC)
# [1] 0.07873188

mean(data$Soil)
# [1] 0.1102174

mean(data$Unknown)
# [1] 41.87645


#### stacked-bar plots #### 

# load data and reformat the data 
# Data set includes species in at least 77 samples with 5 percent relative abundance 
British_samples<-read.csv("British-samples-taxa-5percent-80samples.csv", header = TRUE)
species.df<-as.data.frame(British_samples)
species.long<-melt(species.df, id.vars = "Samples", value.name = "Proportion")

# change the name of the second column to Species 
names(species.long)[2]<-paste("Species")

# create the plot 
ggplot(species.long, aes(fill=Species, y=Proportion, x=Samples)) + 
  geom_bar(position = "fill", stat = "identity", color = "black") +
  theme(axis.ticks.x = element_blank(),
        
        
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(name = "Environment", 
                    values = c("firebrick", "lightskyblue", "thistle", "darkseagreen","wheat", "turquoise",  "blueviolet",
                               "deeppink1", "orange3", "purple4", "royalblue2", "forestgreen", "lightsalmon2", "lightgray",
                               "mistyrose4", "yellow3", "maroon", "sienna1", "cyan", "tomato2", "peachpuff", "gray80"))

#### load data - wAncientCalculus as a source ####

data <- read.csv("ST_British-calculus-wAncient.csv")

##### reformat data #####
data.df<-as.data.frame.matrix(data)

data.long<-melt(data.df, id.vars = c("SampleID"), value.name = "Proportion")
names(data.long)[2]<-paste("Environment")

tiff("British_Sourcetracker_results-wAncientCalculus.tiff", units = "in", width = 6, height = 4, res = 300)
##### stacked-bar plot #####
ggplot() + geom_bar(aes(y = data.long$Proportion, x = data.long$SampleID, fill = data.long$Environment), stat = "identity") +
  labs(x = "Samples", y = "Percentage") +
  theme_bw() +
  scale_fill_manual(name = "Environment", values = c("thistle", "lightskyblue", "firebrick", "darkseagreen", "seashell4", "peachpuff"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)
  )
dev.off()

##### Environment averages #####

mean(data$AncientCalculus)
# [1] 0.9544018

mean(data$Skin)
# [1] 0.03248913

mean(data$EBC)
# [1] 0.001190217

mean(data$Soil)
# [1] 0.001444565

mean(data$Unknown)
# [1] 0.01047428
