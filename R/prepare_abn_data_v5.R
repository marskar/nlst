
# This program prepares person-level data regarding screen-detected abnormalities
# to then be used in analysis_nlst_vX programs (begin 8 Aug 2016)

#### Load data and packages ####
rm(list=ls(all=TRUE))  #Clears console
setwd("/Users/robbinsh_mac/Dropbox/IARC/NCI/NLST/Analysis/Nodule data")
# Read in my personal defined functions
source("/Users/robbinsh_mac/Dropbox/IARC/Resources/R functions/hilary_functions_iarc.R")
# Read in library functions including Stephanie's coxph.risk which I installed from the local tar.gz file
packages <- c("dplyr","ggplot2","sas7bdat","reshape")
lapply(packages, require, c = T)
# Load lung-RADS and abnormalities datasets
lrads <- read.sas7bdat("/Users/robbinsh_mac/Dropbox/IARC/NCI/NLST/Analysis/Nodule data/lungrads1.sas7bdat")
abn <- read.csv("/Users/robbinsh_mac/Dropbox/IARC/NCI/NLST/Analysis/Nodule data/mortalityrisk_abn_080216.csv")

#### Organize Lung-RADS data ####
# First, reshape
lrads.pl <- melt(lrads, id="pid") 
lrads.pl$interval <- 1*(lrads.pl$variable=="slungrad0") +
  2*(lrads.pl$variable=="slungrad1") + 3*(lrads.pl$variable=="slungrad2")
lrads.pl <- subset(lrads.pl, select = -variable)
lrads.pl <- filter(lrads.pl, value!="Not Done")
colnames(lrads.pl) <- c("pid","LRcat","interval")
table(lrads.pl$interval, lrads.pl$LRcat)


#### Organize abnormality data ####
abn.pl <- abn %>% group_by(pid, STUDY_YR) %>%
  summarise(longest.diam = max(SCT_LONG_DIA, na.rm=T),
            #longest.perp.diam = max(SCT_PERP_DIA, na.rm=T),
            any.nodule = as.numeric(any(SCT_AB_DESC==51, na.rm=T)),
            nodule.count = sum(SCT_AB_DESC==51, na.rm=T),
            any.micronodule = as.numeric(any(SCT_AB_DESC==52, na.rm=T)),
            benign.nodule = as.numeric(any(SCT_AB_DESC==53, na.rm=T)),
            atelectasis = as.numeric(any(SCT_AB_DESC==54, na.rm=T)),
            pleur.thick.eff = as.numeric(any(SCT_AB_DESC==55, na.rm=T)),
            adenopathy = as.numeric(any(SCT_AB_DESC==56, na.rm=T)),
            # chest.wall.abn = as.numeric(any(SCT_AB_DESC==57, na.rm=T)),
            consolidation = as.numeric(any(SCT_AB_DESC==58, na.rm=T)),
            emphysema = as.numeric(any(SCT_AB_DESC==59, na.rm=T)),
            # cardio.abn = as.numeric(any(SCT_AB_DESC==60, na.rm=T)),
            opac.fibr = as.numeric(any(SCT_AB_DESC==61, na.rm=T)),
            nod6.not.susp = as.numeric(any(SCT_AB_DESC==62, na.rm=T)),
            other.above = as.numeric(any(SCT_AB_DESC==63, na.rm=T)),
            other.below = as.numeric(any(SCT_AB_DESC==64, na.rm=T)),
            # other.minor = as.numeric(any(SCT_AB_DESC==65, na.rm=T)),
            # any.right.upper = as.numeric(any(SCT_EPI_LOC==1, na.rm=T)),
            any.right.mid = as.numeric(any(SCT_EPI_LOC==2, na.rm=T)),
            # any.right.lower = as.numeric(any(SCT_EPI_LOC==3, na.rm=T)),
            # any.left.upper = as.numeric(any(SCT_EPI_LOC==4, na.rm=T)),
            any.lingula = as.numeric(any(SCT_EPI_LOC==5, na.rm=T)),
            # any.left.lower = as.numeric(any(SCT_EPI_LOC==6, na.rm=T)),
            any.upper = as.numeric(any(SCT_EPI_LOC %in% c(1,4), na.rm=T)),
            any.lower = as.numeric(any(SCT_EPI_LOC %in% c(3,6), na.rm=T)),
            any.spiculation = as.numeric(any(SCT_MARGINS==1, na.rm=T)),
            any.smooth = as.numeric(any(SCT_MARGINS==2, na.rm=T)),
            any.poor.def = as.numeric(any(SCT_MARGINS==3, na.rm=T)),
            any.margin.unab = as.numeric(any(SCT_MARGINS==9, na.rm=T)),
            any.soft.tissue = as.numeric(any(SCT_PRE_ATT==1, na.rm=T)),
            any.GG = as.numeric(any(SCT_PRE_ATT==2, na.rm=T)),
            any.mixed = as.numeric(any(SCT_PRE_ATT==3, na.rm=T)),
            # any.fluid = as.numeric(any(SCT_PRE_ATT==4, na.rm=T)),
            # any.fat = as.numeric(any(SCT_PRE_ATT==6, na.rm=T)),
            any.other.att = as.numeric(any(SCT_PRE_ATT %in% c(4,6,7), na.rm=T)),   # EDITED 21 JUNE 2018 to include fluid + fat + other
            any.unable.att = as.numeric(any(SCT_PRE_ATT==9, na.rm=T)),
            susp.change.att = as.numeric(any(sct_ab_attn==2, na.rm=T)),
            any.growth = as.numeric(any(sct_ab_gwth==2, na.rm=T)),
            any.new.nodule = as.numeric(any(SCT_AB_DESC==51 & sct_ab_preExist==1, na.rm=T)),   # EDITED 12 OCT 2018 to add variable for new nodule
            any.new.nodule.4.7.mm = as.numeric(any(SCT_AB_DESC==51 & SCT_LONG_DIA>=4 & SCT_LONG_DIA<=7 & sct_ab_preExist==1, na.rm=T))  # EDITED 12 OCT 2018 to add variable for new nodules 4-7mm
            )
# abn.pl$lobe.count <- with(abn.pl, any.right.upper+any.right.mid+any.right.lower+
                           #  any.left.upper+any.lingula+any.left.lower)
abn.pl$interval <- 1*(abn.pl$STUDY_YR==0) + 2*(abn.pl$STUDY_YR==1) + 3*(abn.pl$STUDY_YR==2)
abn.pl <- subset(abn.pl, select = -STUDY_YR)
# Set a ceiling of 60 mm to diameter values to prevent influence of extreme outliers
abn.pl$longest.diam <- ifelse(abn.pl$longest.diam>=60, 60, abn.pl$longest.diam)
# abn.pl$longest.perp.diam <- ifelse(abn.pl$longest.perp.diam>=60, 60, abn.pl$longest.perp.diam)
# Add a variable for presence of adenopathy or consolidation
abn.pl$adenop.consol <- as.numeric(abn.pl$adenopathy==1 | abn.pl$consolidation==1)

#### Combine these ####
abn.pl.all <- merge(abn.pl, lrads.pl, by = c("pid", "interval"), all.x=T, all.y=T)
# Make an LRcat variable relevant for negative screen groups
abn.pl.all <- mutate(abn.pl.all, LRcatcol.neg = as.factor(1*(LRcat=="1") + 2*(LRcat=="2") +
                                          3*(LRcat %in% c("3","3 or 4A","3,4A, or 4B","4A","4B","4X","4A or 4B"))))
abn.pl.all$LRcatcol.neg <- ifelse(abn.pl.all$LRcatcol.neg==0, NA, abn.pl.all$LRcatcol.neg)
abn.pl.all$LRcatcol.neg <- factor(abn.pl.all$LRcatcol.neg, levels=c(1, 2, 3), labels = c("1", "2", "3 or higher"))
# Make an LRcat variable relevant for positive screen groups - this assigns 1 to missing!
abn.pl.all <- mutate(abn.pl.all, LRcatcol.pos = as.factor(2*(LRcat=="2") + 3*(LRcat=="3") +
                      4*(LRcat=="4A") + 5*(LRcat=="4B") + 6*(LRcat=="4X") +
                        7*(LRcat %in% c("3 or 4A", "3,4A, or 4B", "4A or 4B"))))
abn.pl.all$LRcatcol.pos <- ifelse(abn.pl.all$LRcatcol.pos==0, NA, abn.pl.all$LRcatcol.pos)
abn.pl.all$LRcatcol.pos <- factor(abn.pl.all$LRcatcol.pos, levels=c(2,3,4,5,6,7), 
                      labels = c("2", "3", "4A", "4B", "4X", "3,4A,or4B"))
# Note: the Lung-RADS category is missing for 485 screens, but 470 are at T2 (15 at T1).
table(abn.pl.all[is.na(abn.pl.all$LRcat),]$interval)

# Save to disk
save(abn.pl.all, file="/Users/robbinsh_mac/Dropbox/IARC/NCI/NLST/Analysis/Nodule data/abn.spl.20181126.rdata")
