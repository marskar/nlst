#### Load packages ####
pkgs = c("here", "dplyr","ggplot2","haven","reshape", "purrr")
lapply(pkgs, require, character.only  = TRUE)

# Read in my personal defined functions
source(here("hilary_functions.R"))

#### Load data  ####
# Load abnormalities dataset https://biometry.nci.nih.gov/plcosub/lung/2009-00516-201204-0017-lung-mortality-risk/nlst/mortalityrisk_abn_080216.csv.zip/@@display-file from https://biometry.nci.nih.gov/plcosub/lung/2009-00516-201204-0017-lung-mortality-risk/nlst
abn <- read.csv(here("mortalityrisk_abn_080216.csv"))

#### Load Lung-RADS data ####

lrads.pl = readRDS('lungrads.rds')

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

abn_lrads_merged <- abn.pl %>%
    # Replace study year with interval
    mutate(interval = case_when(
    STUDY_YR==0 ~ 1,
    STUDY_YR==1 ~ 2, 
    STUDY_YR==2 ~ 3
    )) %>% 
    select(-STUDY_YR) %>% 
    # Set a ceiling of 60 mm to diameter values to prevent influence of extreme outliers
    mutate(longest.diam = ifelse(longest.diam >= 60, 60, longest.diam)) %>% 
    # Add a variable for presence of adenopathy or consolidation
    mutate(adenop.consol = as.numeric(adenopathy == 1 | consolidation == 1)) %>% 
    # Combine the lrads and abn datasets
    merge(lrads.pl, by = c("pid", "interval"), all.x = TRUE, all.y = TRUE) %>% 
    mutate(LRcat = as.factor(LRcat))

saveRDS(abn_lrads_merged, file = "abn_lrads_merged.rds")





    # Make an LRcat variable relevant for negative screen groups
    mutate(LRcatcol.neg = case_when(
        LRcat=="1" ~ 1, 
        LRcat=="2" ~ 2, 
        LRcat %in% c("3","3 or 4A","3,4A, or 4B","4A","4B","4X","4A or 4B") ~ 3
        )) %>% 
    # Make an LRcat variable relevant for positive screen groups - this assigns 1 to missing!)
    mutate(LRcatcol.pos = case_when(
        LRcat=="2" ~ 2,
        LRcat=="3" ~ 3,
        LRcat=="4A" ~ 4, 
        LRcat=="4B" ~ 5,
        LRcat=="4X" ~ 6,
        LRcat %in% c("3 or 4A", "3,4A, or 4B", "4A or 4B") ~ 7))

# Note: the Lung-RADS category is missing for 485 screens, but 470 are at T2 (15 at T1).
table(abn_lrads_merged[is.na(abn_lrads_merged$LRcat),]$interval)
table(abn_lrads_merged[is.na(abn_lrads_merged$LRcat),]$interval)
abn_lrads_merged$LRcat
# Save to disk
save(abn.pl.all, file=here("abn.spl.20181126.rdata"))
