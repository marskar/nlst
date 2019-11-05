#### Load packages ####
library(dplyr)
library(here)
library(readr)


#### Load data  ####
# Load abnormalities dataset https://biometry.nci.nih.gov/plcosub/lung/2009-00516-201204-0017-lung-mortality-risk/nlst/mortalityrisk_abn_080216.csv.zip/@@display-file from https://biometry.nci.nih.gov/plcosub/lung/2009-00516-201204-0017-lung-mortality-risk/nlst
abn <- read_csv(here("data/mortalityrisk_abn_080216.csv"))

#### Load Lung-RADS data ####
# Created by script 01
lrads = read_rds(here('data/lungrads.rds'))

#### Organize abnormality data ####
abn_lrads_merged <- abn %>%
    group_by(pid, STUDY_YR) %>%
    summarise(
        longest_diam = max(SCT_LONG_DIA, na.rm = T),
        longest_perp_diam = max(SCT_PERP_DIA, na.rm = T),
        any_nodule = as.numeric(any(SCT_AB_DESC == 51, na.rm = T)),
        nodule_count = sum(SCT_AB_DESC == 51, na.rm = T),
        any_micronodule = as.numeric(any(SCT_AB_DESC == 52, na.rm = T)),
        benign_nodule = as.numeric(any(SCT_AB_DESC == 53, na.rm = T)),
        atelectasis = as.numeric(any(SCT_AB_DESC == 54, na.rm = T)),
        pleur_thick_eff = as.numeric(any(SCT_AB_DESC == 55, na.rm = T)),
        adenopathy = as.numeric(any(SCT_AB_DESC == 56, na.rm = T)),
        chest_wall_abn = as.numeric(any(SCT_AB_DESC == 57, na.rm = T)),
        consolidation = as.numeric(any(SCT_AB_DESC == 58, na.rm = T)),
        emphysema = as.numeric(any(SCT_AB_DESC == 59, na.rm = T)),
        cardio_abn = as.numeric(any(SCT_AB_DESC == 60, na.rm = T)),
        opac_fibr = as.numeric(any(SCT_AB_DESC == 61, na.rm = T)),
        nod6_not_susp = as.numeric(any(SCT_AB_DESC == 62, na.rm = T)),
        other_above = as.numeric(any(SCT_AB_DESC == 63, na.rm = T)),
        other_below = as.numeric(any(SCT_AB_DESC == 64, na.rm = T)),
        other_minor = as.numeric(any(SCT_AB_DESC == 65, na.rm = T)),
        any_right_upper = as.numeric(any(SCT_EPI_LOC == 1, na.rm = T)),
        any_right_mid = as.numeric(any(SCT_EPI_LOC == 2, na.rm = T)),
        any_right_lower = as.numeric(any(SCT_EPI_LOC == 3, na.rm = T)),
        any_left_upper = as.numeric(any(SCT_EPI_LOC == 4, na.rm = T)),
        any_lingula = as.numeric(any(SCT_EPI_LOC == 5, na.rm = T)),
        any_left_lower = as.numeric(any(SCT_EPI_LOC == 6, na.rm = T)),
        any_upper = as.numeric(any(SCT_EPI_LOC %in% c(1, 4), na.rm = T)),
        any_lower = as.numeric(any(SCT_EPI_LOC %in% c(3, 6), na.rm = T)),
        any_spiculation = as.numeric(any(SCT_MARGINS == 1, na.rm = T)),
        any_smooth = as.numeric(any(SCT_MARGINS == 2, na.rm = T)),
        any_poor_def = as.numeric(any(SCT_MARGINS == 3, na.rm = T)),
        any_margin_unab = as.numeric(any(SCT_MARGINS == 9, na.rm = T)),
        any_soft_tissue = as.numeric(any(SCT_PRE_ATT == 1, na.rm = T)),
        any_GG = as.numeric(any(SCT_PRE_ATT == 2, na.rm = T)),
        any_mixed = as.numeric(any(SCT_PRE_ATT == 3, na.rm = T)),
        any_fluid = as.numeric(any(SCT_PRE_ATT == 4, na.rm = T)),
        any_fat = as.numeric(any(SCT_PRE_ATT == 6, na.rm = T)),
        any_other_att = as.numeric(any(SCT_PRE_ATT %in% c(4, 6, 7), na.rm = T)),
        # EDITED 21 JUNE 2018 to include fluid + fat + other
        any_unable_att = as.numeric(any(SCT_PRE_ATT == 9, na.rm = T)),
        susp_change_att = as.numeric(any(sct_ab_attn == 2, na.rm = T)),
        any_growth = as.numeric(any(sct_ab_gwth == 2, na.rm = T)),
        any_new_nodule = as.numeric(any(
            SCT_AB_DESC == 51 &
                sct_ab_preExist == 1, na.rm = T
        )),
        # EDITED 12 OCT 2018 to add variable for new nodule
        any_new_nodule_4_7mm = as.numeric(
            any(
                SCT_AB_DESC == 51 &
                    SCT_LONG_DIA >= 4 &
                    SCT_LONG_DIA <= 7 &
                    sct_ab_preExist == 1,
                na.rm = T
            )
        )  # EDITED 12 OCT 2018 to add variable for new nodules 4-7mm
    ) %>%
    # Replace study year with interval
    rename(interval = STUDY_YR) %>%
    mutate(
        diam_cat = case_when(
            longest_diam == 0 ~ 1,
            longest_diam > 0 & longest_diam <= 5 ~ 2,
            longest_diam > 5 & longest_diam <= 7 ~ 3,
            longest_diam > 7 & longest_diam <= 10 ~ 4,
            longest_diam > 10 & longest_diam <= 13 ~ 5,
            longest_diam > 13 & longest_diam < 100 ~ 6
        )
    ) %>%
    mutate(diam_cat = as.factor(diam_cat)) %>%
    # Add a variable for presence of adenopathy or consolidation
    mutate(adenop_consol = as.numeric(adenopathy == 1 |
                                          consolidation == 1)) %>%
    # Combine the lrads and abn datasets
    merge(
        lrads,
        by = c("pid", "interval"),
        all.x = TRUE,
        all.y = TRUE
    ) %>%
    mutate(LRcat = as.factor(LRcat)) %>%
    # Make an LRcat variable relevant for negative screen groups
    mutate(LRcatcol.neg = case_when(
        LRcat == "1" ~ 1,
        LRcat == "2" ~ 2,
        LRcat %in% c("3", "3 or 4A", "3,4A, or 4B", "4A", "4B", "4X", "4A or 4B") ~ 3
    )) %>%
    # Make an LRcat variable relevant for positive screen groups - this assigns 1 to missing!)
    mutate(
        LRcatcol.pos = case_when(
            LRcat == "2" ~ 2,
            LRcat == "3" ~ 3,
            LRcat == "4A" ~ 4,
            LRcat == "4B" ~ 5,
            LRcat == "4X" ~ 6,
            LRcat %in% c("3 or 4A", "3,4A, or 4B", "4A or 4B") ~ 7
        )
    ) %>%
    write_rds(here("data/abn_lrads_merged.rds"))
