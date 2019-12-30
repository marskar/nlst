library(h2o) 
h2o.init()
csv <- "data/nlst_abn_lcp_pre.csv"
h2o_data <- h2o.importFile(here(csv))
dim(h2o_data)
names(data)
names(h2o_data)

# prerisks <- rep("log1yrisk", length(x)-1)
# interact_pairs <- mapply(c, prerisks, x[-length(x)], SIMPLIFY=FALSE)
h2o_data$case <- as.factor(h2o_data$case)  #encode the binary response as a factor
h2o.levels(h2o_data$case)  # show the factor levels

splits <- h2o.splitFrame(data = h2o_data, 
                         ratios = c(0.7),  # partition data into 70%, 15% chunks
                         destination_frames = c("train", "test"), # frame ID (not required)
                         seed = 1)  # setting a seed will guarantee reproducibility
train <- splits[[1]]
test <- splits[[2]]

y <- "case"
# Remove the columns that are irrelevant or correlated with the outcome
# x <- setdiff(names(h2o_data), c(
#     y,
#     "candx_days",
#     "cancyr",
#     "canc_rpt_link",
#     "incidence.years",
#     "years.followed",
#     "lung.cancer.death",
#     "pid",
#     "screen_group",
#     "rndgroup",
#     "truefalse_scrnres_ly0",
#     "truefalse_scrnres_ly1",
#     "truefalse_scrnres_ly2",
#     "center",
#     "lss",
#     "scr_res0",
#     "scr_res1",
#     "scr_res2",
#     "scr_iso0",
#     "scr_iso1",
#     "scr_iso2",
#     "other.cause.death"
# ))
x <- c(
    "log1yrisk",
    "emphysema",
    "adenopathy",
    "consolidation"
)

print(x)
length(x)
dim(h2o_data)

h2o_summarize_model <- function(model) {
    print(paste("AIC =", h2o.aic(model)))
    print(paste("AUC =", h2o.auc(model)))
    h2o.varimp_plot(model)
}
h2o_glm_lcrat_only <- h2o.glm(
    x = "log1yrisk",
    y = y,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o_summarize_model(h2o_glm_lcrat_only)

h2o_glm_lcrat_ct <- h2o.glm(
    x = x,
    y = y,
    interaction_pairs = list(
        c("log1yrisk", "emphysema"),
        c("log1yrisk", "consolidation"),
        c("log1yrisk", "adenopathy")
    ),
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o_summarize_model(h2o_glm_lcrat_ct)

x = c(
    "log1yrisk",
    "emphysema",
    "consolidation",
    "adenopathy",
    "max_lcp_score"
    )

h2o_glm_lcrat_ct_lcp <- h2o.glm(
    x = x,
    y = y,
    interaction_pairs = list(
        c("log1yrisk", "emphysema"),
        c("log1yrisk", "consolidation"),
        c("log1yrisk", "adenopathy"),
        c("log1yrisk", "max_lcp_score")
    ),
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o_summarize_model(h2o_glm_lcrat_ct_lcp)


x = c(
    "log1yrisk",
    "emphysema",
    "consolidation",
    "adenopathy",
    "max_lcp_score"
    )

h2o_glm_lcrat_ct_all <- h2o.glm(
    x = x,
    y = y,
    interactions = x,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o_summarize_model(h2o_glm_lcrat_ct_all)

x = c(
    "log1yrisk",
    "emphysema",
    "consolidation",
    "adenopathy",
    "max_lcp_score"
    )

h2o_glm_lcrat_ct_all <- h2o.glm(
    x = x,
    y = y,
    interactions = x,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)
h2o_summarize_model(h2o_glm_lcrat_ct_all)

# 2. LCRAT + CT + LCP_SCORE
# 3. LCRAT + CT + EMPH_SCORE
# 4. LCRAT + CT + LCP_SCORE + EMPH_SCORE
# Describe these models
# Methods: LCRAT + CT (Hilary) + LCP_SCORE (Optellum)
# First, prescreen_risk interactions with all binary features
# Second, prescreen_risk interactions with all discrete (categorical) features
glm_fit1 <- h2o.glm(
    x = x,
    y = y,
    training_frame = train,
    family = "binomial",
    lambda_search = TRUE
)  #Like glm() and glmnet(), h2o.glm() has the family argument
