# ----------------DATA PREPARATION----------------
# Load required package for data handling and analysis
library(tidyverse)

# Import the COVID-19 dataset
data<- read.csv("covidpredict-1.csv")

# Convert categorical variables to factors
# so they are treated correctly in the logistic regression model(for sex is not necessary)
data$comorbidity <- factor(data$comorbidity)
data$female <- factor(data$female)

 
# Split the dataset into development cohort for model training
dev<- subset(data, set == "dev")

# Split the dataset into validation cohort for model evaluation
val<-subset(data, set == "val")

# ---------------- MODEL DEVELOPMENT ----------------

# Fit a multivariable logistic regression model to predict
# in-hospital mortality using the development dataset
# Continuous predictors are assumed to have linear effects
model <- glm(mortality ~ rr + oxygen_sat + urea + crp +
               gcs+ age + female+
               comorbidity  ,
             data = dev,
             family = "binomial"
)

# Display model coefficients, standard errors,
# and statistical significance (p-values)
summary(model)

# Calculate adjusted odds ratios (ORs) and
# corresponding 95% confidence intervals
round(exp(cbind(OR = coef(model), confint(model))), 3)

# ----------reporting the results ( OR and 95% CI) in the table-----------
# Load packages for regression table formatting
library(gt)
library(finalfit)

# Define the dependent (outcome) variable and
# independent (predictor) variables for the model
dependent <- "mortality"
independent <- c("age", "female", "comorbidity", "rr",
                 "oxygen_sat", "gcs", "urea", "crp")


# Generate a formatted summary table of the multivariable
# logistic regression model, including adjusted odds ratios
# and corresponding confidence intervals
glmmulti_full <- dev |>
  glmmulti(dependent, independent) |>
  fit2df(
    explanatory_name = "Variables",
    estimate_name = "Adjusted OR - full model",
  )

# Display the regression results as a formatted table
glmmulti_full |>
  gt()


# ---------------- MODEL VALIDATION ----------------

# Use the developed logistic regression model to predict
# probability of in-hospital mortality in the validation dataset

preds <- predict(model,
                 newdata = val,
                 type = "response"
)


# Inspect predicted probabilities
preds

# Classify patients based on a probability threshold (0.5):
# patients with predicted probability < 0.5 are classified as survival (0),
# and ≥ 0.5 as mortality (1)
preds_outcome <- ifelse(preds < 0.5, 0, 1)


# Create a matrix comparing observed vs predicted outcomes
tab <- table(val$mortality, preds_outcome,
             dnn = c("observed", "predicted"))

# Display matrix
tab

# Calculate overall classification accuracy
# proportion of correctly classified patients
accuracy <- sum(diag(tab)) / sum(tab)
accuracy


# ----------------sensitivity---------------
# Sensitivity (True Positive Rate):
# ability of the model to correctly identify patients who died
# sensitivity = (True Positives)/(True Positives+False Negatives)

sensitivity <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
sensitivity

# ---------------- specificity---------------
# Specificity (True Negative Rate):
# ability of the model to correctly identify survivors
# specificity = (True Negatives)/(True Negatives+False Positives)

specificity <- tab[1, 1] / (tab[1, 1] + tab[1, 2])
specificity


# ----------------AUC and ROC curve---------------
# Load package for ROC analysis
library(pROC)

# Construct ROC curve using predicted probabilities
res <- roc(mortality ~ preds,
           data = val)

# Compute Area Under the Curve (AUC)
auc(res)

# Extract AUC value
res$auc

# plot ROC curve with AUC in title

ggroc(res, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = paste0("AUC = ", round(res$auc, 2)))

# Interpretation:
# The closer the ROC curve is to the top-left corner,
# the better the model discrimination ability.
# AUC closer to 1 indicates better performance..


# ----------------PPV---------------
# Positive Predictive Value (PPV):
# probability that predicted positives are truly positive
#PPV=(TP)/(FP+TP)
PPV <- tab[2,2] / (tab[2,2] + tab[1,2])
PPV

# ----------------NPV---------------
# Negative Predictive Value (NPV):
# probability that predicted negatives are truly negative
#NPV=(TN)/(FN+TN)
NPV <- tab[1,1] / (tab[1,1] + tab[2,1])
NPV

# ---------------- CALIBRATION PLOT ----------------

# Add predicted probabilities to validation dataset
# This is required for calibration assessment
val$preds <- preds

# Define a custom calibration plot function
# This function compares predicted probabilities vs observed outcomes
# and visualises model calibration
val.prob.ci.2 <- function(p, y,
                          xlab = "Predicted probability",
                          ylab = "Observed outcome") {
  
# Scatter plot of predicted probabilities vs observed outcomes
  plot(p, y,
       xlim = c(0, 1),
       ylim = c(0, 1),
       xlab = xlab,
       ylab = ylab,
       main = "Calibration Plot")
  
# Add reference line representing perfect calibration
# (ideal model where predicted = observed)
  abline(0, 1, col = "gray", lwd = 1)
  
# Add smoothed calibration curve using LOESS
# to assess systematic over- or underestimation
  lines(lowess(p, y), col = "blue", lwd = 2)
  
# Add legend for interpretation
  legend("topleft",
         legend = c("Perfect calibration (ideal model)",
                    "Model calibration trend"),
         col = c("gray", "blue"),
         lwd = c(1, 2),
         bty = "n")
}

# Run calibration plot function on validation data
val.prob.ci.2(p = val$preds, y = val$mortality)

# ---------------- DECISION CURVE ANALYSIS ----------------
# Load package for decision curve analysis
library(dcurves)

# Perform decision curve analysis to evaluate clinical utility
# across a range of threshold probabilities
dca(mortality ~ preds,
    data = val,
    time = 1.5,
    thresholds = seq(0, 1, 0.01),
    label = list(pr_failure18 = "Prediction Model")
) %>%
  plot(smooth = TRUE)


# ---------------- EXPECTED VALUE OF PERFECT INFORMATION (EVPI) ----------------
library(MASS)
#Step 1:add the predicted probabilities to validation data
val$preds <- preds

# Number of observations in validation dataset
n <- nrow(val)

# Threshold probability for clinical decision-making
z <- 0.1

#Step 2:The bootstrap Monte Carlo simulation
# Number of bootstrap iterations for Monte Carlo simulation
N <- 10000

# Initialize vectors to store results
NBmodel <- NBall <- NBmax <-rep(0, N)

# Monte Carlo bootstrap simulation to estimate uncertainty
for(i in 1:N){
# Resample validation dataset with replacement
  bsdata <- val[sample(1:n, n, replace = TRUE),]
# Net benefit of treating all patients  
  NBall[i] <- mean(bsdata$mortality - (1 - bsdata$mortality)*z/(1-z))
# Net benefit of prediction model-based decision 
  NBmodel[i] <- mean((bsdata$preds > z) * 
                       (bsdata$mortality - (1 - bsdata$mortality)*z/(1-z)))
}
#Step 3:EVPI calculation
# Calculate Expected Value of Perfect Information (EVPI)
# Difference between perfect decision-making and current model uncertainty
EVPI <- mean(pmax(0, NBmodel, NBall)) - 
  max(0, mean(NBmodel), mean(NBall))

EVPI






