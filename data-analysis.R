library(tidyverse)
library(fixest)
library(vcd)
library(psych)
library(modelr)
library(boot)
library(lme4)
library(kableExtra)
library(texreg)

# Loading in the data ----

# The data on social exchange, or sharing one's resources with others
sex <- read_delim(file = "social-exchange-data.csv", delim = ";")

# The data from the post-experimental questionnaire, concerning one's feelings towards in-group and out-group
sid <- read_delim(file = "social-identity-data.csv", delim = ";") # the social identity data

# Table 1 in the main text ----

tab_1 <- sex |>
  filter(round == 1) |>
  count(visible, pair_type, gave)

tab_1 <- xtabs(n ~ pair_type + gave + visible, data = tab_1)

# Hidden differences

tab_1a <- tab_1[2:1, 2:1, 1]

# Visible differences

tab_1b <- tab_1[2:1, 2:1, 2]

# The odds that resources are shared in round 1 of the social exchange process when the pair is intra- rather than intergroup

summary(oddsratio(tab_1a)) # Hidden differences
summary(oddsratio(tab_1b)) # Visible differences
