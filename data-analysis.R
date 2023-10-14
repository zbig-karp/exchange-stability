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
  count(visible, pair_type, gave) |>
  group_by(visible) |>
  nest()

tab_1 <- tab_1 |>
  mutate(tables = map(data, ~xtabs(n ~ pair_type + gave, data = .)),
         ors = map(tables, loddsratio),
         or_coef = map_dbl(ors, ~.$coefficients),
         or_serr = map_dbl(ors, ~sqrt(.$vcov[1, 1])),
         or_ztest = or_coef/or_serr,
         or_pvalue = 2 * pnorm(abs(or_ztest), lower.tail = FALSE))

# Values in columns or_coef, or_serr, and or_pvalue are reported in the main text, on p. 10. The column tables is a list with 2 x 2 tables of counts for each experimental condition.

# Figure 1: Probability of tie presence by round and the type of pair ----

tab_2 <- sex |>
  group_by(sessionid, visible, split, pair_type, round) |>
  summarise(m = mean(gave)) |>
  group_by(visible, pair_type, round) |>
  summarise(m = mean(m), .groups = "drop")

tab_2 |>
  ggplot(mapping = aes(x = round, y = m, linetype = pair_type)) + 
  geom_line() + 
  facet_wrap(~visible, 
             labeller = labeller(visible = c("No" = "Hidden social differences", 
                                             "Yes" = "Visible social differences"))) + 
  theme_bw() + 
  theme(legend.position = "top") + 
  labs(x = "Round", y = "Probability of tie presence", linetype = "Pair type") + 
  scale_linetype_manual(values = c("dotted", "solid"))

# Figure 2: Proportions of null, asymmetric, and mutual pairs by round, visibility of social differences, and whether the pair is intra-group or inter-group ----

p_pair <- sex |>
  mutate(pair = paste(pmin(subject, alter), pmax(subject, alter), sep = "-")) |>
  group_by(sessionid, visible, split, round, pair_type, pair) |>
  summarise(state = sum(gave), .groups = "drop")

p_pair <- p_pair |>
  count(visible, round, pair_type, state) |>
  group_by(visible, round) |>
  mutate(prop = prop.table(n)) |>
  ungroup()

p_pair <- p_pair |>
  mutate(visible = factor(visible, levels = c("Yes", "No"), labels = c("Visible differences", "Hidden differences")),
         state = factor(state, levels = 0:2, labels = c("Null", "Asymmetric", "Mutual"))) |>
  select(-n)

p_pair |>
  ggplot(mapping = aes(x = round, y = prop, linetype = pair_type)) + 
  geom_line() + 
  facet_grid(visible ~ state) + 
  theme_bw() + 
  theme(legend.position = "top") + 
  labs(x = "Round", y = "Proportion", linetype = "Pair type") + 
  scale_linetype_manual(values = c("dotted", "solid"))

# Figure 3: Average differences in the subjective assessments of the in-group and the out-group in the post-experimental questionnaire ----

# A difference in feelings towards in-group and out-group
sid <- sid |>
  mutate(ident = eval_in - eval_out)

# Converting the data to wide format and adding a new variable -- `score` -- representing average in-group preference across the four domains (dimensions)
sid_wide <- sid  |>
  select(1:5, domain, ident) |>
  pivot_wider(names_from = domain, values_from = ident) |>
  mutate(score = (belong + common + close + liking)/4)

# The principal component analysis on the 4 dimensions of in-group preference (see footnote 8 in the paper)
m.pca1 <-  princomp(formula = ~ belong + close + common + liking, data = sid_wide)

# Cronbach alpha for the index of in-group preference (see p. 11 in the paper)
m.alpha <- alpha(x = select(sid_wide, belong:close))

# Difference in feelings towards the in-group/out-group by experimental condition (see pp. 11-12 in the paper)
m.ttest <- t.test(score ~ visible, data = sid_wide)

# The plot
sid |>
  group_by(domain) |>
  summarise(m = mean(ident), 
            s = sd(ident),
            n = n()) |>
  ggplot(mapping = aes(x = domain, y = m)) + 
  geom_ref_line(h = 0, colour = "#bdbdbd", size = 0.5) +
  geom_errorbar(mapping = aes(ymin = m + qt(0.025, df = n - 1) * s/sqrt(n), 
                              ymax = m + qt(0.975, df = n - 1) * s/sqrt(n)),
                width = 0.1) + 
  geom_point(size = 2) + 
  theme_bw() +
  labs(x = "Domain", y = "Average difference", colour = NULL) + 
  scale_x_discrete(labels = c("Belongingness", "Closeness", "Commonality", "Liking"))

# Table 2: Estimates of the theoretical parameters in the attraction model ----

# A function used to obtain bootstrap estimates of the parameters of interest. It has just two arguments: the dataset and the index of rows defining a bootstrap sample. The dataset in question is sex, the one which contains data on social exchange. The function fits a logistic regression model as specified in Section 1 of the Online Supplement. Then, it extracts the regression coefficients and calculates the values of the parameters in the attraction model using formulas shown in Table 1 of the Online Supplement.

mparams <- function(d, i) {
  d <- d[i, ]
  m1 <- glm(gave ~ received_last * gave_last * pair_type, data = d, family = binomial())
  nd <- expand_grid(gave_last = 0:1, pair_type = c("Inter-group", "Intra-group"), received_last = 0:1)
  nd$pred <- predict(object = m1, newdata = nd, type = "response")
  with(nd,
       c(
         "df" = pred[1],
         "pifb" = 1 - (1 - pred[2])/(1 - pred[1]),
         "ef" = 1 - (1 - pred[3])/(1 - pred[1]),
         "pifw" = 1 - (1 - pred[4])/(1 - pred[3]),
         "dr" = pred[5],
         "pirb" = 1 - (1 - pred[6])/(1 - pred[5]),
         "er" = 1 - (1 - pred[7])/(1 - pred[5]),
         "pirw" = 1 - (1 - pred[8])/(1 - pred[7]),
         "drecf" = (1 - pred[2])/(1 - pred[1]) - (1 - pred[4])/(1 - pred[3]),
         "drecr" = (1 - pred[6])/(1 - pred[5]) - (1 - pred[8])/(1 - pred[7]),
         "quanc" = pred[7] + pred[4] - 2 * pred[7] * pred[4] - pred[5] - pred[2] + 
           2 * pred[5] * pred[2]
       ))
}

# Splitting the data by the experimental condition (i.e., visibility of social differences)

models <- sex |>
  group_by(visible) |>
  nest() 

# Using boot::boot() on the function `mparams` defined above to obtain the parameters estimates along with bootstrap standard errors

models <- models |>
  mutate(mpars = map(data, ~boot(data = ., statistic = mparams, R = 1000)),
         regtabs = map(mpars, ~createTexreg(
           coef.names = c("$d_f$", "$\\pi_{fb}$", "$\\epsilon_f$", "$\\pi_{fw}$", 
                          "$d_r$", "$\\pi_{rb}$", "$\\epsilon_r$", "$\\pi_{rw}$", 
                          "$\\delta_1$", "$\\delta_2$", "quanc"),
           coef = .$t0,
           se = apply(X = .$t, MARGIN = 2, FUN = sd)
         )))

# This prints a table with results on the screen

screenreg(l = models$regtabs, 
          digits = 3, 
          omit.coef = "(delta)|(quanc)", 
          reorder.coef = c(1, 2, 4, 3, 5, 6, 8, 7),
          custom.model.names = c("Visible", "Hidden"),
          custom.header = list("Social differences" = 1:2),
          custom.note = "NB: Bootstrap standard errors in parentheses.",
          groups = list("Tie formation" = 1:4, "Tie renewal" = 5:8))

# This command creates a table with results in the LaTeX format

texreg(l = models$regtabs, 
       digits = 3, 
       omit.coef = "(delta)|(quanc)", 
       reorder.coef = c(1, 2, 4, 3, 5, 6, 8, 7),
       custom.model.names = c("Visible", "Hidden"),
       custom.header = list("Social differences" = 1:2),
       custom.note = "NB: Bootstrap standard errors in parentheses.",
       groups = list("Tie formation" = 1:4, "Tie renewal" = 5:8),
       booktabs = TRUE,
       dcolumn = TRUE, use.packages = FALSE,
       caption = "Estimates of the theoretical parameters in the attraction model",
       caption.above = TRUE,
       label = "tab:params"
)

# Table 3: Estimates of the hypothesised effects and their p-values ----

# Converting the original data into a new dataset where each unordered pair of actors (i, j) in round t is classified as null, asymmetric, or mutual. Specifically, the code below defines 3 variables for each unordered pair of actors: (a) its state in round t, where state refers to whether the pair is null, asymmetric, or mutual; (b) its state in round t - 1; (c) a dummy variable for whether the pair is stable, i.e., maintains its state from round t - 1 to t

sex2 <- sex |>
  mutate(pair = paste(pmin(subject, alter), pmax(subject, alter), sep = " - "),
         pair = paste(sessionid, pair, sep = ":")) |>
  group_by(sessionid, session_date, split, visible, pair_type, pair, round) |>
  summarise(state = sum(gave), 
            .groups = "drop_last") |>
  mutate(is_stable = lag(state) == state,
         state_last = lag(state)) |>
  ungroup()

sex2 <- sex2 |>
  mutate(state = factor(state, levels = 0:2, labels = c("N", "A", "M")),
         state_last = factor(state_last, levels = 0:2, labels = c("N", "A", "M")))

# Splitting the dataset by experimental condition (i.e., whether social differences are visible) and fitting the logistic regression model (10) (see page 13 in the paper) for each condition separately

sex2 <- sex2 |>
  group_by(visible) |>
  nest() |>
  mutate(m1 = map(data, ~glm(is_stable ~ state_last * pair_type, data = ., family = binomial())))

## A function shuffle() which will be used to create permuted datasets. It returns returns the original dataset `sex` with permuted rows. The permutation here means that (a) actors are randomly reassigned to the Klee or Kandinsky groups and (b) rounds are randomly reordered

shuffle <- function(data) {
  
  data <- data |>
    group_by(sessionid, visible, split) |>
    nest()
  
  data <- data |>
    mutate(tab_1 = map(data, ~distinct(., subject, s_cat)),
           tab_1 = map(tab_1, ~mutate(.data = ., subject = sample(x = 1:10, size = 10))),
           tab_1 = map(tab_1, ~arrange(.data = ., subject)),
           tab_1 = map(tab_1, ~filter(.data = ., s_cat == 1)))
  
  data <- data |>
    mutate(data = map2(.x = data, .y = tab_1, .f = ~mutate(
      .data = .x, 
      s_cat_perm = ifelse(subject %in% .y$subject, 1, 2),
      a_cat_perm = ifelse(alter %in% .y$subject, 1, 2),
      pair_type_perm = ifelse(s_cat_perm == a_cat_perm, "Intra-group", "Inter-group"))))
  
  data <- data |>
    mutate(data = map(data,
                      ~group_by(., subject, alter) |> 
                        mutate(round_perm = sample(x = 1:20, size = 20)) |>
                        arrange(subject, alter, round_perm) |>
                        mutate(gave_last_perm = lag(gave),
                               received_last_perm = lag(received)) |>
                        ungroup()))
  
  data <- data |>
    dplyr::select(-tab_1) |>
    unnest(cols = data) |>
    ungroup()
  
  data
  
}

# A list to store permuted datasets
s2_perm <- vector(mode = "list", length = 1000L)

# The permutaed datasets
set.seed(19781011)
for (i in 1:1000) {
  s2_perm[[i]] <- shuffle(data = sex)
  print(i)
}

# Converting the list into a tibble with a list-column
s2_perm <- tibble(data = s2_perm) 

# Transforming each permuted dataset into a new dataset where each pair of actors (i, j) in round t classified as null, asymmetric, or mutual in a similat way as the original dataset earlier on
s2_perm <- s2_perm |>
  mutate(d2 = map(data, ~mutate(.data = ., 
                                pair = paste(pmin(subject, alter), pmax(subject, alter), sep = " - "),
                                pair = paste(sessionid, pair, sep = ":"))))

s2_perm <- s2_perm |>
  mutate(d2 = map(d2, ~{
    group_by(.data = ., sessionid, session_date, split, visible, pair, pair_type_perm, round_perm) |>
      summarise(state = sum(gave), .groups = "drop_last") |>
      mutate(is_stable = lag(state) == state) |>
      ungroup()}))

s2_perm <- s2_perm |>
  mutate(d2 = map(d2, ~mutate(.data = ., state = factor(state, levels = 0:2, labels = c("N", "A", "M")))))

# Fitting logistic regression (10) (see page 13 in the paper) to each permuted dataset, separately for the visible and hidden differences conditions

s2_perm <- s2_perm |>
  mutate(m1 = map(d2, ~glm(is_stable ~ lag(state) * pair_type_perm, 
                           data = ., subset = visible == "Yes" & round_perm > 1,
                           family = binomial(link = "logit"))),
         m2 = map(d2, ~glm(is_stable ~ lag(state) * pair_type_perm, 
                           data = ., subset = visible == "No" & round_perm > 1,
                           family = binomial(link = "logit"))))

# Extracting regression coefficients from the models
s2_perm <- s2_perm |>
  mutate_at(.vars = vars(m1:m2), .funs = list(beta = ~map(., coef)))

# A tibble with regression coefficients from visible social differences

tab_1 <- s2_perm$m1_beta |>
  bind_rows()

# A tibble with regression coefficients from hidden social differences

tab_2 <- s2_perm$m2_beta |>
  bind_rows()

# Binding the two tibbles together and renaming the columns

t2 <- bind_rows(tab_1, tab_2, .id = "treatment") |>
  mutate(treatment = factor(treatment, levels = 1:2, labels = c("Visible", "Hidden"))) |>
  rename("Intercept" = 2, "A" = 3, "M" = 4, "W" = 5, "AW" = 6, "MW" = 7)

# P-values for the visible condition
pv_vis <- t2 |>
  filter(treatment == "Visible") |>
  summarise(p1 = mean(W > abs(coef(sex2$m1[[1]])[4])),
            p2 = mean(W + AW > sum(coef(sex2$m1[[1]])[4:5])),
            p3 = mean(W + MW > sum(coef(sex2$m1[[1]])[c(4, 6)])))

# P-values for the hidden condition
pv_hid <- t2 |>
  filter(treatment != "Visible") |>
  summarise(p1 = mean(W > abs(coef(sex2$m1[[2]])[4])),
            p2 = mean(W + AW > sum(coef(sex2$m1[[2]])[4:5])),
            p3 = mean(W + MW > sum(coef(sex2$m1[[2]])[c(4, 6)])))

# Binding the two tables together
out <- bind_rows(pv_vis, pv_hid, .id = "treatment") |>
  pivot_longer(cols = -1, names_to = c("Hypothesis"), values_to = c("$p$-value"),
               names_transform = list(Hypothesis = ~str_replace(string = ., pattern = "p", replacement = "H")))

# Estimates of the hypothesised effects
est <- sex2 |>
  mutate(m1_coef = map(m1, coef)) |>
  unnest(cols = m1_coef) |>
  select(1, 4) |>
  summarise(H1 = m1_coef[4],
            H2 = sum(m1_coef[4:5]),
            H3 = sum(m1_coef[c(4, 6)])) |>
  pivot_longer(cols = -1, names_to = "Hypothesis", values_to = "Estimate") |>
  mutate(visible = factor(visible, levels = c("Yes", "No"), labels = c("1", "2"))) |>
  rename(treatment = visible)

# Building a table with final results
out <- left_join(est, out) |>
  arrange(Hypothesis, treatment) |>
  mutate(Prediction = rep(c("$a_3 < 0$", "$a_3 + a_4 > 0$ since $c > 0$", "$a_3 + a_5 > 0$"), each = 2),
         treatment = factor(treatment, levels = 1:2, labels = c("Visible", "Hidden"))) |>
  select(2, 5, 1, 3, 4)

out |>
  kable(dec = c(NA, NA, NA, 3, 3), caption = "Estimates of the hypothesised effects and their $p$-values") |>
  kable_styling() |>
  collapse_rows(columns = 1:2, valign = "top")