pkgs <- c(
  "here",
  "mgcv",
  "gratia",
  "readxl",
  "scales",
  "ggplot2",
  "forcats",
  "dplyr",
  "colorspace"
)
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE)

# Growth data
# data are in sheet "growth analysis"

fn <- "data/Sarraude et al. - THs elevation in Japanese quails - dataset.xlsx"
excel_sheets(fn)

q_growth <- read_xlsx(fn, sheet = "growth analysis", na = "NA") |>
  mutate(
    mother_id = factor(motherID),
    egg_id = factor(eggID),
    group = factor(group),
    sex = factor(sex)
  ) |>
  # make the label prettier for graphs
  mutate(
    sex = fct_recode(
      sex,
      "Female" = "F",
      "Male" = "M"
    ),
    group = fct_recode(
      group,
      "T[3]" = "T3",
      "T[4]" = "T4",
      "T[3]~T[4]" = "T3T4"
    )
  )

# plotting code from
ggplot(q_growth, aes(day, mass, color = group)) +
  geom_line(aes(group = eggID)) +
  facet_grid(sex ~ group, labeller = label_parsed) +
  labs(x = "Days since hatching", y = "Body mass (g)", color = "Treatment") +
  scale_color_okabe_ito(
    labels = str2expression,
    limits = c("CO", "T[3]", "T[4]", "T[3]~T[4]"),
    order = c(1, 2, 3, 7)
  ) +
  theme(legend.text = element_text(hjust = 0))

ggplot(q_growth, aes(day, mass, color = group)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom = "pointrange", fatten = 2) +
  facet_wrap(~sex) +
  labs(
    x = "Days since hatching",
    y = "Body mass (g)",
    color = "Treatment",
    tag = "B"
  ) +
  scale_color_okabe_ito(
    labels = str2expression,
    limits = c("CO", "T[3]", "T[4]", "T[3]~T[4]"),
    order = c(1, 2, 3, 7)
  ) +
  theme(legend.text.align = 0)

# model growth curves
ctrl <- gam.control(nthreads = 8, trace = TRUE)

tic()
quail_m1 <- gam(
  mass ~
    s(day) +
      s(day, sex, bs = "sz") +
      s(day, egg_id, bs = "fs", k = 6) +
      s(mother_id, bs = "re"),
  data = q_growth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl
)
toc()

tic()
quail_m2 <- gam(
  mass ~
    s(day, k = 10) +
      s(day, group, bs = "sz") +
      s(day, sex, bs = "sz") +
      s(day, egg_id, bs = "fs", k = 6) +
      s(mother_id, bs = "re"),
  data = q_growth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl
)
toc()

tic()
quail_m3 <- gam(
  mass ~
    s(day, k = 10) +
      s(day, group, bs = "sz") +
      s(day, sex, bs = "sz") +
      s(day, group, sex, bs = "sz") +
      s(day, egg_id, bs = "fs", k = 6) +
      s(mother_id, bs = "re"),
  data = q_growth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl
)
toc()

tic()
quail_m4 <- gam(
  mass ~
    s(day, k = 10) +
      s(day, group, bs = "sz") +
      s(day, sex, bs = "sz") +
      s(day, group, sex, bs = "sz") +
      s(egg_id, bs = "re") +
      s(mother_id, bs = "re"),
  data = q_growth,
  family = Gamma(link = "log"),
  method = "REML",
  control = ctrl
)
toc()

quail_m1_ml <- update(quail_m1, method = "ML")
quail_m2_ml <- update(quail_m2, method = "ML")
quail_m3_ml <- update(quail_m3, method = "ML")
quail_m4_ml <- update(quail_m4, method = "ML")

AIC(quail_m1_ml, quail_m2_ml, quail_m3_ml, quail_m4_ml)

# bam() fits
tic()
quail_b1 <- bam(
  mass ~
    s(day) +
      s(day, sex, bs = "sz") +
      s(day, egg_id, bs = "fs", k = 6) +
      s(mother_id, bs = "re"),
  data = q_growth,
  family = Gamma(link = "log"),
  method = "fREML",
  control = ctrl,
  discrete = TRUE,
  nthreads = 8
)
toc()

tic()
quail_b2 <- bam(
  mass ~
    s(day, k = 10) +
      s(day, group, bs = "sz") +
      s(day, sex, bs = "sz") +
      s(day, egg_id, bs = "fs", k = 6) +
      s(mother_id, bs = "re"),
  data = q_growth,
  family = Gamma(link = "log"),
  method = "fREML",
  control = ctrl,
  discrete = TRUE,
  nthreads = 8
)
toc()

tic()
quail_b3 <- bam(
  mass ~
    s(day, k = 10) +
      s(day, group, bs = "sz") +
      s(day, sex, bs = "sz") +
      s(day, group, sex, bs = "sz") +
      s(day, egg_id, bs = "fs", k = 6) +
      s(mother_id, bs = "re"),
  data = q_growth,
  family = Gamma(link = "log"),
  method = "fREML",
  control = ctrl,
  discrete = TRUE,
  nthreads = 8
)
toc()

tic()
quail_b4 <- bam(
  mass ~
    s(day, k = 10) +
      s(day, group, bs = "sz") +
      s(day, sex, bs = "sz") +
      s(day, group, sex, bs = "sz") +
      s(egg_id, bs = "re") +
      s(mother_id, bs = "re"),
  data = q_growth,
  family = Gamma(link = "log"),
  method = "fREML",
  control = ctrl,
  discrete = TRUE,
  nthreads = 8
)
toc()

quail_b3 |>
  smooths()

quail_b3 |>
  conditional_values(
    condition = c("day", "group", "sex"),
    exclude = c("s(day,egg_id)", "s(mother_id)")
  ) |>
  draw()

quail_ds <- quail_b3 |>
  data_slice(
    day = evenly(day),
    group = evenly(group),
    sex = evenly(sex)
  )

quail_tw1 <- update(quail_b1, family = tw())
quail_tw2 <- update(quail_b2, family = tw())
quail_tw3 <- update(quail_b3, family = tw())
quail_tw4 <- update(quail_b4, family = tw())

AIC(
  quail_b1,
  quail_b2,
  quail_b3,
  quail_b4,
  quail_tw1,
  quail_tw2,
  quail_tw3,
  quail_tw4
)

quail_tw3 |>
  conditional_values(
    condition = c("day", "group", "sex"),
    exclude = c("s(day,egg_id)", "s(mother_id)")
  ) |>
  draw()

quail_tw2 |>
  conditional_values(
    condition = c("day", "group", "sex"),
    exclude = c("s(day,egg_id)", "s(mother_id)")
  ) |>
  ggplot(
    aes(
      x = day,
      y = .fitted,
      colour = group,
      linetype = sex,
      group = interaction(group, sex)
    )
  ) +
  geom_line() +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci, fill = group, colour = NULL),
    alpha = 0.2
  ) +
  facet_wrap(~group, labeller = label_parsed) +
  scale_color_okabe_ito(
    labels = str2expression,
    limits = c("CO", "T[3]", "T[4]", "T[3]~T[4]"),
    order = c(1, 2, 3, 7)
  ) +
  scale_fill_okabe_ito(
    labels = str2expression,
    limits = c("CO", "T[3]", "T[4]", "T[3]~T[4]"),
    order = c(1, 2, 3, 7)
  ) +
  labs(
    x = "Days since hatching",
    y = "Body mass (g)",
    color = "Treatment",
    fill = "Treatment",
    linetype = "Sex",
  ) +
  theme(legend.text = element_text(hjust = 0))

quail_lookup <- q_growth |>
  distinct(egg_id, mother_id, sex, group)

quail_ds <- quail_tw2 |>
  data_slice(
    day = evenly(day),
    egg_id = evenly(egg_id)
  ) |>
  select(-mother_id, -sex, -group) |>
  left_join(
    quail_lookup,
    by = join_by(egg_id == egg_id)
  )

quail_fv <- quail_tw2 |>
  fitted_values(
    data = quail_ds
  )

quail_fv |>
  ggplot(
    aes(
      x = day,
      y = .fitted,
      colour = group,
      group = egg_id
    )
  ) +
  geom_point(
    data = q_growth,
    aes(
      x = day,
      y = mass,
      colour = group
    ),
    size = 1
  ) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci, fill = group, colour = NULL),
    alpha = 0.2
  ) +
  geom_line() +
  facet_grid(sex ~ group, labeller = label_parsed) +
  scale_color_okabe_ito(
    labels = str2expression,
    limits = c("CO", "T[3]", "T[4]", "T[3]~T[4]"),
    order = c(1, 2, 3, 7)
  ) +
  scale_fill_okabe_ito(
    labels = str2expression,
    limits = c("CO", "T[3]", "T[4]", "T[3]~T[4]"),
    order = c(1, 2, 3, 7)
  ) +
  labs(
    x = "Days since hatching",
    y = "Body mass (g)",
    color = "Treatment",
    fill = "Treatment",
  ) +
  theme(legend.text = element_text(hjust = 0))
