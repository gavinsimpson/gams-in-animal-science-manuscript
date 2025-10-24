# Analyse pig growth data

# packages
#library("tictoc")
#library("ggforce")
library("readr")
#library("lubridate")
#library("janitor")
library("mgcv")
library("dplyr")
library("ggplot2")
library("mgcv")
library("gratia")
library("ggdist")
library("patchwork")

pw <- read_csv("data/pig-weight-data.csv") |>
  mutate(
    animal = factor(animal)
  )

pw |>
ggplot(
  aes(
    x = date,
    y = weight_estimate,
    group = animal
  )
) +
geom_line()


ctrl <- gam.control(nthreads = 8, trace = TRUE)

m1 <- gam(
  weight_estimate ~
    animal +
    s(day_of_year, by = animal),
  family = Gamma(link = "log"),
  data = pw,
  method = "REML",
  control = ctrl
)

m2 <- gam(
  weight_estimate ~
    s(day_of_year) +
    s(day_of_year, animal, bs = "sz"),
  family = Gamma(link = "log"),
  data = pw,
  method = "REML",
  control = ctrl
)

m3 <- gam(
  weight_estimate ~
    s(day_of_year, animal, bs = "fs"),
  family = Gamma(link = "log"),
  data = pw,
  method = "REML",
  control = ctrl
)

m4 <- gam(
  weight_estimate ~
    s(day_of_year) +
    s(day_of_year, animal, bs = "fs"),
  family = Gamma(link = "log"),
  data = pw,
  method = "REML",
  control = ctrl
)

AIC(m1, m2, m3, m4)

appraise(m2, method = "simulate")

ds <- pw |>
  select(animal, date) |>
  data_slice(
    animal = evenly(animal),
    date = evenly(date, by = 1) |> as.Date()
  ) |>
  mutate(
    day_of_year = yday(date)
  )

fv_m2 <- m2 |>
  fitted_values(
    data = ds
  )

fv_m2 |>
  ggplot(
    aes(
      x = date,
      y = .fitted,
      group = animal
    )
  ) +
  geom_point(
    data = pw,
    aes(x = date, y = weight_estimate),
    size = 0.8,
    colour = "#56B4E9"
  ) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),
    alpha = 0.2
  ) +
  geom_line() +
  facet_wrap(~ animal) +
  labs(
    x = NULL,
    y = "Weight (kg)"
  )

fv_m2 |>
  ggplot(
    aes(
      x = date,
      y = .fitted,
      group = animal
    )
  ) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),
    alpha = 0.1
  ) +
  geom_line() +
  labs(
    x = NULL,
    y = "Weight (kg)"
  )

# estimates of growth rate on November 15th, 2021
nov15_ds <- m2 |>
  data_slice(
    day_of_year = format(as.Date("2021-11-15"), "%j") |> as.numeric(),
    animal = evenly(animal)
  )

nov15_deriv <- m2 |>
  response_derivatives(
    data = nov15_ds,
    n_sim = 10000,
    eps = 1,
    type = "central",
    focal = "day_of_year"
  )

nov15_plt1 <- nov15_deriv |>
  ggplot(
    aes(x = animal, y = .derivative)
  ) +
  geom_pointrange(
    aes(ymin = .lower_ci, ymax = .upper_ci)
  ) +
  labs(
    x = "Pig", y = expression(Growth ~ rate ~ (day^{-1}))
  )

pigs_deriv <- m2 |>
  derivative_samples(
    data = nov15_ds |> filter(animal %in% c("2", "13", "17")),
    n_sim = 10000,
    eps = 1,
    type = "central",
    focal = "day_of_year"
  )

nov15_plt2 <- pigs_deriv |>
  ggplot(
    aes(x = animal, y = .derivative)
  ) +
  stat_halfeye() +
  labs(
    x = "Pig", y = expression(Growth ~ rate ~ (day^{-1}))
  )

nov15_plt1 + nov15_plt2 +
  plot_layout(ncol = 2, widths = c(0.5, 0.5)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")
