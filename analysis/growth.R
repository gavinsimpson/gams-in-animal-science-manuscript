# Analyse pig growth data from Mona

# 1. be sure to cite the two papers
# 2. only a small subset can be made open, so use 1 pen's worth of pigs

# packages
library("tictoc")
library("ggforce")
library("readr")
library("lubridate")
library("janitor")
library("mgcv")
library("dplyr")
library("ggplot2")
library("mgcv")
library("gratia")

# note data are in a Danish CSV format
base_fn <- "/Users/gavin/work/data/animal-science/"
fn <- paste0(base_fn, "pig-growth-automated-camera-mona/aggregatedWeights-2.csv")

pw_raw <- read_csv2(
    fn,
    col_types = "cDcdc",
    locale = locale()
  ) |>
  janitor::clean_names() |>
  mutate(
    weight_estimate = as.numeric(weight_estimate)
  )

pw_raw |>
ggplot(
  aes(
    x = date,
    y = weight_estimate,
    group = id
  )
) +
geom_line() +
facet_wrap(~ pen)

# pen 803

## lookup for id
pw_lookup <- pw_raw |>
  filter(pen == "803") |>
  distinct(id) |>
  mutate(
    animal = row_number()
  )

pw_data <- pw_raw |>
  filter(pen == "803") |>
  mutate(
    day_of_year = yday(date)
  ) |>
  left_join(
    pw_lookup, by = join_by(id == id)
  ) |>
  select(-id) |>
  mutate(
    animal = factor(animal)
  )

write_csv(pw_data |> select(-pen), file = "data/pig-weight-data.csv")

read_csv("data/pig-weight-data.csv")

pw_data |>
ggplot(
  aes(
    x = date,
    y = weight_estimate,
    group = animal
  )
) +
geom_line()


ctrl <- gam.control(nthreads = 6, trace = TRUE)

m1 <- gam(
  weight_estimate ~
    animal +
    s(day_of_year, by = animal),
    family = Gamma(link = "log"),
    data = pw_data,
    method = "REML",
    control = ctrl
  )

m2 <- gam(
  weight_estimate ~
    s(day_of_year) +
    s(day_of_year, animal, bs = "sz"),
    family = Gamma(link = "log"),
    data = pw_data,
    method = "REML",
    control = ctrl
)

m3 <- gam(
  weight_estimate ~
    s(day_of_year, animal, bs = "fs"),
    family = Gamma(link = "log"),
    data = pw_data,
    method = "REML",
    control = ctrl
)

m4 <- gam(
  weight_estimate ~
    s(day_of_year) +
    s(day_of_year, animal, bs = "fs"),
    family = Gamma(link = "log"),
    data = pw_data,
    method = "REML",
    control = ctrl
)

AIC(m1, m2, m3, m4)

appraise(m2, method = "simulate")

ds <- pw_data |>
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
    data = pw_data,
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


## location scale model
pw_data <- pw_data |>
  mutate(
    rt_rt_weight_count = weight_count^(1/4)
  )

tic()
m2_lss <- gam(
  list(
    weight_estimate ~
      s(day_of_year) +
      s(day_of_year, animal, bs = "sz"),
    ~ s(animal, bs = "re") +
      s(rt_rt_weight_count)
  ),
  family = gammals(),
  data = pw_data,
  method = "REML",
  control = ctrl
)
toc()

summary(m2_lss)
AIC(m2, m2_lss)
draw(m2_lss)

ds_lss <- pw_data |>
  select(animal, date, weight_count) |>
  data_slice(
    animal = evenly(animal),
    date = evenly(date, by = 1) |> as.Date(),
    weight_count = median(weight_count)
  ) |>
  mutate(
    day_of_year = yday(date),
    rt_rt_weight_count = weight_count^(1/4),
    .row = row_number()
  )

fv_m2_lss <- m2_lss |>
  fitted_values(
    data = ds_lss
  ) |>
  filter(.parameter == "location") |>
  left_join(
    ds_lss,
    by = join_by(.row == .row)
  )

fv_m2_lss |>
  ggplot(
    aes(
      x = date,
      y = .fitted,
      group = animal
    )
  ) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),
    alpha = 0.2
  ) +
  geom_line() +
  facet_wrap(~ animal) +
  geom_point(
    data = pw_data,
    aes(x = date, y = weight_estimate)
  )
