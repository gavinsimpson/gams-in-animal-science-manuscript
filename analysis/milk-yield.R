# Data are from Henderson and McCulloch (1990) Transform or Link?
# https://ecommons.cornell.edu/server/api/core/bitstreams/285223da-0e1d-44c9-bb60-2122e3135b55/content
lactation <- data.frame(
  yield = c(0.31, 0.39, 0.50, 0.58, 0.59, 0.64, 0.68, 0.66,
            0.67, 0.70, 0.72, 0.68, 0.65, 0.64, 0.57, 0.48,
            0.46, 0.45, 0.31, 0.33, 0.36, 0.30, 0.26, 0.34,
            0.29, 0.31, 0.29, 0.20, 0.15, 0.18, 0.11, 0.07,
            0.06, 0.01, 0.01),
  week = seq_len(35)
)

# Packages
library("nlstools")
library("mgcv")
library("gratia")
library("ggplot2")

# a reference GLM - not a good fit
m_glm <- glm(
  yield ~ log(week) + week,
  data = lactation,
  family = Gamma("log")
)

glm_fit <- geom_smooth(
  data = lactation,
  aes(
    x = week,
    y = yield
  ),
  method = "glm",
  se = FALSE,
  colour = "#e69f00",
  linewidth = 1.5,
  formula = y ~ log(x) + x,
  method.args = list(family = Gamma("log"))
)

# reference NLS model
m_nls <- nls(
  yield ~ a * week^b * exp(c * week),
  data = lactation,
  start = list(a = 1, b = 1, c = .01)
)

set.seed(1)
m_nls_boot <- nlsBoot(m_nls, niter = 999)

# Gamma model
m_gamma <- gam(
  yield ~ s(week),
  data = lactation,
  method = "REML",
  family = Gamma("log")
)

m_gamma |>
  draw()

m_gamma |>
  conditional_values(
    condition = "week"
  ) |>
  draw() +
  geom_point(
    data = lactation,
    aes(
      x = week,
      y = yield
    )
  ) +
  glm_fit

summary(m_gamma)

m_gamma |>
  appraise(method = "simulate")

# Tweedie model
m_tw <- update(
  m_gamma,
  family = tw()
)

m_tw |>
  draw()

m_tw |>
  conditional_values(
    condition = "week"
  ) |>
  draw() +
  geom_point(
    data = lactation,
    aes(
      x = week,
      y = yield
    )
  ) +
  glm_fit

summary(m_gamma)

AIC(m_gamma, m_tw)

compare_smooths(
  m_gamma, m_tw
) |>
  draw()

m_tw |>
  appraise(method = "simulate")

# Gamma location-scale fit
m_gammals <- gam(
  list(
    yield ~ s(week),
    ~ s(week)
  ),
  data = lactation,
  method = "REML",
  family = gammals()
)

m_gammals |>
  draw()

summary(m_gammals)

# Tweedie location-scale fit
m_twlss <- gam(
  list(
    yield ~ s(week),
    ~ s(week),
    ~ s(week)
  ),
  data = lactation,
  method = "REML",
  family = twlss()
)

m_twlss |>
  draw()

summary(m_twlss)

AIC(m_gammals, m_twlss)
