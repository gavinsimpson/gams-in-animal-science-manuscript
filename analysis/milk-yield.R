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
library("dplyr")
library("patchwork")

# reference NLS model
m_nls <- nls(
  yield ~ a * week^b * exp(c * week),
  data = lactation,
  start = list(a = 1, b = 1, c = .01)
)

# GLM fit

m_glm_gamma <- gam(
  yield ~ log(week) + week,
  data = lactation,
  family = tw(link = "log"),
  method = "ML"
)

m_glm_tw <- gam(
  yield ~ log(week) + week,
  data = lactation,
  family = tw(link = "log"),
  method = "ML"
)

# we want a comparable Tweedie GLM to compare with the Tweedie GAM
m_glm <- m_glm_tw

# bootstrap the NLS model for uncertainty estimates
set.seed(1)
n_boot <- 10000
m_nls_boot <- nlsBoot(m_nls, niter = n_boot)

# Gamma GAM
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
  )#  +
  # glm_fit

summary(m_gamma)

m_gamma |>
  appraise(method = "simulate")

# Tweedie GAM
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
  )#  +
  # glm_fit

m_gamma |> summary()
m_tw |> summary()

AIC(m_gamma, m_tw)

compare_smooths(
  m_gamma, m_tw
) |>
  draw()

m_tw |>
  appraise(method = "simulate")

# Now let's compare some values of interest
#
# * the day of peak yield,
# * the yield at the peak,
# * a rate of decline from the peak to the end of lactation
#
# For Wood's model, these are
#
# * b / c
# * a(b/c)^b * e^-b, and
# * a derivative of the curve

# use the bootstrap stats to get bootstrap distribution of each of these
coefs <- m_nls_boot$coefboot |> data.frame()
# note I fitted exp(c*week), so we need to negate c here
nls_peak_boot <- with(coefs, b / -c)
# max yield at peak week
nls_max_yield <- with(coefs, a * median(nls_peak_boot)^b * exp(-b))

# find peak for the Tweedie GLM - as I fitted c*week, need to negate c
glm_peak <- coef(m_glm)[2] / -coef(m_glm)[3]

# now find it for the glm/gam
n_sim <- 10000 # how many posterior draws
n_grid <- 1000 # how many points on lactation curve shall we estimate at

# create data for prediction
ds <- m_tw |>
  data_slice(
    week = evenly(week, n_grid)
  ) |>
  mutate(
    .row = row_number()
  )

# posterior expectations of tweedie gam
fs_gam <- m_tw |>
  fitted_samples(
    data = ds,
    n = n_sim,
    seed = 2
  ) |> # add on a row identifier for linking with data
  left_join(
    ds, by = join_by(.row)
  )

# posterior expectations of tweedie glm
fs_glm <- m_glm |>
  fitted_samples(
    data = ds,
    n = n_sim,
    seed = 2
  ) |> # add on a row identifier for linking with data
  left_join(
    ds, by = join_by(.row)
  )

# find the peak of each posterior draw
peak_post_gam <- fs_gam |>
  group_by(.draw) |>
  slice_max(.fitted)

peak_post_glm <- fs_glm |>
  group_by(.draw) |>
  slice_max(.fitted)

# create a data frame of results
peak_week_tbl <- data.frame(
  model = factor(c("Wood", "GLM", "GAM"), levels = c("Wood", "GLM", "GAM")),
  .estimate = c(
    nls_peak_boot |> quantile(prob = 0.5),
    peak_post_glm |> pull(week) |> quantile(prob = 0.5),
    peak_post_gam |> pull(week) |> quantile(prob = 0.5)
  ),
  .lower = c(
    nls_peak_boot |> quantile(prob = 0.025),
    peak_post_glm |> pull(week) |> quantile(prob = 0.025),
    peak_post_gam |> pull(week) |> quantile(prob = 0.025)
  ),
  .upper = c(
    nls_peak_boot |> quantile(prob = 0.975),
    peak_post_glm |> pull(week) |> quantile(prob = 0.975),
    peak_post_gam |> pull(week) |> quantile(prob = 0.975)
  )
)

# plot these results, but store plot for later
fat_comp1 <- peak_week_tbl |>
  ggplot(
    aes(
      x = model,
      y = .estimate,
    )
  ) +
  geom_pointrange(
    aes(
      ymin = .lower, ymax = .upper
    )
  ) +
  labs(
    x = NULL,
    y = "Week of peak fat content"
  )

# generate predicted values the NLS model for each of the bootstrap samples
# here, prediction is for the peak
fv_nls <- nlstools::nlsBootPredict(
  m_nls_boot,
  newdata = data.frame(
    week = nls_peak_boot |>
      quantile(prob = 0.5) |>
      rep(2)
  ),
  interval = "confidence"
)

# combine fat content at peak - note we already have the estimate for the peak
# in the posterior draw objects
peak_fat_tbl <- data.frame(
  model = factor(c("Wood", "GLM", "GAM"), levels = c("Wood", "GLM", "GAM")),
  .estimate = c(
    fv_nls[1, 1],
    peak_post_glm |> pull(.fitted) |> quantile(prob = 0.5),
    peak_post_gam |> pull(.fitted) |> quantile(prob = 0.5)
  ),
  .lower = c(
    fv_nls[1, 2],
    peak_post_glm |> pull(.fitted) |> quantile(prob = 0.025),
    peak_post_gam |> pull(.fitted) |> quantile(prob = 0.025)
  ),
  .upper = c(
    fv_nls[1, 3],
    peak_post_glm |> pull(.fitted) |> quantile(prob = 0.975),
    peak_post_gam |> pull(.fitted) |> quantile(prob = 0.975)
  )
)

# ...and plot these data, saving for later
fat_comp2 <- peak_fat_tbl |>
  ggplot(
    aes(
      x = model,
      y = .estimate,
    )
  ) +
  geom_pointrange(
    aes(
      ymin = .lower, ymax = .upper
    )
  ) +
  labs(
    x = NULL,
    y = "Fat content at peak week"
  )

## derivatives...

# Better to just estimate the central finite-difference based first derivative
# at this time point mid way between peak and end of lactation

# nlstools doesn't return the predicted values for each bootstrap sample's
# coefs, so here I hack a function together from nlstools::nlsBootPredict
my_nlsBootPredict <- function(nlsBoot, newdata) {
  nlsformula <- formula(nlsBoot$nls)
  nlsresid <- resid(nlsBoot$nls)
  param <- nlsBoot$coefboot
  bootparam <- nlsBoot$coefboot
  niter <- length(nlsBoot$rse)
  "formula2function" <- function(formu) {
    arg1 <- all.vars(formu)
    arg2 <- vector("list", length(arg1))
    names(arg2) <- arg1
    Args <- do.call("alist", arg2)
    fmodele <- as.function(c(Args, formu))
    return(fmodele)
  }
  f1 <- formula2function(formula(nlsformula)[[3]])
  vardep <- all.vars(nlsformula[[2]])
  varindep <- intersect(all.vars(nlsformula[[3]]), colnames(newdata))
  one.mean.pred <- function(i) {
    do.call(f1, as.list(c(param[i, ], newdata[varindep])))
  }
  sapply(1:niter, one.mean.pred)
}

# function to compute a point in lactation that if halfway between the peak and
# the end of lactation
half_way_pt <- function(peak, tf = 35) {
  peak + ((tf - peak) / 2)
}

# for this then we need predictions at +/- 0.5 days at the mid-way point
# between the peak and the end of lactation

# first, convert 1 week to a half day - h is 1/2 the delta in the finite
# difference
h <- (1 / 7) / 0.5

# compute values at mid point +/- h (half a day)
finite_diff_nls_2 <- my_nlsBootPredict(
  m_nls_boot,
  newdata = data.frame(
    week = half_way_pt(peak_week_tbl[1, 2]) + c(-h, h)
  )
)

# function to compute the central finite difference
central_deriv <- function(x, h) {
  (x[2] - x[1]) / h
}

# compute finite difference 1st deriv; h is doubled here as the total delta is
# 1 day and h if 0.5 days
persistence_nls <- apply(finite_diff_nls_2, 2, central_deriv, h = h * 2)

# compute the derivative for the GLM
persistence_glm <- m_glm |>
  response_derivatives(
    data = data.frame(week = half_way_pt(peak_week_tbl[2, 2])),
    focal = "week",
    order = 1,
    type = "central",
    eps = h,
    seed = 42,
    n_sim = 10000
  )

# compute the derivative for the GAM
persistence_gam <- m_tw |>
  response_derivatives(
    data = data.frame(week = half_way_pt(peak_week_tbl[3, 2])),
    focal = "week",
    order = 1,
    type = "central",
    eps = h,
    seed = 42,
    n_sim = 10000
  )

# combine persistence estimates into a table
# NOTE: response_derivatives is already summarising the posterior distribution
#       so no need for a summary here
persistence_tbl <- data.frame(
  model = factor(c("Wood", "GLM", "GAM"), levels = c("Wood", "GLM", "GAM")),
  .estimate = c(
    persistence_nls |> quantile(prob = 0.5),
    persistence_glm |> pull(.derivative),
    persistence_gam |> pull(.derivative)
  ),
  .lower = c(
    persistence_nls |> quantile(prob = 0.025),
    persistence_glm |> pull(.lower_ci),
    persistence_gam |> pull(.lower_ci)
  ),
  .upper = c(
    persistence_nls |> quantile(prob = 0.975),
    persistence_glm |> pull(.upper_ci),
    persistence_gam |> pull(.upper_ci)
  )
)

fat_comp3 <- persistence_tbl |>
  ggplot(
    aes(
      x = model,
      y = .estimate,
    )
  ) +
  geom_pointrange(
    aes(
      ymin = .lower, ymax = .upper
    )
  ) +
  labs(
    x = NULL,
    y = "Rate of change in fat content"
  )

fat_comp1 + fat_comp2 + fat_comp3 +
  plot_annotation(tag_levels = "a", tag_suffix = ")")
