#' Convert Gaze Coordinates
#'
#' Converts raw gaze coordinates from pixels into degrees of visual angle.
#'
#' @param x Numeric vector of raw coordinates for each gaze sample.
#' @param res Screen resolution (in px).
#' @param dim Screen dimensions (in mm).
#' @param dist Distance between subject and screen (in mm).
#'
#' @return A numeric vector with coordiantes in degrees of visual angle centered at (0,0) for each gaze sample.
#' @export
px2va <- function(x, res, dim, dist) {

  x_centered <- x - res/2 # transform to centered coordinate system (0,0)

  radian <- atan(x_centered/(dist*(res/dim))) # visual angle in radians

  degree <- radian*(180/pi) # visual angle in degrees

  return(degree)
}


# Function to compute the relative angle between subsequent samples

#' Calculate Sample-to-sample Angle
#'
#' Calculates the relative angle between subsequent gaze samples.
#' Uses the backward difference to compute the absolute angle between two samples and
#' the forward difference for the change in angle at each sample.
#'
#' @param x Numeric vector of x coordinates for each gaze sample.
#' @param y Numeric vector of y coordinates for each gaze sample.
#'
#' @return A vector containing the sample-to-sample angle (in radians) for each gaze sample.
#' @importFrom dplyr lead lag
#' @export
calc_theta <- function(x, y) { # x-pos, y-pos

  diff_x <- x - lag(x) # x_t - x_t-1
  diff_y <- y - lag(y) # y_t - y_t-1

  angle <- atan2(diff_y, diff_x) # absolute angle of vector xy_t - xy_t-1
  theta <- lead(angle) - angle # relative angle of two vectors at sample t
  theta.mirrored <- ifelse(theta < 0, theta + 2*pi, theta) # mirror negative relative angles to positive range

  return(theta.mirrored)
}


# Function to compute a windowed difference for the x and y signal

#' Calculate Windowed Diff
#'
#' Calculates the windowed difference for x and y.
#'
#' @param x Numeric vector of x coordinates for each gaze sample.
#' @param y Numeric vector of y coordinates for each gaze sample.
#' @param window_size Window size in millisecond to calculate the direction feature.
#' @param sampling_rate Sampling rate of the eye-tracker (in Hz).
#'
#' @return A list containing the x and y difference vectors.
#' @importFrom dplyr lead lag
#' @export
calc_windowed_diff <- function(x, y, window_size, sampling_rate=1000) {
  # Adjust window size to sampling rate
  window_size <- as.integer(window_size * sampling_rate / 1000)
  l = window_size - 1
  # Calculate direction
  x_diff <- dplyr::lead(x, l) - x
  y_diff <- dplyr::lead(y, l) - y
  for (i in seq(0, l-1)) {
    j <- length(x)-i
    x_diff[j] <- dplyr::lead(x[j:length(x)], i) - x[j]
    y_diff[j] <- dplyr::lead(y[j:length(y)], i) - y[j]
  }
  return(list(x_diff, y_diff))
}


# Function to compute the direction feature according to Startsev et al. (2018)

#' Calculate Direction Feature
#'
#' Calculates the direction feature.
#'
#' @param x Numeric vector of x coordinates for each gaze sample.
#' @param y Numeric vector of y coordinates for each gaze sample.
#' @param window_size Window size in millisecond to calculate the direction feature.
#' @param sampling_rate Sampling rate of the eye-tracker (in Hz).
#'
#' @return A vector containing direction feature for each sample.
#' @importFrom dplyr lead lag
#' @export
calc_direction <- function(x, y, window_size=16, sampling_rate=1000) {
  # Calcualtes the direction feature according to Startsev et al. (2018)
  diffs <- calc_windowed_diff(x, y, window_size=window_size, sampling_rate=sampling_rate)
  x_diff <- diffs[[1]]
  y_diff <- diffs[[2]]  
  direction <- atan2(y_diff, x_diff)
  return(direction)
}


# Function to compute the direction deviation

#' Calculate Direction Deviation
#'
#' Calculates the direction deviation.
#'
#' @param x Numeric vector of x coordinates for each gaze sample.
#' @param y Numeric vector of y coordinates for each gaze sample.
#' @param t Numeric vector of time for each gaze sample.
#' @param window_size Window size in millisecond to calculate the direction feature.
#' @param outer_window_size Window size in millisecond to calculate the outer direction.
#' @param outer_step_size Step size in millisecond to calculate the outer direction.
#' @param sampling_rate Sampling rate of the eye-tracker (in Hz).
#'
#' @return A vector containing the standard deviation of the direction deviation for each sample (NA for first t < outer_step_size) .
#' @importFrom dplyr lead lag
#' @export
calc_direction_deviation <- function(x, y, t, window_size=16, outer_window_size=200, outer_step_size=10, sampling_rate=1000) {
  # Calculates the std of the deviation between the startsev direction and the direction of
  # overlapping outer windows

  # Calculate direction
  direction <- calc_direction(x, y, window_size=window_size, sampling_rate=sampling_rate)
  # Adjust window size to sampling rate
  outer_window_size <- as.integer(outer_window_size * sampling_rate / 1000)
  outer_step_size <- as.integer(outer_step_size * sampling_rate / 1000)
  l = outer_window_size - 1
  # Calculate outer window direction
  outer_window_direction <- calc_direction(x, y, window_size=outer_window_size, sampling_rate=sampling_rate) 
  # Calculate direction deviation
  direction_deviation <- c()
  time = c()
  for (i in seq(1, length(direction), outer_step_size)) {
    direction_deviation <- c(direction_deviation, outer_window_direction[i] - direction[i:(i+l)])
    time <- c(time, t[i:(i+l)])
  }
  # calculate standard deviation of direction deviation per time point
  temp <- aggregate(direction_deviation ~ time, data.frame(time=time, direction_deviation=direction_deviation), sd)

  direction_deviation <- temp$direction_deviation
  return(direction_deviation)
}
