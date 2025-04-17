#' Find optimal designs for benchmark dose analysis
#'
#' @param grad_fun Dose-response gradient function.
#' @param dr_fun Dose-response function. Used to account for non-constant variance.
#' @param cvec c vector for finding c-optimal designs.
#' @param obj Design objective. Can either be "c" or "D".
#' @param theta Vector of local model parameter values.
#' @param d0 Doses in prior experiment.
#' @param n0 Sample size allocations in prior experiment.
#' @param N1 Sample size allowance for current experiment. Used when finding augmented designs.
#' @param num_doses Number of doses in the design.
#' @param max_dose Maximum allowed dose.
#' @param swarm Swarm size for swarm-based optimization algorithm.
#' @param iter Maximum number of iterations for optimization algorithm.
#' @param alg Algorithm from metaheuristicOpt
#' @param ignore_stage1 If TRUE, finds the locally optimal design.
#' @param show_progress If TRUE, show progress bar for optimization.
#'
#' @return
#' @export
#'
#' @examples
#' # find designs for logistic model ###########################################
#' # parameter estimates obtained using ToxicR on deguelin data from drc package
#' theta_logistic = c(-1.5627368, 0.1258373)
#'
#' # D-optimal design
#' set.seed(1136)
#' logistic_Dopt = BMDdesign(
#'   grad_fun = logistic.grad,
#'   dr_fun = logistic.fun,
#'   obj = 'D',
#'   theta = theta_logistic,
#'   num_doses = 2,
#'   max_dose = 50.11872,
#'   swarm = 20,
#'   iter = 100,
#'   alg = 'DE'
#' )
#' logistic_Dopt$x
#' logistic_Dopt$w
#' plot_design(logistic_Dopt)
#' check_eq(logistic_Dopt)
#'
#' # c-optimal design for BMD at BMR=0.1
#' cvec = logistic.bmdgrad(0.1, theta_logistic)
#' set.seed(1140)
#' logistic_copt = BMDdesign(
#'   grad_fun = logistic.grad,
#'   dr_fun = logistic.fun,
#'   obj = 'c',
#'   cvec = cvec,
#'   theta = theta_logistic,
#'   num_doses = 2,
#'   max_dose = 50.11872,
#'   swarm = 30,
#'   iter = 500,
#'   alg = 'DE'
#' )
#' logistic_copt$x
#' logistic_copt$w
#' plot_design(logistic_copt)
#' check_eq(logistic_copt)
#'
#' # 3-parameter Weibull model #################################################
#' # parameter estimates obtained using ToxicR on deguelin data from drc package
#' theta_weibull = c(-1.145036439, 1.804710767, 0.004269448) # note g is reparameterized
#'
#' # D-optimal design
#' set.seed(1152)
#' weibull_Dopt = BMDdesign(
#'   grad_fun = weibull.grad,
#'   dr_fun = weibull.fun,
#'   obj = 'D',
#'   theta = theta_weibull,
#'   num_doses = 3,
#'   max_dose = 50.11872,
#'   swarm = 20,
#'   iter = 500,
#'   alg = 'DE'
#' )
#' weibull_Dopt$x
#' weibull_Dopt$w
#' plot_design(weibull_Dopt)
#' check_eq(weibull_Dopt)
#'
#' # c-optimal design
#' cvec = weibull.bmdgrad(0.1, theta_weibull)
#' set.seed(1156)
#' weibull_copt = BMDdesign(
#'   grad_fun = weibull.grad,
#'   dr_fun = weibull.fun,
#'   obj = 'c',
#'   cvec = cvec,
#'   theta = theta_weibull,
#'   num_doses = 3,
#'   max_dose = 50.11872,
#'   swarm = 20,
#'   iter = 500,
#'   alg = 'DE'
#' )
#' weibull_copt$x
#' weibull_copt$w
#' plot_design(logistic_copt)
#' check_eq(weibull_copt)
BMDdesign = function(grad_fun, dr_fun, cvec = NULL, obj, theta, d0 = NULL, n0 = NULL,
                        N1=NULL, num_doses, max_dose, swarm, iter, alg,
                        ignore_stage1 = T, show_progress = T) {


  if (!ignore_stage1) {
    # compute information matrix from initial design
    F0 = sapply(d0, grad_fun, theta)
    phat0 = sapply(d0, dr_fun, theta)
    v0 = 1/(phat0 * (1 - phat0))
    M0 = n0[1] * F0[, 1] %*% t(F0[, 1]) * v0[1]
    for (i in 2:length(d0)) {
      M0 = M0 + n0[i] * F0[, i] %*% t(F0[, i]) * v0[i]
    }
  }
  else {
    M0 = 0
  }


  # define objective function
  obj_fun = function(vars, ...) {
    pts = length(vars)/2
    x = vars[1:pts]
    w = vars[(pts+1):(2*pts)]

    w = w/sum(w)

    M1 = 0
    for (i in 1:pts) {
      phat_i = dr_fun(x[i], theta)
      v_i = 1 / (phat_i * (1 - phat_i))
      M1_i = w[i] * v_i * grad_fun(x[i], theta) %*% t(grad_fun(x[i], theta))
      M1 = M1 + M1_i
    }

    # augmented information matrix
    if (ignore_stage1) {
      alpha = 1
    }
    else {
      N0 = sum(n0)
      alpha = N1 / (N0 + N1)
    }

    M = alpha*M1 + (1-alpha)*M0

    if (!checkMinv(M))
      return(Inf)
    else {
      if (obj == 'D') {
        obj_val = suppressWarnings(-log(det(M)))
      }
      else if (obj == 'c') {
        cvec = cvec
        Minv = solve(M)
        obj_val = t(cvec) %*% Minv %*% cvec
      }
      else
        stop('objective not defined')

      if (is.na(obj_val) | is.nan(obj_val))
        return(Inf)
      else
        return(obj_val)
    }
  }

  # set up variable bounds and control variables
  pts = num_doses
  rangeVar =  matrix(c(rep(c(0, max_dose), pts), rep(c(0,1), pts)), nrow = 2)
  control = list(numPopulation = swarm, maxIter = iter)

  # call optimizer
  if (show_progress) {
    result = metaheuristicOpt::metaOpt(
      obj_fun,
      optimType = "MIN",
      algorithm = alg,
      numVar = 2*pts,
      rangeVar,
      control,
      seed = NULL
    )
  }
  else {
    # silenced version for use in simulations
    invisible(capture.output({result <- metaheuristicOpt::metaOpt(
      obj_fun,
      optimType = "MIN",
      algorithm = alg,
      numVar = 2*pts,
      rangeVar,
      control,
      seed = NULL
    )}))
  }


  # extract results and process
  vars = result$result
  x = vars[1:pts]
  w = vars[(pts+1):(2*pts)]
  w = w/sum(w)
  # removing the point collapsing
  #x = x[w > 1e-5]
  #w = w[w > 1e-5]
  w = w[order(x)]
  x = x[order(x)]


  # compute 2nd stage information matrix
  M1 = 0
  for (i in 1:length(x)) {
    phat2_i = dr_fun(x[i], theta)
    v2_i = 1 / (phat2_i * (1-phat2_i))
    M1_i =  w[i] * grad_fun(x[i], theta) %*% t(grad_fun(x[i],theta)) * v2_i
    M1 = M1 + M1_i
  }

  # compute full information matrix
  if (ignore_stage1) {
    alpha = 1
  }
  else {
    N0 = sum(n0)
    alpha = N1 / (N0 + N1)
  }
  M = alpha*M1 + (1-alpha)*M0

  # compute objective value
  if (!checkMinv(M)) {
    cvec = NULL
    obj_val = NA
  }
  else {

    if (obj == 'c') {
      cvec = cvec
      obj_val = t(cvec) %*% solve(M) %*% cvec
    }
    else if (obj == 'D') {
      cvec = NULL
      obj_val = -log(det(M))
    }
  }


  return(list(
    x = x,
    w = w,
    obj = obj,
    obj_val = obj_val,
    max_dose = max_dose,
    grad_fun = grad_fun,
    dr_fun = dr_fun,
    theta = theta,
    a = alpha,
    cvec = cvec,
    M0 = M0,
    M1 = M1,
    M = M
  ))
}
