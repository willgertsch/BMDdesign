BMDdesign = function(grad_fun, dr_fun, bmd_grad = NULL, obj, theta, d0, n0,
                        N1, num_doses, max_dose, swarm, iter, alg,
                        ignore_stage1 = F, show_progress = T) {

  # compute cvec if finding c-optimal design

  # compute information matrix from initial design
  F0 = sapply(d0, grad_fun, theta)
  phat0 = sapply(d0, dr_fun, theta)
  v0 = 1/(phat0 * (1 - phat0))
  M0 = n0[1] * F0[, 1] %*% t(F0[, 1]) * v0[1]
  for (i in 2:length(d0)) {
    M0 = M0 + n0[i] * F0[, i] %*% t(F0[, i]) * v0[i]
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
        cvec = bmd_grad(0.1, theta)
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
      cvec = bmd_grad(0.1, theta)
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
