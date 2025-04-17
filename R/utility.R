# utility functions for finding optimal designs
# checks if information matrix is invertible
# returns 1 if invertible and 0 if not
# can optimize easily for 2 dim
checkMinv = function(M) {

  if (class(try(solve(M),silent=T))[1]!="matrix")
    return(0)
  else
    return(1)
}

# look up relevant model information for a model type
mod_info_lookup = function(model_type) {

  if (model_type == 'qlinear') {
    grad_fun = qlinear.grad
    dr_fun = qlinear.fun
    bmd_grad = qlinear.bmdgrad
    pts = 2
  }
  else if (model_type == 'hill') {
    grad_fun = hill.grad
    dr_fun = hill.fun
    bmd_grad = hill.bmdgrad
    pts = 4
  }
  else if (model_type == 'log-logistic') {
    grad_fun = loglogistic.grad
    dr_fun = loglogistic.fun
    bmd_grad = loglogistic.bmdgrad
    pts = 3
  }
  else
    stop('invalid model type')

  return(list(grad_fun = grad_fun, dr_fun = dr_fun, bmd_grad = bmd_grad, pts = pts))
}
