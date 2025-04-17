# equivalence theorem
# checks the optimality of a design
# input is list returned from design finding function
check_eq = function(opt_design) {

  # extract objects from design output
  max_dose = opt_design$max_dose
  grad_fun = opt_design$grad_fun
  dr_fun = opt_design$dr_fun
  theta = opt_design$theta
  a = opt_design$a
  cvec = opt_design$cvec
  M0 = opt_design$M0
  M1 = opt_design$M1
  obj = opt_design$obj

  dose = seq(.0001, max_dose)
  y = sapply(dose, sens, grad_fun, obj, M0, M1, dr_fun, theta, a, cvec)

  plot_dat = data.frame(
    gen_var = y,
    dose = dose
  )

  ggplot2::ggplot(plot_dat, ggplot2::aes(x = dose, y = gen_var)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = 'dose',
      y = 'directional derivative'
    )


}

sens = function(z, grad, obj, M0, M1, dr_fun, theta, a, cvec) {
  dg = grad(z, theta)
  M = (1-a)*M0 + a*M1

  if (obj == 'c') {
    Minv = solve(M)
    dM = Minv %*% cvec %*% t(cvec) %*% Minv
    p = dr_fun(z, theta)
    v = 1/(p*(1-p))
    y = v * t(dg) %*% dM %*% dg - sum(diag(M1 %*% dM))
  }
  else if (obj == 'D') {
    dM = solve(M)
    p = dr_fun(z, theta)
    v = 1/(p*(1-p))
    y = v * t(dg) %*% dM %*% dg - sum(diag(M1 %*% dM))
  }

  return(y)
}
