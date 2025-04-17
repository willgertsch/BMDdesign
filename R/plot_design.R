plot_design = function(opt_design) {

  x = opt_design$x
  w = opt_design$w

  plot_dat = data.frame(
    dose = seq(1e-4, opt_design$max_dose)
  )
  plot_dat$y = opt_design$dr_fun(plot_dat$dose, opt_design$theta)

  ggplot2::ggplot(plot_dat, ggplot2::aes(x = dose, y = y)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::geom_col(
      data = data.frame(
        dose = x,
        weight = w
      ),
      ggplot2::aes(x = dose, y = weight),
      fill = 'red', alpha = 0.6, width = opt_design$max_dose/40
    ) +
    ggplot2::scale_y_continuous(
      name = 'allocation weight',
      sec.axis = ggplot2::dup_axis(name = 'P(Y=1)')
    )
}
