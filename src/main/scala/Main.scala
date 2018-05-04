package io.github.mikerippe.rungekutta

import io.github.mikerippe.rungekutta.RungeKutta

object Main extends App {

  val rungeKutta = new RungeKutta()
  rungeKutta.compute(0, 40)

}
