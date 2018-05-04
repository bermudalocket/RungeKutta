package io.github.mikerippe.rungekutta

import co.theasi.plotly.{Plot, draw}

class RungeKutta {
 
  // parameters
  val sigma: Double = -1
  val alpha: Double = 20 
  val beta: Float = 8f/3

  // step size
  val h = 0.01
  
  // system of equations modeling the Lorenz attractor, can be edited
  def f(t: Double, x: Double, y: Double, z: Double): Double = sigma*(y-x)
  def g(t: Double, x: Double, y: Double, z: Double): Double = x*alpha - x*z - y
  def r(t: Double, x: Double, y: Double, z: Double): Double = x*y - beta*z

  // initial conditions
  val fi = 0
  val gi = 20
  val ri = 20
 
  def compute(a: Double, b: Double) {
    
    // initialize timing
    val timeInit = System.currentTimeMillis

    // calculate iterations
    val N = (b-a)/h
    
    // initialize system of equations
    var t: Double = 0
    var x = scala.collection.mutable.Map[Int, Double]()
    var y = scala.collection.mutable.Map[Int, Double]()
    var z = scala.collection.mutable.Map[Int, Double]()
   
    // apply initial conditions
    x(0) = fi
    y(0) = gi
    z(0) = ri

    // proceed with iterations loop
    for (i <- 0 until (N.toInt - 1)) {

      val k1 = h * f(t, x(i), y(i), z(i))
      val l1 = h * g(t, x(i), y(i), z(i))
      val m1 = h * r(t, x(i), y(i), z(i))

      val k2 = h * f(t + h*0.5, x(i) + 0.5 * k1, y(i) + 0.5 * l1, z(i) + 0.5 * m1)
      val l2 = h * g(t + h*0.5, x(i) + 0.5 * k1, y(i) + 0.5 * l1, z(i) + 0.5 * m1)
      val m2 = h * r(t + h*0.5, x(i) + 0.5 * k1, y(i) + 0.5 * l1, z(i) + 0.5 * m1)

      val k3 = h * f(t + h*0.5, x(i) + 0.5 * k2, y(i) + 0.5 * l2, z(i) + 0.5 * m2)
      val l3 = h * g(t + h*0.5, x(i) + 0.5 * k2, y(i) + 0.5 * l2, z(i) + 0.5 * m2)
      val m3 = h * r(t + h*0.5, x(i) + 0.5 * k2, y(i) + 0.5 * l2, z(i) + 0.5 * m2)

      val k4 = h * f(t + h, x(i) + k3, y(i) + l3, z(i) + m3)
      val l4 = h * g(t + h, x(i) + k3, y(i) + l3, z(i) + m3)
      val m4 = h * r(t + h, x(i) + k3, y(i) + l3, z(i) + m3)
      
      val newx = x(i) + (1f/6) * (k1 + 2*k2 + 2*k3 + k4)
      val newy = y(i) + (1f/6) * (l1 + 2*l2 + 2*l3 + l4)
      val newz = z(i) + (1f/6) * (m1 + 2*m2 + 2*m3 + m4)

      x += (i+1 -> newx)
      y += (i+1 -> newy)
      z += (i+1 -> newz)
      t += h

    } // for

    // end timing
    val timeCompute = System.currentTimeMillis - timeInit
    println("Calculations completed in " + timeCompute + " ms.")
    println("Plotting...")

    // initialize timing for plotting
    val timePlotInit = System.currentTimeMillis

    // create axes
    val xaxis = x.values
    val yaxis = y.values

    // uses plot.ly api
    val plot = Plot().withScatter(xaxis, yaxis)
    draw(plot, "rungekutta")
    
    // calculate plot time
    val timePlotting = System.currentTimeMillis - timeCompute
    println("Plot created and request sent to plot.ly in " + timePlotting + " ms.")

  } // def compute

} // class RungeKutta
