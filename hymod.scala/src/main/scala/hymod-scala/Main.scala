import org.saddle._
// import org.breeze._
import scala.math._

class Hymod {

  private def power(x: Double, y: Double): Double = {
    val a = abs(x)
    return pow(a, y)
  }

  def excess(x_loss: Double, cmax: Double, bexp: Double, prc: Double, pet: Double): Map[String, Double] = {
    var mut_prc = prc
    val xn_prev: Double = x_loss
    val x: Double = (1 - ((bexp + 1) * (xn_prev) / cmax))
    val y: Double = (1 / (bexp + 1))
    val p = power(x, y)
    val ct_prev: Double = cmax * (1 - p)
    val er1: Double = max((mut_prc - cmax + ct_prev), 0.0)
    mut_prc = mut_prc - er1
    //
    val dummy: Double = min(((ct_prev + mut_prc) / cmax), 1)
    var xn: Double = (cmax / (bexp + 1)) * (1 - power((1 - dummy), (bexp + 1)))

    // Calculate Effective rainfall 2
    val er2: Double = max(mut_prc - (xn - xn_prev), 0)

    // Alternative approach
    val evap: Double = (1 - (((cmax / (bexp + 1)) - xn) / (cmax / (bexp + 1)))) * pet // actual ET is linearly related to the soil moisture state
    xn = max(xn - evap, 0) // update state

    val state = Map(
      "excess_runoff_1" -> er1,
      "excess_runoff_2" -> er2,
      "soil_moist" -> xn)

    return state
  }

  def random_params(): collection.mutable.Map[String, Double] = {
    val param_bounds = Map(
      "cmax" -> Map[String, Double]("lower" -> 1.0, "upper" -> 100),
      "bexp" -> Map[String, Double]("lower" -> 0.0, "upper" -> 2.0),
      "alpha" -> Map[String, Double]("lower" -> 0.2, "upper" -> 0.99),
      "ks" -> Map[String, Double]("lower" -> 0.01, "upper" -> 0.5),
      "kq" -> Map[String, Double]("lower" -> 0.5, "upper" -> 1.2))

    val params = collection.mutable.Map("x" -> -999.0)

    for ((k, v) <- param_bounds) {
      val minV: Double = param_bounds(k)("lower")
      val maxV: Double = param_bounds(k)("upper")
      val sampler = breeze.stats.distributions.Uniform(minV, maxV)
      val x = sampler.sample(1)
      params += (k -> x(0))
    }

    return params
  }

  // def hargreaves(forcing: Frame[], tmaxCol:String, tminCol:String, dtCol:String) {
  //   val Gsc = 367
  //   val lhov = 2.257
  //
  //   val dts = forcings[dtCol]
  //   val tmin = forcings[tminCol]
  //   val tmax = forcings[tmaxCol]
  //   val n = len(tmax)
  //   val doy = [x.timetuple().tm_yday for x in dts]
  //
  //   tavg = pd.concat([tmin,tmax],axis=1).mean(axis=1).rename("tavg")
  //
  //   eto = np.zeros(n)
  //
  //   for i,t in enumerate(doy){
  //
  //     b = 2 * pi * (t/365)
  //     Rav = 1.00011 + 0.034221*cos(b) + 0.00128*sin(b) + 0.000719*cos(2*b) + 0.000077*sin(2*b)
  //     Ho = ((Gsc * Rav) * 86400)/1e6
  //
  //     val tmax = focings[tmaxCol]
  //     val tmin = focing[tminCol]
  //     val tavg =
  //
  //     eto[i] = (0.0023 * Ho * (tmax[i]-tmin[i])**0.5 * (tavg[i]+17.8))
  //   }
  //   return eto
  // }

  def simulate(forcing: Frame[Double, Double, Double], cmax: Float, bexp: Float, alpha: Float, ks: Float, kq: Float) {
    val lt_to_m: Double = 0.001

    // HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
    var n: Int = forcing.numRows
    val x_loss: Double = 0.0
    // Initialize slow tank state
    var x_slow: Double = 2.3503 / (ks * 22.5)
    // var x_slow: Double = 0.0 // --> works ok if calibration data starts with low discharge
    // Initialize state(s) of quick tank(s)
    var x_quick = vec.zeros(3)
    val t: Int = 0
    var outflow = vec.zeros(n)
    var output = vec.zeros(n)

    // START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS
    // for (i <- 0 to n) {
    //   var p = forcing.row(precipCol,i)
    //   var pet = forcing.row(petCol,i)
    // }
  }
}

object Main extends App {
  val file = io.CsvFile("/Users/kel/Documents/github/hymod/data/Q9H018_forcing.csv")
  var df = io.CsvParser.parse(List(0, 1, 2, 3, 4, 5))(file).withRowIndex(0).withColIndex(0)
  val idx = Index.make(time.RRule(time.DAILY),time.datetime(1981,1,1), time.datetime(2000,1,1))
  df.setRowIndex(idx.sorted)
  println(idx)
  val model = new Hymod
  val params = model.random_params()
  println(params)
  // val pet =
  println(df.col("tmin"))
  // val q = model.simulation()
  // println("Hello world")
}
