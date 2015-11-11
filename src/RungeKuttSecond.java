import java.util.function.BiFunction;
import java.util.function.Function;

public class RungeKuttSecond implements Method {

	/**
	 * O(h^4)
	 */
	@Override
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, Function<Double, Double> fi,
			Function<Double, Double> fix, double m, double y0, int N){
		double h = (b - a) / (N - 1);
		double[] y = new double[N];
		y[0] = y0;
		for (int i = 0; i < N - 1; i++) {
			double xi = a + i * h;
			double k1 = h * f.apply(xi, y[i]);
			double k2 = h * f.apply(xi + h / 2., y[i] + k1 / 2.);
			double k3 = h * f.apply(xi + h / 2., y[i] + k2 / 2.);
			double k4 = h * f.apply(xi + h, y[i] + k3);
			y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		}
		return y;
	}

	@Override
	public String getName() {
		return "Runge-Kutt's method of 4th order";
	}

	@Override
	public int getP() {
		return 4;
	}
}
