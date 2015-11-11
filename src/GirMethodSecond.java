import java.util.function.BiFunction;
import java.util.function.Function;

public class GirMethodSecond implements Method{

	/**
	 * O(h^4)
	 */
	@Override
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, Function<Double, Double> fi,
			Function<Double, Double> fix, double m, double y0, int N) {
		double h = (b - a) / (N - 1);
		double[] y = new double[N];
		y[0] = y0;

		// Runge-Kutt (second)
		for (int i = 0; i < 3; i++) {
			double xi = a + i * h;
			double k1 = h * f.apply(xi, y[i]);
			double k2 = h * f.apply(xi + h / 2., y[i] + k1 / 2.);
			double k3 = h * f.apply(xi + h / 2., y[i] + k2 / 2.);
			double k4 = h * f.apply(xi + h, y[i] + k3);
			y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		}

		for (int i = 3; i < N - 1; i++) {
			double xi = a + i * h;

			y[i + 1] = ((fix.apply(xi + h) - m * fi.apply(xi + h)) * 12. * h + 48. * y[i] - 36. * y[i - 1] + 16. * y[i - 2]
					- 3. * y[i - 3]) / (25. - 12. * m * h);
		}
		return y;
	}

	@Override
	public String getName() {
		return "Gir's method of 4th order";
	}

	@Override
	public int getP() {
		return 4;
	}
}
