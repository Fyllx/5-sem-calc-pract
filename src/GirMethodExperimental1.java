import java.util.function.BiFunction;
import java.util.function.Function;

public class GirMethodExperimental1 implements Method{

	/**
	 * O(h^3)
	 */
	@Override
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, Function<Double, Double> fi,
			Function<Double, Double> fix, double m, double y0, int N) {
		double h = (b - a) / (N - 1);
		double[] y = new double[N];
		y[0] = y0;

		// Euler
		for (int i = 0; i < 2; i++) {
			double xi = a + i * h;
			y[i+1] = y[i] + h * f.apply(xi, y[i]); 
		}

		for (int i = 2; i < N - 1; i++) {
			double xi = a + i * h;

			y[i + 1] = ((fix.apply(xi + h) - m * fi.apply(xi + h)) * 6 * h + 18 * y[i] - 9 * y[i - 1] + 2 * y[i - 2])
					/ (11 - 6 * m * h);
		}
		return y;
	}

	@Override
	public String getName() {
		return "Gir's method of 3rd order with Euler's first points";
	}
	
	@Override
	public int getP() {
		return 3;
	}
}

