import java.util.function.BiFunction;
import java.util.function.Function;

public class Clarified implements Method {

	private Method meth;

	public Clarified(Method meth) {
		this.meth = meth;
	}

	@Override
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, Function<Double, Double> fi,
			Function<Double, Double> fix, double m, double y0, int N) {
		double[] ydoubled = meth.solve(a, b, f, fi, fix, m, y0, 2 * N - 1);
		double[] yunary = meth.solve(a, b, f, fi, fix, m, y0, N);
		double[] y = new double[N];
		double eps = 0;
		double realeps = 0;
		for (int q = 0; q < N; q++) {
			y[q] = ydoubled[2 * q] + (ydoubled[2 * q] - yunary[q]) / (Math.pow(2, meth.getP()) - 1);
			eps = Math.max(eps, Math.abs(ydoubled[2*q] - yunary[q]) / (Math.pow(2, meth.getP()) - 1));
			realeps = Math.max(realeps, Math.abs(y[q] - fi.apply(a + (b - a) * q * 1. / (N - 1))));
		}
		System.out.println("N = " + N + "; eps = " + eps+"; real eps = "+realeps);
		return y;
	}

	@Override
	public String getName() {
		return "Clarified: " + meth.getName();
	}

	@Override
	public int getP() {
		return meth.getP() + 1;
	}
}
