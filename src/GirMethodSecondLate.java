import java.util.function.Function;

public class GirMethodSecondLate implements MethodLate{

	/**
	 * O(h^4)
	 */
	@Override
	public double[] solve(double tau, double a, double b, Function<Double, Double> fi, int M) {
		double h = (a-tau) / M;
		int N = (int) ((b-tau) / h);
		double[] y = new double[N];
		for(int q=0;q<M+1;q++)
		{
			y[q] = fi.apply(tau+h*q); //tau + h*(M-1) = a
		}

		for (int i = M; i < N - 1; i++) {
			y[i + 1] = (48.*y[i] - 36.*y[i-1] + 16.*y[i-2] - 3.*y[i-3] - 12.*h*y[i+1-M])/25.;
		}
		return y;
	}

	@Override
	public String getName() {
		return "Gir's method of 4th order";
	}

	@Override
	public int getP() {
		return 1;
	}
}
