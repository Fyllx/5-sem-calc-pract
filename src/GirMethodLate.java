import java.util.function.Function;

public class GirMethodLate implements MethodLate{

	/**
	 * O(h^3)
	 */
	@Override
	public double[] solve(double tau, double a, double b, Function<Double, Double> fi, int M) {		
		double h = (a-tau) / (M - 1);
		int N = (int) ((b-tau) / h);
		double[] y = new double[N];
		for(int q=0;q<M;q++)
		{
			y[q] = fi.apply(tau+h*q); //tau + h*(M-1) = a
		}
		
		for (int i = M-1; i < N - 1; i++) {
			double xi = tau + i * h;

			y[i + 1] = (18.*y[i] - 9.*y[i-1] + 2.*y[i-2] - 6*h*y[i+1-M+1])/11.;
		}
		return y;
	}

	@Override
	public String getName() {
		return "Gir's method of 3rd order";
	}
	
	@Override
	public int getP() {
		return 3;
	}
}
