import java.util.function.BiFunction;

public class EilerMethod implements Method {

	/**
	 * O(h)
	 */
	@Override
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, double y0, int N) {
		double h = (b-a) / (N-1); 
		double[] y = new double[N];
		y[0] = y0;
		for (int i = 0; i < N - 1; i++) {
			double xi = a + i * h;
			y[i+1] = y[i] + h * f.apply(xi, y[i]); 
		}
		return y;
	}
	
	@Override
	public String getName()
	{
		return "Eiler's method";
	}

}
