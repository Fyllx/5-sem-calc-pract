import java.util.function.BiFunction;

public class PredCorr implements Method {

	/**
	 * O(h^2)
	 */
	@Override
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, double y0, int N) {
		double h = (b-a) / (N-1); 
		double[] y = new double[N];
		int[] s = new int[N];
		y[0] = y0;
		double ypred, ycorr;
		for (int i = 0; i < N - 1; i++) {
			double xi = a + i * h;
			
			double fi = f.apply(xi, y[i]);
			ypred = y[i] + h * fi;
			while(true)
			{
				ycorr = y[i] + h* (fi + f.apply(xi + h, ypred)) /2.;
				s[i+1]++;
				if(Math.abs(ycorr - ypred) < 1e-3) break;
				ypred = ycorr;
			}
			y[i+1] = ycorr; 
		}
		return y;
	}
	
	@Override
	public String getName()
	{
		return "Predictor-corrector method";
	}

}
