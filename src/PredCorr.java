import java.util.function.BiFunction;
import java.util.function.Function;

public class PredCorr implements Method {

	private double delta;
	private int maxiter = -1;
	public PredCorr(double delta) {
		this.delta = delta;
	}
	
	/**
	 * O(h^2)
	 */
	@Override
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, Function<Double, Double> fi,
			Function<Double, Double> fix, double m, double y0, int N){
		double h = (b-a) / (N-1); 
		double[] y = new double[N];
		int[] s = new int[N];
		y[0] = y0;
		double ypred, ycorr;
		int itermax = 0;
		for (int i = 0; i < N - 1; i++) {
			double xi = a + i * h;
			
			double fati = f.apply(xi, y[i]);
			ypred = y[i] + h * fati;
			int iter = 0;
			while(true)
			{
				ycorr = y[i] + h* (fati + f.apply(xi + h, ypred)) /2.;
				s[i+1]++;
				iter++;
				if(Math.abs((ycorr - ypred)/(ypred)) < delta) break;
				ypred = ycorr;
			}
			itermax = Math.max(itermax, iter);
			y[i+1] = ycorr; 
		}
		maxiter = itermax;
		System.out.println(getName()+": max iter = "+itermax);
		return y;
	}
	
	@Override
	public String getName()
	{
		return "Predictor-corrector method";
	}

	@Override
	public int getP() {
		return 2;
	}
	
	public int getMaxIter()
	{
		return maxiter;
	}
}
