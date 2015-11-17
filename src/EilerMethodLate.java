import java.util.function.Function;

public class EilerMethodLate implements MethodLate {

	/**
	 * O(h)
	 */
	@Override
	public double[] solve(double tau, double a, double b, Function<Double, Double> fi, int M){
		double h = (a-tau) / M;
		int N = (int) ((b-tau) / h);
		double[] y = new double[N];		
		for(int q=0;q<M+1;q++)
		{
			y[q] = fi.apply(tau+h*q); //tau + h*(M-1) = a
		}
		
		for (int i = M; i < N - 1; i++) {
			y[i+1] = y[i] - h * y[i-M]; 
		}
		return y;
	}
	
	@Override
	public String getName()
	{
		return "Euler's method";
	}

	@Override
	public int getP() {
		return 1;
	}
}
