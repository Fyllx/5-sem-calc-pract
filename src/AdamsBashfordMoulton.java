import java.util.function.BiFunction;
import java.util.function.Function;

public class AdamsBashfordMoulton implements Method {

	@Override
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, Function<Double, Double> fi,
			Function<Double, Double> fix, double m, double y0, int N){
		double h = (b-a) / (N-1); 
		double[] y = new double[N];
		int[] s = new int[N];
		y[0] = y0;
		
		//Runge-Kutt's method
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
			
			double f3 = f.apply(xi - 3 * h, y[i-3]);
			double f2 = f.apply(xi - 2 * h, y[i-2]);
			double f1 = f.apply(xi - 1 * h, y[i-1]);
			double f0 = f.apply(xi, y[i]);
			
			double yvol = y[i] + (55. * f0 - 59. * f1 + 37. * f2 - 9. * f3) / 24. * h;
			y[i+1] = y[i] + (9. * f.apply(xi + h, yvol) + 19. * f0 - 5. * f1 + f2) / 24. * h;
		}
		return y;
	}
	
	@Override
	public String getName()
	{
		return "Adams-Bashford-Moulton's method";
	}

	@Override
	public int getP()
	{
		return 4;
	}
}
