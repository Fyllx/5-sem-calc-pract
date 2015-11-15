import java.util.function.BiFunction;
import java.util.function.Function;

public interface MethodLate {
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @param f
	 * @param y0 = fi(a)
	 * @param N
	 * @return
	 */
//	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, double y0, int N);
	public double[] solve(double tau, double a, double b, Function<Double, Double> fi, int M);

	public String getName();
	public int getP();
}
