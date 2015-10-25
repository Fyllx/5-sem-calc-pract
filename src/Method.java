import java.util.function.BiFunction;
import java.util.function.Function;

public interface Method {
	
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
	public double[] solve(double a, double b, BiFunction<Double, Double, Double> f, Function<Double, Double> fi,
			Function<Double, Double> fix, double m, double y0, int N);

	public String getName();
	public int getP();
}
