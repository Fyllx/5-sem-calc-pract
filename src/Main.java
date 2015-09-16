import java.util.function.Function;

public class Main {

	public static void main(String[] args) 
	{		
		double a = -6, b = 6;
		
		Function<Double, Double> f = new Function<Double, Double>() {
			@Override
			public Double apply(Double x) {
				return x + 3 * Math.sin(x) + Math.cos(x);
			}			
		};		
		
//		Result r = EilerMethod.solve(a, b, N, f);		

	}

}