import java.awt.Color;
import java.util.function.BiFunction;
import java.util.function.Function;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.DefaultXYDataset;

public class Main {

	private static Function<Double, Double> fi;
	private static Function<Double, Double> fix;
	private static BiFunction<Double, Double, Double> f;
	private static double m;
	private static double a, b;

	private static int realN = 5000;
	private static double[][] realFunc;
	private static double realmin = Double.MAX_VALUE, realmax = -Double.MAX_VALUE;

	public static void main(String[] args) {
		a = -6;
		b = 6;
		fi = (x) -> (x * x + 3 * Math.sin(x) + Math.cos(x));
		fix = (x) -> (2 * x + 3 * Math.cos(x) - Math.sin(x));
		m = 0.9;
		f = (x, y) -> (fix.apply(x) + m * (y - fi.apply(x)));
		realFunc = new double[2][realN];
		for (int i = 0; i < realN; i++) {
			realFunc[0][i] = a + i * (b - a) / (realN + 1);
			realFunc[1][i] = fi.apply(realFunc[0][i]);
			realmin = realmin > realFunc[1][i] ? realFunc[1][i] : realmin;
			realmax = realmax < realFunc[1][i] ? realFunc[1][i] : realmax;
		}
		
//		int[] N = new int[] {160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960, 81920};
		int[] N = new int[] {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960, 81920};//, 163840};//, 327680};
//		analyze(new EilerMethod(), N);
		
		analyze(new PredCorr(1e-8), N);
//		predCorrEpsilonAnylize(N);
		
//		analyze(new RungeKutt(), N);
		
//		analyze(new GirMethod(), N);
		
//		analyze(new RungeKuttSecond(), N);
		
//		N = new int[] {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240};
		analyze(new AdamsBashfordMoulton(), N);
		
		analyze(new GirMethodSecond(), N);
	}

	private static void analyze(Method m, int[] Narr) {
		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
		chart.removeLegend();
				
		double[] xaxis = new double[]{a, b};
		data.addSeries("X axis", new double[][] { xaxis, new double[2] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		addRealFunction(data, renderer);
		
		double eps[] = new double[Narr.length];		
		double min = realmin;
		double max = realmax;
		for(int q=0;q<Narr.length;q++)
		{
			int N = Narr[q];
			double[] x = new double[N];
			for (int i = 0; i < N; i++) {
				x[i] = a + i * (b - a) / (N - 1);
			}
			double[] y = m.solve(a, b, f, fi, fix, Main.m, fi.apply(a), N);
			
			for(int i=0;i<N;i++)
			{
				eps[q] = Math.max(eps[q], Math.abs(fi.apply(x[i]) - y[i]));
//				min = min > y[i] ? y[i] : min;
//				max = max < y[i] ? y[i] : max;
			}
			
			float fcol = 0.7f * (Narr.length - q) / (Narr.length + 1);
			Color col = new Color(fcol, fcol, fcol);
			data.addSeries("Approximate chart. " + m.getName() + ". N = "+N, new double[][] { x, y });
			renderer.setSeriesPaint(data.getSeriesCount() - 1, col);
		}
		
		System.out.println(m.getName());
		System.out.println("N\teps\tk\tl");
		for(int q=0;q<Narr.length;q++)
		{
			String k = "";
			if(q>0)
			{
				k = ""+Math.log(Math.abs(eps[q-1]/eps[q]))/Math.log(2);
			}
			String l = "";
			if(q>1)
			{
				l = ""+Math.log(Math.abs((eps[q-2]-eps[q-1])/(eps[q-1]-eps[q])))/Math.log(2);
			}
			System.out.print(Narr[q]+" "+eps[q]+" "+k+" "+l+"\n");
		}
		System.out.println("");
		
		double[] Narrdouble = new double[Narr.length];
		for(int q=0;q<Narr.length;q++)
		{
			Narrdouble[q] = 1.0 * Narr[q];
		}
		epsChart(new double[][]{Narrdouble, eps}, m.getP(), "Errors. "+m.getName());
//		min = min - 0.1;
//		max = max + 0.1;		
		
		frame = createFrame(chart, m.getName(), a, b, min, max);
	}
	private static void epsChart(double[][] eps, int p, String title) {
		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		plot.setDomainAxis(new LogarithmicAxis("N"));
		plot.setRangeAxis(new LogarithmicAxis("Eps"));
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
				
		data.addSeries("Eps", eps);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);
		
		double begin = eps[0][0];
		double end = eps[0][eps[0].length-1];
		double[][] ch = new double[2][5000];
		for(int q=0;q<ch[0].length;q++)
		{
			ch[0][q] = begin + q * (end - begin) / (ch[0].length - 1);
			ch[1][q] = Math.pow(ch[0][q], -p);
		}
		data.addSeries("N^(-"+p+")", ch);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
//		data.addSeries("X axis", new double[][] { xaxis, new double[2] });

		frame = createFrame(chart, title);//, 0, end, 0, );
	}

	private static void addRealFunction(DefaultXYDataset data, XYLineAndShapeRenderer renderer) {
		data.addSeries("Real function", realFunc);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
	}
	private static void predCorrEpsilonAnylize(int[] Narr)
	{
		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		plot.setDomainAxis(new LogarithmicAxis("N"));
		plot.setRangeAxis(new LogarithmicAxis("Eps"));
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();

		double[] Narrdouble = new double[Narr.length];
		for(int q=0;q<Narr.length;q++)
		{
			Narrdouble[q] = 1.0 * Narr[q];
		}
		int from = 2;
		int to = 12;
		for(int q=from; q < to; q++)
		{
			double delta = Math.pow((0.1), q);
			PredCorr m = new PredCorr(delta);
			double eps[] = new double[Narr.length];		
			double min = realmin;
			double max = realmax;
			for(int w=0;w<Narr.length;w++)
			{
				int N = Narr[w];
				double[] x = new double[N];
				for (int i = 0; i < N; i++) {
					x[i] = a + i * (b - a) / (N - 1);
				}
				double[] y = m.solve(a, b, f, fi, fix, Main.m, fi.apply(a), N);
				
				for(int i=0;i<N;i++)
				{
					eps[w] = Math.max(eps[w], Math.abs(fi.apply(x[i]) - y[i]));
				}
			}
			float fcol = 0.9f * (to - q) / (to + 1);
			Color col = new Color(fcol, fcol, fcol);
			data.addSeries("Eps, delta = "+delta, new double[][] {Narrdouble, eps} );			
			renderer.setSeriesPaint(data.getSeriesCount() - 1, col);
		}	
		
		double begin = Narrdouble[0];
		double end = Narrdouble[Narrdouble.length-1];
		double[][] ch = new double[2][5000];
		for(int q=0;q<ch[0].length;q++)
		{
			ch[0][q] = begin + q * (end - begin) / (ch[0].length - 1);
			ch[1][q] = Math.pow(ch[0][q], -2);
		}
		data.addSeries("N^(-"+2+")", ch);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
//		data.addSeries("X axis", new double[][] { xaxis, new double[2] });		
		frame = createFrame(chart, "PredCorr error analyze");
	}

	private static ChartFrame createFrame(JFreeChart chart, String name, double xb, double xe, double yb, double ye) {
		XYPlot xyplot = (XYPlot) chart.getPlot();
		xyplot.setBackgroundPaint(Color.white);
		xyplot.setRangeGridlinePaint(Color.gray);
		xyplot.setDomainGridlinePaint(Color.gray);
		xyplot.getRangeAxis().setRange(yb, ye);
		xyplot.getDomainAxis().setRange(xb, xe);

		XYLineAndShapeRenderer lineandshaperenderer = (XYLineAndShapeRenderer) xyplot.getRenderer();
		lineandshaperenderer.setDrawOutlines(true);
		lineandshaperenderer.setUseFillPaint(true);
		lineandshaperenderer.setBaseFillPaint(Color.white);

		ChartFrame frame = new ChartFrame(name, chart);

		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		return frame;
	}
	private static ChartFrame createFrame(JFreeChart chart, String name) {
		XYPlot xyplot = (XYPlot) chart.getPlot();
		xyplot.setBackgroundPaint(Color.white);
		xyplot.setRangeGridlinePaint(Color.gray);
		xyplot.setDomainGridlinePaint(Color.gray);
//		xyplot.getRangeAxis().setRange(yb, ye);
//		xyplot.getDomainAxis().setRange(xb, xe);

		XYLineAndShapeRenderer lineandshaperenderer = (XYLineAndShapeRenderer) xyplot.getRenderer();
		lineandshaperenderer.setDrawOutlines(true);
		lineandshaperenderer.setUseFillPaint(true);
		lineandshaperenderer.setBaseFillPaint(Color.white);

		ChartFrame frame = new ChartFrame(name, chart);

		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.pack();
		frame.setVisible(true);
		return frame;
	}
}
