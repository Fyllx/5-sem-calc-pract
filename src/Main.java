import java.awt.Color;
import java.util.function.BiFunction;
import java.util.function.Function;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
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
		
		doWork(new EilerMethod(), 2560);
		doWork(new PredCorr(), 200);
		doWork(new RungeKutt(), 20);
		doWork(new GirMethod(), 160);
		doWork(new RungeKuttSecond(), 40);		
		doWork(new AdamsBashfordMoulton(), 20);
		doWork(new GirMethodSecond(), 80);
	}

	private static void doWork(Method m, int N) {
		double[] x = new double[N];
		for (int i = 0; i < N; i++) {
			x[i] = a + i * (b - a) / (N - 1);
		}
		double[] y = m.solve(a, b, f, fi.apply(a), N);
		double min = realmin;
		double max = realmax;
		for(int i=0;i<N;i++)
		{
			min = min > y[i] ? y[i] : min;
			max = max < y[i] ? y[i] : max;
		}
		min = min - 0.1;
		max = max + 0.1;

		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();

		data.addSeries("X axis", new double[][] { x, new double[N] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		data.addSeries("Approximate chart. " + m.getName(), new double[][] { x, y });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
		
		addRealFunction(data, renderer);
		
		frame = createFrame(chart, m.getName() + "; N = " + N, a, b, min, max);
	}
	
	private static void doWork(GirMethodAbstract m, int N) {
		double[] x = new double[N];
		for (int i = 0; i < N; i++) {
			x[i] = a + i * (b - a) / (N - 1);
		}
		double[] y = m.solve(a, b, f, fi, fix, Main.m, fi.apply(a), N);
		double min = realmin;
		double max = realmax;
		for(int i=0;i<N;i++)
		{
			min = min > y[i] ? y[i] : min;
			max = max < y[i] ? y[i] : max;
		}
		min = min - 0.1;
		max = max + 0.1;

		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();

		data.addSeries("X axis", new double[][] { x, new double[N] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		data.addSeries("Approximate chart. " + m.getName(), new double[][] { x, y });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
		
		addRealFunction(data, renderer);
		
		frame = createFrame(chart, m.getName() + "; N = " + N, a, b, min, max);
	}

	private static void addRealFunction(DefaultXYDataset data, XYLineAndShapeRenderer renderer) {
		data.addSeries("Real function", realFunc);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);
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
}
