import java.awt.Color;
import java.text.DecimalFormat;
import java.util.function.BiFunction;
import java.util.function.Function;

import javax.swing.JFrame;
import javax.swing.JTable;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.chart.axis.NumberAxis;

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
		a = -3;
		b = 3;
		fi = (x) -> (x * x + 5 * Math.sin(2 * x) + Math.cos(x));
		fix = (x) -> (2 * x + 10 * Math.cos(2 * x) - Math.sin(x));
		m = 0.9;
		f = (x, y) -> (fix.apply(x) + m * (y - fi.apply(x)));
		realFunc = new double[2][realN];
		for (int i = 0; i < realN; i++) {
			realFunc[0][i] = a + i * (b - a) / (realN + 1);
			realFunc[1][i] = fi.apply(realFunc[0][i]);
			realmin = realmin > realFunc[1][i] ? realFunc[1][i] : realmin;
			realmax = realmax < realFunc[1][i] ? realFunc[1][i] : realmax;
		}

		// analyze(new EilerMethod(), new int[] {10, 20, 40, 80, 160, 320, 640,
		// 1280, 2560, 5120, 10240, 20480, 40960, 81920, 163840, 327680,
		// 655360, 1310720, 2621440, 5242880, 10485760, 20971520}, false);
		// analyze(new EilerMethod(), new int[] {10, 20, 40, 80, 160, 320, 640,
		// 1280, 2560, 5120, 10240, 20480, 40960, 81920}, true);

//		int[] N = new int[] { 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960, 81920, 163840,
//				327680, 655360, 1310720, 2621440};
		// analyze(new PredCorr(1e-15), N, false);
		// predCorrEpsilonAnylize(N);
		//
		// N = new int[] {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120,
		// 10240, 20480, 40960, 81920, 163840, 327680};
		// analyze(new RungeKutt(), N, false);

		 int[] N = new int[] {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480};
		// analyze(new GirMethod(), N, false);

		// N = new int[] {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120,
		// 10240, 20480, 40960, 81920, 163840};
		//
		// analyze(new RungeKuttSecond(), N, false);
		//
		//// N = new int[] {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120,
		// 10240};
		// analyze(new AdamsBashfordMoulton(), N, false);
		//
		// analyze(new GirMethodSecond(), N, false);

		Method[] allMethods = new Method[] { new EilerMethod(), new PredCorr(1e-10), new RungeKutt(), new GirMethod(),
				new RungeKuttSecond(), new AdamsBashfordMoulton(), new GirMethodSecond() };
		// compare(80, allMethods);
//		compareEps(allMethods, N);


//		Method[] girsVariations = new Method[] { new GirMethodExperimental1(), new GirMethodExperimental2(),
//				new GirMethod(), new GirMethodExperimental4() };
//		compareEps(girsVariations, N);
		
//		Method[] rungeclarified = new Method[]{new RungeKuttSecond(), new Clarified(new Clarified(new RungeKuttSecond()))};
//		compare(20, rungeclarified);
//		compareEps(rungeclarified, N);
//		analyze(rungeclarified[0], N, false);
//		analyze(rungeclarified[1], N, false);
		
		analyze(new EilerMethod(), new int[] {3710}, true);
	}

	private static void compareEps(Method[] meth, int[] Narr) {
		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
		// chart.removeLegend();

		double[] xaxis = new double[] { a, b };
		data.addSeries("X axis", new double[][] { xaxis, new double[2] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		addRealFunction(data, renderer);

		double min = realmin;
		double max = realmax;

		double[][] epsMethod = new double[meth.length][Narr.length];

		for (int q = 0; q < meth.length; q++) {
			Method m = meth[q];
			String[] columnNames = { "N", "Eps" };
			Object[][] epsTableData = new Object[Narr.length+1][2];
			epsTableData[0] = columnNames;
			for (int w = 0; w < Narr.length; w++) {
				double[] x = new double[Narr[w]];
				for (int i = 0; i < Narr[w]; i++) {
					x[i] = a + i * (b - a) / (Narr[w] - 1);
				}
				double[] y = m.solve(a, b, f, fi, fix, Main.m, fi.apply(a), Narr[w]);

				for (int e = 0; e < y.length; e++) {
					epsMethod[q][w] = Math.max(epsMethod[q][w], Math.abs(fi.apply(x[e]) - y[e]));
				}
			}
			for (int w = 1; w < Narr.length+1; w++) {
				epsTableData[w][0] = ""+Narr[w-1];
				epsTableData[w][1] = epsMethod[q][w - 1];
			}
			JFrame epstableFrame = new JFrame("Epstable, "+m.getName());
			// tableFrame.setLayout(new );
			epstableFrame.add(new JTable(epsTableData, columnNames));
			epstableFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			epstableFrame.pack();
			epstableFrame.setVisible(true);
		}
		double[] Narrdouble = new double[Narr.length];
		for (int q = 0; q < Narr.length; q++) {
			Narrdouble[q] = 1.0 * Narr[q];
		}
		double[][][] epschart = new double[meth.length][][];
		String[] names = new String[meth.length];
		int[] p = new int[meth.length];
		for (int q = 0; q < meth.length; q++) {
			names[q] = meth[q].getName();
			epschart[q] = new double[][] { Narrdouble, epsMethod[q] };
			p[q] = meth[q].getP();
		}

		epsChart(epschart, names, new int[0], "Errors comparsion");
	}

	private static void compare(int N, Method[] meth) {
		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
		// chart.removeLegend();

		double[] xaxis = new double[] { a, b };
		data.addSeries("X axis", new double[][] { xaxis, new double[2] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		addRealFunction(data, renderer);

		double min = realmin;
		double max = realmax;

		double[] epsMethod = new double[meth.length];
		for (int q = 0; q < meth.length; q++) {
			Method m = meth[q];
			double[] x = new double[N];
			for (int i = 0; i < N; i++) {
				x[i] = a + i * (b - a) / (N - 1);
			}
			double[] y = m.solve(a, b, f, fi, fix, Main.m, fi.apply(a), N);

			for (int w = 0; w < y.length; w++) {
				epsMethod[q] = Math.max(epsMethod[q], Math.abs(fi.apply(x[w]) - y[w]));
			}

			data.addSeries(m.getName(), new double[][] { x, y });
			// renderer.setSeriesPaint(data.getSeriesCount() - 1, col);
			// renderer.setSeries
		}
		frame = createFrame(chart, "Comparsion", a, b, min, max);

		String[] columnNames = { "Method", "Eps" };
		Object[][] tableData = new Object[meth.length + 1][2];
		tableData[0] = columnNames;
		for (int q = 1; q < meth.length + 1; q++) {
			tableData[q][0] = meth[q - 1].getName();
			tableData[q][1] = epsMethod[q - 1];
		}
		JFrame tableFrame = new JFrame("Comparsion table");
		// tableFrame.setLayout(new );
		tableFrame.add(new JTable(tableData, columnNames));
		tableFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		tableFrame.pack();
		tableFrame.setVisible(true);
	}

	private static void analyze(Method m, int[] Narr, boolean drawCharts) {

		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
		chart.removeLegend();

		double[] xaxis = new double[] { a, b };
		data.addSeries("X axis", new double[][] { xaxis, new double[2] });
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

		addRealFunction(data, renderer);

		double eps[] = new double[Narr.length];
		double min = realmin;
		double max = realmax;
		for (int q = 0; q < Narr.length; q++) {
			int N = Narr[q];
			double[] x = new double[N];
			for (int i = 0; i < N; i++) {
				x[i] = a + i * (b - a) / (N - 1);
			}
			double[] y = m.solve(a, b, f, fi, fix, Main.m, fi.apply(a), N);

			for (int i = 0; i < N; i++) {
				eps[q] = Math.max(eps[q], Math.abs(fi.apply(x[i]) - y[i]));
				// min = min > y[i] ? y[i] : min;
				// max = max < y[i] ? y[i] : max;
			}

			float fcol = 0.7f * (Narr.length - q) / (Narr.length + 1);
			Color col = new Color(fcol, fcol, fcol);
			if (drawCharts) {
				data.addSeries("Approximate chart. " + m.getName() + ". N = " + N, new double[][] { x, y });
				renderer.setSeriesPaint(data.getSeriesCount() - 1, col);
			}
		}

		// System.out.println(m.getName());
		// System.out.println("N\teps\tk\tl");
		String[] k = new String[Narr.length];
		String[] l = new String[Narr.length];
		for (int q = 0; q < Narr.length; q++) {
			String kk = "";
			if (q > 0) {
				kk = "" + Math.log(Math.abs(eps[q - 1] / eps[q])) / Math.log(2);
			}
			k[q] = kk;
			String ll = "";
			if (q > 1) {
				ll = "" + Math.log(Math.abs((eps[q - 2] - eps[q - 1]) / (eps[q - 1] - eps[q]))) / Math.log(2);
			}
			l[q] = ll;
			// System.out.print(Narr[q]+" "+eps[q]+" "+k[q]+" "+l[q]+"\n");
		}
		// System.out.println("");

		double[] Narrdouble = new double[Narr.length];
		for (int q = 0; q < Narr.length; q++) {
			Narrdouble[q] = 1.0 * Narr[q];
		}
		epsChart(new double[][][] { { Narrdouble, eps } }, new String[] { "Errors. " + m.getName() },
				new int[] { m.getP() }, "Errors. " + m.getName());
		// min = min - 0.1;
		// max = max + 0.1;
		if (drawCharts)
			frame = createFrame(chart, m.getName(), a, b, min, max);

		///////////////////////////// tables
		String[] columnNames = { "N", "Eps", "K", "L" };
		Object[][] tableData = new Object[Narr.length + 1][4];
		tableData[0] = columnNames;
		for (int q = 1; q < Narr.length + 1; q++) {
			tableData[q][0] = Narr[q - 1];
			tableData[q][1] = eps[q - 1];
			tableData[q][2] = k[q - 1];
			tableData[q][3] = l[q - 1];
		}
		JFrame tableFrame = new JFrame(m.getName() + ". Table");
		// tableFrame.setLayout(new );
		tableFrame.add(new JTable(tableData, columnNames));
		tableFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		tableFrame.pack();
		tableFrame.setVisible(true);
	}

	private static void epsChart(double[][][] eps, String[] names, int[] p, String title) {
		DefaultXYDataset data;
		ChartFrame frame;
		JFreeChart chart;

		data = new DefaultXYDataset();
		chart = ChartFactory.createXYLineChart(" ", "X", "Y", data);
		// chart.getLegend().setVisible(false);
		XYPlot plot = (XYPlot) chart.getPlot();

		LogarithmicAxis domain = new LogarithmicAxis("N");
		domain.setNumberFormatOverride(new DecimalFormat("0E0"));
		LogarithmicAxis range = new LogarithmicAxis("Eps");
		range.setNumberFormatOverride(new DecimalFormat("0E0"));

		plot.setDomainAxis(domain);
		plot.setRangeAxis(range);

		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();

		for (int q = 0; q < eps.length; q++) {
			data.addSeries(names[q], eps[q]);
			if (eps.length == 1)
				renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.black);

			// data.addSeries("X axis", new double[][] { xaxis, new double[2]
			// });
		}
		double begin = eps[0][0][0];
		double end = eps[0][0][eps[0][0].length - 1];
		for (int q = 0; q < p.length; q++) {
			double[][] ch = new double[2][5000];
			for (int w = 0; w < ch[0].length; w++) {
				ch[0][w] = begin + w * (end - begin) / (ch[0].length - 1);
				ch[1][w] = Math.pow(ch[0][w], -p[q]);
			}
			data.addSeries("N^(-" + p[q] + ")", ch);
			renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
		}

		frame = createFrame(chart, title);// , 0, end, 0, );
	}

	private static void addRealFunction(DefaultXYDataset data, XYLineAndShapeRenderer renderer) {
		data.addSeries("Real function", realFunc);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
	}

	private static void predCorrEpsilonAnylize(int[] Narr) {
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
		for (int q = 0; q < Narr.length; q++) {
			Narrdouble[q] = 1.0 * Narr[q];
		}
		int from = 2;
		int to = 16;
		Object[][] iterTableData = new Object[Narr.length + 1][to - from + 1];
		Object[][] epsTableData = new Object[Narr.length + 1][to - from + 1];
		DecimalFormat formatter = new DecimalFormat("0.000E0");

		String[] columnNames = new String[to - from + 1];
		columnNames[0] = "N";
		iterTableData[0][0] = "N";
		epsTableData[0][0] = "N";
		for (int q = from; q < to; q++) {
			double delta = Math.pow((0.1), q);
			iterTableData[0][q - from + 1] = formatter.format(delta);
			epsTableData[0][q - from + 1] = formatter.format(delta);
			columnNames[q - from + 1] = "" + delta;

			PredCorr m = new PredCorr(delta);
			double[] eps = new double[Narr.length];
			int[] iters = new int[Narr.length];
			double min = realmin;
			double max = realmax;
			for (int w = 0; w < Narr.length; w++) {
				int N = Narr[w];
				double[] x = new double[N];
				for (int i = 0; i < N; i++) {
					x[i] = a + i * (b - a) / (N - 1);
				}
				double[] y = m.solve(a, b, f, fi, fix, Main.m, fi.apply(a), N);
				iterTableData[w + 1][q - from + 1] = m.getMaxIter();

				for (int i = 0; i < N; i++) {
					eps[w] = Math.max(eps[w], Math.abs(fi.apply(x[i]) - y[i]));
				}
				epsTableData[w + 1][q - from + 1] = formatter.format(eps[w]);
			}
			float fcol = 0.9f * (to - q) / (to + 1);
			Color col = new Color(fcol, fcol, fcol);
			data.addSeries("Eps, delta = " + delta, new double[][] { Narrdouble, eps });
			renderer.setSeriesPaint(data.getSeriesCount() - 1, col);
		}
		for (int w = 0; w < Narr.length; w++) {
			iterTableData[w + 1][0] = Narr[w];
			epsTableData[w + 1][0] = Narr[w];
		}
		/// iterTable
		JFrame itertableFrame = new JFrame("Pred-corr analyze, iter count");
		// tableFrame.setLayout(new );
		itertableFrame.add(new JTable(iterTableData, columnNames));
		itertableFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		itertableFrame.pack();
		itertableFrame.setVisible(true);
		/// iterTable

		/// epsTable
		JFrame epstableFrame = new JFrame("Pred-corr analyze, eps");
		// tableFrame.setLayout(new );
		epstableFrame.add(new JTable(epsTableData, columnNames));
		epstableFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		epstableFrame.pack();
		epstableFrame.setVisible(true);
		/// epsTable

		double begin = Narrdouble[0];
		double end = Narrdouble[Narrdouble.length - 1];
		double[][] ch = new double[2][5000];
		for (int q = 0; q < ch[0].length; q++) {
			ch[0][q] = begin + q * (end - begin) / (ch[0].length - 1);
			ch[1][q] = Math.pow(ch[0][q], -2);
		}
		data.addSeries("N^(-" + 2 + ")", ch);
		renderer.setSeriesPaint(data.getSeriesCount() - 1, Color.red);
		// data.addSeries("X axis", new double[][] { xaxis, new double[2] });
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
		// xyplot.getRangeAxis().setRange(yb, ye);
		// xyplot.getDomainAxis().setRange(xb, xe);

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
