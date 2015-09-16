
public class Point {
	private double x, y, u, eps;

	public Point() {
		x = y = u = eps = 0;
	}

	public Point(double x, double y, double u, double eps) {
		this.x = x;
		this.y = y;
		this.u = u;
		this.eps = eps;
	}

	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}

	public double getU() {
		return u;
	}

	public void setU(double u) {
		this.u = u;
	}

	public double getEps() {
		return eps;
	}

	public void setEps(double eps) {
		this.eps = eps;
	}
	
}
