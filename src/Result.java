import java.util.ArrayList;

public class Result {	
	ArrayList<Point> set;
	
	public Result() {
		set = new ArrayList<>();
	}
	
	public void add(Point p)
	{
		set.add(p);
	}
	public Point get(int i)
	{
		return set.get(i);
	}

}
