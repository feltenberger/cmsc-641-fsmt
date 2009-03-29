package edu.umbc.algorithms.fmst;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import org.apache.commons.io.FileUtils;

import edu.umbc.algorithms.old.FairSMT;

/**
 * @author dave
 *
 */
public class GraphUtils {
	/**
	 * A random number generator for our utils class.
	 */
	private static Random rand = new Random();
	/**
	 * Generate a random point within the range specified.
	 * The offset will be from any edge in the plane bounded by the range.
	 * e.g. a range of x=10, y=10 and offset=1 means the x and y values will be
	 * between 1 and 9 instead of 0 and 10. This is mostly useful for visualization.
	 * @param xRange
	 * @param yRange
	 * @param offset
	 * @return
	 */
	public static Point randomPoint(int xRange, int yRange, int offset) {
		int x = rand.nextInt(xRange-(2*offset)) + offset;
		int y = rand.nextInt(yRange-(2*offset)) + offset;
		return new Point(x, y);
	}

	/**
	 * Generate a random point and also generate numNeighbors
	 * additional points that are neighbors to the new point.
	 * @param numNeighbors
	 * @return
	 */
	public static List<Point> randomPoint(int xRange, int yRange, int offset, int numNeighbors) {
		List<Point> newPoints = new ArrayList<Point>();
		Point p = randomPoint(xRange, yRange, offset);
		newPoints.add(p);
		for(int i = 0; i < numNeighbors; i++) {
			Point neighbor = randomPoint(xRange, yRange, offset);
			newPoints.add(neighbor);
			//neighbor.addNeighbor(p);
			//p.addNeighbor(neighbor);
		}
		return newPoints;
	}

	/**
	 * @param points
	 */
	public static void sortByWeight(List<Point> points) {
		Collections.sort(points, new Comparator<Point>() {
			public int compare(Point o1, Point o2) {
				if(o1.w < o2.w)
					return -1;
				else if(o1.w > o2.w)
					return 1;
				return 0;
			}
		});
	}
	public static void sortByX(List<Point> points) {
		Collections.sort(points, new Comparator<Point>() {
			public int compare(Point o1, Point o2) {
				if(o1.x < o2.x)
					return -1;
				else if(o1.x > o2.x)
					return 1;
				return 0;
			}
		});
	}

	/**
	 * sort such that te1[0]<te1[1]<...<te1[NNh-1]
	 * @param te1
	 * @param te2
	 * @param te3
	 * @param NNh
	 */
	public static void sort(double te1[], double te2[], double te3[], int NNh) {
		int kk, kks, ii, jj, mm;
		double b1, b2, b3, c1, c2, c3;
		kks = (int) (NNh / 2);
		for (kk = kks; kk >= 1; kk--) {
			ii = kk;
			b1 = te1[ii - 1];
			b2 = te2[ii - 1];
			b3 = te3[ii - 1];
			while (2 * ii <= NNh) {
				jj = 2 * ii;
				if (jj + 1 <= NNh) {
					if (te1[jj - 1] < te1[jj]) {
						jj++;
					}
				}
				if (te1[jj - 1] <= b1) {
					break;
				}
				te1[ii - 1] = te1[jj - 1];
				te2[ii - 1] = te2[jj - 1];
				te3[ii - 1] = te3[jj - 1];
				ii = jj;
			}// wend
			te1[ii - 1] = b1;
			te2[ii - 1] = b2;
			te3[ii - 1] = b3;
		}// next kk
		for (mm = NNh - 1; mm >= 1; mm--) {
			c1 = te1[mm];
			c2 = te2[mm];
			c3 = te3[mm];
			te1[mm] = te1[0];
			te2[mm] = te2[0];
			te3[mm] = te3[0];
			ii = 1;
			while (2 * ii <= mm) {
				kk = 2 * ii;
				if (kk + 1 <= mm) {
					if (te1[kk - 1] <= te1[kk]) {
						kk++;
					}
				}
				if (te1[kk - 1] <= c1) {
					break;
				}
				te1[ii - 1] = te1[kk - 1];
				te2[ii - 1] = te2[kk - 1];
				te3[ii - 1] = te3[kk - 1];
				ii = kk;
			}// wend
			te1[ii - 1] = c1;
			te2[ii - 1] = c2;
			te3[ii - 1] = c3;
		}// next mm
	}

	/**
	 * Dumps a Steiner minimum tree to the file system.
	 * @param file
	 * @param smt
	 * @param includeSteinerNodes
	 * @throws Exception
	 */
	public static void dumpSMT(String file, FairSMT smt, boolean includeSteinerNodes) throws Exception {
		File out = new File(file);
		StringBuffer sb = new StringBuffer();
		sb.append(smt.height).append(" ").append(smt.width).append("\n");
		for(Point p : smt.getPoints()) {
			// if we should skip steiner nodes and this is one, then skip it.
			if(!includeSteinerNodes && p.steiner > 0)
				continue;

			sb.append(p.x).append(" ").
				append(p.y).append(" ").
				append(p.w).append(" ").
				append(p.steiner).append("\n");
		}
		FileUtils.writeStringToFile(out, sb.toString());
	}

	/**
	 * @param file
	 * @return
	 * @throws Exception
	 */
	@SuppressWarnings("unchecked")
	public static FairSMT readSMT(String file) throws Exception {
		FairSMT fsmt = null;

		File in = new File(file);
		if(!in.exists())
			return null;

		List<String> lines = FileUtils.readLines(in);
		List<Point> allPoints = new ArrayList<Point>();
		int height = 0;
		int width = 0;
		for(String line : lines) {
			if(line == null)
				continue;
			String[] parts = line.split(" ");
			// this should be the first line...
			if(parts.length == 2) {
				height = Integer.parseInt(parts[0]);
				width = Integer.parseInt(parts[1]);
			}
			else if(parts.length == 4) {
				double x = Double.parseDouble(parts[0]);
				double y = Double.parseDouble(parts[1]);
				double w = Double.parseDouble(parts[2]);
				int steiner = Integer.parseInt(parts[3]);
				Point p = new Point(x, y);
				p.steiner = steiner;

				if(w <= 0.0)
					p.generateW();
				else
					p.w = w;
				allPoints.add(p);
			}
		}

		fsmt = new FairSMT(width, height, allPoints.size());
		fsmt.getPoints().clear();
		fsmt.getPoints().addAll(allPoints);
		return fsmt;
	}
	/**
	 * The euclidean distance between two points.
	 * 
	 * @param d1
	 * @param d2
	 * @param d3
	 * @param d4
	 * @return
	 */
	public static double euclideanDistance(double d1, double d2, double d3, double d4) {
		double dw;
		dw = Math.pow(Math.pow(d3 - d1, 2.0) + Math.pow(d4 - d2, 2.0), 0.5);
		return dw;
	}

	
}
