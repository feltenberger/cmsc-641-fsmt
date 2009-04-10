package edu.umbc.algorithms.fmst.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.awt.geom.Point2D;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.SerializationUtils;
import org.apache.log4j.Logger;

import edu.umbc.algorithms.fmst.FairSMT;
import edu.umbc.algorithms.fmst.Point;


/**
 * @author dave
 */
public class GraphUtils {
    /**
     * the logger....
     */
    private static final Logger log = Logger.getLogger(GraphUtils.class);
    /**
     * A random number generator for our utils class.
     */
    private static Random rand = new Random();

    /**
     * Generate a random point within the range specified.
     * The offset will be from any edge in the plane bounded by the range.
     * e.g. a range of x=10, y=10 and offset=1 means the x and y values will be
     * between 1 and 9 instead of 0 and 10. This is mostly useful for visualization.
     *
     * @param xRange
     * @param yRange
     * @param offset
     * @return
     */
    public static Point randomPoint(int xRange, int yRange, int offset) {
        int x = rand.nextInt(xRange - (2 * offset)) + offset;
        int y = rand.nextInt(yRange - (2 * offset)) + offset;
        return new Point(x, y);
    }

    /**
     * Generate a random point and also generate numNeighbors
     * additional points that are neighbors to the new point.
     *
     * @param numNeighbors
     * @return
     */
    public static List<Point> randomPoint(int xRange, int yRange, int offset, int numNeighbors) {
        List<Point> newPoints = new ArrayList<Point>();
        Point p = randomPoint(xRange, yRange, offset);
        newPoints.add(p);
        for (int i = 0; i < numNeighbors; i++) {
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
        sortByX(points, true);
    }

    /**
     * sort the points by their weight.
     *
     * @param points
     */
    public static void sortByWeight(List<Point> points, final boolean includeSteinerNodes) {
        Collections.sort(points, new Comparator<Point>() {
            public int compare(Point o1, Point o2) {
                // if we shouldn't include steiner nodes, then always push them to the right/end of the array.
                if (!includeSteinerNodes) {
                    if (o1.isSteiner() && !o2.isSteiner())
                        return 1;
                    else if (o1.isSteiner() && o2.isSteiner())
                        return 0;
                    else if (!o1.isSteiner() && o2.isSteiner())
                        return -1;
                }

                if (o1.w < o2.w)
                    return -1;
                else if (o1.w > o2.w)
                    return 1;
                return 0;
            }
        });
    }

    /**
     * sort the points by the x value.
     *
     * @param points
     * @param includeSteinerNodes
     */
    public static void sortByX(List<Point> points, final boolean includeSteinerNodes) {
        Collections.sort(points, new Comparator<Point>() {
            public int compare(Point o1, Point o2) {
                if (!includeSteinerNodes) {
                    if (o1.isSteiner() && !o2.isSteiner())
                        return 1;
                    else if (o1.isSteiner() && o2.isSteiner())
                        return 0;
                    else if (!o1.isSteiner() && o2.isSteiner())
                        return -1;
                }
                if (o1.x < o2.x)
                    return -1;
                else if (o1.x > o2.x)
                    return 1;
                return 0;
            }
        });
    }

    /**
     * sort such that te1[0]<te1[1]<...<te1[NNh-1]
     * this is the crap ass one that was written in the original and is too
     * ugly to bother modifying.
     *
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
     *
     * @param file
     * @param smt
     * @param includeSteinerNodes
     * @throws Exception
     */
    public static void saveSMTPlainText(String file, FairSMT smt, boolean includeSteinerNodes) throws Exception {
        File out = new File(file);
        if (!out.getName().endsWith("smt"))
            out = new File(out.getAbsolutePath() + ".smt");

        StringBuffer sb = new StringBuffer();
        sb.append(smt.height).append(" ").append(smt.width).append("\n");
        for (Point p : smt.getPoints()) {
            // if we should skip steiner nodes and this is one, then skip it.
            if (!includeSteinerNodes && p.steiner > 0)
                continue;

            sb.append(p.x).append(" ").
                    append(p.y).append(" ").
                    append(p.w).append(" ").
                    append(p.steiner).append("\n");
        }
        FileUtils.writeStringToFile(out, sb.toString());
    }

    /**
     * Serializes the Steiner Minimum Tree and saves it in the specified file.
     *
     * @param file
     * @param smt
     * @throws Exception
     */
    public static void writeSerializedSMT(String file, FairSMT smt) throws Exception {
        // serialize the FairSMT object into a byte[]
        byte[] msSerialized = SerializationUtils.serialize(smt);
        // store the bytes to a file.
        FileUtils.writeByteArrayToFile(new File(file), msSerialized);
    }

    /**
     * Loads the reference SMT object.
     *
     * @param file
     * @return
     * @throws Exception
     */
    public static FairSMT readSerializedSMT(String file) throws Exception {
        FairSMT referenceSMT = null;
        File reference = new File(file);
        if (reference.exists()) {
            byte[] tmp = FileUtils.readFileToByteArray(reference);
            try {
                log.info("Trying to load " + file);
                referenceSMT = (FairSMT) SerializationUtils.deserialize(tmp);
                log.info("Successfully loaded " + file);
            } catch (Exception e) {
                log.error("Error loading the SMT!", e);
                throw e;
            }
        }
        return referenceSMT;
    }

    /**
     * @param file
     * @return
     * @throws Exception
     */
    @SuppressWarnings("unchecked")
    public static FairSMT readPlainTextSMT(String file) throws Exception {
        FairSMT fsmt = null;

        File in = new File(file);
        if (!in.exists())
            return null;

        List<String> lines = FileUtils.readLines(in);
        List<Point> allPoints = new ArrayList<Point>();
        int height = 0;
        int width = 0;
        for (String line : lines) {
            if (line == null)
                continue;
            String[] parts = line.split(" ");
            // this should be the first line...
            if (parts.length == 2) {
                height = Integer.parseInt(parts[0]);
                width = Integer.parseInt(parts[1]);
            } else if (parts.length == 4) {
                double x = Double.parseDouble(parts[0]);
                double y = Double.parseDouble(parts[1]);
                double w = Double.parseDouble(parts[2]);
                int steiner = Integer.parseInt(parts[3]);
                Point p = new Point(x, y);
                p.steiner = steiner;

                if (w <= 0.0)
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
     * @param x1 x-coordinate of first point
     * @param y1 y-coordinate of first point
     * @param x2 x-coordinate of second point
     * @param y2 y-coordinate of second point
     * @return distance between first and second point
     */
    public static double euclideanDistance(double x1, double y1, double x2, double y2) {
        return Math.pow(Math.pow(x2 - x1, 2.0) + Math.pow(y2 - y1, 2.0), 0.5);
    }

    /**
     * Assume a line whose starting coordinate is (x1,y1) and final coordinate is (x2, y2)
     * This function returns the coordinates of new point which sits along the line and the distance
     * between initial point and new point = distance
     * @param x1 x-coordinate of initial point
     * @param y1 y-coordinate of initial point
     * @param x2 x-coordinate of final point
     * @param y2 y-coordinate of final point
     * @param distance distance from initial point
     * @return coordinate of new point
     */

    public static Point2D getCoordinates(double x1, double y1, double x2, double y2, double distance) {
        double nx, ny;
        Point2D result = new Point2D.Double();
        double D = euclideanDistance(x1, y1, x2, y2);
        if (distance >= D) {
            result.setLocation(x2, y2);
            return result;
        }

        double abs_diff_y = Math.abs(y2 - y1) * distance / D;
        double abs_diff_x = Math.abs(x2 - x1) * distance / D;

        if (y2 > y1) {
            ny = y1 + abs_diff_y;
        } else if (y2 < y1) {
            ny = y1 - abs_diff_y;
        } else {
            ny = y1;
        }
        if (x2 > x1) {
            nx = x1 + abs_diff_x;
        } else if (x2 < x1) {
            nx = x1 - abs_diff_x;
        } else {
            nx = x1;
        }

        result.setLocation(nx, ny);
        return result;
    }
    
    /**
     * calculates the slope of the line going through two points
     */
    public static double slope(Point a, Point b)
	{
		return (a.y -b.y)/(a.x-b.x);
	}
    
    /**
     * calculates the center point of the circle which has points a,b and c on its circumference
     */
    public static Point center(Point a, Point b, Point c)
	{
		double ma = slope(b,a);
		double mb = slope(c,b);

		Point center = new Point(0.0,0.0);
		center.x = (ma * mb *(a.y - c.y) + mb * (a.x + b.x) - ma * (b.x + c.x)) / (2 * (mb - ma));
		center.y = -1 * (1/ma)*(center.x - (a.x + b.x) / 2) + (a.y + b.y) / 2;

		return center;
	}

	/**
	 * Gets the standard deviation of the power consumption rate.
	 * 
	 * @param minPoints
	 * @param numNodesInSMT
	 * @return standard deviation of power consumption rate in graph.
	 */
	public static double getStandardDevOfPCR(List<Point>minPoints, int numNodesInSMT) {
		double mean = getMeanPCR(minPoints, numNodesInSMT);
		double variance = 0;
		int counter = 0;
		for (int i = 0; i < numNodesInSMT; i++) {
			Point p1 = minPoints.get(i);
			//if(p1.isSteiner())
			{
				variance += Math.pow((p1.getPCR()-mean),2);
				counter++;
			}
		}
		return Math.sqrt(variance/counter);
	}

	/**
	 * Gets the average power consumption rate of all points in the graph.
	 * 
	 * @param minPoints
	 * @param numNodesInSMT
	 * @return average power consumption rate of nodes in the graph
	 */
	public static double getMeanPCR(List<Point>minPoints, int numNodesInSMT) {
		double avg = 0;
		int counter = 0;
		for (int i = 0; i < numNodesInSMT; i++) {
			Point p1 = minPoints.get(i);
			//if(p1.isSteiner())
			{
				avg += p1.getPCR();
				counter++;
			}
		}
		return avg/counter;
	}

	/**
	 * Dumps statistics.
	 */
	public static void printStatistics(FairSMT smt) {
		List<Point>minPoints = smt.getMinPoints();
		int numNodesInSMT = smt.getNumNodesInSMT();
		double max = Double.NEGATIVE_INFINITY;
		double min = Double.POSITIVE_INFINITY;
		double total = 0;
		int numPoints = 0;
		for (Point p : minPoints )
		{
			double temp = p.getPCR();
			if (temp > max) max = temp;
			if (temp < min) min = temp;

			total += temp;
			numPoints++;				
			//System.out.println("PCR = " + p.getPCR());
		}
		System.out.println("STDEV     = " + GraphUtils.getStandardDevOfPCR(minPoints, numNodesInSMT));		
		System.out.println("MAX POWER = " + max);
		System.out.println("MIN POWER = " + min);
		System.out.println("AVG POWER = " + total / numPoints);
		System.out.println("TOTAL POW = " + total);
		System.out.println("K VALUE   = " + smt.getMaxRelayNodes());
		System.out.println("TRGT STDD = " + smt.getTargetPCR());
		System.out.println("TOTAL POW = " + total);
		System.out.println("------------------------------------");
	}

	/**
	 * appends the results of a make fair run to the provided string.
	 * 
	 * @param results
	 * @param smt
	 * @param runNumber the iteration
	 * @return
	 */
	public static String appendResults(String results, FairSMT smt, int runNumber) {
		StringBuffer sb = new StringBuffer(results).append("\n");

		List<Point>minPoints = smt.getMinPoints();
		int numNodesInSMT = smt.getNumNodesInSMT();
		double max = Double.NEGATIVE_INFINITY;
		double min = Double.POSITIVE_INFINITY;
		double total = 0;
		int numPoints = 0;
		for (Point p : minPoints )
		{
			double temp = p.getPCR();
			if (temp > max) max = temp;
			if (temp < min) min = temp;

			total += temp;
			numPoints++;				
			//System.out.println("PCR = " + p.getPCR());
		}

		double stdDev = GraphUtils.getStandardDevOfPCR(minPoints, numNodesInSMT);
		// tab delimited list of values
		// run num	std dev	max pcr	min pcr	avg pcr	total pcr	k	target pcr

		// the run number
		sb.append(runNumber).append("\t")
		// standard deviation
		.append(stdDev).append("\t")
		// maximum consuming node
		.append(max).append("\t")
		// minimum consuming node
		.append(min).append("\t")
		// average power consumption
		.append(total / numPoints).append("\t")
		// total power consumption
		.append(total).append("\t")
		// k value
		.append(smt.getMaxRelayNodes()).append("\t")
		// target power consumption rate.
		.append(smt.getTargetPCR()).append("\t")
		;

		return sb.toString();
	}
}
