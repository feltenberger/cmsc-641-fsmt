package edu.umbc.algorithms.fmst;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.apache.log4j.Logger;
import edu.umbc.algorithms.fmst.util.GraphUtils;
import javax.swing.JOptionPane;

import edu.umbc.algorithms.fmst.adt.MaxHeap;


/**
 * @author dave
 *
 */
public class FairSMT extends JPanel implements Runnable {
	private static final long serialVersionUID = 1268750618801149582L;
	private static final transient Logger log = Logger.getLogger(FairSMT.class);
	/**
	 * the thread that calculates all the steiner points
	 */
	private transient Thread computationThread;
	private transient Thread optimizationThread;

	/**
	 * the colors associated with drawing the canvas
	 */
	private Color bgColor;
	private Color textColor;
	private Color edgeColor;
	private Color nodeColor;
	private Color steinerNodeColor;

	/**
	 * the number of iterations where we tried to minimize the SMT
	 */
	private int numIterations;
	/**
	 * the number of times we found a better SMT than in the previous
	 * iteration.
	 */
	private int numBetterTreesFound;
	/**
	 * the size (in number of nodes) of the Steiner minimum tree,
	 * which includes Steiner nodes + non-Steiner nodes.
	 */
	private int numNodesInSMT;
	/**
	 * the number of "real", i.e. non-Steiner, nodes in the tree
	 */
	private int numNonSteinerNodes = 20;
	/**
	 * the height of the applet
	 */
	public int height = 300;
	/**
	 * the width of the applet
	 */
	public int width = 300;
	/**
	 * total number of nodes, including steiner nodes.
	 */
	//private int this.points.size();
	/**
	 * part of the number of convex hull nodes
	 */
	private int coch;
	/**
	 * the other part of the convex hull nodes
	 */
	private int coch2;

	private double steinerX, steinerY;
	/**
	 * the current minimum tree length
	 */
	private double minTreeLen;
	/**
	 * is this SMT running?
	 */
	private boolean isRunning = false;
	/**
	 * has this SMT run already?
	 */
	private boolean hasRun = false;

	/**
	 * the original points from the graph -- added once in the run method.
	 */
	private List<Point> originalPoints = new ArrayList<Point>();
	/**
	 * all the points in the graph
	 */
	private List<Point> points = Collections.synchronizedList(new ArrayList<Point>());
	/**
	 * a secondary array for temporary storage
	 */
	private List<Point> temporaryPoints = Collections.synchronizedList(new ArrayList<Point>());
	/**
	 * all the minimum points in the graph
	 */
	private List<Point> minPoints = Collections.synchronizedList(new ArrayList<Point>());

	// I have no idea wtf these are for...
	private double aa[] = new double[200];
	private double bb[] = new double[200];
	private double aa2[] = new double[200];
	private double bb2[] = new double[200];

	/**
	 * something to do with connections or something?
	 */
	private int cn[] = new int[2000];
	private int cn2[] = new int[200];
	private int connect[][] = new int[200][10];

	// add any new variables below this line so serialization isn't broken!
	// ********************************************************************
	/**
	 * the list of edges in the tree
	 */
	private List<Edge> edges = Collections.synchronizedList(new ArrayList<Edge>());
	/**
	 * the list of edges in the min tree
	 */
	private List<Edge> minEdges = Collections.synchronizedList(new ArrayList<Edge>());

	private Boolean runOnce = true;
	private Boolean runMakeFair = false;
	private int sleepTime = 100;

	/**
	 *
	 */
	public FairSMT() {
		init();
	}

	/**
	 * @param width
	 * @param height
	 * @param numNonSteinerNodes
	 */
	public FairSMT(int width, int height, int numNonSteinerNodes) {
		this.width = width;
		this.height = height;
		this.numNonSteinerNodes = numNonSteinerNodes;
		init();
	}

	/**
	 * Initialize the object
	 */
	private void init() {
		bgColor = Color.black;
		nodeColor = Color.white;
		textColor = Color.green;
		edgeColor = Color.green;
		steinerNodeColor = Color.red;
		setPreferredSize(new Dimension(width, height));
		setSize(getPreferredSize());

		generateRandomNonSteinerNodes();
	}

	/**
	 * Generates the number of non-nodes that this object is configured to have.
	 */
	private void generateRandomNonSteinerNodes() {
		// initialize the tree by adding numNodes nodes randomly on the canvas.
		// make the steinerArray all 0, i.e. no steiner nodes.
		for (int k = 0; k < numNonSteinerNodes; k++) {
			double x = Math.random() * (width - 30) + 15;
			double y = Math.random() * (width - 30) + 15;
			Point p = new Point(x, y);
			p.generateW();
			this.points.add(p);
		}
	}

	/**
	 * sets the x1 value.
	 * @param index
	 * @param x
	 */
	private void x(int index, double x) {
		this.points.get(index).x = x;
	}

	/**
	 * gets the x1 value
	 * @param index
	 * @return
	 */
	private double x(int index) {
		return this.points.get(index).x;
	}
	private int e1(int index) {
		if(index >= this.edges.size())
			return 0;
		return this.edges.get(index).index1;
		//return e1[index];
	}
	private void e1(int index, int value) {
		if(index < this.edges.size())
			this.edges.get(index).index1 = value;
		else
			getClass();
		//this.e1[index] = value;
	}
	private int e2(int index) {
		if(index >= this.edges.size())
			return 0;
		return this.edges.get(index).index2;
		//return e2[index];
	}
	private void e2(int index, int value) {
		if(index < this.edges.size())
			this.edges.get(index).index2 = value;
		else
			getClass();
		//this.e2[index] = value;
	}

	/**
	 * @param index
	 */
	private void copyPointToSecondary(int index) {
		Point p = this.points.get(index);
		if(this.temporaryPoints.size() == index)
			this.temporaryPoints.add(p.copy());
		else
			this.temporaryPoints.set(index, p.copy());
	}

	/**
	 * gets the x2 value
	 * @param index
	 * @return
	 */
	private double x2(int index) {
		return this.temporaryPoints.get(index).x;
	}

	/**
	 * gets the y1 value
	 * @param index
	 * @return
	 */
	private double y(int index) {
		return points.get(index).y;
	}
	/**
	 * sets the y1 value
	 * @param index
	 * @param y
	 */
	private void y(int index, double y) {
		points.get(index).y = y;
	}
	/**
	 * sets the weight
	 * @param index index
	 * @param w new weight
	 */
	private void w(int index, double w) {
		points.get(index).w = w;
	}
	/**
	 * returns the weight
	 * @param index index
	 * @return weight
	 */
	private double w(int index) {
		return points.get(index).w;
	}

	/**
	 * gets the y2 value
	 * @param index index
	 * @return y-value
	 */
	private double y2(int index) {
		return this.temporaryPoints.get(index).y;
	}

	/**
	 * gets whether this is a steiner node or not
	 * @param index index
	 * @return 1 if the node is a steiner point 0 otherwise
	 */
	private int isSteiner(int index) {
		return this.points.get(index).steiner;
	}
	/**
	 * sets the value of isSteinerNode
	 * @param index
	 * @param val
	 */
	private void setSteiner(int index, int val) {
		this.points.get(index).steiner = val;
	}

	/**
	 * compute convex hull
	 */
	public void convexHull() {
		double xmax = x(numNonSteinerNodes - 1);
		int firstNodeIsMaxX = 0;
		coch = 0;
		coch2 = 0;
		double lto = 0.0;
		double xx, yy;

		double half_pi = Math.PI / 2;

		while (lto != x2(0)) {
			xx = x(0);
			yy = y(0);
			for (int i = 1; i < numNonSteinerNodes; i++) {
				w(i, Math.asin((y(i) - yy) / GraphUtils.euclideanDistance(x(i), y(i), xx, yy)) + half_pi);
				if (x(i) > xx) {
					w(i, 3 * Math.PI / 2 - Math.asin((y(i) - yy) / GraphUtils.euclideanDistance(x(i), y(i), xx, yy)));
				}
				w(i, w(i) + Math.PI);
				if (w(i) > 2 * Math.PI) {
					w(i, w(i) - 2 * Math.PI);
				}
				if (firstNodeIsMaxX == 1) {
					w(i, Math.asin((y(i) - yy)
							/ GraphUtils.euclideanDistance(x(i), y(i), xx, yy))
							+ half_pi);
					if (x(i) > xx) {
						w(i, 3 * Math.PI / 2 - Math.asin((y(i) - yy) / GraphUtils.euclideanDistance(x(i), y(i), xx, yy)));
					}
					if (w(i) > 2 * Math.PI) {
						w(i, w(i) - 2 * Math.PI);
					}
				}
			} // next i

			x(0, xx);
			y(0, yy);
			w(0, 100.0);

			// sort by the w-values, whatever they are...
			//GraphUtils.sort(w1, x1, y1, numNonSteinerNodes);
			GraphUtils.sortByWeight(points, false);

			lto = x(0);
			if (firstNodeIsMaxX == 0) {
				aa[coch] = (y(0) - yy) / (x(0) - xx);
				bb[coch] = y(0) - aa[coch] * x(0);
				coch++;
			} else {
				aa2[coch2] = (y(0) - yy) / (x(0) - xx);
				bb2[coch2] = y(0) - aa2[coch2] * x(0);
				coch2++;
			}
			if (x(0) == xmax) {
				firstNodeIsMaxX = 1;
			}
		}// wend

		for (int k = 0; k < numNonSteinerNodes; k++) {
			x(k, x2(k));
			y(k, y2(k));
		}
	}// ch

	/**
	 * @param index
	 * @param edge
	 */
	private void insertOrAddEdge(int index, Edge edge) {
		if(this.edges.size() == index)
			this.edges.add(edge);
		else
			this.edges.set(index, edge);
	}

	/**
	 * compute min spanning tree
	 */
	private void mst() {
		double ei[] = new double[20000];
		double ej[] = new double[20000];
		double edgeDistances[] = new double[20000];

		int ncv[] = new int[200];
		int cv[][] = new int[200][200];
		int edgeIndex = 0;
		int count = 0;
		for (int imst = 0; imst < this.points.size() - 1; imst++) {
			if (isSteiner(imst) < 2) {
				for (int jmst = imst + 1; jmst < this.points.size(); jmst++) {
					if (isSteiner(jmst) < 2) {
						ei[count] = imst + 0.0;
						ej[count] = jmst + 0.0;
						edgeDistances[count] = GraphUtils.euclideanDistance(x(imst), y(imst), x(jmst),
								y(jmst));
						count++;
					}// sl[j]
				}// next j
			}// sl[i]
		}// next i

		// sort by the edge distances (i think?)
		GraphUtils.sort(edgeDistances, ei, ej, count);

		// initialize all the cn array to 0.
		for (int i = 0; i < this.points.size(); i++) {
			cn[i] = 0;
		}
		for (int i = 0; i < count; i++) {
			int label = 0;
			for (int j = 0; j < ncv[(int) ei[i]]; j++) {
				if (cv[(int) ei[i]][j] == (int) ej[i]) {
					label = 1;
				}
			}// j
			if (label == 0) {
				synchronized(this) {
					Edge edge = new Edge(null, null);
					edge.index1 = (int) ei[i];
					edge.index2 = (int) ej[i];
					this.insertOrAddEdge(edgeIndex, edge);

					e1(edgeIndex, (int) ei[i]);
					e2(edgeIndex, (int) ej[i]);

					edgeIndex++;
					//log.info(Thread.currentThread().getName());

					//connect[edge.index1][cn[edge.index1]] = edge.index2;
					connect[e1(edgeIndex - 1)][cn[e1(edgeIndex - 1)]] = e2(edgeIndex - 1);
					//cn[edge.index1]++;
					cn[e1(edgeIndex - 1)]++;

					//connect[edge.index2][cn[edge.index2]] = edge.index1;
					connect[e2(edgeIndex - 1)][cn[e2(edgeIndex - 1)]] = e1(edgeIndex - 1);
					//cn[edge.index2]++;
					cn[e2(edgeIndex - 1)]++;
				}

				int ni = ncv[(int) ei[i]];
				int nj = ncv[(int) ej[i]];

				cv[(int) ei[i]][ncv[(int) ei[i]]] = (int) ej[i];
				ncv[(int) ei[i]]++;
				for (int j = 0; j < nj; j++) {
					cv[(int) ei[i]][ncv[(int) ei[i]]] = cv[(int) ej[i]][j];
					ncv[(int) ei[i]]++;
					cv[cv[(int) ej[i]][j]][ncv[cv[(int) ej[i]][j]]] = (int) ei[i];
					ncv[cv[(int) ej[i]][j]]++;
					for (int k = 0; k < ni; k++) {
						cv[cv[(int) ej[i]][j]][ncv[cv[(int) ej[i]][j]]] = cv[(int) ei[i]][k];
						ncv[cv[(int) ej[i]][j]]++;
					}
				}// j

				cv[(int) ej[i]][ncv[(int) ej[i]]] = (int) ei[i];
				ncv[(int) ej[i]]++;
				for (int j = 0; j < ni; j++) {
					cv[(int) ej[i]][ncv[(int) ej[i]]] = cv[(int) ei[i]][j];
					ncv[(int) ej[i]]++;
					cv[cv[(int) ei[i]][j]][ncv[cv[(int) ei[i]][j]]] = (int) ej[i];
					ncv[cv[(int) ei[i]][j]]++;
					for (int k = 0; k < nj; k++) {
						cv[cv[(int) ei[i]][j]][ncv[cv[(int) ei[i]][j]]] = cv[(int) ej[i]][k];
						ncv[cv[(int) ei[i]][j]]++;
					}
				}// j
			}// label=0
		}// i
	}// mst

	/**
	 * compute value of the object function, in this case, the length of the
	 * total tree
	 * @return
	 */
	private double totalTreeLength() {
		double total = 0.0;
		for (int i = 0; i < this.points.size() - 1; i++) {
			total = total + GraphUtils.euclideanDistance(x(e1(i)), y(e1(i)), x(e2(i)), y(e2(i)));
			//total = total + euclideanDistance(x(e1(i)), y(e1[i]), x(e2[i]), y(e2[i]));
		}
		return total;
	}

	//
	private void bohe() {
		int numTotalNodesTemp;

		synchronized(points) {
			// copy arrays from x1 to x2
			for(int i = numNonSteinerNodes; i < this.points.size(); i++) {
				copyPointToSecondary(i);
				cn[i] = cn[i];
			}

			// i thinks this is saying that if this is a steiner node and it's on the
			// convex hull, make the steiner value 2 so we can remove it below
			for(int i = numNonSteinerNodes; i < this.points.size(); i++) {
				if (isSteiner(i) == 1 && cn[i] != 3) {
					setSteiner(i, 2);
					//isSteinerNode[hi] = 2;
				}// sl[i]==1 && cn[i]<3
			}// wend i

			// here we remove any steiner nodes with values higher than 1
			numTotalNodesTemp = numNonSteinerNodes;
			for(int i = numNonSteinerNodes; i < this.points.size(); i++) {
				if (isSteiner(i) == 1) {
					//x1[numTotalNodesTemp] = x2(hi);
					x(numTotalNodesTemp, x2(i));
					//y1[numTotalNodesTemp] = y2(hi);
					y(numTotalNodesTemp, y2(i));
					cn[numTotalNodesTemp] = cn2[i];
					numTotalNodesTemp++;
				}// sl[hi]==1
			}// i

			if(numTotalNodesTemp < this.points.size()) {
				while(numTotalNodesTemp < this.points.size()) {
					// keep plucking values off the end
					this.points.remove(this.points.size()-1);
				}
				//log.info("gotta remove some nodes!");
			}

			//this.points.size() = numTotalNodesTemp;
			for(int hi = numNonSteinerNodes; hi < this.points.size(); hi++) {
				setSteiner(hi, 1);
			}
		}

	}// bohe

	/**
	 * Compute the Steiner point (steinerX,steinerY) of points at indexes ij1, ij2, ij3
	 * @param ij1
	 * @param ij2
	 * @param ij3
	 */
	private void computeSteinerPoint(int ij1, int ij2, int ij3) {
		double ka, se, xz, yz;
		final double one_third_pi = Math.PI / 3.0;

		ka = (y(ij2) - y(ij1)) / (x(ij2) - x(ij1));
		se = y(ij1) - ka * x(ij1);
		xz = x(ij2) - x(ij1);
		yz = y(ij2) - y(ij1);
		if (y(ij3) < ka * x(ij3) + se && x(ij1) > x(ij2)) {
			steinerX = x(ij1) + Math.cos(one_third_pi) * xz + Math.sin(one_third_pi) * yz;
			steinerY = y(ij1) - Math.sin(one_third_pi) * xz + Math.cos(one_third_pi) * yz;
		}
		if (y(ij3) > ka * x(ij3) + se && x(ij1) < x(ij2)) {
			steinerX = x(ij1) + Math.cos(one_third_pi) * xz + Math.sin(one_third_pi) * yz;
			steinerY = y(ij1) - Math.sin(one_third_pi) * xz + Math.cos(one_third_pi) * yz;
		}
		if (y(ij3) < ka * x(ij3) + se && x(ij1) < x(ij2)) {
			steinerX = x(ij1) + Math.cos(-one_third_pi) * xz + Math.sin(-one_third_pi) * yz;
			steinerY = y(ij1) - Math.sin(-one_third_pi) * xz + Math.cos(-one_third_pi) * yz;
		}
		if (y(ij3) > ka * x(ij3) + se && x(ij1) > x(ij2)) {
			steinerX = x(ij1) + Math.cos(-one_third_pi) * xz + Math.sin(-one_third_pi) * yz;
			steinerY = y(ij1) - Math.sin(-one_third_pi) * xz + Math.cos(-one_third_pi) * yz;
		}
	}

	/**
	 * adds a new steiner node to the graph based on the following criteria:
	 * 1.) the steinerX,steinerY point that was calculated in computeSteinerPoint()
	 *      which averages the three points' locations and makes steinerX,steinerY that average.
	 * 2.) some other complicated crap involving averaging the distances and
	 *      adjusting the placement of points by the averages and stuff...?
	 *
	 * @param mij1
	 * @param mij2
	 * @param mij3
	 */
	private void addSteinerNodeToGraph(int mij1, int mij2, int mij3) {
		addSteinerNodeToGraph(this.points.size(), mij1, mij2, mij3);
		//this.points.size()++;
	}// addSteinerNodeToGraph

	/**
	 * adds a new steiner node to the graph based on the following criteria:
	 * 1.) the steinerX,steinerY point that was calculated in computeSteinerPoint()
	 *      which averages the three points' locations and makes steinerX,steinerY that average.
	 * 2.) some other complicated crap involving averaging the distances and
	 *      adjusting the placement of points by the averages and stuff...?
	 *
	 * @param steinerIndex
	 * @param mij1
	 * @param mij2
	 * @param mij3
	 */
	private void addSteinerNodeToGraph(int steinerIndex, int mij1, int mij2, int mij3) {
		double avgX, avgY, an, bn, distToMidpoint;

		// average the first, second, and steiner points' x-coord
		avgX = (x(mij1) + x(mij2) + steinerX) / 3.0;
		// average the first, second, and steiner points' y-coord
		avgY = (y(mij1) + y(mij2) + steinerY) / 3.0;
		log.debug("avgx = " + avgX + ", avgy = " + avgY);

		// calculate the distance from the first point to the midpoint
		distToMidpoint = GraphUtils.euclideanDistance(x(mij1), y(mij1), avgX, avgY);
		log.debug("dist to midpoint: " + distToMidpoint);

		// adjust all the points' coordinates by the average
		x(mij1, x(mij1) - avgX);
		y(mij1, y(mij1) - avgY);
		x(mij2, x(mij2) - avgX);
		y(mij2, y(mij2) - avgY);
		x(mij3, x(mij3) - avgX);
		y(mij3, y(mij3) - avgY);

		log.debug("steinerX = " + steinerX + ", y = " + steinerY);
		// same for the steiner point...
		steinerX = steinerX - avgX;
		steinerY = steinerY - avgY;
		log.debug("steinerX = " + steinerX + ", y = " + steinerY);

		// calculate the slope of the line from the third point to the steiner point
		an = (y(mij3) - steinerY) / (x(mij3) - steinerX);
		// dunno what this is...
		bn = steinerY - an * steinerX;
		log.debug("an = " + an + ", bn = " + bn);

		// if the steiner point's x coord is to the left of the third point's coord
		// do {whatever the hell it's doing} to add a new steiner node at the end of the array of points
		if (steinerX < x(mij3)) {
			double x =  (-an * bn + Math.pow(an * an * bn * bn
					- (bn * bn - distToMidpoint * distToMidpoint) * (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());

			// add a steiner point if we're at the end of the array.
			if(steinerIndex == this.points.size())
				this.points.add(new Point(0.0, 0.0));

			x(steinerIndex, x);

			double y = an * x(steinerIndex) + bn + 0.0001 * Math.sin(Math.PI * 2.0 * Math.random());
			y(steinerIndex, y);
			log.debug("x = " + x + ", y = " + y);

			//isSteinerNode[steinerIndex] = 1;
			setSteiner(steinerIndex, 1);
		} else {
			// otherwise, do something different to add the steiner node at the end of the array of points
			double x = (-an * bn - Math.pow(an * an * bn * bn
					- (bn * bn - distToMidpoint * distToMidpoint) * (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());

			// add a steiner point if we're at the end of the array.
			if(steinerIndex == this.points.size())
				this.points.add(new Point(0.0, 0.0));

			x(steinerIndex, x);

			double y = an * x(steinerIndex) + bn + 0.0001 * Math.sin(Math.PI * 2.0 * Math.random());
			y(steinerIndex, y);

			log.debug("x = " + x + ", y = " + y);

			setSteiner(steinerIndex, 1);
		}
		x(mij1, x(mij1) + avgX);
		y(mij1, y(mij1) + avgY);
		x(mij2, x(mij2) + avgX);
		y(mij2, y(mij2) + avgY);
		x(mij3, x(mij3) + avgX);
		y(mij3, y(mij3) + avgY);

		steinerX = steinerX + avgX;
		steinerY = steinerY + avgY;

		x(steinerIndex, x(steinerIndex) + avgX);
		y(steinerIndex, y(steinerIndex) + avgY);
		//if(num == 20)
		//	System.exit(20);
		//num++;
		log.debug("-----------------------------------");
	}
	//private static int num = 0;

	private void tase1(int it1) {
		int qint;
		double x, x2, y, y2, qa, qb, qc, qw, minq;

		x = x(connect[it1][0]);
		y = y(connect[it1][0]);
		x2 = x(connect[it1][1]);
		y2 = y(connect[it1][1]);
		qa = GraphUtils.euclideanDistance(x, y, x2, y2);
		qb = GraphUtils.euclideanDistance(x, y, x(it1), y(it1));
		qc = GraphUtils.euclideanDistance(x2, y2, x(it1), y(it1));
		qw = Math.acos((qb * qb + qc * qc - qa * qa) / (2.0 * qb * qc));
		minq = qw;
		qint = 1;
		x = x(connect[it1][1]);
		y = y(connect[it1][1]);
		x2 = x(connect[it1][2]);
		y2 = y(connect[it1][2]);
		qa = GraphUtils.euclideanDistance(x, y, x2, y2);
		qb = GraphUtils.euclideanDistance(x, y, x(it1), y(it1));
		qc = GraphUtils.euclideanDistance(x2, y2, x(it1), y(it1));
		qw = Math.acos((qb * qb + qc * qc - qa * qa) / (2.0 * qb * qc));
		if (qw < minq) {
			minq = qw;
			qint = 2;
		}
		x = x(connect[it1][2]);
		y = y(connect[it1][2]);
		x2 = x(connect[it1][0]);
		y2 = y(connect[it1][0]);
		qa = GraphUtils.euclideanDistance(x, y, x2, y2);
		qb = GraphUtils.euclideanDistance(x, y, x(it1), y(it1));
		qc = GraphUtils.euclideanDistance(x2, y2, x(it1), y(it1));
		qw = Math.acos((qb * qb + qc * qc - qa * qa) / (2.0 * qb * qc));
		if (qw < minq) {
			minq = qw;
			qint = 3;
		}
		if (qint == 1) {
			computeSteinerPoint(it1, connect[it1][0], connect[it1][1]);
			addSteinerNodeToGraph(it1, connect[it1][0], connect[it1][1]);
		}
		else if (qint == 2) {
			computeSteinerPoint(it1, connect[it1][1], connect[it1][2]);
			addSteinerNodeToGraph(it1, connect[it1][1], connect[it1][2]);
		}
		else if (qint == 3) {
			computeSteinerPoint(it1, connect[it1][2], connect[it1][0]);
			addSteinerNodeToGraph(it1, connect[it1][2], connect[it1][0]);
		}
	}// tase1

	private void tase2(int it2) {
		computeSteinerPoint(it2, connect[it2][0], connect[it2][1]);
		addSteinerNodeToGraph(it2, connect[it2][0], connect[it2][1]);
	}// tase2

	private void tuika() {
		int tlabel=0, tra=0, tra2=0;
		double one_third_pi = Math.PI / 3.0;
		while (tlabel == 0 && tra < 10) {
			tra++;
			int ter = 0;
			for (int it = 0; it < numNonSteinerNodes; it++) {
				if (cn[it] > 2) {
					tra2++;
					ter = 1;
					tase1(it);
					mst();
					hantei();
					if (this.points.size() > 2 * numNonSteinerNodes - 2) {
						break;
					}
				}// cn[it]==3
			}// it
			for (int i = 0; i < numNonSteinerNodes; i++) {
				if (cn[i] == 2) {
					double xx, xx2, yy, yy2, qa, qb, qc, qw;
					xx = x(connect[i][0]);
					yy = y(connect[i][0]);
					xx2 = x(connect[i][1]);
					yy2 = y(connect[i][1]);
					qa = GraphUtils.euclideanDistance(xx, yy, xx2, yy2);
					qb = GraphUtils.euclideanDistance(xx, yy, x(i), y(i));
					qc = GraphUtils.euclideanDistance(xx2, yy2, x(i), y(i));
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					if (qw < 2.0 * one_third_pi) {
						tra2++;
						ter = 2;
						tase2(i);
						mst();
						hantei();
						if (this.points.size() > 2 * numNonSteinerNodes - 2) {
							break;
						}
					}// qw>2.0*pi/3.0
				}// cn[it]==2
			}// it
			if (ter == 0) {
				tlabel = 1;
			}
		}// tlabel
		if (tra2 > 0) {
			bohe();
			mst();
			hantei();
		}
	}// tuika

	/**
	 * this loops through all the steiner nodes in the base tree
	 * and does stuff with them...  i don't know exactly what :-)
	 */
	private void mrlonely() {
		double xx, xx2, yy, yy2, qa, qb, qc, qw, maxq, qwe;
		int tlabel=0, tra=0;
		double one_third_pi = Math.PI / 3.0;
		while (tlabel == 0) {
			tra++;
			int ter = 0;
			for (int i = numNonSteinerNodes; i < this.points.size(); i++) {
				if (cn[i] == 3) {
					xx = x(connect[i][0]);
					yy = y(connect[i][0]);
					xx2 = x(connect[i][1]);
					yy2 = y(connect[i][1]);
					qa = GraphUtils.euclideanDistance(xx, yy, xx2, yy2);
					qb = GraphUtils.euclideanDistance(xx, yy, x(i), y(i));
					qc = GraphUtils.euclideanDistance(xx2, yy2, x(i), y(i));
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = GraphUtils.euclideanDistance(qw, 0.0, 2.0 * one_third_pi, 0.0);
					maxq = qwe;
					xx = x(connect[i][1]);
					yy = y(connect[i][1]);
					xx2 = x(connect[i][2]);
					yy2 = y(connect[i][2]);
					qa = GraphUtils.euclideanDistance(xx, yy, xx2, yy2);
					qb = GraphUtils.euclideanDistance(xx, yy, x(i), y(i));
					qc = GraphUtils.euclideanDistance(xx2, yy2, x(i), y(i));
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = GraphUtils.euclideanDistance(qw, 0.0, 2.0 * one_third_pi, 0.0);
					if (qwe > maxq) {
						maxq = qwe;
					}
					xx = x(connect[i][2]);
					yy = y(connect[i][2]);
					xx2 = x(connect[i][0]);
					yy2 = y(connect[i][0]);
					qa = GraphUtils.euclideanDistance(xx, yy, xx2, yy2);
					qb = GraphUtils.euclideanDistance(xx, yy, x(i), y(i));
					qc = GraphUtils.euclideanDistance(xx2, yy2, x(i), y(i));
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = GraphUtils.euclideanDistance(qw, 0.0, 2.0 * one_third_pi, 0.0);
					if (qwe > maxq) {
						maxq = qwe;
					}
					if (maxq > Math.PI / 180.0) {
						ter++;
						computeSteinerPoint(connect[i][0], connect[i][1],
								connect[i][2]);
						addSteinerNodeToGraph(i, connect[i][0], connect[i][1],
								connect[i][2]);
						mst();
						hantei();
					}// maxq>pi/180.0
				}// cn[it1]==2
			}// it1
			if (ter == 0) {
				tlabel = 1;
			}
		}// tlabel
	}// mrlonely

	/**
	 *
	 */
	private void hantei() {
		double treeLen = totalTreeLength();
		// if the current tree length is less than the previous minimum tree,
		// then copy the current tree to the min tree
		if (treeLen < minTreeLen) {
			// the current tree is now the minimum tree
			numBetterTreesFound++;
			minTreeLen = treeLen;
			synchronized(this) {
				// copy all the x,y coords and whether they're steiner min nodes
				this.minPoints = Collections.synchronizedList(new ArrayList<Point>());
				for(Point p : this.points) {
					this.minPoints.add(p.copy());
				}

				// update the number of nodes in SMT
				numNodesInSMT = this.points.size();

				// copy over all the edges
				this.minEdges = Collections.synchronizedList(new ArrayList<Edge>());
				for (int i = 0; i < this.points.size() - 1; i++) {
					this.minEdges.add(this.edges.get(i).copy());
				}
			}
			// repaint this em effin' tree
			repaint();
		}
	}// hantei

	public void run() {

		if (runMakeFair == true)
		{
			makeFair();
		}
		else
		{

			numBetterTreesFound = 0;

			// sort the points by the x-coordinate.
			//GraphUtils.sort(x1, y1, w1, numNonSteinerNodes);
			if(this.originalPoints.size() > 0) {
				this.points.clear();
				this.points.addAll(originalPoints);
			}
			else {
				GraphUtils.sortByX(points, false);
				originalPoints.addAll(points);
			}

			// copy all the points to the secondary points
			for (int k = 0; k < numNonSteinerNodes; k++) {
				this.copyPointToSecondary(k);
				/*
			x2(k, x(k));
			y2(k, y(k));
				 */
			}

			convexHull();

			minTreeLen = Double.MAX_VALUE;
			numIterations = 0;
			isRunning = true;
			while(isRunning) {
				numIterations++;

				repaintIfNecessary();

				//int steinerCount = 0;
				// this.points.size() = numNonSteinerNodes;
				points.clear();
				points.addAll(originalPoints);
				for (int k = 0; k < numNonSteinerNodes - 2; k++) {
					if (Math.random() < 0.5) {
						int chrule = 0;
						int br2;
						double randX = 0.0, randY = 0.0;
						while (chrule == 0) {
							randX = Math.random() * (width - 30) + 15;
							randY = Math.random() * (height - 30) + 15;
							br2 = 0;
							for (int h = 0; h < coch; h++) {
								if (randY > aa[h] * randX + bb[h]) {
									br2 = 1;
								}
							}
							for (int h = 0; h < coch2; h++) {
								if (randY < aa2[h] * randX + bb2[h]) {
									br2 = 1;
								}
							}
							if (br2 == 0) {
								chrule = 1;
							}
						}// chrule

						Point p = new Point(randX, randY);
						p.steiner = 1;
						this.points.add(p);
						//x(numNonSteinerNodes + steinerCount, x1kouho);
						//y(numNonSteinerNodes + steinerCount, y1kouho);
						//setSteiner(numNonSteinerNodes + steinerCount, 1);

						//steinerCount++;
						//this.points.size() = numNonSteinerNodes + steinerCount;
					}
				}

				mst();
				hantei();
				bohe();
				mst();
				bohe();
				mst();
				hantei();
				mrlonely();
				tuika();
			}
		}// next;;
	}// run

	public void update(Graphics g) {
		paint(g);
	}

	/* Given a set of points V, a number of maximum relay nodes k and a desired fairness measure alpha 
	 * Find a connected graph that achieves fairness >= alpha with additional nodes<=k
	 */
	public void makeFair()
	{
		//NOTE: adjust the constants (MAX_RELAY_NODES and TARGET_STDEV in Constants.java)
		
		//overhead of our data structure - construct the list of neighbors
		if (runOnce)
		{
			makeNeighborList();
			runOnce = false;
		}
		boolean showPrompt = true;  //show a prompt when target StDev is reached
		
		//do this loop while we still have resources to add
		//for (int i = (numNodesInSMT - numNonSteinerNodes); i < Constants.MAX_RELAY_NODES; i++)
		for (int i = 0; i < Constants.MAX_RELAY_NODES; i++)
		{			
			// do the while loop until the system has stabilized
			double stdev = Double.POSITIVE_INFINITY;
			while (stdev > getStandardDevOfPCR())
			{
				stdev = getStandardDevOfPCR();
				moveSteinerNodes_geometric();			
				try {
					repaint();
					Thread.sleep(sleepTime * 2); 
				} catch (Exception ignored) { }
			}
			log.info("Done with moveSteinerNodes_geometric.");

			addFairSteinerNode();
			log.info("Done with addFairSteinerNode. Added "+i+"/"+Constants.MAX_RELAY_NODES + " relay nodes.");

			moveSteinerNodes();
			log.info("Done with moveSteinerNodes.");
			
			//check if the target StDev is reached
			if (getStandardDevOfPCR()< Constants.TARGET_STDEV && showPrompt)
			{			
				int n = JOptionPane.showConfirmDialog(
						null,
						"Target Standard Deviation of "+ Constants.TARGET_STDEV +" reached.\n" +
						"The actual StDev is " + getStandardDevOfPCR() +" using " + i + " additional relay nodes\n\n"+
						"Do you want to continue and ignore future messages?\n",
						"Target STDEV reached",
						JOptionPane.YES_NO_OPTION);
				if (n == 1) //if no is pressed;
				{
					return;
				}
				showPrompt = false;			
			}
			
			//print some stats - this probably should be replace by writing some
			//stats to a file (i.e. for each k value)
			printStatistics();
		}
		
		repaint();
		printStatistics();
		log.info("Done making the graph fair. I did my best with the constraints given.");
	}

	/**
	 * Attempts to make the graph more fair
	 */
	public void makeFair_old()
	{			
		log.info("starting makeFair()");
		boolean showPrompt = true;

		if (runOnce) //only construct the neighbor list once (should be moved somewhere else)
		{
			makeNeighborList();
			runOnce = false;
		}

		int newNodeCount = 0;		
		while (newNodeCount < Constants.MAX_RELAY_NODES)
		{		
			double stdev = getStandardDevOfPCR();
			int converge = 0;

			while (stdev > Constants.TARGET_STDEV && converge <= Constants.CONVERGENCE_CUTOFF)
			{
				moveSteinerNodes_geometric();
				double temp = getStandardDevOfPCR();
				if (Math.abs(temp-stdev) < Constants.CONVERGENCE_THRESHOLD)
				{
					converge++;
				}
				else
				{
					//converge = 0;
				}
				stdev = temp;							
			}
			printStatistics();
			if(stdev <= Constants.TARGET_STDEV)
			{
				System.out.println("Target STDEV (<"+ Constants.TARGET_STDEV +") reached! The STDEV = " + stdev);

				if (showPrompt)
				{
					int n = JOptionPane.showConfirmDialog(
							null,
							"Target Standard Deviation Reached.\n" +
							"Do you want to continue and ignore future messages?\n" +
							"Note: make sure that all edges are in range (i.e. are green)",
							"Target STDEV reached",
							JOptionPane.YES_NO_OPTION);
					if (n == 1) //if no is pressed;
					{
						return;
					}
					showPrompt = false;
				}

			}	
			else
				System.out.println("Convergence occured  with a STDEV = " + stdev);

			log.info("Done moving Nodes");

			// Add new nodes
			if (newNodeCount < Constants.MAX_RELAY_NODES)
			{
				newNodeCount++;
				log.info("Adding new Node " + newNodeCount + "/"+Constants.MAX_RELAY_NODES);
				addFairSteinerNode();				
			}
			printStatistics();
			converge = 0;
			while (stdev > Constants.TARGET_STDEV && converge <= Constants.CONVERGENCE_CUTOFF)
			{				
				moveSteinerNodes_geometric();
				double temp = getStandardDevOfPCR();
				if (Math.abs(temp-stdev) < Constants.CONVERGENCE_THRESHOLD)
				{
					converge++;
				}
				else
				{
					//converge = 0;
				}
				stdev = temp;							
			}
			repaint();
		}
		log.info("Done with makeFair()");

		repaint();		
	}



	/**
	 * only add one steiner node to the graph and see if it
	 * maximizes the fairness
	 * WOW, this is probably the ugliest code I've ever written.
	 */
	private void addFairSteinerNode()
	{
		double max = Double.NEGATIVE_INFINITY;		
		Point maxP = null;
		int maxPIndex = -1;
		for (int i = 0; i < minPoints.size(); i++)
		{
			double temp = minPoints.get(i).getPCR();
			if (temp > max)
			{
				max = temp;
				maxP = minPoints.get(i);
				maxPIndex = i;
			}
		}		
		Point neighbor = maxP.getFarthestNeighbor();

		//find the edge to replace
		int removeIndex = -1;
		int neighborIndex = -1;
		for (int i = 0; i < minEdges.size(); i++)
		{
			if (maxP.isSameAs(this.minPoints.get(minEdges.get(i).index1)) && 
					neighbor.isSameAs(this.minPoints.get(minEdges.get(i).index2)))
			{
				removeIndex = i;
				neighborIndex = minEdges.get(i).index2;
			}
			if (maxP.isSameAs(this.minPoints.get(minEdges.get(i).index2)) && 
					neighbor.isSameAs(this.minPoints.get(minEdges.get(i).index1)))
			{
				removeIndex = i;
				neighborIndex = minEdges.get(i).index1;
			}			
		}
		if (removeIndex == -1) log.error("Remove Index problem");


		minEdges.remove(removeIndex);

		Point newP = new Point((maxP.x + neighbor.x)/2, (maxP.y + neighbor.y)/2);
		//newP.generateW();
		newP.steiner = 1;
		newP.neighbors.add(maxP);
		newP.neighbors.add(neighbor);

		minPoints.add(newP);
		numNodesInSMT++;

		int newPIndex = -1;
		for (int i = 0; i < minPoints.size(); i++)
		{
			if (newP.isSameAs(minPoints.get(i))) newPIndex = i;
		}

		if (newPIndex == -1) log.error("Index problem");

		Edge e1 = new Edge(maxP,newP);
		e1.index1 = maxPIndex;
		e1.index2 = newPIndex;
		Edge e2 = new Edge(neighbor,newP);
		e2.index1 = neighborIndex;
		e2.index2 = newPIndex;
		minEdges.add(e1);
		minEdges.add(e2);

		// remove the old neighbors
		maxP.neighbors.clear();
		neighbor.neighbors.clear();
		for (int i = 0 ; i < minEdges.size(); i++)
		{
			if (maxP.isSameAs(this.minPoints.get(minEdges.get(i).index1)))
			{
				maxP.neighbors.add(this.minPoints.get(minEdges.get(i).index2));
			}
			if (maxP.isSameAs(this.minPoints.get(minEdges.get(i).index2)))
			{
				maxP.neighbors.add(this.minPoints.get(minEdges.get(i).index1));
			}
			if (neighbor.isSameAs(this.minPoints.get(minEdges.get(i).index1)))
			{
				neighbor.neighbors.add(this.minPoints.get(minEdges.get(i).index2));
			}
			if (neighbor.isSameAs(this.minPoints.get(minEdges.get(i).index2)))
			{
				neighbor.neighbors.add(this.minPoints.get(minEdges.get(i).index1));
			}			
		}

	}

	private void printStatistics()
	{
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
		System.out.println("STDEV     = " + getStandardDevOfPCR());		
		System.out.println("MAX POWER = " + max);
		System.out.println("MIN POWER = " + min);
		System.out.println("AVG POWER = " + total / numPoints);
		System.out.println("TOTAL POW = " + total);
		System.out.println("------------------------------------");
	}


	/**
	 *  generates the neighbor list for each point
	 */
	private void makeNeighborList()
	{
		for (int i = 0; i < this.minEdges.size(); i++)
		{
			Point p1 = this.minPoints.get(this.minEdges.get(i).index1);
			Point p2 = this.minPoints.get(this.minEdges.get(i).index2);
			p1.neighbors.add(p2);
			p2.neighbors.add(p1);
		}
	}

	/**
	 * Moves true Steiner Nodes (i.e. movable nodes that have exactly 3 neighbors) by assuming that
	 * the 3 neighbors are points on a circle. The Steiner Node is then moved to the center of this
	 * circle. Since the 3 points can be organized in 3! different ways (i.e. p1,p2,p3; p1,p3,p2; 
	 * ....p3,p2,p1), the arrangement that minimizes the Steiner node's PCR (or the radius of the
	 * circle) is chosen. 
	 */
	private void moveSteinerNodes_geometric()
	{
		List<Point> list = new ArrayList<Point>();

		for (int i = 0; i < numNodesInSMT; i++)
		{
			Point p1 = this.minPoints.get(i);
			if(p1.isSteiner())
				list.add(p1);
		}

		for (int i = 0; i < list.size(); i++)
		{
			Point steiner = (Point) list.get(i);

			//only do this on real Steiner Nodes (i.e. movable nodes with 3 neighbors)
			if (steiner.neighbors.size() != 3)
			{
				//log.info("Cannot use moveSteinerNodes_geometric - not enough or too many neighbors");
			}
			else
			{
				Point temp = new Point(steiner.x,steiner.y);
				Point p1 = steiner.neighbors.get(0);
				Point p2 = steiner.neighbors.get(1);
				Point p3 = steiner.neighbors.get(2);

				double currentPCR = steiner.getPCR(); //use the lowest PCR			
				//##########################################				
				Point center = GraphUtils.center(p1,p2,p3);
				steiner.x = center.x;
				steiner.y = center.y;
				if (currentPCR > steiner.getPCR())
				{
					temp.x = center.x;
					temp.y = center.y;
					currentPCR = steiner.getPCR();
				}				
				//##########################################
				center = GraphUtils.center(p1,p3,p2);
				steiner.x = center.x;
				steiner.y = center.y;
				if (currentPCR > steiner.getPCR())
				{
					temp.x = center.x;
					temp.y = center.y;
					currentPCR = steiner.getPCR();
				}
				//##########################################
				center = GraphUtils.center(p2,p1,p3);
				steiner.x = center.x;
				steiner.y = center.y;
				if (currentPCR > steiner.getPCR())
				{
					temp.x = center.x;
					temp.y = center.y;
					currentPCR = steiner.getPCR();
				}
				//##########################################
				center = GraphUtils.center(p2,p3,p1);
				steiner.x = center.x;
				steiner.y = center.y;
				if (currentPCR > steiner.getPCR())
				{
					temp.x = center.x;
					temp.y = center.y;
					currentPCR = steiner.getPCR();
				}
				//##########################################
				center = GraphUtils.center(p3,p1,p2);
				steiner.x = center.x;
				steiner.y = center.y;
				if (currentPCR > steiner.getPCR())
				{
					temp.x = center.x;
					temp.y = center.y;
					currentPCR = steiner.getPCR();
				}
				//##########################################
				center = GraphUtils.center(p3,p2,p1);
				steiner.x = center.x;
				steiner.y = center.y;
				if (currentPCR > steiner.getPCR())
				{
					temp.x = center.x;
					temp.y = center.y;
					currentPCR = steiner.getPCR();
				}

				// set the best temporary point that minimizes PCR as
				//the new steiner point
				steiner.x = temp.x;
				steiner.y = temp.y;				
			}
		}
	}




	private void moveSteinerNodes_old()
	{
		MaxHeap heap = new MaxHeap();
		for (int i = 0; i < numNodesInSMT; i++) {
			Point p1 = this.minPoints.get(i);
			if(p1.isSteiner())
				heap.add(p1);
		}

		double prevSTDEV = Double.POSITIVE_INFINITY;
		double currentSTDEV = getStandardDevOfPCR();

		while (prevSTDEV > currentSTDEV)
		{
			log.info(prevSTDEV + "::::::" + currentSTDEV);
			prevSTDEV = currentSTDEV;

			Point source = (Point) heap.extractMax();			

			double threshold = 5;
			boolean repeat = true;
			log.info("outside");

			int count = 0;
			while (repeat)
			{
				if (count > 9)
				{
					count++;
				}
				repeat = false;
				Point target = source.getFarthestNeighbor();

				log.info("moving");
				moveTowards(source, target);
				Point lowestPCRNeighbor = source.getLowestPCRNeighbor();

				moveAway(source, lowestPCRNeighbor);			
				try
				{
					repaint();
					Thread.sleep(sleepTime);
				}
				catch (Exception ignored) { }

				List<Double> distances = new ArrayList<Double>();
				for (int i = 0; i < source.neighbors.size(); i++)
				{
					Point n = source.neighbors.get(i);
					distances.add(GraphUtils.euclideanDistance(n.x, n.y, source.x, source.y));
				}

				for (int i = 0; i < distances.size(); i++)
				{
					for (int j = 0; j < distances.size(); j++)
					{
						if (i != j)
						{
							if (Math.abs(distances.get(i) - distances.get(j)) > threshold)
							{
								log.info("abs: "+ Math.abs(distances.get(i) - distances.get(j)));
								repeat = true;
							}							
						}
					}					
				}				
			}

			heap.add(source);
			heap.rebuildHeap();

			currentSTDEV = getStandardDevOfPCR();
		}

	}


	/**
	 *  Added by Fatih Senel
	 */
	private void moveSteinerNodes()
	{
		MaxHeap heap = new MaxHeap();
		for (int i = 0; i < numNodesInSMT; i++)
		{
			Point p1 = this.minPoints.get(i);
			if(p1.isSteiner()) heap.add(p1);
		}

		for (int i = 0; i < heap.size(); i++) {
			Point source = (Point) heap.extractMax();
			//if(source.getPCR() > Constants.MAX_PCR_ALLOWED)
			{
				Point target = source.getFarthestNeighbor();                
				moveTowards(source, target);
				
				//Point lowestPCRNeighbor = source.getLowestPCRNeighbor();
				//moveAway(source, lowestPCRNeighbor);
				try {
					repaint();
					Thread.sleep(sleepTime);
				} catch (Exception ignored) { }
			}
			heap.rebuildHeap();
		}
	}

	/**
	 *  moves a Steiner Point away from the neighbor with the lowest PCR
	 */
	private void moveAway(Point source, Point target) {

		Point2D s = new Point2D.Double();
		Point2D t = new Point2D.Double();
		Point2D previousLocation = new Point2D.Double();
		s.setLocation(source.x, source.y);
		previousLocation.setLocation(source.x, source.y);
		t.setLocation(target.x, target.y);

		double step_size = 1;

		int check = 0;
		//while (target.getPCR() < source.getAVGPCR())
		{
			check ++;
			//if (check > 100) break;
			Point2D currentLocation = GraphUtils.getCoordinates(source.x, source.y, target.x, target.y, -step_size);
			source.x = currentLocation.getX();
			source.y = currentLocation.getY();
		}
	}



	/**
	 * moves a Steiner Point towards the furthest neighbor
	 */
	private void moveTowards(Point source, Point target) {

		Point2D s = new Point2D.Double();
		Point2D t = new Point2D.Double();
		Point2D previousLocation = new Point2D.Double();
		s.setLocation(source.x, source.y);
		previousLocation.setLocation(source.x, source.y);
		t.setLocation(target.x, target.y);

		// in each step move 1 unit forward towards the target
		double step_size = 1;
		//previous PCR : Stop moving when previousPCR is better than current one.
		double previousPCR = source.getPCR();

		int distance = (int) Math.floor(GraphUtils.euclideanDistance(source.x, source.y, target.x, target.y));
		for (int i = 0; i < distance; i += step_size)
		{			
			Point2D currentLocation = GraphUtils.getCoordinates(source.x, source.y, target.x, target.y, step_size);
			source.x = currentLocation.getX();
			source.y = currentLocation.getY();
			double currentPCR = source.getPCR();

			if (currentPCR < previousPCR) {
				previousPCR = currentPCR;
				previousLocation = currentLocation;
			} else {
				//move backwards
				source.x = previousLocation.getX();
				source.y = previousLocation.getY();
			}
		}
	}

	private double getStandardDevOfPCR() {
		double mean = getMeanPCR();
		double variance = 0;
		int counter = 0;
		for (int i = 0; i < numNodesInSMT; i++) {
			Point p1 = this.minPoints.get(i);
			//if(p1.isSteiner())
			{
				variance += Math.pow((p1.getPCR()-mean),2);
				counter++;
			}
		}
		return Math.sqrt(variance/counter);
	}

	private double getMeanPCR(){
		double avg = 0;
		int counter = 0;
		for (int i = 0; i < numNodesInSMT; i++) {
			Point p1 = this.minPoints.get(i);
			//if(p1.isSteiner())
			{
				avg += p1.getPCR();
				counter++;
			}
		}
		return avg/counter;
	}



	/**
	 * based on some stupid logic that only marginally makes sense.
	 */
	private void repaintIfNecessary() {
		int timlabel = 0;
		if (numNonSteinerNodes > 20) {
			if (numIterations % 3 == 0) {
				timlabel = 1;
			}
		}
		if (numNonSteinerNodes > 9 && numNonSteinerNodes < 21) {
			if (numIterations % 30 == 0) {
				timlabel = 1;
			}
		}
		if (numNonSteinerNodes < 10) {
			if (numIterations % 100 == 0) {
				timlabel = 1;
			}
		}
		if (timlabel == 1) {
			repaint();
		}
	}

	public void paint(java.awt.Graphics g) {
		// set and fill in the background color.
		g.setColor(bgColor);
		g.fillRect(1, 1, width, height);

		// had concurrency issues, so synchronize the entire shebang here
		synchronized(this) {
			/*
			// draw the lines of each edge
			// note that this will not actually show up on the canvas since the bgcolor is still the active color.
			for (int i = 0; i < this.points.size() - 1; i++) {
				g.drawLine((int) x(e1(i)), (int) y(e1(i)), (int) x(e2(i)),
						(int) y(e2(i)));
			}
			 */

			// draw all of the nodes.  if it's a steiner node, make it a square; otherwise, a circle.
			g.setColor(nodeColor);
			for (int i = 0; i < numNodesInSMT; i++) {
				Point p = this.minPoints.get(i);
				if (p.steiner == 0) {
					g.setColor(nodeColor);
					g.fillOval(p.xInt() - 3, p.yInt() - 3, 6, 6);
				} else {
					g.setColor(steinerNodeColor);
					//g.drawRect(p.xInt() - 5, p.yInt() - 5, 10, 10);
					g.fillOval(p.xInt() - 3, p.yInt() - 3, 6, 6);
				}
			}

			// draw the edges connecting all the minimum nodes.
			g.setColor(edgeColor);
			for (int i = 0; i < minEdges.size(); i++) {
				//Point p1 = this.minPoints.get(e1min[i]);
				Point p1 = this.minPoints.get(this.minEdges.get(i).index1);
				Point p2 = this.minPoints.get(this.minEdges.get(i).index2);
				if (GraphUtils.euclideanDistance(p1.x, p1.y, p2.x, p2.y) <= Constants.TRANSMISSION_RANGE)
				{
					g.setColor(edgeColor);
					g.drawLine(p1.xInt(), p1.yInt(), p2.xInt(), p2.yInt());
				}
				else
				{
					g.setColor(Color.gray);
					g.drawLine(p1.xInt(), p1.yInt(), p2.xInt(), p2.yInt());
				}
			}
		}

		drawKey(g);
	}// paint

	/**
	 * draw the key at the top of the canvas.
	 * @param g graphics
	 */
	private void drawKey(Graphics g) {
		g.setColor(textColor);

		g.drawString("Nodes", 15, 15); g.drawString("=", 65, 15); g.drawString("" + numNonSteinerNodes, 80, 15);
		g.drawString("Steiner", 15, 30);g.drawString("=", 65, 30); g.drawString("" + (numNodesInSMT - numNonSteinerNodes), 80, 30);
		g.drawString("Total", 15, 45);g.drawString("=", 65, 45); g.drawString("" + numNodesInSMT, 80, 45);


		Double max = Double.NEGATIVE_INFINITY;
		Double min = Double.POSITIVE_INFINITY;
		Double total = 0.0;
		int numPoints = 0;
		for (Point p : minPoints )
		{
			Double temp = p.getPCR();
			if (temp > max) max = temp;
			if (temp < min) min = temp;

			total += temp;
			numPoints++;
		}
		Double stdev = getStandardDevOfPCR();

		g.drawString("Power Usage", 115, 15);
		g.drawString("Stdev", 115, 30); g.drawString("=", 160, 30); g.drawString("" + stdev.intValue(), 175, 30);
		g.drawString("Max", 115, 45);g.drawString("=", 160, 45); g.drawString("" + max.intValue(), 175, 45);
		g.drawString("Min", 115, 60);g.drawString("=", 160, 60); g.drawString("" + min.intValue(), 175, 60);
		g.drawString("Avg", 115, 75);g.drawString("=", 160, 75); g.drawString("" + ((Double)(total / numPoints)).intValue(), 175, 75);
		g.drawString("Total", 115, 90);g.drawString("=", 160, 90); g.drawString("" + total.intValue(), 175, 90);


		g.drawString("Max Stdev", 225, 15); g.drawString("=", 290, 15); g.drawString("" + Constants.TARGET_STDEV.intValue(), 305, 15);
		g.drawString("Max PCR", 225, 30);g.drawString("=", 290, 30); g.drawString("" + Constants.MAX_PCR_ALLOWED.intValue(), 305, 30);
		g.drawString("Max New", 225, 45);g.drawString("=", 290, 45); g.drawString("" + Constants.MAX_RELAY_NODES, 305, 45);
		g.drawString("T Range", 225, 60);g.drawString("=", 290, 60); g.drawString("" + Constants.TRANSMISSION_RANGE.intValue(), 305, 60);
		g.drawString("Conv. Cut", 225, 75);g.drawString("=", 290, 75); g.drawString("" + Constants.CONVERGENCE_CUTOFF, 305, 75);
		g.drawString("Conv. Thr", 225, 90);g.drawString("=", 290, 90); g.drawString("" + Constants.CONVERGENCE_THRESHOLD, 305, 90);

		g.drawString("Min Tree Length =" + minTreeLen, 350, 15);
		g.drawString("Iteration = " + numIterations, 350, 30);
		g.drawString("Better Trees = " + numBetterTreesFound, 350, 45);

		g.drawLine(width -120, height - 20, width - 20, height - 20);
		g.drawLine(width -120, height - 20, width -120, height - 25);
		g.drawLine(width -70, height - 20, width -70, height - 25);
		g.drawLine( width - 20, height - 20,  width - 20, height - 25);
		g.drawString("0", width -130, height - 20);
		g.drawString("100", width - 19, height - 20);
	}

	/**
	 * @return
	 */
	public List<Point> getPoints() {
		return this.points;
	}

	/**
	 * start the execution
	 */
	public void start() {
		if(!this.hasRun) {
			if(this.isRunning || this.computationThread != null)
				throw new IllegalStateException("Can't start an SMT that's already running!");
			this.computationThread = new Thread(this);
			this.computationThread.start();
		}
		else {			
			log.warn("This FairSMT has already run.  We won't run it again. " +
					"There's a bug in the logic when it restarts that causes index out " +
					"of bound exceptions because the Japanese dude who originally wrote it did a crappy job. " +
			"I don't want to figure out how to fix it. :-)");
		}
	}

	public void startMakeFair()
	{
		runMakeFair = true;
		this.optimizationThread = new Thread(this);
		this.optimizationThread.start();		
	}

	/**
	 * stop the execution
	 */
	public void stop() {
		if(!this.isRunning || this.computationThread == null)
			throw new IllegalStateException("Can't stop an SMT that isn't running!");

		// reset counter variables
		this.coch = 0;
		this.coch2 = 0;
		// make this stop running.
		this.setRunning(false);
		this.hasRun = true;
		//this.secondaryPoints.clear();
		//this.edges.clear();
		this.computationThread = null;		
	}

	/**
	 * is the run loop going?
	 * @return
	 */
	public boolean isRunning() {
		return this.isRunning;
	}

	/**
	 * stop this SMT from running.
	 * @param isRunning
	 */
	private void setRunning(boolean isRunning) {
		this.isRunning = false;
	}

	/**
	 * @param args
	 */
	public static void main(String...args) {
		FairSMT ms = new FairSMT(800, 800, 100);
		JFrame f = new JFrame();
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		f.add(ms);
		f.pack();
		f.setVisible(true);

		ms.start();
	}
}
