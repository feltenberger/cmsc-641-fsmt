package edu.umbc.algorithms.old;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import edu.umbc.algorithms.fmst.Edge;
import edu.umbc.algorithms.fmst.GraphUtils;
import edu.umbc.algorithms.fmst.Point;

/**
 * @author dave
 *
 */
public class FairSMT extends JPanel implements Runnable {
	private static final long serialVersionUID = 1268750618801149582L;
	//private static final transient Logger log = Logger.getLogger(FairSMT.class);
	/**
	 * the thread that calculates all the steiner points
	 */
	private transient Thread computationThread;

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
	private int numNodesIncludingSteiner;
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
	 * all the points in the graph
	 */
	private List<Point> points = Collections.synchronizedList(new ArrayList<Point>());
	/**
	 * a secondary array for temporary storage
	 */
	private List<Point> secondaryPoints = Collections.synchronizedList(new ArrayList<Point>());
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
		steinerNodeColor = Color.cyan;
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
		//this.e2[index] = value;
	}

	/**
	 * @param index
	 * @param p
	 */
	private void copyPointToSecondary(int index) {
		Point p = this.points.get(index);
		if(this.secondaryPoints.size() == index)
			this.secondaryPoints.add(p.copy());
		else
			this.secondaryPoints.set(index, p.copy());
	}
	/**
	 * gets the x2 value
	 * @param index
	 * @return
	 */
	private double x2(int index) {
		return this.secondaryPoints.get(index).x;
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
	 * @param index
	 * @param w
	 */
	private void w(int index, double w) {
		points.get(index).w = w;
	}
	/**
	 * returns the weight
	 * @param index
	 * @return
	 */
	private double w(int index) {
		return points.get(index).w;
	}

	/**
	 * gets the y2 value
	 * @param index
	 * @return
	 */
	private double y2(int index) {
		return this.secondaryPoints.get(index).y;
	}

	/**
	 * gets whether this is a steiner node or not
	 * @param index
	 * @return
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
			GraphUtils.sortByWeight(points);

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
		for (int imst = 0; imst < numNodesIncludingSteiner - 1; imst++) {
			if (isSteiner(imst) < 2) {
				for (int jmst = imst + 1; jmst < numNodesIncludingSteiner; jmst++) {
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
		for (int i = 0; i < numNodesIncludingSteiner; i++) {
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
		for (int i = 0; i < numNodesIncludingSteiner - 1; i++) {
			total = total + GraphUtils.euclideanDistance(x(e1(i)), y(e1(i)), x(e2(i)), y(e2(i)));
			//total = total + euclideanDistance(x(e1(i)), y(e1[i]), x(e2[i]), y(e2[i]));
		}
		return total;
	}

	//
	private void bohe() {
		int numTotalNodesTemp;

		// copy arrays from x1 to x2
		for(int hi = numNonSteinerNodes; hi < numNodesIncludingSteiner; hi++) {
			copyPointToSecondary(hi);
			/*
			x2(hi, x(hi));
			y2(hi, y(hi));
			*/
			cn[hi] = cn[hi];
		}

		// 
		for(int hi = numNonSteinerNodes; hi < numNodesIncludingSteiner; hi++) {
			if (isSteiner(hi) == 1 && cn[hi] != 3) {
				setSteiner(hi, 2);
				//isSteinerNode[hi] = 2;
			}// sl[i]==1 && cn[i]<3
		}// wend hi

		numTotalNodesTemp = numNonSteinerNodes;
		for(int hi = numNonSteinerNodes; hi < numNodesIncludingSteiner; hi++) {
			if (isSteiner(hi) == 1) {
				//x1[numTotalNodesTemp] = x2(hi);
				x(numTotalNodesTemp, x2(hi));
				//y1[numTotalNodesTemp] = y2(hi);
				y(numTotalNodesTemp, y2(hi));
				cn[numTotalNodesTemp] = cn2[hi];
				numTotalNodesTemp++;
			}// sl[hi]==1
		}// hi
		numNodesIncludingSteiner = numTotalNodesTemp;
		for(int hi = numNonSteinerNodes; hi < numNodesIncludingSteiner; hi++) {
			setSteiner(hi, 1);
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
		double one_third_pi = Math.PI / 3.0;

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
		addSteinerNodeToGraph(numNodesIncludingSteiner, mij1, mij2, mij3);
		numNodesIncludingSteiner++;
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

		// calculate the distance from the first point to the midpoint
		distToMidpoint = GraphUtils.euclideanDistance(x(mij1), y(mij1), avgX, avgY);

		// adjust all the points' coordinates by the average
		x(mij1, x(mij1) - avgX);
		y(mij1, y(mij1) - avgY);
		x(mij2, x(mij2) - avgX);
		y(mij2, y(mij2) - avgY);
		x(mij3, x(mij3) - avgX);
		y(mij3, y(mij3) - avgY);

		// same for the steiner point...
		steinerX = steinerX - avgX;
		steinerY = steinerY - avgY;

		// calculate the slope of the line from the third point to the steiner point
		an = (y(mij3) - steinerY) / (x(mij3) - steinerX);
		// dunno what this is...
		bn = steinerY - an * steinerX;

		// if the steiner point's x coord is to the left of the third point's coord
		// do {whatever the hell it's doing} to add a new steiner node at the end of the array of points
		if (steinerX < x(mij3)) {
			double tmp =  (-an * bn + Math.pow(an * an * bn * bn
					- (bn * bn - distToMidpoint * distToMidpoint) * (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());
			x(steinerIndex, tmp);

			tmp = an * x(steinerIndex) + bn + 0.0001
			* Math.sin(Math.PI * 2.0 * Math.random());
			y(steinerIndex, tmp);

			//isSteinerNode[steinerIndex] = 1;
			setSteiner(steinerIndex, 1);
		} else {
			// otherwise, do something different to add the steiner node at the end of the array of points
			double tmp = (-an * bn - Math.pow(an * an * bn * bn
					- (bn * bn - distToMidpoint * distToMidpoint) * (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());
			x(steinerIndex, tmp);

			tmp = an * x(steinerIndex) + bn + 0.0001
			* Math.sin(Math.PI * 2.0 * Math.random());
			y(steinerIndex, tmp);

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
	}

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
					if (numNodesIncludingSteiner > 2 * numNonSteinerNodes - 2) {
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
						if (numNodesIncludingSteiner > 2 * numNonSteinerNodes - 2) {
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
			for (int i = numNonSteinerNodes; i < numNodesIncludingSteiner; i++) {
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
				numNodesInSMT = numNodesIncludingSteiner;
	
				// copy over all the edges
				this.minEdges = Collections.synchronizedList(new ArrayList<Edge>());
				for (int i = 0; i < numNodesIncludingSteiner - 1; i++) {
					this.minEdges.add(this.edges.get(i).copy());
				}
			}
			// repaint this em effin' tree
			repaint();
		}
	}// hantei

	public void run() {
		numBetterTreesFound = 0;

		// sort the points by the x-coordinate.
		//GraphUtils.sort(x1, y1, w1, numNonSteinerNodes);
		GraphUtils.sortByX(points);

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
		for (;;) {
			numIterations++;

			repaintIfNecessary();

			int steinerCount = 0;
			numNodesIncludingSteiner = numNonSteinerNodes;
			for (int k = 0; k < numNonSteinerNodes - 2; k++) {
				if (Math.random() < 0.5) {
					int chrule = 0;
					int br2;
					double x1kouho = 0.0, y1kouho = 0.0;
					while (chrule == 0) {
						x1kouho = Math.random() * (width - 30) + 15;
						y1kouho = Math.random() * (height - 30) + 15;
						br2 = 0;
						for (int h = 0; h < coch; h++) {
							if (y1kouho > aa[h] * x1kouho + bb[h]) {
								br2 = 1;
							}
						}
						for (int h = 0; h < coch2; h++) {
							if (y1kouho < aa2[h] * x1kouho + bb2[h]) {
								br2 = 1;
							}
						}
						if (br2 == 0) {
							chrule = 1;
						}
					}// chrule

					Point p = new Point(x1kouho, y1kouho);
					p.steiner = 1;
					this.points.add(p);
					//x(numNonSteinerNodes + steinerCount, x1kouho);
					//y(numNonSteinerNodes + steinerCount, y1kouho);
					//setSteiner(numNonSteinerNodes + steinerCount, 1);

					steinerCount++;
					numNodesIncludingSteiner = numNonSteinerNodes + steinerCount;
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
		}// next;;
	}// run

	public void update(Graphics g) {
		paint(g);
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

		// draw the lines of each edge
		// note that this will not actually show up on the canvas since the bgcolor is still the active color.
		for (int i = 0; i < numNodesIncludingSteiner - 1; i++) {
			g.drawLine((int) x(e1(i)), (int) y(e1(i)), (int) x(e2(i)),
					(int) y(e2(i)));
		}

		// had concurrency issues, so synchronize the entire shebang here
		synchronized(this) {
			// draw all of the nodes.  if it's a steiner node, make it a square; otherwise, a circle.
			g.setColor(nodeColor);
			for (int i = 0; i < numNodesInSMT; i++) {
				Point p = this.minPoints.get(i);
				if (p.steiner == 0) {
					g.setColor(nodeColor);
					g.fillOval(p.xInt() - 3, p.yInt() - 3, 6, 6);
				} else {
					g.setColor(steinerNodeColor);
					g.drawRect(p.xInt() - 5, p.yInt() - 5, 10, 10);
				}
			}
	
			// draw the edges connecting all the minimum nodes.
			g.setColor(edgeColor);
			for (int i = 0; i < numNodesInSMT - 1; i++) {
				//Point p1 = this.minPoints.get(e1min[i]);
				Point p1 = this.minPoints.get(this.minEdges.get(i).index1);
				Point p2 = this.minPoints.get(this.minEdges.get(i).index2);
				g.drawLine(p1.xInt(), p1.yInt(), p2.xInt(), p2.yInt());
			}
		}

		drawKey(g);
	}// paint

	/**
	 * draw the key at the top of the canvas.
	 * @param g
	 */
	private void drawKey(Graphics g) {
		g.setColor(textColor);
		g.drawString("N=" + numNonSteinerNodes, 15, 15);
		g.drawString("M=" + (numNodesInSMT - numNonSteinerNodes), 15, 30);
		g.drawString("N+M=" + numNodesInSMT, 65, 15);
		g.drawString("counter=" + numIterations, 65, 30);
		g.drawString("d=" + minTreeLen, 145, 15);
		g.drawString("counter2=" + numBetterTreesFound, 145, 30);
		g.drawLine(380, 15, 480, 15);
		g.drawLine(380, 15, 380, 10);
		g.drawLine(430, 15, 430, 10);
		g.drawLine(480, 15, 480, 10);
		g.drawString("0", 370, 15);
		g.drawString("100", 481, 15);
	}

	/**
	 * start the background thread
	 */
	public void start() {
		if (computationThread == null) {
			computationThread = new Thread(this);
			computationThread.start();
		}
	}

	/**
	 * stop the background thread
	 */
	@SuppressWarnings("all")
	public void stop() {
		if (computationThread != null) {
			computationThread.stop();
			computationThread = null;
		}
	}

	/**
	 * @param args
	 */
	public static void main(String...args) {
		FairSMT ms = new FairSMT(800, 800, 100);
		JFrame f = new JFrame();
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		ms.start();

		f.add(ms);
		f.pack();
		f.setVisible(true);
	}

	/**
	 * @return
	 */
	public List<Point> getPoints() {
		return this.points;
	}
}
