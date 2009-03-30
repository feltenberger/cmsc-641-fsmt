package edu.umbc.algorithms.old;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

import javax.swing.JPanel;

/**
 * @author dave
 *
 */
public class MinSteiner extends JPanel {
	//private static transient final Logger log = Logger.getLogger(MinSteiner.class);
	private static final long serialVersionUID = 1268750618801149582L;
	/**
	 * the thread that calculates all the steiner points
	private transient Thread computationThread;
	 */

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
	private int height = 300;
	/**
	 * the width of the applet
	 */
	private int width = 300;
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
	 * the x-coordinate of each point at item i
	 * corresponding index in y1 for y-coord.
	 */
	private double x1[] = new double[200];
	/**
	 * the y-coord of each point at item i;
	 * corresponding index in x1 for x-coord.
	 */
	private double y1[] = new double[200];

	/**
	 * the x coordinate of the min tree's points
	 */
	private double xmin[] = new double[200];
	/**
	 * the y coordinate of the min tree's points
	 */
	private double ymin[] = new double[200];

	private double w1[] = new double[200];
	private double x2[] = new double[200];
	private double y2[] = new double[200];
	private double aa[] = new double[200];
	private double bb[] = new double[200];
	private double aa2[] = new double[200];
	private double bb2[] = new double[200];

	/**
	 * is the point at index i a Steiner point?  1=yes, 0=no
	 */
	private int isSteinerNode[] = new int[200];
	/**
	 * is the point at index i a Steiner _minimum_ point?  1=yes, 0=no.
	 */
	private int isSteinerMinNode[] = new int[200];

	/**
	 * contains the indices of the first node in the edge.
	 * e.g.,
	 *     element e1[10]=2
	 * means that the first point in edge 10 is at the 
	 * x-coord  x1[2]  and the y-coord of  y1[2]
	 */
	private int e1[] = new int[200];
	/**
	 * the same as e1, only for the second node attached
	 * to the edge.
	 */
	private int e2[] = new int[200];

	/**
	 * index of first element in edges of minimum tree
	 */
	private int e1min[] = new int[200];
	/**
	 * index of second element in edges of minimum tree
	 */
	private int e2min[] = new int[200];

	private int cn[] = new int[2000];
	private int cn2[] = new int[200];
	private int connect[][] = new int[200][10];

	/**
	 * 
	 */
	public MinSteiner() {
		init();
	}
	/**
	 * @param width
	 * @param height
	 * @param numNonSteinerNodes
	 */
	public MinSteiner(int width, int height, int numNonSteinerNodes) {
		//this.width = width;
		//this.height = height;
		//this.numNonSteinerNodes = numNonSteinerNodes;
		this.width = 300;
		this.height = 300;
		this.numNonSteinerNodes = 3;
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
		// generate random steiner nodes
		double half_pi = Math.PI / 2;
		// initialize the tree by adding numNodes nodes randomly on the canvas.
		// make the steinerArray all 0, i.e. no steiner nodes.
		/*
		for (int k = 0; k < numNonSteinerNodes; k++) {
			x1[k] = Math.random() * (width - 30) + 15;
			y1[k] = Math.random() * (height - 45) + 30;
			isSteinerNode[k] = 0;
			w1[k] = Math.asin((y1[k] - 0.0) / euclideanDistance(x1[k], y1[k], 0.0, 0.0)) + half_pi;
		}
		*/
		x1[0] = 100.0;
		y1[0] = 50.0;
		isSteinerNode[0] = 0;
		w1[0] = Math.asin((y1[0] - 0.0) / euclideanDistance(x1[0], y1[0], 0.0, 0.0)) + half_pi;

		x1[1] = 150.0;
		y1[1] = 50.0;
		isSteinerNode[1] = 0;
		w1[1] = Math.asin((y1[1] - 0.0) / euclideanDistance(x1[1], y1[1], 0.0, 0.0)) + half_pi;

		x1[2] = 125.0;
		y1[2] = 290.0;
		isSteinerNode[2] = 0;
		w1[2] = Math.asin((y1[2] - 0.0) / euclideanDistance(x1[2], y1[2], 0.0, 0.0)) + half_pi;

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
	private double euclideanDistance(double d1, double d2, double d3, double d4) {
		double dw;
		dw = Math.pow(Math.pow(d3 - d1, 2.0) + Math.pow(d4 - d2, 2.0), 0.5);
		return dw;
	}

	/**
	 * sort such that te1[0]<te1[1]<...<te1[NNh-1]
	 * @param te1
	 * @param te2
	 * @param te3
	 * @param NNh
	 */
	private void sort(double te1[], double te2[], double te3[], int NNh) {
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
	 * compute convex hull
	 */
	public void convexHull() {
		double xmax = x1[numNonSteinerNodes - 1];
		int orikaeshi = 0;
		coch = 0;
		coch2 = 0;
		double lto = 0.0;
		double xx, yy;

		double half_pi = Math.PI / 2;

		while (lto != x2[0]) {
			xx = x1[0];
			yy = y1[0];
			for (int i = 1; i < numNonSteinerNodes; i++) {
				w1[i] = Math
						.asin((y1[i] - yy) / euclideanDistance(x1[i], y1[i], xx, yy))
						+ half_pi;
				if (x1[i] > xx) {
					w1[i] = 3
							* Math.PI
							/ 2
							- Math.asin((y1[i] - yy)
									/ euclideanDistance(x1[i], y1[i], xx, yy));
				}
				w1[i] = w1[i] + Math.PI;
				if (w1[i] > 2 * Math.PI) {
					w1[i] = w1[i] - 2 * Math.PI;
				}
				if (orikaeshi == 1) {
					w1[i] = Math.asin((y1[i] - yy)
							/ euclideanDistance(x1[i], y1[i], xx, yy))
							+ half_pi;
					if (x1[i] > xx) {
						w1[i] = 3
								* Math.PI
								/ 2
								- Math.asin((y1[i] - yy)
										/ euclideanDistance(x1[i], y1[i], xx, yy));
					}
					if (w1[i] > 2 * Math.PI) {
						w1[i] = w1[i] - 2 * Math.PI;
					}
				}
			} // next i

			x1[0] = xx;
			y1[0] = yy;
			w1[0] = 100.0;

			// sort by the w-values, whatever they are...
			sort(w1, x1, y1, numNonSteinerNodes);
			lto = x1[0];
			if (orikaeshi == 0) {
				aa[coch] = (y1[0] - yy) / (x1[0] - xx);
				bb[coch] = y1[0] - aa[coch] * x1[0];
				coch++;
			} else {
				aa2[coch2] = (y1[0] - yy) / (x1[0] - xx);
				bb2[coch2] = y1[0] - aa2[coch2] * x1[0];
				coch2++;
			}
			if (x1[0] == xmax) {
				orikaeshi = 1;
			}
		}// wend

		for (int k = 0; k < numNonSteinerNodes; k++) {
			x1[k] = x2[k];
			y1[k] = y2[k];
		}
	}// ch

	/**
	 * compute min spanning tree
	 */
	void mst() {
		double ei[] = new double[20000];
		double ej[] = new double[20000];
		double edgeDistances[] = new double[20000];

		int ncv[] = new int[200];
		int cv[][] = new int[200][200];
		int koko = 0;
		int count = 0;
		for (int imst = 0; imst < numNodesIncludingSteiner - 1; imst++) {
			if (isSteinerNode[imst] < 2) {
				for (int jmst = imst + 1; jmst < numNodesIncludingSteiner; jmst++) {
					if (isSteinerNode[jmst] < 2) {
						ei[count] = imst + 0.0;
						ej[count] = jmst + 0.0;
						edgeDistances[count] = euclideanDistance(x1[imst], y1[imst], x1[jmst],
								y1[jmst]);
						count++;
					}// sl[j]
				}// next j
			}// sl[i]
		}// next i

		// sort by the edge distances (i think?)
		sort(edgeDistances, ei, ej, count);

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
				e1[koko] = (int) ei[i];
				e2[koko] = (int) ej[i];
				koko++;
				connect[e1[koko - 1]][cn[e1[koko - 1]]] = e2[koko - 1];
				cn[e1[koko - 1]]++;
				connect[e2[koko - 1]][cn[e2[koko - 1]]] = e1[koko - 1];
				cn[e2[koko - 1]]++;

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
			total = total
					+ euclideanDistance(x1[e1[i]], y1[e1[i]], x1[e2[i]],
							y1[e2[i]]);
		}
		return total;
	}

	//
	private void bohe() {
		int numTotalNodesTemp;

		// copy arrays from x1 to x2
		for(int hi = numNonSteinerNodes; hi < numNodesIncludingSteiner; hi++) {
			x2[hi] = x1[hi];
			y2[hi] = y1[hi];
			cn[hi] = cn[hi];
		}

		// 
		for(int hi = numNonSteinerNodes; hi < numNodesIncludingSteiner; hi++) {
			if (isSteinerNode[hi] == 1 && cn[hi] != 3) {
				isSteinerNode[hi] = 2;
			}// sl[i]==1 && cn[i]<3
		}// wend hi

		numTotalNodesTemp = numNonSteinerNodes;
		for(int hi = numNonSteinerNodes; hi < numNodesIncludingSteiner; hi++) {
			if (isSteinerNode[hi] == 1) {
				x1[numTotalNodesTemp] = x2[hi];
				y1[numTotalNodesTemp] = y2[hi];
				cn[numTotalNodesTemp] = cn2[hi];
				numTotalNodesTemp++;
			}// sl[hi]==1
		}// hi
		numNodesIncludingSteiner = numTotalNodesTemp;
		for(int hi = numNonSteinerNodes; hi < numNodesIncludingSteiner; hi++) {
			isSteinerNode[hi] = 1;
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

		ka = (y1[ij2] - y1[ij1]) / (x1[ij2] - x1[ij1]);
		se = y1[ij1] - ka * x1[ij1];
		xz = x1[ij2] - x1[ij1];
		yz = y1[ij2] - y1[ij1];
		if (y1[ij3] < ka * x1[ij3] + se && x1[ij1] > x1[ij2]) {
			steinerX = x1[ij1] + Math.cos(one_third_pi) * xz + Math.sin(one_third_pi) * yz;
			steinerY = y1[ij1] - Math.sin(one_third_pi) * xz + Math.cos(one_third_pi) * yz;
		}
		if (y1[ij3] > ka * x1[ij3] + se && x1[ij1] < x1[ij2]) {
			steinerX = x1[ij1] + Math.cos(one_third_pi) * xz + Math.sin(one_third_pi) * yz;
			steinerY = y1[ij1] - Math.sin(one_third_pi) * xz + Math.cos(one_third_pi) * yz;
		}
		if (y1[ij3] < ka * x1[ij3] + se && x1[ij1] < x1[ij2]) {
			steinerX = x1[ij1] + Math.cos(-one_third_pi) * xz + Math.sin(-one_third_pi) * yz;
			steinerY = y1[ij1] - Math.sin(-one_third_pi) * xz + Math.cos(-one_third_pi) * yz;
		}
		if (y1[ij3] > ka * x1[ij3] + se && x1[ij1] > x1[ij2]) {
			steinerX = x1[ij1] + Math.cos(-one_third_pi) * xz + Math.sin(-one_third_pi) * yz;
			steinerY = y1[ij1] - Math.sin(-one_third_pi) * xz + Math.cos(-one_third_pi) * yz;
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
		avgX = (x1[mij1] + x1[mij2] + steinerX) / 3.0;
		// average the first, second, and steiner points' y-coord
		avgY = (y1[mij1] + y1[mij2] + steinerY) / 3.0;
		//log.info("avgx = " + avgX + ", avgy = " + avgY);

		// calculate the distance from the first point to the midpoint
		distToMidpoint = euclideanDistance(x1[mij1], y1[mij1], avgX, avgY);
		//log.info("dist to midpoint: " + distToMidpoint);

		// adjust all the points' coordinates by the average
		x1[mij1] = x1[mij1] - avgX;
		y1[mij1] = y1[mij1] - avgY;
		x1[mij2] = x1[mij2] - avgX;
		y1[mij2] = y1[mij2] - avgY;
		x1[mij3] = x1[mij3] - avgX;
		y1[mij3] = y1[mij3] - avgY;

		//log.info("steinerX = " + steinerX + ", y = " + steinerY);
		// same for the steiner point...
		steinerX = steinerX - avgX;
		steinerY = steinerY - avgY;
		//log.info("steinerX = " + steinerX + ", y = " + steinerY);

		// calculate the slope of the line from the third point to the steiner point
		an = (y1[mij3] - steinerY) / (x1[mij3] - steinerX);
		// dunno what this is...
		bn = steinerY - an * steinerX;
		//log.info("an = " + an + ", bn = " + bn);

		// if the steiner point's x coord is to the left of the third point's coord
		// do {whatever the hell it's doing} to add a new steiner node at the end of the array of points
		if (steinerX < x1[mij3]) {
			x1[steinerIndex] = (-an * bn + Math.pow(an * an * bn * bn
					- (bn * bn - distToMidpoint * distToMidpoint) * (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());
			y1[steinerIndex] = an * x1[steinerIndex] + bn + 0.0001
					* Math.sin(Math.PI * 2.0 * Math.random());
			isSteinerNode[steinerIndex] = 1;
		} else {
			// otherwise, do something different to add the steiner node at the end of the array of points
			x1[steinerIndex] = (-an * bn - Math.pow(an * an * bn * bn
					- (bn * bn - distToMidpoint * distToMidpoint) * (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());
			if(x1[steinerIndex] > 999999.9999)
				System.out.println(x1[steinerIndex]);
			y1[steinerIndex] = an * x1[steinerIndex] + bn + 0.0001
					* Math.sin(Math.PI * 2.0 * Math.random());
			isSteinerNode[steinerIndex] = 1;
		}
		x1[mij1] = x1[mij1] + avgX;
		y1[mij1] = y1[mij1] + avgY;
		x1[mij2] = x1[mij2] + avgX;
		y1[mij2] = y1[mij2] + avgY;
		x1[mij3] = x1[mij3] + avgX;
		y1[mij3] = y1[mij3] + avgY;
		steinerX = steinerX + avgX;
		steinerY = steinerY + avgY;
		x1[steinerIndex] = x1[steinerIndex] + avgX;
		y1[steinerIndex] = y1[steinerIndex] + avgY;
		//if(num == 20)
		//	System.exit(0);
		//num++;
		//log.info("-----------------------------------");
	}// msbn2
	//private static int num = 0;

	private void tase1(int it1) {
		int qint;
		double x, x2, y, y2, qa, qb, qc, qw, minq;

		x = x1[connect[it1][0]];
		y = y1[connect[it1][0]];
		x2 = x1[connect[it1][1]];
		y2 = y1[connect[it1][1]];
		qa = euclideanDistance(x, y, x2, y2);
		qb = euclideanDistance(x, y, x1[it1], y1[it1]);
		qc = euclideanDistance(x2, y2, x1[it1], y1[it1]);
		qw = Math.acos((qb * qb + qc * qc - qa * qa) / (2.0 * qb * qc));
		minq = qw;
		qint = 1;
		x = x1[connect[it1][1]];
		y = y1[connect[it1][1]];
		x2 = x1[connect[it1][2]];
		y2 = y1[connect[it1][2]];
		qa = euclideanDistance(x, y, x2, y2);
		qb = euclideanDistance(x, y, x1[it1], y1[it1]);
		qc = euclideanDistance(x2, y2, x1[it1], y1[it1]);
		qw = Math.acos((qb * qb + qc * qc - qa * qa) / (2.0 * qb * qc));
		if (qw < minq) {
			minq = qw;
			qint = 2;
		}
		x = x1[connect[it1][2]];
		y = y1[connect[it1][2]];
		x2 = x1[connect[it1][0]];
		y2 = y1[connect[it1][0]];
		qa = euclideanDistance(x, y, x2, y2);
		qb = euclideanDistance(x, y, x1[it1], y1[it1]);
		qc = euclideanDistance(x2, y2, x1[it1], y1[it1]);
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
					xx = x1[connect[i][0]];
					yy = y1[connect[i][0]];
					xx2 = x1[connect[i][1]];
					yy2 = y1[connect[i][1]];
					qa = euclideanDistance(xx, yy, xx2, yy2);
					qb = euclideanDistance(xx, yy, x1[i], y1[i]);
					qc = euclideanDistance(xx2, yy2, x1[i], y1[i]);
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
					xx = x1[connect[i][0]];
					yy = y1[connect[i][0]];
					xx2 = x1[connect[i][1]];
					yy2 = y1[connect[i][1]];
					qa = euclideanDistance(xx, yy, xx2, yy2);
					qb = euclideanDistance(xx, yy, x1[i], y1[i]);
					qc = euclideanDistance(xx2, yy2, x1[i], y1[i]);
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = euclideanDistance(qw, 0.0, 2.0 * one_third_pi, 0.0);
					maxq = qwe;
					xx = x1[connect[i][1]];
					yy = y1[connect[i][1]];
					xx2 = x1[connect[i][2]];
					yy2 = y1[connect[i][2]];
					qa = euclideanDistance(xx, yy, xx2, yy2);
					qb = euclideanDistance(xx, yy, x1[i], y1[i]);
					qc = euclideanDistance(xx2, yy2, x1[i], y1[i]);
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = euclideanDistance(qw, 0.0, 2.0 * one_third_pi, 0.0);
					if (qwe > maxq) {
						maxq = qwe;
					}
					xx = x1[connect[i][2]];
					yy = y1[connect[i][2]];
					xx2 = x1[connect[i][0]];
					yy2 = y1[connect[i][0]];
					qa = euclideanDistance(xx, yy, xx2, yy2);
					qb = euclideanDistance(xx, yy, x1[i], y1[i]);
					qc = euclideanDistance(xx2, yy2, x1[i], y1[i]);
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = euclideanDistance(qw, 0.0, 2.0 * one_third_pi, 0.0);
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
			// copy all the x,y coords and whether they're steiner min nodes
			for (int i = 0; i < numNodesIncludingSteiner; i++) {
				xmin[i] = x1[i];
				ymin[i] = y1[i];
				isSteinerMinNode[i] = isSteinerNode[i];
			}

			// update the number of nodes in SMT
			numNodesInSMT = numNodesIncludingSteiner;

			// copy over all the edges
			for (int i = 0; i < numNodesIncludingSteiner - 1; i++) {
				e1min[i] = e1[i];
				e2min[i] = e2[i];
			}
			// repaint this em effin' tree
			repaint();
		}
	}// hantei

	public void run() {
		numBetterTreesFound = 0;

		// sort the points by the x-coordinate.
		sort(x1, y1, w1, numNonSteinerNodes);

		// 
		for (int k = 0; k < numNonSteinerNodes; k++) {
			x2[k] = x1[k];
			y2[k] = y1[k];
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
					x1[numNonSteinerNodes + steinerCount] = x1kouho;
					y1[numNonSteinerNodes + steinerCount] = y1kouho;
					isSteinerNode[numNonSteinerNodes + steinerCount] = 1;
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
			g.drawLine((int) x1[e1[i]], (int) y1[e1[i]], (int) x1[e2[i]],
					(int) y1[e2[i]]);
		}

		// draw all of the nodes.  if it's a steiner node, make it a square; otherwise, a circle.
		g.setColor(nodeColor);
		for (int i = 0; i < numNodesInSMT; i++) {
			if (isSteinerMinNode[i] == 0) {
				g.setColor(nodeColor);
				g.fillOval((int) xmin[i] - 3, (int) ymin[i] - 3, 6, 6);
			} else {
				g.setColor(steinerNodeColor);
				g.drawRect((int) xmin[i] - 5, (int) ymin[i] - 5, 10, 10);
			}
		}

		// draw the edges connecting all the minimum nodes.
		g.setColor(edgeColor);
		for (int i = 0; i < numNodesInSMT - 1; i++) {
			g.drawLine((int) xmin[e1min[i]], (int) ymin[e1min[i]],
					(int) xmin[e2min[i]], (int) ymin[e2min[i]]);
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
		/*
		if (computationThread == null) {
			computationThread = new Thread(this);
			computationThread.start();
		}
		*/
		run();
	}

	/**
	 * stop the background thread
	 */
	@SuppressWarnings("all")
	public void stop() {
		/*
		if (computationThread != null) {
			computationThread.stop();
			computationThread = null;
		}
		*/
	}
}
