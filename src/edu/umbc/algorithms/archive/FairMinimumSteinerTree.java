package edu.umbc.algorithms.archive;

import java.awt.Color;
import java.awt.Event;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JFrame;

/**
 * Code lifted from: http://www.nirarebakun.com/graph/emsteinercli.html
 * 
 * @author Dave
 * 
 */
public class FairMinimumSteinerTree extends JFrame {
	private static final long serialVersionUID = 1L;
	private boolean isRunning = true;
	private Color bgColor;
	private int tim, mouselabel, counter2;
	private int numPoints, height, width, ni, nj, minNN;
	private int orikaeshi, sco, NN, koko, coch, coch2;
	private double lto, xmax, dd2;
	private double tx, ty;
	private double xx, yy, mind;
	private Color colors[] = new Color[13];

	private Point p[] = new Point[2000];
	private int ncv[] = new int[200];
	private int cv[][] = new int[200][200];
	private double ei[] = new double[20000];
	private double ej[] = new double[20000];
	private double ed[] = new double[20000];

	// the x and y coordinates of the points
	private double x1[] = new double[200];
	private double y1[] = new double[200];
	
	private double xmin[] = new double[200];
	private double ymin[] = new double[200];
	private double w1[] = new double[200];
	private double x2[] = new double[200];
	private double y2[] = new double[200];
	private double aa[] = new double[200];
	private double bb[] = new double[200];
	private double aa2[] = new double[200];
	private double bb2[] = new double[200];

	// something to do with steiner points
	private int sl[] = new int[200];
	private int sl2[] = new int[200];

	// a steiner minimum point -- 1 for yes, 0 for no.
	private int slmin[] = new int[200];

	// edge 1 and 2  (these are the indexes to items in other arrays?)
	private int e1[] = new int[200];
	private int e2[] = new int[200];

	// the min edges
	private int e1min[] = new int[200];
	private int e2min[] = new int[200];

	private int cn[] = new int[2000];
	private int cn2[] = new int[200];
	private int connect[][] = new int[200][10];

	public static void main(String...args) {
		
	}

	/**
	 * Initialize the Applet.
	 * 
	 * @see java.applet.Applet#init()
	 */
	public void init() {
		bgColor = Color.black;
		colors[0] = Color.white;
		colors[1] = Color.green;
		colors[2] = Color.cyan;
		height = 350;
		width = 500;
		numPoints = 0;
		this.addMouseListener(new MouseListener() {

			public void mouseClicked(MouseEvent e) {
				handleMouseClick(e.getX(), e.getY());
			}
			public void mouseEntered(MouseEvent e) {
			}

			public void mouseExited(MouseEvent e) {
			}

			public void mousePressed(MouseEvent e) {
			}

			public void mouseReleased(MouseEvent e) {
			}
			
		}
		);
	}

	/**
	 * 
	 */
	public FairMinimumSteinerTree() {
		init();
	}

	/**
	 * Calculate the Eclidean distance between the points.
	 * 
	 * @param d1
	 * @param d2
	 * @param d3
	 * @param d4
	 * @return
	 */
	private double euclDist(double d1, double d2, double d3, double d4) {
		double dw;
		dw = Math.pow(Math.pow(d3 - d1, 2.0) + Math.pow(d4 - d2, 2.0), 0.5);
		return dw;
	}

	/**
	 * Sort such that te1[0] < te1[1] < ... te1[NNh-1]
	 * 
	 * @param te1
	 * @param te2
	 * @param te3
	 * @param NNh
	 */
	private void heapv(double te1[], double te2[], double te3[], int NNh) {
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
	}// heapv

	/**
	 * Compute the Convex Hull
	 */
	private void convexHull() {
		int i;
		xmax = x1[numPoints - 1];
		orikaeshi = 0;
		coch = 0;
		coch2 = 0;
		lto = 0.0;

		while (lto != x2[0]) {
			xx = x1[0];
			yy = y1[0];
			for (i = 1; i < numPoints; i++) {
				w1[i] = Math
						.asin((y1[i] - yy) / euclDist(x1[i], y1[i], xx, yy))
						+ Math.PI / 2;
				if (x1[i] > xx) {
					w1[i] = 3
							* Math.PI
							/ 2
							- Math.asin((y1[i] - yy)
									/ euclDist(x1[i], y1[i], xx, yy));
				}
				w1[i] = w1[i] + Math.PI;
				if (w1[i] > 2 * Math.PI) {
					w1[i] = w1[i] - 2 * Math.PI;
				}
				if (orikaeshi == 1) {
					w1[i] = Math.asin((y1[i] - yy)
							/ euclDist(x1[i], y1[i], xx, yy))
							+ Math.PI / 2;
					if (x1[i] > xx) {
						w1[i] = 3
								* Math.PI
								/ 2
								- Math.asin((y1[i] - yy)
										/ euclDist(x1[i], y1[i], xx, yy));
					}
					if (w1[i] > 2 * Math.PI) {
						w1[i] = w1[i] - 2 * Math.PI;
					}
				}
			}// next i
			x1[0] = xx;
			y1[0] = yy;
			w1[0] = 100.0;
			heapv(w1, x1, y1, numPoints);
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
	}// ch

	/**
	 * Compute the minimum spanning tree.
	 */
	private void mst() {
		int imst, jmst, co, i, j, k;
		int tyao, mstcou;
		tyao = 0;
		mstcou = 0;
		while (tyao == 0 && mstcou < 2) {
			tyao = 1;
			co = 0;
			for (imst = 0; imst < NN - 1; imst++) {
				if (sl[imst] < 2) {
					for (jmst = imst + 1; jmst < NN; jmst++) {
						if (sl[jmst] < 2) {
							ei[co] = imst + 0.0;
							ej[co] = jmst + 0.0;
							ed[co] = euclDist(x1[imst], y1[imst], x1[jmst],
									y1[jmst]);
							co++;
						}// sl[j]
					}// next j
				}// sl[i]
			}// next i

			heapv(ed, ei, ej, co);
			for (i = 0; i < NN; i++) {
				ncv[i] = 0;
			}
			for (i = 0; i < NN; i++) {
				cn[i] = 0;
			}
			koko = 0;
			for (i = 0; i < co; i++) {
				if (mouselabel == 1 && mstcou == 0) {
					tyao = 0;
					mstcou++;
					break;
				}
				if (mouselabel == 1 && mstcou == 1) {
					tyao = 0;
					mstcou++;
				}
				int label = 0;
				for (j = 0; j < ncv[(int) ei[i]]; j++) {
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

					ni = ncv[(int) ei[i]];
					nj = ncv[(int) ej[i]];

					cv[(int) ei[i]][ncv[(int) ei[i]]] = (int) ej[i];
					ncv[(int) ei[i]]++;
					for (j = 0; j < nj; j++) {
						cv[(int) ei[i]][ncv[(int) ei[i]]] = cv[(int) ej[i]][j];
						ncv[(int) ei[i]]++;
						cv[cv[(int) ej[i]][j]][ncv[cv[(int) ej[i]][j]]] = (int) ei[i];
						ncv[cv[(int) ej[i]][j]]++;
						for (k = 0; k < ni; k++) {
							cv[cv[(int) ej[i]][j]][ncv[cv[(int) ej[i]][j]]] = cv[(int) ei[i]][k];
							ncv[cv[(int) ej[i]][j]]++;
						}
					}// j

					cv[(int) ej[i]][ncv[(int) ej[i]]] = (int) ei[i];
					ncv[(int) ej[i]]++;
					for (j = 0; j < ni; j++) {
						cv[(int) ej[i]][ncv[(int) ej[i]]] = cv[(int) ei[i]][j];
						ncv[(int) ej[i]]++;
						cv[cv[(int) ei[i]][j]][ncv[cv[(int) ei[i]][j]]] = (int) ej[i];
						ncv[cv[(int) ei[i]][j]]++;
						for (k = 0; k < nj; k++) {
							cv[cv[(int) ei[i]][j]][ncv[cv[(int) ei[i]][j]]] = cv[(int) ej[i]][k];
							ncv[cv[(int) ei[i]][j]]++;
						}
					}// j
				}// label=0
			}// i
		}// wend tyao
	}// mst

	/**
	 * Calculate the length of the the total tree.
	 * 
	 * @return
	 */
	private double treeLength() {
		double distance;
		int i;
		distance = 0.0;
		for (i = 0; i < NN - 1; i++) {
			distance = distance
					+ euclDist(x1[e1[i]], y1[e1[i]], x1[e2[i]], y1[e2[i]]);
		}
		return distance;
	}

	private void bohe() {
		if (numPoints > 2) {
			int hi, shinNN;
			for (hi = numPoints; hi < NN; hi++) {
				x2[hi] = x1[hi];
				y2[hi] = y1[hi];
				sl2[hi] = sl[hi];
				cn[hi] = cn[hi];
			}
			for (hi = numPoints; hi < NN; hi++) {
				if (sl[hi] == 1 && cn[hi] != 3) {
					sl[hi] = 2;
				}// sl[i]==1 && cn[i]<3
			}// wend hi
			shinNN = numPoints;
			for (hi = numPoints; hi < NN; hi++) {
				if (sl[hi] == 1) {
					x1[shinNN] = x2[hi];
					y1[shinNN] = y2[hi];
					cn[shinNN] = cn2[hi];
					shinNN++;
				}// sl[hi]==1
			}// hi
			NN = shinNN;
			for (hi = numPoints; hi < NN; hi++) {
				sl[hi] = 1;
			}
		}
	}// bohe

	/**
	 * Compute (tx,ty) as Steiner point os ij1,ij2,ij3
	 * 
	 * @param ij1
	 * @param ij2
	 * @param ij3
	 */
	private void computeSteinerPoint(int ij1, int ij2, int ij3) {
		double ka, se, xz, yz;
		ka = (y1[ij2] - y1[ij1]) / (x1[ij2] - x1[ij1]);
		se = y1[ij1] - ka * x1[ij1];
		xz = x1[ij2] - x1[ij1];
		yz = y1[ij2] - y1[ij1];
		if (y1[ij3] < ka * x1[ij3] + se && x1[ij1] > x1[ij2]) {
			tx = x1[ij1] + Math.cos(Math.PI) * xz + Math.sin(Math.PI / 3.0) * yz;
			ty = y1[ij1] - Math.sin(Math.PI / 3.0) * xz + Math.cos(Math.PI / 3.0) * yz;
		}
		if (y1[ij3] > ka * x1[ij3] + se && x1[ij1] < x1[ij2]) {
			tx = x1[ij1] + Math.cos(Math.PI / 3.0) * xz + Math.sin(Math.PI / 3.0) * yz;
			ty = y1[ij1] - Math.sin(Math.PI / 3.0) * xz + Math.cos(Math.PI / 3.0) * yz;
		}
		if (y1[ij3] < ka * x1[ij3] + se && x1[ij1] < x1[ij2]) {
			tx = x1[ij1] + Math.cos(-Math.PI / 3.0) * xz + Math.sin(-Math.PI / 3.0) * yz;
			ty = y1[ij1] - Math.sin(-Math.PI / 3.0) * xz + Math.cos(-Math.PI / 3.0) * yz;
		}
		if (y1[ij3] > ka * x1[ij3] + se && x1[ij1] > x1[ij2]) {
			tx = x1[ij1] + Math.cos(-Math.PI / 3.0) * xz + Math.sin(-Math.PI / 3.0) * yz;
			ty = y1[ij1] - Math.sin(-Math.PI / 3.0) * xz + Math.cos(-Math.PI / 3.0) * yz;
		}
	}

	private void musoubana(int mij1, int mij2, int mij3) {
		double xc, yc, an, bn, rr;
		xc = (x1[mij1] + x1[mij2] + tx) / 3.0;
		yc = (y1[mij1] + y1[mij2] + ty) / 3.0;
		rr = euclDist(x1[mij1], y1[mij1], xc, yc);
		x1[mij1] = x1[mij1] - xc;
		y1[mij1] = y1[mij1] - yc;
		x1[mij2] = x1[mij2] - xc;
		y1[mij2] = y1[mij2] - yc;
		x1[mij3] = x1[mij3] - xc;
		y1[mij3] = y1[mij3] - yc;
		tx = tx - xc;
		ty = ty - yc;
		an = (y1[mij3] - ty) / (x1[mij3] - tx);
		bn = ty - an * tx;
		if (tx < x1[mij3]) {
			x1[NN] = (-an * bn + Math.pow(an * an * bn * bn - (bn * bn - rr * rr)
					* (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());
			y1[NN] = an * x1[NN] + bn + 0.0001
					* Math.sin(Math.PI * 2.0 * Math.random());
			sl[NN] = 1;
			NN++;
		} else {
			x1[NN] = (-an * bn - Math.pow(an * an * bn * bn - (bn * bn - rr * rr)
					* (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());
			y1[NN] = an * x1[NN] + bn + 0.0001
					* Math.sin(Math.PI * 2.0 * Math.random());
			sl[NN] = 1;
			NN++;
		}
		x1[mij1] = x1[mij1] + xc;
		y1[mij1] = y1[mij1] + yc;
		x1[mij2] = x1[mij2] + xc;
		y1[mij2] = y1[mij2] + yc;
		x1[mij3] = x1[mij3] + xc;
		y1[mij3] = y1[mij3] + yc;
		tx = tx + xc;
		ty = ty + yc;
		x1[NN - 1] = x1[NN - 1] + xc;
		y1[NN - 1] = y1[NN - 1] + yc;
	}// musoubana

	private void msbn2(int msbn, int mij1, int mij2, int mij3) {
		double xc, yc, an, bn, rr;
		xc = (x1[mij1] + x1[mij2] + tx) / 3.0;
		yc = (y1[mij1] + y1[mij2] + ty) / 3.0;
		rr = euclDist(x1[mij1], y1[mij1], xc, yc);
		x1[mij1] = x1[mij1] - xc;
		y1[mij1] = y1[mij1] - yc;
		x1[mij2] = x1[mij2] - xc;
		y1[mij2] = y1[mij2] - yc;
		x1[mij3] = x1[mij3] - xc;
		y1[mij3] = y1[mij3] - yc;
		tx = tx - xc;
		ty = ty - yc;
		an = (y1[mij3] - ty) / (x1[mij3] - tx);
		bn = ty - an * tx;
		if (tx < x1[mij3]) {
			x1[msbn] = (-an * bn + Math.pow(an * an * bn * bn - (bn * bn - rr * rr)
					* (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());
			y1[msbn] = an * x1[msbn] + bn + 0.0001
					* Math.sin(Math.PI * 2.0 * Math.random());
			sl[msbn] = 1;
		} else {
			x1[msbn] = (-an * bn - Math.pow(an * an * bn * bn - (bn * bn - rr * rr)
					* (1.0 + an * an), 0.5))
					/ (1.0 + an * an)
					+ 0.0001
					* Math.cos(Math.PI * 2.0 * Math.random());
			y1[msbn] = an * x1[msbn] + bn + 0.0001
					* Math.sin(Math.PI * 2.0 * Math.random());
			sl[msbn] = 1;
		}
		x1[mij1] = x1[mij1] + xc;
		y1[mij1] = y1[mij1] + yc;
		x1[mij2] = x1[mij2] + xc;
		y1[mij2] = y1[mij2] + yc;
		x1[mij3] = x1[mij3] + xc;
		y1[mij3] = y1[mij3] + yc;
		tx = tx + xc;
		ty = ty + yc;
		x1[msbn] = x1[msbn] + xc;
		y1[msbn] = y1[msbn] + yc;
	}// msbn2

	private void tase1(int it1) {
		int qint;
		double xx2, yy2, qa, qb, qc, qw, minq;
		xx = x1[connect[it1][0]];
		yy = y1[connect[it1][0]];
		xx2 = x1[connect[it1][1]];
		yy2 = y1[connect[it1][1]];
		qa = euclDist(xx, yy, xx2, yy2);
		qb = euclDist(xx, yy, x1[it1], y1[it1]);
		qc = euclDist(xx2, yy2, x1[it1], y1[it1]);
		qw = Math.acos((qb * qb + qc * qc - qa * qa) / (2.0 * qb * qc));
		minq = qw;
		qint = 1;
		xx = x1[connect[it1][1]];
		yy = y1[connect[it1][1]];
		xx2 = x1[connect[it1][2]];
		yy2 = y1[connect[it1][2]];
		qa = euclDist(xx, yy, xx2, yy2);
		qb = euclDist(xx, yy, x1[it1], y1[it1]);
		qc = euclDist(xx2, yy2, x1[it1], y1[it1]);
		qw = Math.acos((qb * qb + qc * qc - qa * qa) / (2.0 * qb * qc));
		if (qw < minq) {
			minq = qw;
			qint = 2;
		}
		xx = x1[connect[it1][2]];
		yy = y1[connect[it1][2]];
		xx2 = x1[connect[it1][0]];
		yy2 = y1[connect[it1][0]];
		qa = euclDist(xx, yy, xx2, yy2);
		qb = euclDist(xx, yy, x1[it1], y1[it1]);
		qc = euclDist(xx2, yy2, x1[it1], y1[it1]);
		qw = Math.acos((qb * qb + qc * qc - qa * qa) / (2.0 * qb * qc));
		if (qw < minq) {
			minq = qw;
			qint = 3;
		}
		if (qint == 1) {
			computeSteinerPoint(it1, connect[it1][0], connect[it1][1]);
			musoubana(it1, connect[it1][0], connect[it1][1]);
		}
		if (qint == 2) {
			computeSteinerPoint(it1, connect[it1][1], connect[it1][2]);
			musoubana(it1, connect[it1][1], connect[it1][2]);
		}
		if (qint == 3) {
			computeSteinerPoint(it1, connect[it1][2], connect[it1][0]);
			musoubana(it1, connect[it1][2], connect[it1][0]);
		}
	}// tase1

	private void tase2(int it2) {
		computeSteinerPoint(it2, connect[it2][0], connect[it2][1]);
		musoubana(it2, connect[it2][0], connect[it2][1]);
	}// tase2

	private void tuika() {
		int tlabel, tra, tra2;
		tlabel = 0;
		tra = 0;
		tra2 = 0;
		while (tlabel == 0 && tra < 10) {
			tra++;
			int ter = 0;
			for (int it = 0; it < numPoints; it++) {
				if (cn[it] > 2) {
					tra2++;
					ter = 1;
					tase1(it);
					mst();
					hantei();
					if (NN > 2 * numPoints - 2) {
						break;
					}
				}// cn[it]==3
			}// it
			for (int it = 0; it < numPoints; it++) {
				if (cn[it] == 2) {
					double xx2, yy2, qa, qb, qc, qw;
					xx = x1[connect[it][0]];
					yy = y1[connect[it][0]];
					xx2 = x1[connect[it][1]];
					yy2 = y1[connect[it][1]];
					qa = euclDist(xx, yy, xx2, yy2);
					qb = euclDist(xx, yy, x1[it], y1[it]);
					qc = euclDist(xx2, yy2, x1[it], y1[it]);
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					if (qw < 2.0 * Math.PI / 3.0) {
						tra2++;
						ter = 2;
						tase2(it);
						mst();
						hantei();
						if (NN > 2 * numPoints - 2) {
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

	private void mrlonely() {
		double xx2, yy2, qa, qb, qc, qw, maxq, qwe;
		int tlabel, tra, tra2;
		tlabel = 0;
		tra = 0;
		tra2 = 0;
		while (tlabel == 0) {
			tra++;
			int ter = 0;
			for (int it1 = numPoints; it1 < NN; it1++) {
				if (cn[it1] == 3) {
					xx = x1[connect[it1][0]];
					yy = y1[connect[it1][0]];
					xx2 = x1[connect[it1][1]];
					yy2 = y1[connect[it1][1]];
					qa = euclDist(xx, yy, xx2, yy2);
					qb = euclDist(xx, yy, x1[it1], y1[it1]);
					qc = euclDist(xx2, yy2, x1[it1], y1[it1]);
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = euclDist(qw, 0.0, 2.0 * Math.PI / 3.0, 0.0);
					maxq = qwe;
					xx = x1[connect[it1][1]];
					yy = y1[connect[it1][1]];
					xx2 = x1[connect[it1][2]];
					yy2 = y1[connect[it1][2]];
					qa = euclDist(xx, yy, xx2, yy2);
					qb = euclDist(xx, yy, x1[it1], y1[it1]);
					qc = euclDist(xx2, yy2, x1[it1], y1[it1]);
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = euclDist(qw, 0.0, 2.0 * Math.PI / 3.0, 0.0);
					if (qwe > maxq) {
						maxq = qwe;
					}
					xx = x1[connect[it1][2]];
					yy = y1[connect[it1][2]];
					xx2 = x1[connect[it1][0]];
					yy2 = y1[connect[it1][0]];
					qa = euclDist(xx, yy, xx2, yy2);
					qb = euclDist(xx, yy, x1[it1], y1[it1]);
					qc = euclDist(xx2, yy2, x1[it1], y1[it1]);
					qw = Math.acos((qb * qb + qc * qc - qa * qa)
							/ (2.0 * qb * qc));
					qwe = euclDist(qw, 0.0, 2.0 * Math.PI / 3.0, 0.0);
					if (qwe > maxq) {
						maxq = qwe;
					}
					if (maxq > Math.PI / 180.0) {
						ter++;
						computeSteinerPoint(connect[it1][0], connect[it1][1],
								connect[it1][2]);
						msbn2(it1, connect[it1][0], connect[it1][1],
								connect[it1][2]);
						mst();
						hantei();
					}// maxq>pi/180.0
				}// cn[it1]==2
			}// it1
			if (ter == 0) {
				tlabel = 1;
			}
		}// tlabel
		if (tra2 > 0) {
			bohe();
			mst();
			hantei();
		}
	}// mrlonely

	private void hantei() {
		if (numPoints > 2) {
			int ih;
			dd2 = treeLength();
			if (dd2 < mind) {
				counter2++;
				mind = dd2;
				for (ih = 0; ih < NN; ih++) {
					xmin[ih] = x1[ih];
					ymin[ih] = y1[ih];
					slmin[ih] = sl[ih];
				}
				minNN = NN;
				for (ih = 0; ih < NN - 1; ih++) {
					e1min[ih] = e1[ih];
					e2min[ih] = e2[ih];
				}
				repaint();
			}
		}
	}// hantei

	private void hontai() {
		mouselabel = 0;
		int timlabel = 0;
		while (mouselabel == 0) {
			tim++;
			timlabel = 0;
			if (numPoints > 20) {
				if (tim % 3 == 0) {
					timlabel = 1;
				}
			}
			if (numPoints > 9 && numPoints < 21) {
				if (tim % 30 == 0) {
					timlabel = 1;
				}
			}
			if (numPoints < 10) {
				if (tim % 100 == 0) {
					timlabel = 1;
				}
			}
			if (timlabel == 1) {
				repaint();
			}
			sco = 0;
			NN = numPoints;
			for (int k = 0; k < numPoints - 2; k++) {
				if (Math.random() < 0.5) {
					int chrule, br2;
					double x1kouho = 100.0, y1kouho = 100.0;
					chrule = 0;
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
					x1[numPoints + sco] = x1kouho;
					y1[numPoints + sco] = y1kouho;
					sl[numPoints + sco] = 1;
					sco++;
					NN = numPoints + sco;
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
		}// mouselabel
	}// hontai

	public void run() {
		counter2 = 0;
		xx = 0.0;
		yy = 0.0;
		//width = 500;
		//height = 350;

		while(isRunning) {
			for (int k = 0; k < numPoints; k++) {
				x1[k] = p[k].x + (Math.cos(k * 1.1) + k) * 0.0001;
				y1[k] = p[k].y + (Math.sin(k * 1.1) + k) * 0.001;
				sl[k] = 0;
				cn[k] = 0;
				w1[k] = Math
						.asin((y1[k] - yy) / euclDist(x1[k], y1[k], xx, yy))
						+ Math.PI / 2;
			}
			heapv(x1, y1, w1, numPoints);

			for (int k = 0; k < numPoints; k++) {
				x2[k] = x1[k];
				y2[k] = y1[k];
				xmin[k] = x1[k];
				ymin[k] = y1[k];
				slmin[k] = 0;
			}
			minNN = numPoints;
			NN = numPoints;
			if (numPoints > 2) {
				convexHull();
				mind = 99999.9;
				tim = 0;
				mst();
				hantei();
				hontai();
			}// N>1
		}// next;;
	}// run

	public void update(Graphics g) {
		paint(g);
	}

	public void paint(java.awt.Graphics g) {
		g.setColor(bgColor);
		g.fillRect(1, 1, width, height);
		g.setColor(Color.pink);
		for (int i = 0; i < NN - 1; i++) {
			g.drawLine((int) x1[e1[i]], (int) y1[e1[i]], (int) x1[e2[i]],
					(int) y1[e2[i]]);
		}

		g.setColor(colors[0]);
		for (int i = 0; i < minNN; i++) {
			if (slmin[i] == 0) {
				g.setColor(colors[0]);
				g.fillOval((int) xmin[i] - 3, (int) ymin[i] - 3, 6, 6);
			}
			// else{
			if (slmin[i] == 1) {
				g.setColor(colors[2]);
				g.drawRect((int) xmin[i] - 5, (int) ymin[i] - 5, 10, 10);
			}
		}

		g.setColor(colors[1]);
		if (numPoints == 2) {
			g.drawLine((int) x1[0], (int) y1[0], (int) x1[1], (int) y1[1]);
		} else {
			for (int i = 0; i < minNN - 1; i++) {
				g.drawLine((int) xmin[e1min[i]], (int) ymin[e1min[i]],
						(int) xmin[e2min[i]], (int) ymin[e2min[i]]);
			}
		}

		paintLegend(g);
	}// paint

	private void paintLegend(Graphics g) {
		g.drawString("N=" + numPoints, 15, 15);
		g.drawString("M=" + (minNN - numPoints), 15, 30);
		g.drawString("N+M=" + minNN, 65, 15);
		g.drawString("counter=" + tim, 65, 30);
		if (numPoints < 2) {
			g.drawString("d= 0.0", 145, 15);
		}
		if (numPoints == 2) {
			g.drawString("d=" + euclDist(x1[0], y1[0], x1[1], y1[1]), 145, 15);
		}
		if (numPoints > 2) {
			g.drawString("d=" + mind, 145, 15);
		}
		g.drawString("counter2=" + counter2, 145, 30);
		g.drawLine(380, 15, 480, 15);
		g.drawLine(380, 15, 380, 10);
		g.drawLine(430, 15, 430, 10);
		g.drawLine(480, 15, 480, 10);
		g.drawString("0", 370, 15);
		g.drawString("100", 481, 15);
	}

	/**
	 * Handles the mouse click.
	 * @param x
	 * @param y
	 * @return
	 */
	public boolean handleMouseClick(int x, int y) {
		if (numPoints < 50) {
			tim = 0;
			mind = 99999.9;
			p[numPoints] = new Point(x, y);
			x1[numPoints] = p[numPoints].x
					+ (Math.cos(numPoints * 1.1) + numPoints) * 0.0001;
			y1[numPoints] = p[numPoints].y
					+ (Math.sin(numPoints * 1.1) + numPoints) * 0.001;
			xmin[numPoints] = p[numPoints].x
					+ (Math.cos(numPoints * 1.1) + numPoints) * 0.0001;
			ymin[numPoints] = p[numPoints].y
					+ (Math.sin(numPoints * 1.1) + numPoints) * 0.001;
			sl[numPoints] = 0;
			slmin[numPoints] = 0;
			numPoints++;
			minNN = numPoints;
			NN = numPoints;
			counter2 = 0;
		} else {
			numPoints = 0;
			minNN = 0;
			NN = 0;
		}
		mouselabel = 1;
		repaint();
		return true;
	}
}
