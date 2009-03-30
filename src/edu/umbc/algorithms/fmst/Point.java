package edu.umbc.algorithms.fmst;

import java.awt.Graphics;
import java.io.Serializable;

/**
 * @author dave
 *
 */
@SuppressWarnings("serial")
public class Point implements Serializable {
	public static final int POINT_SIZE = 6;
	public double x;
	public double y;
	public double w;
	/**
	 * 0 = no, >=1 is yes
	 */
	public int steiner = 0;

	/**
	 * @param x
	 * @param y
	 */
	public Point(int x, int y) {
		this.x = x;
		this.y = y;
	}
	/**
	 * @param x
	 * @param y
	 */
	public Point(double x, double y) {
		this.x = x;
		this.y = y;
	}
	/**
	 * @return
	 */
	public int xInt() {
		return (int)x;
	}
	/**
	 * @return
	 */
	public int yInt() {
		return (int)y;
	}

	/**
	 * This point draws itself, taking into consideration 
	 * any offsets to make graphics look prettier, etc.
	 * @param g
	 */
	public void draw(Graphics g) {
		int offset = POINT_SIZE / 2;
		g.fillOval(this.xInt()-offset, this.yInt()-offset, POINT_SIZE, POINT_SIZE);
	}

	/**
	 * make a copy of yourself
	 * @return
	 */
	public Point copy() {
		Point p = new Point(x, y);
		p.steiner = steiner;
		p.w = w;
		return p;
	}

	/**
	 * Generates the weight value.
	 */
	public void generateW() {
		double half_pi = Math.PI / 2;
		this.w = Math.asin((y - 0.0) / GraphUtils.euclideanDistance(x, y, 0.0, 0.0)) + half_pi;
	}
	/**
	 * @return
	 */
	public boolean isSteiner() {
		return this.steiner > 0;
	}
}
