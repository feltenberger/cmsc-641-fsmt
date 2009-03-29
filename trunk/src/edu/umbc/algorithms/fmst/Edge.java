package edu.umbc.algorithms.fmst;

import java.awt.Graphics;

/**
 * @author dave
 *
 */
public class Edge {
	/**
	 * this is here for backwards compatibility with the japanese guy's steiner code
	 */
	public int index1;
	public int index2;
	/**
	 * The first point attached to the edge.
	 */
	public Point p1;
	/**
	 * The other point attached to the edge.
	 */
	public Point p2;

	/**
	 * @param p1
	 * @param p2
	 */
	public Edge(Point p1, Point p2) {
		this.p1 = p1;
		this.p2 = p2;
	}

	/**
	 * Draw this edge.
	 * @param g
	 */
	public void draw(Graphics g) {
		g.drawLine(p1.xInt(), p1.yInt(), p2.xInt(), p2.yInt());
	}

	/**
	 * copies itself and returns the copy
	 * @return
	 */
	public Edge copy() {
		Edge edge = new Edge(p1, p2);
		edge.index1 = this.index1;
		edge.index2 = this.index2;
		return edge;
	}
}
