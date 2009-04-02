package edu.umbc.algorithms.fmst;

import edu.umbc.algorithms.fmst.util.GraphUtils;

import java.awt.Graphics;

/**
 * @author dave
 */
public class Edge implements Comparable {
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
     * @param p1 first point
     * @param p2 second point
     */
    public Edge(Point p1, Point p2) {
        this.p1 = p1;
        this.p2 = p2;
    }

    /**
     * Draw this edge.
     *
     * @param g graphics
     */
    public void draw(Graphics g) {
        g.drawLine(p1.xInt(), p1.yInt(), p2.xInt(), p2.yInt());
    }

    /**
     * @return copies itself and returns the copy
     */
    public Edge copy() {
        Edge edge = new Edge(p1, p2);
        edge.index1 = this.index1;
        edge.index2 = this.index2;
        return edge;
    }

    public double getDistance() {
        return GraphUtils.euclideanDistance(p1.x, p1.y, p2.x, p2.y);
    }

    public int compareTo(Object o) {
        Edge e = (Edge) o;
        if (this.getDistance() == e.getDistance()) {
            return 0;
        } else if (this.getDistance() > e.getDistance()) {
            return 1;
        } else {
            return -1;
        }
    }
}
