package edu.umbc.algorithms.fmst;

import java.awt.Graphics;
import java.io.Serializable;
import java.util.ArrayList;

import edu.umbc.algorithms.fmst.util.GraphUtils;

/**
 * @author dave
 */
@SuppressWarnings("serial")
public class Point implements Serializable, Comparable {
    public static final int POINT_SIZE = 6;
    public double x;
    public double y;
    public double w;

    /**
     * fatih
     * this list is used to determine the neighboring information
     */
    public ArrayList<Point> neighbors = new ArrayList<Point>();


    /**
     * 0 = no, >=1 is yes
     */
    public int steiner = 0;

    /**
     * @param x coordinate
     * @param y coordinate
     */
    public Point(int x, int y) {
        this.x = x;
        this.y = y;
    }

    /**
     * @param x coordinate
     * @param y coordinate
     */
    public Point(double x, double y) {
        this.x = x;
        this.y = y;
    }

    /**
     * @return int value of x coordinate
     */
    public int xInt() {
        return (int) x;
    }

    /**
     * @return integer value of y coordinate
     */
    public int yInt() {
        return (int) y;
    }

    /**
     * This point draws itself, taking into consideration
     * any offsets to make graphics look prettier, etc.
     *
     * @param g graphics
     */
    public void draw(Graphics g) {
        int offset = POINT_SIZE / 2;
        g.fillOval(this.xInt() - offset, this.yInt() - offset, POINT_SIZE, POINT_SIZE);
    }

    /**
     * make a copy of yourself
     *
     * @return copied object
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
     * @return true if the node is a steiner point
     */
    public boolean isSteiner() {
        return this.steiner > 0;
    }


    /**
     * Power Consumption for a given node
     *
     * @return maximum edge weight to the power 2
     */
    public double getPCR() {
    	double max = -1;
        for (int i = 0; i < neighbors.size(); i++) {
            Point neighbor = neighbors.get(i);
            double dist = GraphUtils.euclideanDistance(x, y, neighbor.x, neighbor.y);
            if (dist > max) {
                max = dist;
            }
        }
        return Math.pow(max,2);
    }
    
    public double getAVGPCR() {
    	double avg = 0;
        for (int i = 0; i < neighbors.size(); i++) {
            Point neighbor = neighbors.get(i);
            double dist = GraphUtils.euclideanDistance(x, y, neighbor.x, neighbor.y);
            
            avg += dist;            
        }
        return Math.pow(avg/neighbors.size(),2);
    }
       
    public Point getLowestPCRNeighbor()
    {
        double minPCR = Double.POSITIVE_INFINITY;
        Point lowest = null;
        for (int i = 0; i < neighbors.size(); i++) {
            Point neighbor = neighbors.get(i);
           
            if (neighbor.getPCR() < minPCR)
            {
            	minPCR = neighbor.getPCR();
            	lowest = neighbor;            	
            }                   
        }
        return lowest;
    }
    

    public Point getFarthestNeighbor(){
        double max = -1;
        int index = -1;
        for (int i = 0; i < neighbors.size(); i++) {
            Point neighbor = neighbors.get(i);
            double dist = GraphUtils.euclideanDistance(x, y, neighbor.x, neighbor.y);
            if (dist > max) {
                max = dist;
                index = i;
            }
        }
        return neighbors.get(index);

    }

    public boolean isSameAs(Point p)
    {
    	if (p.w == this.w && p.x == this.x && p.y == this.y)
    		return true;
    	return false;
    }
    
    public int compareTo(Object o) {
        Point e = (Point) o;
        if (this.getPCR() == e.getPCR()) {
            return 0;
        } else if (this.getPCR() > e.getPCR()) {
            return 1;
        } else {
            return -1;
        }
    }
    
}
