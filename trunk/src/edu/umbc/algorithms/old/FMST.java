package edu.umbc.algorithms.old;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.apache.log4j.Logger;

import edu.umbc.algorithms.fmst.Edge;
import edu.umbc.algorithms.fmst.Point;
import edu.umbc.algorithms.fmst.util.GraphUtils;

/**
 * A Fair Minimum Steiner Tree.
 * 
 * @author dave
 * 
 */
public class FMST extends JPanel {
	private static final long serialVersionUID = 2854736955451806258L;
	/**
	 * The logger...
	 */
	private static final Logger log = Logger.getLogger(FMST.class);
	/**
	 * The points in the graph.
	 */
	private List<Point> points = new ArrayList<Point>();
	/**
	 * All the edges in the graph.
	 */
	private List<Edge> edges = new ArrayList<Edge>();
	/**
	 * The width of the window.
	 */
	private int width = 500;
	/**
	 * The height of the window.
	 */
	private int height = 500;
	/**
	 * The maximum number of nodes the tree can have.
	 */
	private int maxNodes = 10;
	/**
	 * The maximun number of neighbors each node can have. (Is this
	 * needed/relevant?)
	 */
	private int maxNeighbors = 3;

	/**
	 * Create a new Fair Minimum Steiner Tree.
	 */
	public FMST() {
		log.info("Creating new FMST.");
		for (int i = 0; i < maxNodes; i++) {
			this.points.addAll(GraphUtils.randomPoint(width, height, 5,
					maxNeighbors));
		}

		
		this.setPreferredSize(new Dimension(this.width + 20, this.height + 20));
		this.setSize(width + 20, height + 20);
	}

	/**
	 * Update the FMST.
	 * 
	 * @see javax.swing.JComponent#update(java.awt.Graphics)
	 */
	public void update(Graphics g) {
		paint(g);
	}

	/**
	 * Draw the Steiner tree.
	 * 
	 * @see java.awt.Container#paint(java.awt.Graphics)
	 */
	public void paint(Graphics g) {
		// Outline the active area.
		g.setColor(Color.black);
		g.drawLine(0, 0, width, 0);
		g.drawLine(0, 0, 0, height);
		g.drawLine(0, height, width, height);
		g.drawLine(width, 0, width, height);

		// draw each point
		for (Point point : points) {
			point.draw(g);
			// g.fillOval(point.x, point.y, 10, 10);
		}

		// draw each edge
		for (Edge edge : edges) {
			edge.draw(g);
			// g.drawLine(point.x, point.y, neighbor.x, neighbor.y);
		}

	}

	public static void main(String... args) {
		FMST fmst = new FMST();
		JFrame outerFrame = new JFrame();
		outerFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		outerFrame.add(fmst);
		outerFrame.pack();
		outerFrame.setVisible(true);
	}
}
