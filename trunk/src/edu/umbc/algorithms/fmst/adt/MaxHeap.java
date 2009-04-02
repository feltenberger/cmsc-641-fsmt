package edu.umbc.algorithms.fmst.adt;

/**
 * @author : Fatih Senel
 * Date: Apr 2, 2009
 * Time: 3:42:23 AM
 */
import java.util.*;



/**
 * Implementation of a binary max heap.
 *
 * @author Ron Weiss (ronw@ee.columbia.edu)
 */
public class MaxHeap extends Heap
{
    public MaxHeap()
    {
        super();
    }

    /**
     *  Use given Comparator for all comparisons between elements in
     *  this Heap.  Otherwise rely on compareTo methods and Comparable
     *  Objects.
     * @param c comparator of nodes
     */
    public MaxHeap(Comparator c)
    {
        super(c);
    }

    public MaxHeap(int capacity)
    {
        super(capacity);
    }

    public MaxHeap(Collection c)
    {
        super(c);
    }

    public Object extractMax()
    {
        return remove(0);
    }

    /**
     * Compare two Objects in this heap - wrapper around
     * compareTo/Comparator.compare
     */
    protected int cmp(int node1, int node2)
    {
        return -super.cmp(node1, node2);
    }


}
