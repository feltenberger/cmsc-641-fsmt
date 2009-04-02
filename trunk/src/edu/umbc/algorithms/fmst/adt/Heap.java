package edu.umbc.algorithms.fmst.adt;

/**
 * @author : Fatih Senel
 * Date: Apr 2, 2009
 * Time: 3:27:42 AM
 */
import java.util.*;

// I can't believe that the huge Java API doesn't already include such
// a basic data structure

/**
 * Abstract implementation of the basic functions needed for a binary
 * qheap using java.util.Vector as a backend.  Unlike
 * java.util.TreeSet, this data structure can handle duplicate
 * entries.
 *
 * @author Ron Weiss (ronw@ee.columbia.edu)
 */
@SuppressWarnings({"unchecked"})
public abstract class Heap extends Vector
{
    private Comparator comp = null;

    public Heap()
    {
        super();
    }

    /**
     *  Use given Comparator for all comparisons between elements in
     *  this Heap.  Otherwise rely on compareTo methods and Comparable
     *  Objects.
     * @param c comparator of nodes
     */
    public Heap(Comparator c)
    {
        super();
        comp = c;
    }

    public Heap(int capacity)
    {
        super(capacity);
    }

    public Heap(Collection c)
    {
        super();
        addAll(c);
    }

    public Object remove(int index)
    {
        Object o = get(index);

        set(index, get(size()-1));
        removeElementAt(size()-1);
        heapify(index);

        return o;
    }

    public boolean remove(Object o)
    {
        boolean found = false;
        for(int i = 0; i < size(); i++)
        {
            if(o == null ? get(i) == null : o.equals(get(i)))
            {
                found = true;
                remove(i);

                break;
            }
        }

        return found;
    }

    public boolean add(Object o)
    {
        boolean b = super.add(o);

        for(int node = size()-1; node > 0;)
        {
            int parent = (node-1)/2;

            if(cmp(node, parent) < 0)
            {
                // swap them and reheapify
                Object tmp = get(node);
                set(node, get(parent));
                set(parent, tmp);
            }

            node = parent;
        }

        //System.out.print("\nContents: ");
        //for(int x = 0; x < size(); x++)
        //    System.out.print(get(x) + " ");
        //System.out.println();

        return b;
    }

    public boolean addAll(Collection c)
    {
        boolean b = super.addAll(c);
        rebuildHeap();

        return(b);
    }

    /**
     * Ensure that every element in this heap obeys the heap property.
     * Runs in linear time.
     *
     * This is meant to be called if/when the Comparator associated
     * with this object is modified.
     */
    public void rebuildHeap()
    {
        // do the whole linear time build-heap thing
        for(int i = size()/2; i >= 0; i--)
            heapify(i);
    }

    /**
     * Perform an in place heap sort on the data stored in this heap.
     * After calling sort, a call to this objects iterator() method
     * will iterate through the data stored in the heap in ascending
     * sorted order.
     */
    public void sort()
    {
        for(int x = size()-1; x > 0; x--)
        {
            // swap end of heap with the root, then heapify whats
            // left.
            Object tmp = get(x);
            set(x, get(0));
            set(0, tmp);

            heapify(0, x);
        }
    }

    /**
     * Compare two Objects in this heap - wrapper around
     * compareTo/Comparator.compare
     * @param node1 node 1
     * @param node2 second node
     * @return comparison result
     */
    protected int cmp(int node1, int node2)
    {
        int c;
        if(comp != null)
            c = comp.compare(get(node1), get(node2));
        else
            c = ((Comparable)get(node1)).compareTo(get(node2));

        return c;
    }

    /**
     * Ensure that the subtree rooted at node obeys the heap property
     * @param node root fo subtree
     * @param size size
     */
    private void heapify(int node, int size)
    {
        if(node > size)
            return;

        int left = (node+1)*2-1;
        int right = (node+1)*2;

        int minidx = node;
        if(left < size && cmp(left, node) < 0)
            minidx = left;
        if(right < size && cmp(right, node) < 0 && cmp(right, left) < 0)
            minidx = right;

        if(minidx != node)
        {
            // swap them and recurse on the subtree rooted at minidx
            Object tmp = get(node);
            set(node, get(minidx));
            set(minidx, tmp);
            heapify(minidx, size);
        }
    }

    /**
     * Ensure that the subtree rooted at node obeys the heap property
     * @param node root of subtree
     */
    private void heapify(int node)
    {
        heapify(node, size());
    }
}

