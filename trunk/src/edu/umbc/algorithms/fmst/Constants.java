package edu.umbc.algorithms.fmst;

/**
 * @author : Fatih Senel
 *         Date: Apr 2, 2009
 *         Time: 12:48:37 PM
 */
public interface Constants {
    public static final Double DEFAULT_TARGET_STDEV = 100.0;  //to be need to changes
    public static final Double MAX_PCR_ALLOWED = 4000.0;
    public static final int DEFAULT_MAX_RELAY_NODES = 5;
    public static final Double TRANSMISSION_RANGE = Math.sqrt(MAX_PCR_ALLOWED);
    public static final int CONVERGENCE_CUTOFF = 5;  //when we see the same STDEV  CONVERGENCE_CUTOFF times in a row, the we converged
    public static final int CONVERGENCE_THRESHOLD = 10;  //the difference of the standard deviation between loops must be greater
                                     //than this value, otherwise we consider it as "converging"
}
