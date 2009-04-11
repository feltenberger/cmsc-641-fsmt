package edu.umbc.algorithms.fmst.test;

import java.io.File;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;

import edu.umbc.algorithms.fmst.FairSMT;
import edu.umbc.algorithms.fmst.util.GraphUtils;

/**
 * @author Dave
 *
 */
public class FSMTDataCollection {
	private static final Logger log = Logger.getLogger(FSMTDataCollection.class);
	public static final String HEADER = "run num	std dev	max pcr	min pcr	avg pcr	total pcr	num terminal	k	target pcr";
	/**
	 * creates a fair smt with the given parameters.
	 * 
	 * @param heightAndWidth
	 * @param numNodes
	 * @param pcrGoal
	 * @param k
	 * @return
	 */
	private FairSMT generateFSMT(int heightAndWidth, int numNodes, double pcrGoal, int k) {
		FairSMT fsmt = new FairSMT(heightAndWidth, heightAndWidth, numNodes, k);
		configureSMT(fsmt, pcrGoal, k);
		return fsmt;
	}

	/**
	 * configures the fsmt.
	 * @param fsmt
	 * @param pcrGoal
	 * @param k
	 */
	private void configureSMT(FairSMT fsmt, double pcrGoal, int k) {
		fsmt.setTargetPCR(pcrGoal);
		fsmt.setMaxRelayNodes(k);
		fsmt.setRunWithGraphicsEnabled(false);
	}

	/**
	 * runs the FSMT
	 * 
	 * @param fsmt
	 * @throws Exception
	 */
	private void runFSMT(FairSMT fsmt) throws Exception {
		try {
			// generate the steiner tree.
			fsmt.start();
			// sleep for 7 seconds to let it generate the tree.
			Thread.sleep(1000 * 5);
			// stop it.
			fsmt.stop();

			// sleep to make sure the thread stopped.
			Thread.sleep(5);

			// start the "make fair" process.
			fsmt.startMakeFair();
		} catch (Exception e) {
			fsmt.stop();
			throw e;
		}
	}

	/**
	 * @param minHW
	 * @param maxHW
	 * @param hwIncrement
	 * @param minK
	 * @param maxK
	 * @param runsPerSize
	 * @throws Exception
	 */
	private void runTest(int minHW, int maxHW, int hwIncrement, int minNumNodes, int maxNumNodes, int numNodeIncrement, int minK, int maxK, int runsPerSize) throws Exception {

		//FairSMT fsmt = GraphUtils.readPlainTextSMT("data/600x600-10-nodes.smt");
		//dataCollector.configureSMT(fsmt, 0.0, 5);
		//dataCollector.runFSMT(fsmt);

		String results = HEADER;
		for(int hw = minHW; hw <= maxHW; hw += hwIncrement) {
			for(int k = minK; k <= maxK; k++) {
				for(int nn = minNumNodes; nn <= maxNumNodes; nn += numNodeIncrement){
					for(int run = 0; run < runsPerSize; run++) {
						FairSMT fsmt = null;
						try {
							fsmt = this.generateFSMT(hw, nn, 0.0, k);
							this.runFSMT(fsmt);
							results = GraphUtils.appendResults(results, fsmt, run);
						}
						catch(Throwable t) {
							log.error("Error!  Continuing with next test...", t);
						}
					}
				}
			}
		}

		// dump the results of the runs to a file.
		FileUtils.writeStringToFile(new File("data/results_" + System.currentTimeMillis() + ".csv"), results);
	}

	/**
	 * starting point...
	 * 
	 * @param args
	 */
	public static void main(String...args) throws Exception {
		FSMTDataCollection dataCollector = new FSMTDataCollection();
		// try all sizes from 100x100 to 600x600 (increments of 100)
		//	number of nodes from 10 to 50 (increments of 10)
		//		from k=1 to 20
		//			run each configuration 10 times
		//dataCollector.runTest(100, 600, 200, 50, 50, 5, 1, 5, 3);
		dataCollector.runTest(100, 600, 100, 10, 50, 10, 1, 20, 10);
	}
}
