package edu.umbc.algorithms.fmst;

import org.apache.log4j.Logger;

import edu.umbc.algorithms.fmst.util.GraphUtils;

/**
 * this just boot-straps the creation of the visualizations etc.
 * @author dave
 *
 */
public class StartHereFairSMT {
	private static final Logger log = Logger.getLogger(StartHereFairSMT.class);

	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(String...args) throws Exception {
		log.info("Starting now...");
		//FairSMT blah = GraphUtils.readPlainTextSMT("data/dave_test.smt");
		//showSMT(blah);

		//FairSMT ms = GraphUtils.readPlainTextSMT("data/reference_fair.smt");
		FairSMT ms = GraphUtils.readPlainTextSMT("data/3_nodes.smt");
		if(ms == null) {
			// create a new SMT.
			ms = new FairSMT(500, 500, 50);
			GraphUtils.saveSMTPlainText("data/reference_fair.smt", ms, false);
		}

		// visualize the min steiner trees
		showSMT(ms);
	}

	/**
	 * just add the min steiner object to a jframe and show/start it
	 * @param ms
	 */
	@SuppressWarnings("all")
	private static void showSMT(FairSMT ms) {
		FairSMTViz viz = new FairSMTViz(ms, true);
	}

}
