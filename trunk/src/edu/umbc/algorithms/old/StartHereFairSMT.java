package edu.umbc.algorithms.old;

import java.io.File;

import javax.swing.JFrame;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.SerializationUtils;
import org.apache.log4j.Logger;

import edu.umbc.algorithms.fmst.GraphUtils;


public class StartHereFairSMT {
	private static final Logger log = Logger.getLogger(StartHereFairSMT.class);

	public static void main(String...args) throws Exception {
		FairSMT blah = GraphUtils.readSMT("data/dave_test.txt");
		createNewJFrame(blah);

		if(true)
			return;

		FairSMT ms = loadReferenceSMT("data/reference_fair.smt");
		FairSMT msPt = GraphUtils.readSMT("data/reference_fair.txt");

		if(ms == null) {
			ms = new FairSMT(800, 800, 100);
			// serialize the FairSMT object into a byte[]
			byte[] msSerialized = SerializationUtils.serialize(ms);
			FileUtils.writeByteArrayToFile(new File("data/reference_fair.smt"), msSerialized);
			// deserialize the FairSMT
			// FairSMT msDeserialized = (FairSMT)SerializationUtils.deserialize(msSerialized);
		}
		if(msPt == null) {
			GraphUtils.dumpSMT("data/reference_fair.txt", ms, false);
			msPt = GraphUtils.readSMT("data/reference_fair.txt");
		}

		// visualize the min steiner trees
		createNewJFrame(ms);
		createNewJFrame(msPt);
		//createNewJFrame(msDeserialized);

		// sleep for 10 seconds
		Thread.sleep(50*1000);

		// stop the min steiner trees after sleeping
		ms.stop();
		msPt.stop();
		//msDeserialized.stop();
	}

	/**
	 * just add the min steiner object to a jframe and show/start it
	 * @param ms
	 */
	private static void createNewJFrame(FairSMT ms) {
		JFrame f = new JFrame();
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);


		f.add(ms);
		f.pack();
		f.setVisible(true);
		ms.start();
	}

	/**
	 * Loads the reference SMT object.
	 * @param file
	 * @return
	 * @throws Exception
	 */
	private static FairSMT loadReferenceSMT(String file) throws Exception {
		FairSMT referenceSMT = null;
		File reference = new File(file);
		if(reference.exists()) {
			byte[] tmp = FileUtils.readFileToByteArray(reference);
			try {
				log.info("Trying to load " + file);
				referenceSMT = (FairSMT)SerializationUtils.deserialize(tmp);
				log.info("Successfully loaded " + file);
			} catch (Exception e) {
				log.error("Error loading the SMT!");
			}
		}
		return referenceSMT;
	}

}
