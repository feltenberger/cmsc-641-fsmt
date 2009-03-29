package edu.umbc.algorithms.old;

import java.io.File;

import javax.swing.JFrame;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.SerializationUtils;
import org.apache.log4j.Logger;


public class StartHere {
	private static final Logger log = Logger.getLogger(StartHere.class);

	public static void main(String...args) throws Exception {
		MinSteiner ms = loadReferenceSMT("data/reference.smt");
		ms = null;

		if(ms == null) {
			ms = new MinSteiner(800, 800, 100);
			// serialize the MinSteiner object into a byte[]
			//byte[] msSerialized = SerializationUtils.serialize(ms);
			//FileUtils.writeByteArrayToFile(new File("data/reference.smt"), msSerialized);
			// deserialize the MinSteiner
			// MinSteiner msDeserialized = (MinSteiner)SerializationUtils.deserialize(msSerialized);
		}

		// visualize the min steiner trees
		createNewJFrame(ms);
		//createNewJFrame(msDeserialized);

		// sleep for 10 seconds
		Thread.sleep(50*1000);

		// stop the min steiner trees after sleeping
		ms.stop();
		//msDeserialized.stop();
	}

	/**
	 * just add the min steiner object to a jframe and show/start it
	 * @param ms
	 */
	private static void createNewJFrame(MinSteiner ms) {
		JFrame f = new JFrame();
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		ms.start();

		f.add(ms);
		f.pack();
		f.setVisible(true);
	}

	/**
	 * Loads the reference SMT object.
	 * @param file
	 * @return
	 * @throws Exception
	 */
	private static MinSteiner loadReferenceSMT(String file) throws Exception {
		MinSteiner referenceSMT = null;
		File reference = new File(file);
		if(reference.exists()) {
			byte[] tmp = FileUtils.readFileToByteArray(reference);
			try {
				referenceSMT = (MinSteiner)SerializationUtils.deserialize(tmp);
			} catch (Exception e) {
				log.error("Error loading the SMT!");
			}
		}
		return referenceSMT;
	}

}
