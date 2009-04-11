package edu.umbc.algorithms.fmst.test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

import org.apache.commons.io.FileUtils;

/**
 * @author Dave
 *
 */
public class DataCrunchingTest extends TestCase {
	//private static final Logger log = Logger.getLogger(DataCrunchingTest.class);
	private List<Container> values = new ArrayList<Container>();
	//private List<Container> averages = new ArrayList<Container>();
	private String filePrefix = "data/results_1239437137722";

	/* (non-Javadoc)
	 * @see junit.framework.TestCase#setUp()
	 */
	@SuppressWarnings("unchecked")
	protected void setUp() throws Exception {
		super.setUp();
		List<String> lines = FileUtils.readLines(new File(filePrefix + ".csv"));
		boolean firstLine = true;
		for(String line : lines) {
			if(firstLine) {
				firstLine = false;
				continue;
			}
			if(line != null && !"".equals(line.trim())) {
				String[] lineParts = line.split(",");
				Container c = new Container();
				c.run = Integer.parseInt(lineParts[0]);
				c.std_dev = Double.parseDouble(lineParts[1]);
				c.max_pcr = Double.parseDouble(lineParts[2]);
				c.min_pcr = Double.parseDouble(lineParts[3]);
				c.avg_pcr = Double.parseDouble(lineParts[4]);
				c.total_pcr = Double.parseDouble(lineParts[5]);
				c.num_terminals = Integer.parseInt(lineParts[6]);
				c.k = Integer.parseInt(lineParts[7]);
				c.target_pcr = Double.parseDouble(lineParts[8]);
				values.add(c);
			}
		}
	}

	/* (non-Javadoc)
	 * @see junit.framework.TestCase#tearDown()
	 */
	protected void tearDown() throws Exception {
		super.tearDown();
	}

	/**
	 * @throws Exception
	 */
	public void testCalculateAverages() throws Exception {
		Container tmp = new Container();
		StringBuffer results = new StringBuffer(FSMTDataCollection.HEADER.replaceAll("\t", ",")+",canvas size").append("\n");
		int numLines = 0;
		int size = 100;
		int prevNumTerminals = 0;
		int prevRunNum = Integer.MAX_VALUE;
		for(Container c : values) {
			if(c.run < prevRunNum && numLines > 0) {
				appendResults(tmp, results, size);

				// start averages over.
				tmp = new Container();
				if(prevNumTerminals != c.num_terminals && c.num_terminals == 10 && c.k == 1)
					size += 100;
			}
			tmp.run++;
			tmp.std_dev += c.std_dev;
			tmp.max_pcr += c.max_pcr;
			tmp.min_pcr += c.min_pcr;
			tmp.avg_pcr += c.avg_pcr;
			tmp.total_pcr += c.total_pcr;
			tmp.num_terminals += c.num_terminals;
			tmp.k += c.k;
			prevRunNum = c.run;
			prevNumTerminals = c.num_terminals;
			numLines++;
		}
		appendResults(tmp, results, size);
		FileUtils.writeStringToFile(new File(filePrefix + "_avgs_k-15.csv"), results.toString());
	}

	/**
	 * append the averages to the string buffer
	 * @param tmp
	 * @param results
	 * @param size
	 */
	private void appendResults(Container tmp, StringBuffer results, int size) {
		tmp.std_dev = tmp.std_dev / tmp.run;
		tmp.max_pcr = tmp.max_pcr / tmp.run;
		tmp.min_pcr = tmp.min_pcr / tmp.run;
		tmp.avg_pcr = tmp.avg_pcr / tmp.run;
		tmp.total_pcr = tmp.total_pcr / tmp.run;
		tmp.num_terminals = tmp.num_terminals / tmp.run;
		tmp.k = tmp.k / tmp.run;
		tmp.target_pcr = tmp.target_pcr / tmp.run;

		//if(tmp.num_terminals != 20 || size != 500)
		if(tmp.k != 15)
			return;

		results.append(tmp.run).append(",")
		.append(tmp.std_dev).append(",")
		.append(tmp.max_pcr).append(",")
		.append(tmp.min_pcr).append(",")
		.append(tmp.avg_pcr).append(",")
		.append(tmp.total_pcr).append(",")
		.append(tmp.num_terminals).append(",")
		.append(tmp.k).append(",")
		.append(tmp.target_pcr).append(",")
		.append(size).append("\n");
	}

	private class Container {
		int run;
		double std_dev;
		double max_pcr;
		double min_pcr;
		double avg_pcr;
		double total_pcr;
		int num_terminals;
		int k;
		double target_pcr;
	}
}
