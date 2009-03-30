package edu.umbc.algorithms.fmst;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

import org.apache.log4j.Logger;

import edu.umbc.algorithms.fmst.util.GraphUtils;
import edu.umbc.algorithms.fmst.util.SimpleDialog;

/**
 * @author dave
 * 
 */
public class FairSMTViz extends JFrame {
	private static final long serialVersionUID = -6627308653284975372L;
	private static final Logger log = Logger.getLogger(FairSMTViz.class);
	private FairSMT smt = null;

	/**
	 * @param smt
	 */
	public FairSMTViz(FairSMT smt, boolean autoStart) {
		super();
		this.smt = smt;
		initMenus();
		initialize();
		if(autoStart)
			this.smt.start();
	}

	/**
	 * initialize options
	 */
	private void initialize() {
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.add(smt);
		this.pack();
		this.setVisible(true);
	}

	/**
	 * create the menu options
	 */
	private void initMenus() {
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		JMenu fileMenu = new JMenu("File");
		JMenu actionsMenu = new JMenu("Actions");
		menuBar.add(fileMenu);
		menuBar.add(actionsMenu);

		/* File Menu */

		JMenuItem newAction = new JMenuItem("New");
		JMenuItem loadAction = new JMenuItem("Load");
		JMenuItem saveAction = new JMenuItem("Save");
		JMenuItem exitAction = new JMenuItem("Exit");
		fileMenu.add(newAction);
		fileMenu.add(loadAction);
		fileMenu.add(exitAction);

		newAction.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				// get the height and width
				SimpleDialog dialog = new SimpleDialog(FairSMTViz.this, true, "What width / height (nxn)?");
				int heightAndWidth = Integer.parseInt(dialog.getValue());
				
				// get the number of nodes
				dialog = new SimpleDialog(FairSMTViz.this, true, "How many nodes to start with?");
				int numNodes = Integer.parseInt(dialog.getValue());
				
				// construct and load the new smt
				FairSMT smt = new FairSMT(heightAndWidth, heightAndWidth, numNodes);
				FairSMTViz.this.loadNewSMT(smt, true);
			}
		});
		loadAction.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fc = new JFileChooser();
				fc.setCurrentDirectory(new File("data"));
				int returnVal = fc.showOpenDialog(FairSMTViz.this);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					if (file.getName().endsWith("smt")) {
						try {
							FairSMT smt = GraphUtils.readPlainTextSMT(file.getAbsolutePath());
							FairSMTViz.this.loadNewSMT(smt, true);
							log.info("Loaded SMT: " + file.getName() + ".");
						} catch (Exception e1) {
							log.error("Error loading the SMT!", e1);
						}
					} else {
						log.info("This isn't an .smt file...not loading.");
					}
				}
			}
		});
		saveAction.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFileChooser fc = new JFileChooser();
				fc.setCurrentDirectory(new File("data"));
				int returnVal = fc.showSaveDialog(FairSMTViz.this);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					try {
						GraphUtils.saveSMTPlainText(file.getAbsolutePath(),
								FairSMTViz.this.smt, false);
					} catch (Exception e1) {
						log.error("Error saving SMT...", e1);
					}
					log.info("Saving SMT to: " + file.getName() + ".");
				}
			}
		});
		exitAction.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if(FairSMTViz.this.smt.isRunning())
					FairSMTViz.this.smt.stop();
				// sleep for a little to make sure it stops.
				try { Thread.sleep(100); } catch (Exception ignore) { }
				// exit.
				System.exit(0);
			}
		});

		/* Actions menu */
		JMenuItem stopAction = new JMenuItem("Stop");
		JMenuItem startAction = new JMenuItem("Start");
		actionsMenu.add(startAction);
		actionsMenu.add(stopAction);

		/* File Menu Options */
		stopAction.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (smt.isRunning()) {
					smt.stop();
				} else {
					log.info("It's already stopped!");
				}
			}
		});
		startAction.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (!smt.isRunning()) {
					smt.start();
				} else {
					log.info("SMT is already running!");
				}
			}
		});
	}

	/**
	 * stops and unloads the old one if it's running, 
	 * then loads the new one
	 * @param smt
	 * @param autoStart
	 */
	public void loadNewSMT(FairSMT smt, boolean autoStart) {
		if (FairSMTViz.this.smt.isRunning())
			FairSMTViz.this.smt.stop();
		FairSMTViz.this.remove(FairSMTViz.this.smt);

		FairSMTViz.this.smt = smt;
		// FairSMTViz.this.add(FairSMTViz.this.smt);
		FairSMTViz.this.initialize();
		if(autoStart)
			FairSMTViz.this.smt.start();
	}

}
