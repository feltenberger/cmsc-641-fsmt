package edu.umbc.algorithms.fmst.util;

//Fri Oct 25 18:07:43 EST 2004
//
// Written by Sean R. Owens, released to the public
// domain.  Share and enjoy.  http://darksleep.com/player

// A very simple custom dialog that takes a string as a parameter,
// displays it in a JLabel, along with two Jbuttons, one labeled Yes,
// and one labeled No, and waits for the user to click one of them.

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * modified from the above URL.
 * @author dave
 *
 */
public class SimpleDialog extends JDialog implements ActionListener {
	private static final long serialVersionUID = -7724235027148030532L;
	private JPanel myPanel = null;
	private JButton yesButton = null;
	private JButton noButton = null;
	private JTextField textField = null;
	private boolean clickedOk = false;
	private String value = null;

	/**
	 * did the user click ok?
	 * @return
	 */
	public boolean getClickedOk() {
		return this.clickedOk;
	}
	/**
	 * what value did the user click?
	 * @return
	 */
	public String getValue() {
		return this.value;
	}

	/**
	 * @param frame
	 * @param modal
	 * @param myMessage
	 */
	public SimpleDialog(JFrame frame, boolean modal, String myMessage) {
		super(frame, modal);

		// add the message to the panel.
		myPanel = new JPanel();
		getContentPane().add(myPanel);
		myPanel.add(new JLabel(myMessage));

		// add the text to the field.
		textField = new JTextField();
		textField.setPreferredSize(new Dimension(50, 20));
		textField.setSize(50, 20);
		myPanel.add(textField);

		// add yes button
		yesButton = new JButton("Okay");
		yesButton.addActionListener(this);
		myPanel.add(yesButton);

		// add no button
		noButton = new JButton("Canel");
		noButton.addActionListener(this);
		myPanel.add(noButton);

		pack();
		setLocationRelativeTo(frame);
		setVisible(true);
	}

	/* (non-Javadoc)
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e) {
		if (yesButton == e.getSource()) {
			clickedOk = true;
			value = textField.getText();
			setVisible(false);
		} else if (noButton == e.getSource()) {
			clickedOk = false;
			value = textField.getText();
			setVisible(false);
		}
	}

}