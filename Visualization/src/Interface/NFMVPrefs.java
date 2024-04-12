/*
 * NFMVPrefs.java
 * NairnFEAMPMViz Application
 * 
 * Created by John Nairn on 2/17/08
 * Copyright (c) 2008-2010 RSAC Software. All rights reserved
 */

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.util.ArrayList;

import javax.swing.*;

import geditcom.JNFramework.*;
import geditcom.plot2d.*;

public class NFMVPrefs extends JNPreferences implements ActionListener
{	
	private static final long serialVersionUID = 1L;

	// Results window split locations
	public static String ResultsTopSplitKey = "Results Top Split";
	public static double ResultsTopSplitDef = 0.35;
	public static String ResultsSplitKey = "Results Split";
	public static double ResultsSplitDef = 0.50;

	// Commands window split location
	public static String CommandsSplitKey = "Commands Split";
	public static double CommandsSplitDef = 0.60;

	// path to NairnMPM binary
	public static String NairnMPMKey = "NairnMPM Path";
	public static String NairnMPMDef = "$(bundle)";

	// path to NairnMPM dtd
	public static String NairnMPMDTDKey = "NairnMPM DTD File Path";
	public static String NairnMPMDTDDef = "$(bundle)";
	public static String NairnMPMValidateKey = "Validate MPM Calculations";
	public static boolean NairnMPMValidateDef = true;

	// path to NairnFEA binary
	public static String NairnFEAKey = "NairnFEA Path";
	public static String NairnFEADef = "$(bundle)";

	// path to NairnFEA DTD
	public static String NairnFEADTDKey = "NairnFEA DTD File Path";
	public static String NairnFEADTDDef = "$(bundle)";
	public static String NairnFEAValidateKey = "Validate MPM Calculations";
	public static boolean NairnFEAValidateDef = true;

	// path to shell command
	public static String ShellKey = "Shell Path";
	public static String ShellDef = "$(windows)";

	// path to workspace folder
	public static String WorkSpaceKey = "Workspace Folder";
	public static String WorkSpaceDef = "";

	// REMOTE_ACCESS
	// remote server
	public static String RemoteServerKey = "Remote Server";
	public static String RemoteUserKey = "Remote Username";
	public static String RemoteUserPassKey = "Remote Password";
	// code exec location (local or remote)
	public static String RemoteMPMPathKey = "Remote Path to MPM Code";
	public static String RemoteMPMDTDKey = "Remote Path to MPM DTD";
	public static String RemoteFEAPathKey = "Remote Path to FEA Code";
	public static String RemoteFEADTDKey = "Remote Path to FEA DTD";
	
	// Multiple Remove Servers
	public static String RemoteServerListKey = "Remote Server List";
	public static String RemoteServerListDef = "";
	public static String CodeExecLocationKey = "Code Exec Location";
	public static String CodeExecLocationDef = "local";

	// plot spectrum
	public static String SpectrumKey = "PlotSpectrum";
	public static int SpectrumDef = 1;
	public static String NumContoursKey = "NumberOfPlotContours";
	public static int NumContoursDef = 0;
	public static String NumSubelementsKey = "NumberOfSubelements";
	public static int NumSubelementsDef = 4;
	
	private JComboBox<String> families;
	private JComboBox<Integer> sizes;
	private JComboBox<String> ofamilies;
	private JComboBox<Integer> osizes;
	private JComboBox<String> rfamilies;
	private JComboBox<Integer> rsizes;
	
	// window fonts
	public static String ScriptFontKey = "Script Editing Font";
	public static String ScriptFontDef = "Monospaced";
	public static String ScriptFontSizeKey = "Script Editing Font Size";
	public static int ScriptFontSizeDef = 14;
	public static String OutputFontKey = "Output Font";
	public static String OutputFontDef = "Monospaced";
	public static String OutputFontSizeKey = "Output Font Size";
	public static int OutputFontSizeDef = 14;
	public static String ResultsFontKey = "Results Font";
	public static String ResultsFontDef = "Monospaced";
	public static String ResultsFontSizeKey = "Output Font Size";
	public static int ResultsFontSizeDef = 14;

	// plot colors
	public static String backColorKey = "Background Color";
	public static Color backColorDef = new Color(84, 89, 109); // 0.329f,0.349f,0.427f
	public static String meshLineColorKey = "Mesh Line Color";
	public static Color meshLineColorDef = new Color(191, 191, 191); // 0.75f,0.75f,0.75f
	public static String meshNodesColorKey = "Nodal Point Color";
	public static Color meshNodesColorDef = new Color(212, 212, 212); // 0.83f,0.83f,0.83f
	public static String textColorKey = "Text Label Color";
	public static Color textColorDef = new Color(255, 255, 255); // 1.f,1.f,1.f
	
	public static String cplaneColorKey = "Crack Plane Color";
	public static Color cplaneColorDef = new Color(0,0,0);
	public static String caboveColorKey = "Crack Top Surface Color";
	public static Color caboveColorDef = new Color(0,0,255);
	public static String cbelowColorKey = "Crack Top Surface Color";
	public static Color cbelowColorDef = new Color(255,0,0);
	public static String dispMeshLineColorKey = "Displaced Mesh Line Color";
	public static Color dispMeshLineColorDef = new Color(191, 191, 191); // 0.75f,0.75f,0.75f

	// plot colors not yet editable
	public static String tempBCColorKey = "Temperature BC Color";
	public static Color tempBCColorDef = Color.red;
	public static String concBCColorKey = "Temperature BC Color";
	public static Color concBCColorDef = Color.blue;
	
	// plot settings
	public static String maxElongKey = "Max Elongation Plotted";
	public static double maxElongDef = 5.0;									// <1 for no limit

	// plot options
	public static String PlotOptionsFlagsKey = "PlotOptionsFlags";
	public static int PlotOptionsFlagsDef = 0x007F;
	public static String ParticleSizeKey = "ParticleSize";
	public static int ParticleSizeDef = 100;
	
	private JTextField mpmCodePath = new JTextField();
	private JTextField mpmDTDPath = new JTextField();
	private JTextField feaCodePath = new JTextField();
	private JTextField feaDTDPath = new JTextField();
	private JTextField shellPath = new JTextField();
	private JTextField workPath = new JTextField();
	private JTextField numContours = new JTextField();
	private JTextField numSubelements = new JTextField();
	private JFileChooser chooser = new JFileChooser();

	// prefs
	private JComboBox<String> spectrumBox;

	private JCheckBox validateMPM = new JCheckBox("Validate",NairnMPMValidateDef);
	private JCheckBox validateFEA = new JCheckBox("Validate",NairnFEAValidateDef);

	// REMOTE_ACCESS
	private JComboBox<PlotMenuItem> remoteList=new JComboBox<PlotMenuItem>();
	private JTextField remoteServerAddr = new JTextField();
	private JTextField remoteUsername = new JTextField();
	private JPasswordField remoteUserpass = new JPasswordField();
	private JTextField mpmCodePathFld = new JTextField();
	private JTextField mpmDTDPathFld = new JTextField();
	private JTextField feaCodePathFld = new JTextField();
	private JTextField feaDTDPathFld = new JTextField();
	private JRadioButton rdbtnExecLocal = new JRadioButton("Local");
	private JRadioButton rdbtnExecRemote = new JRadioButton("Remote");
	
	private ArrayList<String> servers;
	private ArrayList<String> currentServer = new ArrayList<String>(7);
	public int currentServerIndex;
	public static String serverDelim = ",";
	
	private static boolean currentRemoteMode = false;
	
	public static final int ServerIP = 0;
	public static final int ServerUser = 1;
	public static final int ServerPWD = 2;
	public static final int ServerMPM = 3;
	public static final int ServerMPMDTD = 4;
	public static final int ServerFEA = 5;
	public static final int ServerFEADTD = 6;

	// ----------------------------------------------------------------------------
	// constructor
	// ----------------------------------------------------------------------------

	public NFMVPrefs()
	{	super("Preferences");

		JTabbedPane tabbedPane = new JTabbedPane();
		tabbedPane.setFocusable(false);

		// Remove panel
		tabbedPane.addTab("Code", null, codePanel(),"Settings for running calculations");

		// Remote preferences pane
		tabbedPane.addTab("Remote", null, remotePanel(), "Remote server settings");

		// color preferences pane ----------------------------------------
		tabbedPane.addTab("Colors", null, colorPanel(), "Select colors used in plots");

		// color preferences pane ----------------------------------------
		tabbedPane.addTab("Fonts", null, fontPanel(), "Select fonts for windows");

		// FEA and MPM preferences pane
		// -----------------------------------------
		//JPanel panel4 = new JPanel();
		//tabbedPane.addTab("FEA and MPM", null, panel4,"FEA and MPM plot settings");

		getContentPane().add(tabbedPane);
		setSize(600, 480);
		setLocation(60, 30 + JNApplication.menuBarHeight());
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
	}
	
	// create code preferences panel
	public JPanel codePanel()
	{	JPanel panel1 = new JPanel(new GridLayout(6,1));

		panel1.add(filePathPanel("MPM Code", mpmCodePath, prefs.get(
				NairnMPMKey, NairnMPMDef), false, null));
		panel1.add(filePathPanel("MPM DTD File", mpmDTDPath, prefs.get(
				NairnMPMDTDKey, NairnMPMDTDDef), prefs.getBoolean(
				NairnMPMValidateKey, NairnMPMValidateDef), validateMPM));
		panel1.add(filePathPanel("FEA Code", feaCodePath, prefs.get(
				NairnFEAKey, NairnFEADef), false, null));
		panel1.add(filePathPanel("FEA DTD File", feaDTDPath, prefs.get(
				NairnFEADTDKey, NairnFEADTDDef), prefs.getBoolean(
				NairnFEAValidateKey, NairnFEAValidateDef), validateFEA));
		panel1.add(filePathPanel("Shell Command", shellPath, prefs.get(
				ShellKey, ShellDef), false, null));
		panel1.add(filePathPanel("Work Space Directory", workPath, prefs.get(
				WorkSpaceKey, WorkSpaceDef), false, null));

		return panel1;
	}
	
	// REMOTE_ACCESS
	// build panel for remote connections
	public JPanel remotePanel()
	{	JPanel panel3 = new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
		panel3.setLayout(gridbag);
		
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.insets=new Insets(3,6,3,6);
		c.weighty = 0.;

		// Remote server list -------------------------------------------
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 1;
		JLabel label = new JLabel("Servers:");
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		panel3.add(label);
		
		c.gridx=1;
		c.weightx = 3.0;
		c.gridwidth = 3;
		remoteList.removeAllItems();
		remoteList.addItem(new PlotMenuItem("Pick a server..."));
		servers = getPrefList(RemoteServerListKey,RemoteServerListDef,serverDelim);
		if(servers.size()==1)
		{	// delete it if empty
			if(servers.get(0).length()==0)
				servers.clear();
		}
		
		// loop through servers and add valid ones to the menu
		int i=0;
		while(i<servers.size())
		{	String serverName = servers.get(i);
			ArrayList<String> config = getPrefList(serverName,"",serverDelim);
			if(config.size()<7)
			{	// not valid
				servers.remove(i);
			}
			else
			{	remoteList.addItem(new PlotMenuItem(serverName));
				i++;
			}
		}
		currentServer.clear();
		currentServerIndex = 0;			// item in menu (1 to # of servers)
		
		// finish menu and add to panel
		remoteList.addItem(new PlotMenuItem("Add server..."));
		remoteList.addItem(new PlotMenuItem("Delete current server..."));
		remoteList.setToolTipText("Pick existing server or create a new one");
		gridbag.setConstraints(remoteList,c);
		panel3.add(remoteList);
		
		// when use menu, insert the expression
		remoteList.addItemListener(new ItemListener()
		{	public void itemStateChanged(ItemEvent e)
			{	int numItems = remoteList.getItemCount();
				int newIndex = remoteList.getSelectedIndex();
				if(e.getStateChange()==ItemEvent.SELECTED)
				{	if(newIndex==numItems-2)
					{	// add a server (not added to refs until name is changed
						String addName = "server.ip.or.name";
						int extra = 1;
						while(!uniqueServer(addName))
						{	addName = "server.ip.or.name"+"."+extra;
							extra++;
						}
						currentServer.clear();
						currentServer.add(addName);
						currentServer.add("username");
						currentServer.add("password");
						currentServer.add("/usr/local/bin/NairnMPM");
						currentServer.add("/usr/local/bin/NairnMPM.dtd");
						currentServer.add("/usr/local/bin/NairnFEA");
						currentServer.add("/usr/local/bin/NairnFEA.dtd");
						
						// save in preferences (change list of servers and key for that server)
						servers.add(addName);
						setPrefList(RemoteServerListKey,servers,serverDelim);
						setPrefList(addName,currentServer,serverDelim);
						
						// select it (which will get it inserted)
						remoteList.insertItemAt(new PlotMenuItem(addName),numItems-2);
						remoteList.setSelectedIndex(numItems-2);
					}
					else if(newIndex==numItems-1)
					{	// delete serve data and from server list
						if(currentServerIndex>0)
						{	String delServer = currentServer.get(ServerIP);
							prefs.remove(delServer);
							servers.remove(currentServerIndex-1);
							setPrefList(RemoteServerListKey,servers,serverDelim);
							remoteList.removeItemAt(currentServerIndex);
							currentServerIndex = -1;
							currentServer.clear();
							remoteList.setSelectedIndex(0);
						}
					}
					else if(newIndex==0)
					{	// pick zero item
						currentServerIndex = -1;
						currentServer.clear();
						insertServer();
					}
					else if(newIndex!=currentServerIndex)
					{	// switch to an existing server
						PlotMenuItem pm=(PlotMenuItem)remoteList.getSelectedItem();
						String serverName = pm.toString();
						if(newIndex!=currentServerIndex);
						{	currentServer = getPrefList(serverName,"",serverDelim);
							currentServerIndex = newIndex;
							insertServer();
						}
					}
				}
			}
		});

		// Server name --------------------------
		String toolTip = "Enter server address (e.g. 'mpm.fsl.orst.edu')";
		remoteLine(panel3,gridbag,c,"Server Name:",remoteServerAddr,RemoteServerKey,toolTip);
		
		// User name and password line (empty for now) -------------------------
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 1;
		label = new JLabel("Username:");
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		panel3.add(label);

		c.gridx=1;
		c.weightx = 1.5;
		c.gridwidth = 1;
		remoteUsername.setActionCommand(RemoteUserKey);
		remoteUsername.addActionListener(this);
		remoteUsername.addFocusListener(new PrefFocusListener(RemoteUserKey));
		remoteUsername.setToolTipText("Enter user name for server '"+remoteServerAddr.getText()+"'");
		gridbag.setConstraints(remoteUsername,c);
		panel3.add(remoteUsername);
		
		c.gridx=2;
		c.weightx = 0.;
		c.gridwidth = 1;
		label = new JLabel("Password:");
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		panel3.add(label);
		
		c.gridx=3;
		c.weightx = 1.5;
		c.gridwidth = 1;
		remoteUserpass.setActionCommand(RemoteUserPassKey);
		remoteUserpass.addActionListener(this);
		remoteUserpass.addFocusListener(new PrefFocusListener(RemoteUserPassKey));
		remoteUserpass.setToolTipText("Enter password for user '"+remoteUsername.getText()+"' on server '"+remoteServerAddr.getText()+"'");
		gridbag.setConstraints(remoteUserpass,c);
		panel3.add(remoteUserpass);
		
		// MPM Code -------------------------------------------
		toolTip = "Enter full path to MPM executeable (can use '~' to indicate home directory)";
		remoteLine(panel3,gridbag,c,"MPM Code Path:",mpmCodePathFld,RemoteMPMPathKey,toolTip);
		toolTip = "Enter full path to MPM dtd file (cannot use '~' to indicate home directory)";
		remoteLine(panel3,gridbag,c,"MPM DTD Path:",mpmDTDPathFld,RemoteMPMDTDKey,toolTip);
		
		// FEA Code -------------------------------------------
		toolTip = "Enter full path to FEA executeable (can use '~' to indicate home directory)";
		remoteLine(panel3,gridbag,c,"FEA Code Path:",feaCodePathFld,RemoteFEAPathKey,toolTip);
		toolTip = "Enter full path to FEA dtd file (cannot use '~' to indicate home directory)";
		remoteLine(panel3,gridbag,c,"FEA DTD Path:",feaDTDPathFld,RemoteFEADTDKey,toolTip);
		
		// Remove/local buttons
		ButtonGroup execLocation = new ButtonGroup();
		rdbtnExecRemote.addActionListener(this);
		rdbtnExecRemote.setFocusable(false);
		rdbtnExecLocal.addActionListener(this);
		rdbtnExecLocal.setFocusable(false);
		rdbtnExecLocal.setActionCommand("code exec location");
		execLocation.add(rdbtnExecLocal);
		GridBagConstraints gbc_rdbtnExecLocal = new GridBagConstraints();
		gbc_rdbtnExecLocal.fill = GridBagConstraints.BOTH;
		gbc_rdbtnExecLocal.insets = new Insets(18, 20, 3, 5);
		gbc_rdbtnExecLocal.gridx = 0;
		gbc_rdbtnExecLocal.gridy = 8;
		panel3.add(rdbtnExecLocal, gbc_rdbtnExecLocal);
		rdbtnExecRemote.setActionCommand("code exec location");
		execLocation.add(rdbtnExecRemote);
		GridBagConstraints gbc_rdbtnExecRemote = new GridBagConstraints();
		gbc_rdbtnExecRemote.insets = new Insets(0, 20, 0, 5);
		gbc_rdbtnExecRemote.fill = GridBagConstraints.BOTH;
		gbc_rdbtnExecRemote.gridx = 0;
		gbc_rdbtnExecRemote.gridy = 9;
		panel3.add(rdbtnExecRemote, gbc_rdbtnExecRemote);
		
		if(restoreRemoteMode())
			rdbtnExecRemote.setSelected(true);
		else
			rdbtnExecLocal.setSelected(true);
		
		// empty fill space on the bottom
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 4;
		c.weighty = 10;
		c.fill = GridBagConstraints.VERTICAL;
		label = new JLabel(" ");
		gridbag.setConstraints(label,c);
		panel3.add(label);

		return panel3;
	}
	
	// build panel for remote connections
	public JPanel colorPanel()
	{	JPanel panel2 = new JPanel();

		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints pc = new GridBagConstraints(); // for the pane
		panel2.setLayout(gridbag);

		// spectrum box
		GridBagLayout spectrumBag = new GridBagLayout();
		JNBoxPanel spectrumPan = new JNBoxPanel("Color Plot Spectrum",spectrumBag);
		GridBagConstraints c = new GridBagConstraints(); // for the box panel

		JLabel section = new JLabel("Type: ");
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.WEST;
		c.gridx = 0;
		c.gridwidth = 1;
		c.weightx = 0.0;
		c.insets = new Insets(6, 6, 1, 0); // tlbr
		spectrumBag.setConstraints(section, c);
		spectrumPan.add(section);

		String[] lines = { "Rainbow 1 (Purple to Red)",
				"Rainbow 2 (Blue to Red)", "Rainbow 3 (Purple to Orange)",
				"Grayscale (Black to White)", "Cool Diverging (Blue to Red)",
				"Viridis", "Magma", "Inferno", "Plasma", "Turbo" };
		spectrumBox = new JComboBox<String>(lines);
		spectrumBox.setSelectedIndex(ColorPicker.getSpectrumType());
		spectrumBox.addActionListener(this);
		spectrumBox.setFocusable(false);
		spectrumBox.setToolTipText("Select one of the available spectra for coloring the plots.");
		c.insets = new Insets(6, 3, 1, 6); // tlbr
		c.gridx = 1;
		c.weightx = 1.0;
		spectrumBag.setConstraints(spectrumBox, c);
		spectrumPan.add(spectrumBox);

		JLabel section1 = new JLabel("Contours: ");
		c.insets = new Insets(6, 6, 1, 0); // tlbr
		c.gridx = 0;
		c.weightx = 0.0;
		spectrumBag.setConstraints(section1, c);
		spectrumPan.add(section1);

		numContours.setText("" + ColorPicker.getNumberOfContours());
		numContours.setEditable(true);
		numContours.setColumns(5);
		numContours.setToolTipText("Enter 1 for continuous colors or >1 for number of discrete colors and hit return or enter");
		numContours.addActionListener(this);
		numContours.setActionCommand("change contours");
		numContours.addFocusListener(new PrefFocusListener(NumContoursKey));
		c.insets = new Insets(6, 3, 1, 6); // tlbr
		c.gridx = 1;
		c.weightx = 1.0;
		spectrumBag.setConstraints(numContours, c);
		spectrumPan.add(numContours);

		JLabel section2 = new JLabel("Subelements: ");
		c.insets = new Insets(6, 6, 1, 0); // tlbr
		c.gridx = 0;
		c.weightx = 0.0;
		spectrumBag.setConstraints(section2, c);
		spectrumPan.add(section2);

		numSubelements.setText("" + ElementBase.getSubelementDensity());
		numSubelements.setEditable(true);
		numSubelements.setColumns(5);
		numSubelements.setToolTipText("Enter number of subelements (nXn) to use for mesh plots.");
		numSubelements.addActionListener(this);
		numSubelements.setActionCommand("change subelement density");
		numSubelements.addFocusListener(new PrefFocusListener(NumSubelementsKey));
		c.insets = new Insets(6, 3, 1, 6); // tlbr
		c.gridx = 1;
		c.weightx = 1.0;
		spectrumBag.setConstraints(numSubelements, c);
		spectrumPan.add(numSubelements);

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridwidth = 1;
		c.weightx = 1.0;

		// add spectrum box to pan
		pc.fill = GridBagConstraints.HORIZONTAL;
		pc.gridx = 0;
		pc.gridwidth = 1;
		pc.weightx = 1.0;
		pc.weighty = 0.;
		pc.insets = new Insets(1, 6, 1, 6); // tlbr
		pc.anchor = GridBagConstraints.CENTER;
		gridbag.setConstraints(spectrumPan, pc);
		panel2.add(spectrumPan);

		// plot element color wells box
		GridBagLayout wellsBag = new GridBagLayout();
		JNBoxPanel wellsPan = new JNBoxPanel("Plot Elements", wellsBag);
		GridBagConstraints ec = new GridBagConstraints(); // for the plot

		JLabel elab1 = new JLabel("Background: ");
		ec.fill = GridBagConstraints.NONE;
		ec.anchor = GridBagConstraints.EAST;
		ec.insets = new Insets(10, 0, 0, 0); // tlbr
		ec.gridx = 0;
		ec.weightx = 0.5;
		ec.gridwidth = 1;
		wellsBag.setConstraints(elab1, ec);
		wellsPan.add(elab1);

		JNColorWell bgWell = new JNColorWell(getPrefColor(backColorKey,backColorDef));
		bgWell.addActionListener(this, "bgColor", "Choose Background Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 1;
		ec.weightx = 1.0;
		wellsBag.setConstraints(bgWell, ec);
		wellsPan.add(bgWell);

		JLabel elab2 = new JLabel("Mesh Lines: ");
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx = 2;
		ec.weightx = 0.5;
		wellsBag.setConstraints(elab2, ec);
		wellsPan.add(elab2);

		JNColorWell mlWell = new JNColorWell(getPrefColor(meshLineColorKey,meshLineColorDef));
		mlWell.addActionListener(this, "mlColor", "Choose Mesh Line Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 3;
		ec.weightx = 1.0;
		wellsBag.setConstraints(mlWell, ec);
		wellsPan.add(mlWell);

		JLabel elab3 = new JLabel("Labels: ");
		ec.insets = new Insets(10, 0, 0, 0); // tlbr
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx = 0;
		ec.weightx = 0.5;
		ec.gridwidth = 1;
		wellsBag.setConstraints(elab3, ec);
		wellsPan.add(elab3);

		JNColorWell labWell = new JNColorWell(getPrefColor(textColorKey,textColorDef));
		labWell.addActionListener(this, "labColor","Choose Color for Text Labels");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 1;
		ec.weightx = 1.0;
		wellsBag.setConstraints(labWell, ec);
		wellsPan.add(labWell);

		JLabel elab4 = new JLabel("Nodal Points: ");
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx = 2;
		ec.weightx = 0.5;
		wellsBag.setConstraints(elab4, ec);
		wellsPan.add(elab4);

		JNColorWell npWell = new JNColorWell(getPrefColor(meshNodesColorKey,meshNodesColorDef));
		npWell.addActionListener(this, "npColor", "Choose Nodal Point Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 3;
		ec.weightx = 1.0;
		wellsBag.setConstraints(npWell, ec);
		wellsPan.add(npWell);

		JLabel elab5 = new JLabel("Disp. Mesh Lines: ");
		ec.fill = GridBagConstraints.NONE;
		ec.anchor = GridBagConstraints.EAST;
		ec.insets = new Insets(10, 0, 0, 0); // tlbr
		ec.gridx = 0;
		ec.weightx = 0.5;
		ec.gridwidth = 1;
		wellsBag.setConstraints(elab5, ec);
		wellsPan.add(elab5);

		JNColorWell dmlWell = new JNColorWell(getPrefColor(dispMeshLineColorKey,dispMeshLineColorDef));
		dmlWell.addActionListener(this, "dmlColor", "Choose Displaced Mesh Line Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 1;
		ec.weightx = 1.0;
		wellsBag.setConstraints(dmlWell, ec);
		wellsPan.add(dmlWell);

		JLabel elab6 = new JLabel("Crack Planes: ");
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx = 2;
		ec.weightx = 0.5;
		wellsBag.setConstraints(elab6, ec);
		wellsPan.add(elab6);

		JNColorWell cpWell = new JNColorWell(getPrefColor(cplaneColorKey,cplaneColorDef));
		cpWell.addActionListener(this, "cpColor", "Choose Crack Plane Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 3;
		ec.weightx = 1.0;
		wellsBag.setConstraints(cpWell, ec);
		wellsPan.add(cpWell);

		JLabel elab7 = new JLabel("Above Crack Planes: ");
		ec.fill = GridBagConstraints.NONE;
		ec.anchor = GridBagConstraints.EAST;
		ec.insets = new Insets(10, 0, 6, 0); // tlbr
		ec.gridx = 0;
		ec.weightx = 0.5;
		ec.gridwidth = 1;
		wellsBag.setConstraints(elab7, ec);
		wellsPan.add(elab7);

		JNColorWell caboveWell = new JNColorWell(getPrefColor(caboveColorKey,caboveColorDef));
		caboveWell.addActionListener(this, "caboveColor", "Choose Crack Above Surface Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 1;
		ec.weightx = 1.0;
		wellsBag.setConstraints(caboveWell, ec);
		wellsPan.add(caboveWell);

		JLabel elab8 = new JLabel("Below Crack Planes: ");
		ec.anchor = GridBagConstraints.EAST;
		ec.gridx = 2;
		ec.weightx = 0.5;
		wellsBag.setConstraints(elab8, ec);
		wellsPan.add(elab8);

		JNColorWell cbelowWell = new JNColorWell(getPrefColor(cbelowColorKey,cbelowColorDef));
		cbelowWell.addActionListener(this, "cbelowColor", "Choose Crack Below Surface Color");
		ec.anchor = GridBagConstraints.WEST;
		ec.gridx = 3;
		ec.weightx = 1.0;
		wellsBag.setConstraints(cbelowWell, ec);
		wellsPan.add(cbelowWell);

		// add plot element color wells box to pane
		gridbag.setConstraints(wellsPan, pc);
		panel2.add(wellsPan);
		
		// empty fill on the bottom
		pc.fill = GridBagConstraints.VERTICAL;
		pc.weighty = 10.;
		JLabel label = new JLabel(" ");
		gridbag.setConstraints(label,pc);
		panel2.add(label);
		
		return panel2;
	}
	
	// build panel for remote connections
	public JPanel fontPanel()
	{	JPanel panel4 = new JPanel();

		GridBagLayout panebag = new GridBagLayout();
		GridBagConstraints pc = new GridBagConstraints(); // for the pane
		panel4.setLayout(panebag);

		// font box
		GridBagLayout gridbag = new GridBagLayout();
		JNBoxPanel fontPan = new JNBoxPanel("Window Fonts",gridbag);
		GridBagConstraints c = new GridBagConstraints(); // for the box panel
		
		// Script font
		families = makeFontNameMenu(prefs.get(ScriptFontKey, ScriptFontDef));
		sizes = makeFontSizesMenu(prefs.getInt(ScriptFontSizeKey, ScriptFontSizeDef));
		makeFontFields(c,"Editing:",gridbag,fontPan,families,sizes);

		// Output font
		ofamilies = makeFontNameMenu(prefs.get(OutputFontKey, OutputFontDef));
		osizes = makeFontSizesMenu(prefs.getInt(OutputFontSizeKey, OutputFontSizeDef));
		makeFontFields(c,"Output:",gridbag,fontPan,ofamilies,osizes);

		// Results font
		rfamilies = makeFontNameMenu(prefs.get(ResultsFontKey, ResultsFontDef));
		rsizes = makeFontSizesMenu(prefs.getInt(ResultsFontSizeKey, ResultsFontSizeDef));
		makeFontFields(c,"Results:",gridbag,fontPan,rfamilies,rsizes);

		// add spectrum box to pan
		pc.fill = GridBagConstraints.HORIZONTAL;
		pc.gridx = 0;
		pc.gridwidth = 1;
		pc.weightx = 1.0;
		pc.weighty = 0.;
		pc.insets = new Insets(1, 6, 1, 6); // tlbr
		pc.anchor = GridBagConstraints.CENTER;
		panebag.setConstraints(fontPan, pc);
		panel4.add(fontPan);
		
		// empty fill on the bottom
		pc.fill = GridBagConstraints.VERTICAL;
		pc.weighty = 10.;
		JLabel label = new JLabel(" ");
		panebag.setConstraints(label,pc);
		panel4.add(label);
		
		return panel4;
	}
	
	protected void makeFontFields(GridBagConstraints c,String title,GridBagLayout gridbag,
			JNBoxPanel fontPan,JComboBox<String> fontNames,JComboBox<Integer> fontSizes)
	{
		// label for type of font
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.WEST;
		c.gridx=0;
		c.gridwidth = 1;
		c.weightx = 0.0;
		c.insets=new Insets(6, 6, 1, 0);			// tlbr
		JLabel nameLabel=new JLabel(title);
		gridbag.setConstraints(nameLabel,c);
		fontPan.add(nameLabel);
		
		// menu of font families
		c.gridx++;
		c.weightx = 1.0;
		c.insets = new Insets(6, 3, 1, 6);
		gridbag.setConstraints(fontNames,c);
		fontPan.add(fontNames);
		
		// menu of font sizes
		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.CENTER;
		c.gridx++;
		c.insets=new Insets(6, 3, 1, 6);			// tlbr
		gridbag.setConstraints(fontSizes,c);
		fontPan.add(fontSizes);		
	}
	
	// create combo box for font names
	protected JComboBox<String> makeFontNameMenu(String eFontName)
	{	// font menu in inspector
		String [] names=GraphicsEnvironment.getLocalGraphicsEnvironment().getAvailableFontFamilyNames(null);
		JComboBox<String> fNames=new JComboBox<String>(names);
		fNames.addActionListener(this);
		fNames.setFocusable(false);
		fNames.setActionCommand("NewFontName");
		
		// select font
		int i;
		for(i=0;i<fNames.getItemCount();i++)
		{	if(((String)(fNames.getItemAt(i))).equals(eFontName))
			{	fNames.setSelectedIndex(i);
				break;
			}
		}
		
		return fNames;
	}

	// create combo box
	protected JComboBox<Integer> makeFontSizesMenu(int eFontSize)
	{	// font menu in inspector
		int i;
		JComboBox<Integer> fSizes=new JComboBox<Integer>();
		int [] allowSizes={8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 24, 28, 32, 36 };
		for(i=0;i<allowSizes.length;i++)
			fSizes.addItem(new Integer(allowSizes[i]));
		fSizes.addActionListener(this);
		fSizes.setFocusable(false);
		fSizes.setActionCommand("NewFontSize");
		
		// select size
		for(i=0;i<fSizes.getItemCount();i++)
		{	if(((Integer)(fSizes.getItemAt(i))).intValue()==eFontSize)
			{	fSizes.setSelectedIndex(i);
				break;
			}
		}

		return fSizes;
	}
	
	// line on remote panel with label and text field
	public void remoteLine(JPanel panel3,GridBagLayout gridbag,GridBagConstraints c,String tlab,
			JTextField theFld,String aPrefKey,String toolTip)
	{
		c.gridx=0;
		c.weightx = 0.;
		c.gridwidth = 1;
		JLabel label = new JLabel(tlab);
		label.setHorizontalAlignment(JLabel.RIGHT);
		gridbag.setConstraints(label,c);
		panel3.add(label);
		
		c.gridx=1;
		c.weightx = 3.0;
		c.gridwidth = 3;
		theFld.setActionCommand(aPrefKey);
		theFld.addActionListener(this);
		theFld.addFocusListener(new PrefFocusListener(aPrefKey));
		theFld.setToolTipText(toolTip);
		gridbag.setConstraints(theFld,c);
		panel3.add(theFld);
	}
	
	// file remote fields with current server
	public void insertServer()
	{	if(currentServer.size()>=7)
		{	remoteServerAddr.setText(currentServer.get(ServerIP));
			remoteUsername.setText(currentServer.get(ServerUser));
			remoteUserpass.setText(currentServer.get(ServerPWD));
			mpmCodePathFld.setText(currentServer.get(ServerMPM));
			mpmDTDPathFld.setText(currentServer.get(ServerMPMDTD));
			feaCodePathFld.setText(currentServer.get(ServerFEA));
			feaDTDPathFld.setText(currentServer.get(ServerFEADTD));
		}
		else
		{	remoteServerAddr.setText("");
			remoteUsername.setText("");
			remoteUserpass.setText("");
			mpmCodePathFld.setText("");
			mpmDTDPathFld.setText("");
			feaCodePathFld.setText("");
			feaDTDPathFld.setText("");
		}
	}
	
	// check if new name is unique
	public boolean uniqueServer(String addr)
	{	for(int i=0;i<servers.size();i++)
			if(addr.equals(servers.get(i))) return false;
		return true;
	}

	// Create panel for file path entry
	public JPanel filePathPanel(String pathName, JTextField pathField,
			String pathText, boolean validPref, JCheckBox validate)
	{	JPanel filePanel = new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		filePanel.setLayout(gridbag);

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridwidth = 1;
		c.weightx = 1.0;
		c.insets = new Insets(0, 6, 0, 0); // tlbr
		JLabel nameLabel = new JLabel(pathName + " Path:");
		gridbag.setConstraints(nameLabel, c);
		filePanel.add(nameLabel);

		int rowWidth = 2;
		if (validate != null)
		{	c.gridx++;
			validate.setSelected(validPref);
			validate.setActionCommand(pathName + " Use");
			validate.addActionListener(this);
			validate.setToolTipText("Check to validate input files using the provided DTD file");
			gridbag.setConstraints(validate, c);
			filePanel.add(validate);
			rowWidth++;
		}
		else if(pathName.equals("MPM Code"))
		{	c.gridx++;
			c.fill = GridBagConstraints.NONE;
			JButton bundle = new JButton("Use Bundle");
			bundle.setActionCommand("Use Bundle");
			bundle.addActionListener(this);
			bundle.setToolTipText("Click to use binaries in the bundle folder");
			bundle.setFocusable(false);
			gridbag.setConstraints(bundle, c);
			filePanel.add(bundle);
			rowWidth++;
		}

		c.fill = GridBagConstraints.NONE;
		c.anchor = GridBagConstraints.EAST;
		c.gridx++;
		c.insets = new Insets(0, 0, 0, 6); // tlbr
		JButton change = new JButton("Change...");
		if(NairnFEAMPMViz.isWindowsOS() && pathName.contentEquals("Shell Command"))
		{	change.setEnabled(false);
			change.setToolTipText("Shell command not used in Windows");
		}
		else
		{	change.setActionCommand(pathName);
			change.addActionListener(this);
			change.setToolTipText("Click to change the " + pathName + " path");
		}
		change.setFocusable(false);
		gridbag.setConstraints(change, c);
		filePanel.add(change);

		c.anchor = GridBagConstraints.CENTER;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridwidth = rowWidth;
		c.insets = new Insets(0, 6, 0, 9); // tlbr
		pathField.setText(pathText);
		pathField.setEditable(false);
		pathField.setColumns(45);
		pathField.setToolTipText("The current " + pathName
				+ " path; click 'Change...' to change it.");
		gridbag.setConstraints(pathField, c);
		filePanel.add(pathField);

		return filePanel;
	}

	//----------------------------------------------------------------------------
	// handle application commands
	//----------------------------------------------------------------------------
	
	public void actionPerformed(ActionEvent e)
	{   String theCmd=e.getActionCommand();

		if(theCmd.equals("MPM Code"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(bundle)")>=0)
				newPath = "$(bundle)";
			mpmCodePath.setText(newPath);
			prefs.put(NairnMPMKey,newPath);
		}
		
		else if(theCmd.equals("MPM DTD File Use"))
		{	prefs.putBoolean(NairnMPMValidateKey,validateMPM.isSelected());
		}
		
		else if(theCmd.equals("MPM DTD File"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(bundle)")>=0)
				newPath = "$(bundle)";
			mpmDTDPath.setText(newPath);
			prefs.put(NairnMPMDTDKey,newPath);
		}
		
		else if(theCmd.equals("FEA Code"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(bundle)")>=0)
				newPath = "$(bundle)";
			feaCodePath.setText(newPath);
			prefs.put(NairnFEAKey,newPath);
		}
		
		else if(theCmd.equals("FEA DTD File Use"))
		{	prefs.putBoolean(NairnFEAValidateKey,validateFEA.isSelected());
		}
		
		else if(theCmd.equals("Use Bundle"))
		{	String newPath = "$(bundle)";
			mpmCodePath.setText(newPath);
			prefs.put(NairnMPMKey,newPath);
			mpmDTDPath.setText(newPath);
			prefs.put(NairnMPMDTDKey,newPath);
			feaCodePath.setText(newPath);
			prefs.put(NairnFEAKey,newPath);
			feaDTDPath.setText(newPath);
			prefs.put(NairnFEADTDKey,newPath);
		}
		
		else if(theCmd.equals("bgColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(backColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("labColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(textColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("mlColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(meshLineColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("npColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(meshNodesColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("dmlColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(dispMeshLineColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("cpColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(cplaneColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("caboveColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(caboveColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("cbelowColor"))
		{	JNColorWell cWell=(JNColorWell)e.getSource();
			setPrefColor(cbelowColorKey,cWell.getWellColor());
			return;
		}
		
		else if(theCmd.equals("FEA DTD File"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(bundle)")>=0)
				newPath = "$(bundle)";
			feaDTDPath.setText(newPath);
			prefs.put(NairnFEADTDKey,newPath);
		}
		
		else if(theCmd.equals("Shell Command"))
		{	int result=chooser.showOpenDialog(this);
			if(result==JFileChooser.CANCEL_OPTION) return;
			String newPath=chooser.getSelectedFile().getPath();
			if(newPath.indexOf("$(windows)")>=0)
				newPath = "$(windows)";
			shellPath.setText(newPath);
			prefs.put(ShellKey,newPath);
		}
		
		else if(theCmd.equals("Work Space Directory"))
		{	chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int result=chooser.showOpenDialog(this);
			String newPath=chooser.getSelectedFile().getPath();
			chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			if(result==JFileChooser.CANCEL_OPTION) return;
			workPath.setText(newPath);
			prefs.put(WorkSpaceKey,newPath);
			setWorkspace(NairnFEAMPMViz.GetAppChooser());
		}
		
		else if(theCmd.equals("change contours"))
		{	try
			{	int newContours=Integer.valueOf(numContours.getText());
				if(newContours<1) newContours=1;
				int oldNum=ColorPicker.getNumberOfContours();
				if(newContours!=oldNum)
				{	prefs.putInt(NumContoursKey,newContours);
					ColorPicker.setNumberOfContours();
				}
			}
			catch(Exception ie)
			{	Toolkit.getDefaultToolkit().beep();
			}
		}

		else if(theCmd.equals("change subelement density"))
		{	try
			{	int newDensity=Integer.valueOf(numSubelements.getText());
				if(newDensity<1) newDensity=1;
				int oldNum=ElementBase.getSubelementDensity();
				if(newDensity!=oldNum)
				{	prefs.putInt(NumSubelementsKey,newDensity);
					ElementBase.setSubelementDensity();
				}
			}
			catch(Exception ie)
			{	Toolkit.getDefaultToolkit().beep();
			}
		}

		else if (theCmd.equalsIgnoreCase("code exec location"))
		{	if (this.rdbtnExecLocal.isSelected())
			{	prefs.put(CodeExecLocationKey, "local");
				currentRemoteMode = false;
			}
			else if (this.rdbtnExecRemote.isSelected())
			{	prefs.put(CodeExecLocationKey, "remote");
				currentRemoteMode = true;
			}
		}
		
		else if(theCmd.equalsIgnoreCase("NewFontName"))
		{	JComboBox<?> cb = (JComboBox<?>)e.getSource();
			if(cb==families)
			{	prefs.put(ScriptFontKey,(String)cb.getSelectedItem());
			}
			else if(cb==ofamilies)
			{	prefs.put(OutputFontKey,(String)cb.getSelectedItem());
			}
			else if(cb==rfamilies)
			{	prefs.put(ResultsFontKey,(String)cb.getSelectedItem());
			}
		}
		
		else if(theCmd.equalsIgnoreCase("NewFontSize"))
		{	JComboBox<?> cb = (JComboBox<?>)e.getSource();
			if(cb==sizes)
			{	prefs.putInt(ScriptFontSizeKey,((Integer)cb.getSelectedItem()).intValue());
			}
			else if(cb==osizes)
			{	prefs.putInt(OutputFontSizeKey,((Integer)cb.getSelectedItem()).intValue());
			}
			else if(cb==rsizes)
			{	prefs.putInt(ResultsFontSizeKey,((Integer)cb.getSelectedItem()).intValue());
			}
		}
		
		else
		{	// the rest. Check remote first
			boolean valid = changeRemoteOption(theCmd);
			
			if(!valid)
			{	JComboBox<?> cb = (JComboBox<?>) e.getSource();
				if(cb == spectrumBox)
				{	int oldType = ColorPicker.getSpectrumType();
					int newType = cb.getSelectedIndex();
					if(newType != oldType)
					{	prefs.putInt(SpectrumKey, newType);
						ColorPicker.setSpectrumType();
					}
				}
				else
					System.out.println("Unrecognized preferences commane: "+theCmd);
			}
		}
	}
	
	// action or focus change for remote field
	public boolean changeRemoteOption(String theCmd)
	{
		if (theCmd.equals(RemoteServerKey))
		{	if(currentServerIndex<1) return true;
		
			// if name not changed, do nothing
			String newName = remoteServerAddr.getText();
			if(newName.equals(currentServer.get(ServerIP))) return true;
			
			// make sure unique name
			if(!uniqueServer(newName))
			{	remoteServerAddr.setText(currentServer.get(ServerIP));
				JNUtilities.showMessage(this, "The new server name is not valid. It matches an existing server.");
				return true;
			}
			
			// remove old key (new one set below)
			prefs.remove(currentServer.get(ServerIP));
			
			// change name and menu item (committed below)
			currentServer.set(ServerIP,newName);	
			PlotMenuItem pm=(PlotMenuItem)remoteList.getItemAt(currentServerIndex);
			pm.setString(newName);
			remoteList.repaint();
			
			// update servers list
			servers.set(currentServerIndex-1, newName);
			setPrefList(RemoteServerListKey,servers,serverDelim);
		}

		if(theCmd.equals(RemoteUserKey))
		{	if(currentServerIndex<1) return true;
			currentServer.set(ServerUser, remoteUsername.getText());	
		}

		else if (theCmd.equals(RemoteUserPassKey))
		{	if(currentServerIndex<1) return true;
			currentServer.set(ServerPWD,new String(remoteUserpass.getPassword()));
		}

		else if (theCmd.equals(RemoteMPMPathKey))
		{	if(currentServerIndex<1) return true;
			currentServer.set(ServerMPM,mpmCodePathFld.getText());
		}

		else if (theCmd.equals(RemoteMPMDTDKey))
		{	if(currentServerIndex<1) return true;
			currentServer.set(ServerMPMDTD,mpmDTDPathFld.getText());
		}

		else if (theCmd.equals(RemoteFEAPathKey))
		{	if(currentServerIndex<1) return true;
			currentServer.set(ServerFEA,feaCodePathFld.getText());
		}

		else if (theCmd.equals(RemoteFEADTDKey))
		{	if(currentServerIndex<1) return true;
			currentServer.set(ServerFEADTD,feaDTDPathFld.getText());
		}
		
		else
			return false;
		
		// update now
		setPrefList(currentServer.get(0),currentServer,serverDelim);
		return true;
	}

	// ----------------------------------------------------------------------------
	// class methods
	// ----------------------------------------------------------------------------

	// create preferences storage object
	public static void initializePrefs()
	{	createPrefs("com/geditcom/NairnFEAMPMViz");
		restoreRemoteMode();
	
		// default close action
		NFMVPrefs.setDefaultQuitCloseLastWindow(false);

		// restore open docs
		NFMVPrefs.setDefaultRestoreOpenDocs(true);
		
		// default settings
		ColorPicker.setNumberOfContours();
		ColorPicker.setSpectrumType();
		ElementBase.setSubelementDensity();
	}

	// set a chooser to current work space
	public static void setWorkspace(JFileChooser chooser)
	{	String workPath=prefs.get(WorkSpaceKey,WorkSpaceDef);
		if(workPath.length()>0)
		{	try
			{	File workDirectory=new File(workPath);
				chooser.setCurrentDirectory(workDirectory);
			}
			catch(Exception e) {}
		}
	}
	
	// return true it set to run remotely
	// REMOTE_ACCESS
	public static boolean getRemoteMode() { return currentRemoteMode; }
	public static void setRemoteMode(boolean tempRemote) { currentRemoteMode = tempRemote; }
	public static boolean restoreRemoteMode()
	{	if(prefs!=null)
			currentRemoteMode = prefs.get(CodeExecLocationKey, CodeExecLocationDef).equalsIgnoreCase("remote");
		else
			currentRemoteMode = false;
		return currentRemoteMode;
	}
	public static ArrayList<String> getServerInfo()
	{	NFMVPrefs prefWindow = (NFMVPrefs)NFMVPrefs.appPreferencesWindow();
		if(prefWindow!=null)
		{	if(prefWindow.currentServerIndex>0)
				return prefWindow.currentServer;
		}
		
		// use the first one
		ArrayList<String> availableServers = getPrefList(RemoteServerListKey,RemoteServerListDef,serverDelim);
		if(availableServers.size()<1) return null;
		
		// look up the first one
		String firstServer = availableServers.get(0);
		ArrayList<String> useServer = getPrefList(firstServer,"",serverDelim);
		if(useServer.size()<7) return null;

		return useServer;
	}
	
	// listenere for fields tied to string preferences
	public class PrefFocusListener implements FocusListener
	{
		private String prefKey;
		
		public PrefFocusListener(String aKey)
		{	prefKey = aKey;
		}
		
		public void focusLost(FocusEvent e)
		{	if(prefKey.equals(NumContoursKey))
			{	try
				{	int newContours=Integer.valueOf(numContours.getText());
					if(newContours<1) newContours=1;
					int oldNum=ColorPicker.getNumberOfContours();
					if(newContours!=oldNum)
					{	prefs.putInt(NumContoursKey,newContours);
						ColorPicker.setNumberOfContours();
					}
				}
				catch(Exception ie)
				{	Toolkit.getDefaultToolkit().beep();
				}
			}
			else if(prefKey.equals(NumSubelementsKey))
			{	try
				{	int newDensity=Integer.valueOf(numSubelements.getText());
					if(newDensity<1) newDensity=1;
					int oldNum=ElementBase.getSubelementDensity();
					if(newDensity!=oldNum)
					{	prefs.putInt(NumSubelementsKey,newDensity);
						ElementBase.setSubelementDensity();
					}
				}
				catch(Exception ie)
				{	Toolkit.getDefaultToolkit().beep();
				}
			}
			else
			{	if(currentServerIndex>0)
					changeRemoteOption(prefKey);
			}
		}
		
		public void focusGained(FocusEvent e)
		{
		}
	}

}
