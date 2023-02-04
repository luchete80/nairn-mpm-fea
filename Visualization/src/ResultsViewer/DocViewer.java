/*
 * DocViewer.java
 * NairnFEAMPMViz Application
 * 
 * Created by John Nairn on Fri Mar 05 2004.
 * Copyright (c) 2004 RSAC Software. All rights reserved.
 */

import java.awt.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import javax.swing.*;
import geditcom.JNFramework.*;

public class DocViewer extends JNDocument
{
	private static final long serialVersionUID = 1L;

	// instance variables
	
	public ResultsDocument resDoc=new ResultsDocument();		// one results document in this viewer
	private TextDisplay display;		// split view with sections and text of sections
	public ControlPanel controls;		// controls to set plot options
	private JSplitPane full;			// window split view - TextDisplay on top, ControlPanel on bottom
	private boolean startWithMesh=false;
	private CmdViewer commandsWindow=null;

	// Actions
	private ShowPartnerAction showPartnerCommand = new ShowPartnerAction();
	private ScaleResultsAction scaleResultsCommand = new ScaleResultsAction();
	private StartPlotAction startPlotCommand = new StartPlotAction();
	private ExportVTKAction exportVTKAction = new ExportVTKAction();
	
	public DocViewer(boolean openMesh)
	{   super("DocViewer");
	
		setFramePrefs("Results Window Width",600,"Results Window Height",800);
	
		makeMenuBar();
		
		Container content=getContentPane( );
		resDoc.setDocController(this);
		startWithMesh=openMesh;
		
		// tool bar icons
		CmdViewer.addDefaultToolBar(this);
		
		addToolBarBar();
		Class<?> baseClass=JNApplication.main.getClass();
		ImageIcon newMPM=new ImageIcon(baseClass.getResource("Resources/image-x-generic.png"));
		addToolBarIcon(newMPM,null,"Open plot using the currently selected options.",startPlotCommand);
		ImageIcon showCmds=new ImageIcon(baseClass.getResource("Resources/commands-editor.png"));
		addToolBarIcon(showCmds,null,"Show associated commands window (if available).",showPartnerCommand);

		// top displays sections and text of selected section
		Font theFont = new Font(NFMVPrefs.prefs.get(NFMVPrefs.ResultsFontKey,
				NFMVPrefs.ResultsFontDef),Font.PLAIN,
				NFMVPrefs.prefs.getInt(NFMVPrefs.ResultsFontSizeKey,
						NFMVPrefs.ResultsFontSizeDef));
		display=new TextDisplay(resDoc,theFont);
		
		// control area
		controls=new ControlPanel(this);
		JScrollPane bottom=new JScrollPane(controls);
		
		// full window split pane
		full=new JSplitPane(JSplitPane.VERTICAL_SPLIT,display,bottom);
		
		// add to content
		content.add(full,BorderLayout.CENTER);
		
		// size and location
		finishFrameworkWindow(true);
	}
	
	// make menu bar on launch
	protected void makeMenuBar()
	{
		// Create the menu bar
		JMenuBar menuBar = new JMenuBar();
		if(!JNApplication.isMacLNF())
			menuBar.add(defaultApplicationMenu());		// Application menu
		menuBar.add(CmdViewer.defaultFileMenu(this));				// File menu
		NairnFEAMPMViz.addExamplesMenu(menuBar,"File");
		
		// Edit menu
		JMenu menu = new JMenu("Edit");
		menuBar.add(menu);
		menu.add(JNApplication.main.getOpenFindAction());
		menu.add(JNApplication.main.getFindAgainAction());
		
		// Analyze menu
		menu = new JMenu("Analyze");
		menuBar.add(menu);
		menu.add(startPlotCommand);
		menu.add(scaleResultsCommand);
		menu.add(exportVTKAction);

		// Window
		menu = new JMenu("Window");
		menuBar.add(menu);
		if(JNApplication.isMacLNF())
		{	menu.add(JNApplication.main.getOpenHelpAction());
		}
		menu.add(showPartnerCommand);
		menu.addSeparator();
		setWindowMenu(menu);

		// add the menu bar
		setJMenuBar(menuBar);
	}
	
	// read data, throw Exception error
	public void loadTextFromFile(String fileText) throws Exception
	{	display.readFile(getFile(),fileText);
		controls.fileHasLoaded();
		setTitle(resDoc.name);
	}
	
	// reread the file when it has changed (i.e., reuse window for new calculations)
	public void loadNewTextFromFile()
	{
		// read the file's text (UTF-8 encoding)
		String fileText;
		try
		{	FileInputStream theFile = new FileInputStream(getFile());
			BufferedReader read = new BufferedReader(new InputStreamReader(theFile, "UTF-8"));
		
			StringBuffer buffer = new StringBuffer();
			String str = read.readLine();
			if(str!=null)
			{	buffer.append(str);
				while((str=read.readLine()) != null)
				{	buffer.append((char)0x0A);
					buffer.append(str);
				}
			}
			
			read.close();
			theFile.close();
			
			// convert to a string
			fileText = buffer.toString();
			loadTextFromFile(fileText);
		}
		catch (Exception re)
		{	JNUtilities.showMessage(null,"Error rereading calculation results '"+getFile().getName()+": " + re);
			return;
		}
	}
	
	//----------------------------------------------------------------------------
	// handle application commands
	//----------------------------------------------------------------------------
	
	// only find in this read only document
	public void doFindReplaceAction(int frAction)
	{	display.textPane.doFindReplaceAction(frAction);
	}
	
	// start new plot using current settings
	public void startNewPlot(int newType)
	{	Class<?> currentClass=null;
		switch(newType)
		{	case LoadArchive.PARTICLE_PLOT:
				currentClass=MPMParticlePlotWindow.class;
			case LoadArchive.MESH_PLOT:
				if(newType!=LoadArchive.PARTICLE_PLOT)
				{	currentClass = resDoc.isMPMAnalysis() ?
						MeshPlotWindow.class : FEAMeshPlotWindow.class ;
				}
				MoviePlotWindow movieFrame=getMovieFrame();
				if(!currentMovieIsOK(currentClass,movieFrame))
				{	// current one not available - create a new window
					if(newType==LoadArchive.PARTICLE_PLOT)
						movieFrame=(MoviePlotWindow)(new MPMParticlePlotWindow(resDoc,this));
					else if(resDoc.isMPMAnalysis())
						movieFrame=(MoviePlotWindow)(new MeshPlotWindow(resDoc,this));
					else
						movieFrame=(MoviePlotWindow)(new FEAMeshPlotWindow(resDoc,this));
					
					// manually set to current setting in main window controls
					controls.plotOpened();
					int newComponent=controls.getPlotComponent(-1,true,null);
					movieFrame.setTitle(resDoc.name+" ("+PlotQuantity.plotName(newComponent,resDoc)+")");
					movieFrame.setVisible(true);
					movieFrame.toFront();
				}
				else
				{	// reuse the current movie plot window
					int newComponent=controls.getPlotComponent(-1,true,null);
					movieFrame.setTitle(resDoc.name+" ("+PlotQuantity.plotName(newComponent,resDoc)+")");
					movieFrame.setVisible(true);
					movieFrame.toFront();
					movieFrame.beginNewIndexNewComponent(controls.getArchiveIndex(),newComponent);
					if(newType==controls.getPlotType())
						movieFrame.movieControls.syncPlotQuantityMenus();
				}
				JNNotificationCenter.getInstance().postNotification("PlotUnitsChanged",this,null);
				break;
			
			case LoadArchive.TIME_PLOT:
				TimePlotWindow timeFrame=getTimeFrame();
				if(timeFrame==null)
				{	timeFrame=new TimePlotWindow(this);
					controls.plotOpened();
				}
				
				try
				{	timeFrame.addPlot(controls,null);
				}
				catch(Exception tpe)
				{	JNApplication.appBeep();
					JOptionPane.showMessageDialog(timeFrame,tpe.getLocalizedMessage());
				}
				
				break;
			
			case LoadArchive.MESH2D_PLOT:
				XYPlotWindow xyPlotFrame=getXYPlotFrame();
				if(xyPlotFrame==null)
				{	xyPlotFrame=new XYPlotWindow(this);
					controls.plotOpened();
				}
				
				try
				{	xyPlotFrame.addPlot(controls,null);
				}
				catch(Exception tpe)
				{	JNApplication.appBeep();
					JOptionPane.showMessageDialog(xyPlotFrame,tpe.getLocalizedMessage());
				}
				
				break;
			
			default:
				break;
		}
	}
	
	// start mesh 2D plot in script control
	public ISListType scriptTimeplot(ISDictType settings) throws Exception
	{	// look for global results
		String menutext = (String)settings.gcis_objectForKey("menutext");
		if(menutext==null)
			throw new Exception("timeplot settings does not specity what to plot in 'menutext'");
		
		// trap global data plot
		if(menutext.equals("Global Results"))
		{	ISListType plotResults = new ISListType(null,false);
			resDoc.collectTimePlotData(settings,plotResults);
			if(plotResults.count()<2)
				throw new Exception("The timeplot command found no plot data");
			return plotResults;
		}

		TimePlotWindow timeFrame=getTimeFrame();
		if(timeFrame==null)
		{	timeFrame=new TimePlotWindow(this);
			controls.plotOpened();
		}
		int startPlots = timeFrame.getNumberOfPlots();
		timeFrame.addPlot(controls,settings);
		
		// wait for run to be done
		while(true)
		{	Thread.sleep(250);
			if(timeFrame==null) break;
			if(!timeFrame.isPlotting())
				break;
		}
		
		if(timeFrame==null)
			throw new Exception("Scripted plot failed");
		int endPlots = timeFrame.getNumberOfPlots();
		if(startPlots==endPlots)
			throw new Exception("Scripted plot failed");
		
		// get plot results
		ISListType plotResults = timeFrame.getLastPlot();
		
		// close unless "close" says "no"
		String closePlot = (String)settings.gcis_objectForKey("close");
		if(closePlot!=null)
		{	if(!closePlot.equals("no"))
				timeFrame.dispose();
		}
		else
			timeFrame.dispose();
		
		// return the results
		return plotResults;
	}

	// start mesh 2D plot in script control
	public ISListType scriptXYplot(ISDictType settings) throws Exception
	{	XYPlotWindow xyPlotFrame=getXYPlotFrame();
		if(xyPlotFrame==null)
		{	xyPlotFrame=new XYPlotWindow(this);
			controls.plotOpened();
		}
		int startPlots = xyPlotFrame.getNumberOfPlots();
		xyPlotFrame.addPlot(controls,settings);
		
		// wait for run to be done
		while(true)
		{	Thread.sleep(250);
			if(xyPlotFrame==null) break;
			if(!xyPlotFrame.isPlotting())
				break;
		}
		
		if(xyPlotFrame==null)
			throw new Exception("Scripted plot failed");
		int endPlots = xyPlotFrame.getNumberOfPlots();
		if(startPlots==endPlots)
			throw new Exception("Scripted plot failed");
		
		// get plot results
		ISListType plotResults = xyPlotFrame.getLastPlot();
		
		// close unless "close" says "no"
		String closePlot = (String)settings.gcis_objectForKey("close");
		if(closePlot!=null)
		{	if(!closePlot.equals("no"))
				xyPlotFrame.dispose();
		}
		else
			xyPlotFrame.dispose();
		
		// return the results
		return plotResults;
	}
	
	// see if current movie is same class and requested movie
	private boolean currentMovieIsOK(Class<?> theClass,MoviePlotWindow movieFrame)
	{	if(movieFrame==null) return false;
		if(theClass.equals(movieFrame.getClass())) return true;
		movieFrame.dispose();
		return false;
	}
	
	// draw check on check mesh run
	public void checkMeshNow()
	{	controls.setCheckMeshItem();
		if(resDoc.isMPMAnalysis())
			startNewPlot(LoadArchive.PARTICLE_PLOT);
		else
			startNewPlot(LoadArchive.MESH_PLOT);
	}
	
	// scale mesh dialog
	public void scaleResults()
	{	// run the dialog box
		ScaleResultsDialog sr=new ScaleResultsDialog(this);
		if(sr.getClickedButton()==JNDialog.CANCEL_BUTTON) return;
		
		boolean changeUnits=false;
		if(resDoc.units.setLengthScaleIndex(sr.getLengthIndex()))
		{	changeUnits=true;
		}
		if(resDoc.units.setTimeScaleIndex(sr.getTimeIndex()))
		{	changeUnits=true;
		}
		
		// change of if new units were selected
		if(changeUnits)
		{	// close plot windows
			MoviePlotWindow mp=getMovieFrame();
			if(mp!=null) mp.windowClosing(null);
			TimePlotWindow tp=getTimeFrame();
			if(tp!=null) tp.windowClosing(null);
			XYPlotWindow xyp=getXYPlotFrame();
			if(xyp!=null) xyp.windowClosing(null);
			
			// reread the data
			try
			{	resDoc.clear(false);
				resDoc.DecodeFileSections(getFile());
				controls.updateTimeDisplay();
			}
			catch (Exception e)
			{	String emsg = e.getMessage();
				if(emsg == null)
					emsg = "Scanner error probably due to corrupted or misformatted data.";
				JOptionPane.showMessageDialog(this,"The analysis failed to read for rescaling:\n   "
												+ emsg);
				windowClosing(null);
			}
		}
	}

	//----------------------------------------------------------------------------
	// Accessors
	//----------------------------------------------------------------------------
	
	public MoviePlotWindow getMovieFrame()
	{	return (MoviePlotWindow)getFirstChildOfType("movieFrame");
	}
	
	public TimePlotWindow getTimeFrame()
	{	return (TimePlotWindow)getFirstChildOfType("timeFrame");
	}
	
	public XYPlotWindow getXYPlotFrame()
	{	return (XYPlotWindow)getFirstChildOfType("xyPlotFrame");
	}
	
	// Scripting attributes for internal scripts for results document
	public Object gcis_getObjectAttribute(String attr,CmdViewer server)
	{
		if(attr.contentEquals("archives"))
		{	ISListType archs = new ISListType(null);
			for(int i=0;i<resDoc.mpmArchives.size();i++)
				archs.gcis_addObject(resDoc.mpmArchives.get(i));
			return archs;
		}
		
		return null;
	}

	// Scripting attributes for internal scripts for GEDCOMObject
	public String gcis_getAttribute(String [] atoms,int i,CmdViewer server)
	{
	    String attr = server.grabAtom(atoms,i);
	    
	    if(attr.equals("energy"))
			return resDoc.getEnergy();
		
		else if(attr.equals("section"))
		{	i++;
			if(i >= atoms.length)
				return "ERROR: The section name is missing";
			attr = server.grabAtom(atoms,i);
			
			// check for quoted text
			int clen = attr.length();
			if(clen > 1)
			{	if(attr.charAt(0) == '"' && attr.charAt(clen - 1) == '"')
					attr =  attr.substring(1, clen - 1);
			}

			return resDoc.section(attr);
		}
	    
		else if(attr.equals("name"))
		{	File theFile = getFile();
			if(theFile!=null)
				return theFile.getName();
			return "Untitled";
		}

		else if(attr.equals("path"))
		{	File theFile = getFile();
			if(theFile!=null)
				return theFile.getPath();
			return "Untitled";
		}

		else if(attr.equals("folder"))
		{	File theFile = getFile();
			if(theFile!=null)
				return theFile.getParent();
			return "";
		}
	    
		else if(attr.equals("class"))
			return "ResultsDocument";

		return null;
	}
	
	// set command window when opened by simulation here
	public void setCommandsWindow(CmdViewer theCmds) { commandsWindow=theCmds; }

	//----------------------------------------------------------------------------
	// Window events
	//----------------------------------------------------------------------------
	
	public void windowOpened(WindowEvent e)
	{	double splitLoc=NFMVPrefs.prefs.getDouble(NFMVPrefs.ResultsTopSplitKey,
													NFMVPrefs.ResultsTopSplitDef);
		if(splitLoc<10.)
		{	Dimension winSize = JNPreferences.getFrameSize(this);
			splitLoc = (double)winSize.width/2.;
		}
	    display.setDividerLocation((int)splitLoc);
	    
		splitLoc=NFMVPrefs.prefs.getDouble(NFMVPrefs.ResultsSplitKey,
													NFMVPrefs.ResultsSplitDef);
		if(splitLoc<10.)
		{	Dimension winSize = JNPreferences.getFrameSize(this);
			splitLoc = (double)winSize.height/2.;
		}
		full.setDividerLocation((int)splitLoc);
		
		super.windowOpened(e);
		
		if(startWithMesh)
		{	controls.setCheckMeshItem();
			if(resDoc.isMPMAnalysis())
				startNewPlot(LoadArchive.PARTICLE_PLOT);
			else
				startNewPlot(LoadArchive.MESH_PLOT);
		}
	}
			
	public void windowClosing(WindowEvent e)
	{
		double loc=(double)display.getDividerLocation();
		NFMVPrefs.prefs.putDouble(NFMVPrefs.ResultsTopSplitKey,loc);
	    
	    loc=(double)full.getDividerLocation();
		NFMVPrefs.prefs.putDouble(NFMVPrefs.ResultsSplitKey,loc);

		super.windowClosing(e);
		
		// inform the commandsWindow
		if(commandsWindow!=null) commandsWindow.setLinkedResults(null);
	}
	
	public void childWindowDidClose(JNChildWindow child)
	{	super.childWindowDidClose(child);
		//might need to change button name
		controls.plotOpened();
	}

	// open dialog to extrack VTK files
	public void startExtractVTK()
	{	// dialog to extract VTK
		new ExtractVTK(this);
	}

	//----------------------------------------------------------------------------
	// Actions as inner classes
	//----------------------------------------------------------------------------
	
	// action to shaw partner menu command
	protected class ShowPartnerAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public ShowPartnerAction()
		{	super("Show Commands");
		}
 
		public void actionPerformed(ActionEvent e)
		{	if(commandsWindow!=null)
				commandsWindow.toFront();
			else
				JNApplication.appBeep();
		}
	}

	// action to shaw partner menu command
	protected class ScaleResultsAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public ScaleResultsAction()
		{	super("Scale Results...");
		}
 
		public void actionPerformed(ActionEvent e) { scaleResults(); }
	}
	
	// action to show partner menu command
	protected class StartPlotAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public StartPlotAction()
		{	super("Start Plot");
		}
 
		public void actionPerformed(ActionEvent e) { startNewPlot(controls.getPlotType()); }
	}

	// action to start a plot
	public StartPlotAction getStartPlotAction() { return startPlotCommand; }

		// action to export VTK
	protected class ExportVTKAction extends JNAction
	{	private static final long serialVersionUID = 1L;

		public ExportVTKAction()
		{	super("Export Particle VTKs...");
		}
 
		public void actionPerformed(ActionEvent e) { startExtractVTK(); }
	}
}
