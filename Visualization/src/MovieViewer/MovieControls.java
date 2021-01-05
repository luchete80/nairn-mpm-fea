/*******************************************************************
	MovieControls.java
	NairnFEAMPMViz

	Created by John Nairn on 3/1/07.
	Copyright 2007 RSAC Software. All rights reserved.
*******************************************************************/

import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import java.net.*;
import geditcom.JNFramework.*;

public class MovieControls extends JPanel
{
	static final long serialVersionUID=11L;
	
	//----------------------------------------------------------------------------
	// variables and constants
	//----------------------------------------------------------------------------
	
	public static final int HEIGHT=44;

	private JButton rewindMov=new JButton();
	private JButton playMov=new JButton();
	public TimeSelector selectTime=null;
	public JComboBox<PlotMenuItem> pquant=new JComboBox<PlotMenuItem>();
	public JComboBox<String> pcmpnt=new JComboBox<String>();
	public boolean disableStartPlot=false;
	
	// particle size
	private JLabel sizeSelected=new JLabel("PS: 100",JLabel.LEFT);
	public JSlider mpmParticleSize=new JSlider(JSlider.HORIZONTAL,0,200,5);
	int particleSize=100;
	
	private JLabel elongSelected=new JLabel("Max F: 20",JLabel.LEFT);
	public JSlider maxElongSlider=new JSlider(JSlider.HORIZONTAL,0,20,1);

	// axes
	private String xchar="x";
	private String ychar="y";
	private String zchar="z";
	private String totalchar="magnitude";
	
	// hold variable
	private ResultsDocument resDoc;

	//----------------------------------------------------------------------------
	// initialize
	//----------------------------------------------------------------------------
	
	MovieControls(int width,ResultsDocument gResDoc,MoviePlotWindow movieCtrl,DocViewer gDocView)
	{   super();
		setLayout(null);
		setPreferredSize(new Dimension(width,HEIGHT));
		setBackground(Color.lightGray);
		resDoc=gResDoc;
		
		int hpos=8;
		if(gResDoc.isMPMAnalysis())
		{	// rewind
			setButtonIcon("Rewind.png",rewindMov);
			rewindMov.setSize(rewindMov.getPreferredSize());
			int centerLoc=(HEIGHT-rewindMov.getHeight())>>1;
			rewindMov.setLocation(hpos,centerLoc);
			rewindMov.setActionCommand("Rewind");
			rewindMov.addActionListener(movieCtrl);
			rewindMov.setFocusPainted(false);
			add(rewindMov);
			
			// play button
			setButtonIcon("Play.png",playMov);
			playMov.setSize(playMov.getPreferredSize());
			centerLoc=(HEIGHT-playMov.getHeight())>>1;
			hpos+=rewindMov.getWidth()+6;
			playMov.setLocation(hpos,centerLoc);
			playMov.setActionCommand("Play");
			playMov.addActionListener(movieCtrl);
			playMov.setFocusPainted(false);
			add(playMov);
			
			LookAndFeel laf=UIManager.getLookAndFeel();
			if(laf!=null)
			{	if(laf.isNativeLookAndFeel())
				{	playMov.setBackground(Color.lightGray);
					rewindMov.setBackground(Color.lightGray);
				}
			}

			// time selector
			selectTime=new TimeSelector(movieCtrl,gResDoc.docCtrl);
			selectTime.setSize(new Dimension(140,HEIGHT));
			hpos+=playMov.getWidth()+12;
			selectTime.setLocation(hpos,0);
			add(selectTime);
			
			hpos+=selectTime.getWidth()+20;
		}
		
		if(resDoc.isAxisymmetric())
		{	xchar="r";
			ychar="z";
			zchar="t";
		}
		
		// data panel
		MeshPlotData meshData=new MeshPlotData(movieCtrl.getPlotView(),gResDoc);
		meshData.setSize(new Dimension(225,HEIGHT));
		meshData.setLocation(hpos,0);
		meshData.setBackground(Color.lightGray);
		add(meshData);
		
		hpos+=meshData.getWidth()+3;
		
		// quantity menu
		JComboBox<PlotMenuItem> quant=gResDoc.docCtrl.controls.getQuantityMenu();
		int n=quant.getItemCount();
		int i;
		for(i=0;i<n;i++)
		{	pquant.addItem(quant.getItemAt(i));
		}
		pquant.setSize(pquant.getPreferredSize());
		pquant.setLocation(hpos,(HEIGHT-pquant.getHeight())/2);
		add(pquant);
		hpos+=pquant.getWidth()+3;
		
		// component menu
		JComboBox<String> cmpnt=gResDoc.docCtrl.controls.getComponentMenu();
		n=cmpnt.getItemCount();
		for(i=0;i<n;i++)
		{	pcmpnt.addItem(cmpnt.getItemAt(i));
		}
		pcmpnt.setSize(new Dimension(100,pquant.getHeight()));
		pcmpnt.setLocation(hpos,(HEIGHT-pquant.getHeight())/2);
		add(pcmpnt);
		hpos+=pcmpnt.getWidth()+3;
		
		// select quantity
		pquant.setSelectedIndex(quant.getSelectedIndex());
		setComponentMenu();
		pcmpnt.setSelectedIndex(cmpnt.getSelectedIndex());
		
		// when quantity changes, update component menu, update parent controls, redraw plot
		pquant.addItemListener(new ItemListener()
		{	public void itemStateChanged(ItemEvent e)
			{	setComponentMenu();
				if(e.getStateChange()==ItemEvent.SELECTED)
					JNNotificationCenter.getInstance().postNotification("PlotQuantityChanged",resDoc.docCtrl,pquant);
			}
		});
		
		pcmpnt.addItemListener(new ItemListener()
		{	public void itemStateChanged(ItemEvent e)
			{	if(e.getStateChange()==ItemEvent.SELECTED)
					JNNotificationCenter.getInstance().postNotification("PlotQuantityChanged",resDoc.docCtrl,pcmpnt);
			}
		});
		
		// particle size
		if(gResDoc.isMPMAnalysis())
		{	mpmParticleSize.setSize(new Dimension(100,15));
			mpmParticleSize.setLocation(hpos,6);
			mpmParticleSize.setBackground(Color.lightGray);
			mpmParticleSize.setFocusable(false);
			add(mpmParticleSize);
		
			sizeSelected.setFont(new Font("sanserif",Font.PLAIN,10));
			sizeSelected.setSize(sizeSelected.getPreferredSize());
			sizeSelected.setLocation(hpos+(mpmParticleSize.getWidth()-sizeSelected.getWidth())/2,
						HEIGHT-sizeSelected.getHeight()-6);
			add(sizeSelected);
			hpos+=mpmParticleSize.getWidth()+3;
		
			mpmParticleSize.addChangeListener(new ChangeListener()
			{   public void stateChanged(ChangeEvent e)
				{	particleSize=mpmParticleSize.getValue();
					sizeSelected.setText("PS: "+particleSize);
					if(!mpmParticleSize.getValueIsAdjusting())
						JNNotificationCenter.getInstance().postNotification("ParticleSizeChanged",resDoc.docCtrl,mpmParticleSize);
				}
			});
			particleSize=(int)gResDoc.docCtrl.controls.getParticleSize();
			mpmParticleSize.setValue(particleSize);
			mpmParticleSize.setToolTipText("Scale particle size as percent of cell size (default is 100%)");
			
			maxElongSlider.setSize(new Dimension(80,15));
			maxElongSlider.setLocation(hpos,6);
			maxElongSlider.setBackground(Color.lightGray);
			maxElongSlider.setFocusable(false);
			add(maxElongSlider);
			
			elongSelected.setFont(new Font("sanserif",Font.PLAIN,10));
			elongSelected.setSize(elongSelected.getPreferredSize());
			elongSelected.setLocation(hpos+(maxElongSlider.getWidth()-elongSelected.getWidth())/2,
						HEIGHT-elongSelected.getHeight()-6);
			add(elongSelected);
			hpos+=maxElongSlider.getWidth()+3;
			
			maxElongSlider.addChangeListener(new ChangeListener()
			{   public void stateChanged(ChangeEvent e)
				{	int newMaxElong=maxElongSlider.getValue();
					elongSelected.setText("Max F: "+newMaxElong);
					NFMVPrefs.prefs.putDouble(NFMVPrefs.maxElongKey,(double)newMaxElong);
					if(!maxElongSlider.getValueIsAdjusting())
						JNNotificationCenter.getInstance().postNotification("MaxElongChanged",resDoc.docCtrl,mpmParticleSize);
				}
			});
			double maxElong = NFMVPrefs.prefs.getDouble(NFMVPrefs.maxElongKey,
					NFMVPrefs.maxElongDef);
			maxElongSlider.setValue((int)(maxElong+0.5));
			elongSelected.setText("Max F: "+(int)(maxElong+0.5));
			maxElongSlider.setToolTipText("Set maximum elongation plot (=0 to for no limit)");
			
		}

	}
	
	//----------------------------------------------------------------------------
	// accessors
	//----------------------------------------------------------------------------
	
	// help to set an image of a button
	private void setButtonIcon(String icon,JButton theBtn)
	{	URL btnImage=NairnFEAMPMViz.class.getResource("Resources/"+icon);
		theBtn.setIcon(new ImageIcon(btnImage));
	}
	
	// reset selected  index in mesh view when about to replot and true or false if changed
	public boolean setArchiveIndex(int newIndex) { return selectTime!=null ? selectTime.setArchiveIndex(newIndex) : false; }
	
	public int getArchiveIndex() { return selectTime!=null ? selectTime.getArchiveIndex() : 0; }
	public boolean incrementArchiveIndex() { return selectTime!=null ? selectTime.incrementArchiveIndex() : false; }
	public boolean decrementArchiveIndex() { return selectTime!=null ? selectTime.decrementArchiveIndex() : false; }
	
	// called when movie started or stopped
	public void setPlaying(boolean playing)
	{	if(playing)
			setButtonIcon("Pause.png",playMov);
		else
			setButtonIcon("Play.png",playMov);
	}
	
	// called when state changed in quantity menu
	public void setComponentMenu()
	{
		PlotMenuItem pm=(PlotMenuItem)pquant.getSelectedItem();
		if(pm==null) return;
		
		int numItems = pcmpnt.getItemCount();
		switch(pm.getTag())
		{   case PlotQuantity.MPMSIGMAX:
			case PlotQuantity.MPMEPSX:
			case PlotQuantity.MPMPLEPSX:
			case PlotQuantity.MPMEPSTOTX:
			case PlotQuantity.MESHSIGMAX:
			case PlotQuantity.MESHSTRAINX:
			case PlotQuantity.MESHELEMSIGMAX:
				if(numItems!=4 || !pcmpnt.getItemAt(0).equals(xchar+xchar))
				{	pcmpnt.removeAllItems();
					pcmpnt.addItem(xchar+xchar);
					pcmpnt.addItem(ychar+ychar);
					pcmpnt.addItem(xchar+ychar);
					pcmpnt.addItem(zchar+zchar);
				}
				pcmpnt.setEnabled(true);
				break;
			
			
			case PlotQuantity.MPMVELX:
			case PlotQuantity.MPMDISPX:
				if(numItems!=3 || !pcmpnt.getItemAt(0).equals(xchar))
				{	pcmpnt.removeAllItems();
					pcmpnt.addItem(xchar);
					pcmpnt.addItem(ychar);
					pcmpnt.addItem(totalchar);
				}
				pcmpnt.setEnabled(true);
				break;
				
			case PlotQuantity.MESHDISPX:
			case PlotQuantity.MESHFORCEX:
				if(numItems!=2 || !pcmpnt.getItemAt(0).equals(xchar))
				{	pcmpnt.removeAllItems();
					pcmpnt.addItem(xchar);
					pcmpnt.addItem(ychar);
				}
				pcmpnt.setEnabled(true);
				break;
			
			case PlotQuantity.MPMSPINVELOCITYX:
			case PlotQuantity.MPMSPINMOMENTUMX:
				if(resDoc.is3D())
				{	if(numItems!=3 || !pcmpnt.getItemAt(0).equals(xchar))
					{	pcmpnt.removeAllItems();
						pcmpnt.addItem(xchar);
						pcmpnt.addItem(ychar);
						pcmpnt.addItem(zchar);
					}
				}
				else
				{	if(numItems!=1 || !pcmpnt.getItemAt(0).equals(zchar))
					{	pcmpnt.removeAllItems();
						pcmpnt.addItem(zchar);
					}
				}
				pcmpnt.setEnabled(true);
				break;
				
			case PlotQuantity.MPMDCDX:
				if(numItems!=2 || !pcmpnt.getItemAt(0).equals("dc/d"+xchar))
				{	pcmpnt.removeAllItems();
					if(resDoc.hasPorePressure)
					{	pcmpnt.addItem("dp/d"+xchar);
						pcmpnt.addItem("dp/d"+ychar);
					}
					else
					{	pcmpnt.addItem("dc/d"+xchar);
						pcmpnt.addItem("dc/d"+ychar);
					}
				}
				pcmpnt.setEnabled(true);
				break;
				
			case PlotQuantity.INTERFACETRACTION_N:
				if(numItems!=2 || !pcmpnt.getItemAt(0).equals("normal"))
				{	pcmpnt.removeAllItems();
					pcmpnt.addItem("normal");
					pcmpnt.addItem("tangential");
				}
				pcmpnt.setEnabled(true);
				break;
				
			default:
				pcmpnt.setEnabled(false);
				break;
		}
	}
	
	// get plot component from current selection
	public int getPlotComponent()
	{	PlotMenuItem pm=(PlotMenuItem)pquant.getSelectedItem();
		if(pm==null) return -1;
		int plotComponent=pm.getTag();
		
		// adjust component menus - add component selected in component menu
		int extra;
		switch(plotComponent)
		{   
			case PlotQuantity.MPMVELX:
			{
				extra=pcmpnt.getSelectedIndex();
				int nextVal = resDoc.is3D()?3:2;
				if(extra == nextVal) extra=PlotQuantity.SHIFT4MAGNITUDE;
				break;
			}
			case PlotQuantity.MPMDISPX:
			{
				extra=pcmpnt.getSelectedIndex();
				int nextVal = resDoc.is3D()?3:2;
				if(extra == nextVal) extra=PlotQuantity.SHIFT4MAGNITUDE;
				break;
			}
				
			case PlotQuantity.MPMSIGMAX:
			case PlotQuantity.MPMEPSX:
			case PlotQuantity.MPMPLEPSX:
			case PlotQuantity.MPMEPSTOTX:
			case PlotQuantity.MPMDCDX:
			case PlotQuantity.MESHSIGMAX:
			case PlotQuantity.MESHDISPX:
			case PlotQuantity.MESHSTRAINX:
			case PlotQuantity.MESHDVDX:
			case PlotQuantity.MESHELEMSIGMAX:
			case PlotQuantity.MESHFORCEX:
			case PlotQuantity.INTERFACETRACTION_N:
				extra=pcmpnt.getSelectedIndex();
				if(extra>=0) plotComponent+=extra;
				break;
			case PlotQuantity.MPMSPINVELOCITYX:
			case PlotQuantity.MPMSPINMOMENTUMX:
				if(resDoc.is3D())
				{	extra=pcmpnt.getSelectedIndex();
					if(extra>=0) plotComponent+=extra;
				}
				else
					plotComponent+=2;	// only z in 2D
				break;
			default:
				break;
		}
		return plotComponent;
	}
	
	// synch quantity and component menu when replot
	public void syncPlotQuantityMenus()
	{	disableStartPlot=true;
		JComboBox<PlotMenuItem> plotQuant=resDoc.docCtrl.controls.getQuantityMenu();
		if(plotQuant.getSelectedIndex()!=pquant.getSelectedIndex())
			pquant.setSelectedIndex(plotQuant.getSelectedIndex());
		JComboBox<String> plotCmpnt=resDoc.docCtrl.controls.getComponentMenu();
		if(plotCmpnt.getSelectedIndex()!=pcmpnt.getSelectedIndex())
			pcmpnt.setSelectedIndex(plotCmpnt.getSelectedIndex());
		JSlider partSize=resDoc.docCtrl.controls.getParticleSizeSlider();
		if(partSize.getValue()!=mpmParticleSize.getValue())
			mpmParticleSize.setValue(partSize.getValue());
		disableStartPlot=false;
	}

}
