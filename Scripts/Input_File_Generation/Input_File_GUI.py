import re, wx, os, sys
from wx.lib.pubsub import pub

#********************************************************************************
#   Cassandra - An open source atomistic Monte Carlo software package
#   developed at the University of Notre Dame.
#   http://cassandra.nd.edu
#   Prof. Edward Maginn <ed@nd.edu>
#   Copyright (2013) University of Notre Dame du Lac
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#********************************************************************************
#
#   Created by Brian Keene on February 4, 2015
#
#********************************************************************************

# Below, each panel is described as a separate class, placed in appropriate parent panels,
# all of which are children of the main Frame.  Each class contains separate sizers and
# storage variables, with event handling for each type of widget used on that panel.

# The fundamental variables for the GUI are the ensemble selected (which determines the number
# of boxes used in the simulation) and the number of species in simulation.

# Panel One - Basic information
class p1a(wx.Panel):
  ## Declare the class variables
  runname_store = ""
  ensembletype_store = ""
  ensembleindex_store = 0
  nbrspecies_store = 1
  nbr_boxes = 1
  simulation_dirname = ""
  # user input stored for later use
  # note - for variables designed to have default values, or only one option,
  # set initial values for the class variables.  If changed by the user in the interface,
  # event handling will change the value for that instance.
  ## initialize the panel

  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

    ##### Create a sizer unique to this panel
    grid1a = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid1a)
    self.dirname = ""

    # Static text labels
    self.staticlbls = []
    alertstring = "Please provide the following information before continuing to other panels."
    labelstr = [alertstring, "Run Name: ","Simulation Directory:", "Ensemble:", "Number of Species:"]
    for i in range(5):
      self.staticlbls.append(wx.StaticText(self,label=labelstr[i]))
      if i>0:
        grid1a.Add(self.staticlbls[i],pos=(i+2,1))
    grid1a.Add(self.staticlbls[0],pos=(1,1),span=(1,4)) #requires span to maintain normal column widths

    # text control widget
    self.runnamequery=wx.TextCtrl(self,value="",name="run_name")
    grid1a.Add(self.runnamequery,pos=(3,2),flag=wx.EXPAND)
    self.Bind(wx.EVT_TEXT , self.textfunc1a , self.runnamequery)
    self.runnamequery.SetToolTipString("The simulation run name.  Do not use file or directory level restricted characters.  Do not use extensions.")
    # choice control widgets
    nbrspec_options = ["","1","2","3","4","5","6"]
    self.specieschoice = wx.Choice(self,-1,choices=nbrspec_options,name="specieschoice")
    self.specieschoice.SetToolTipString("Select the number of species in simulation.")
    grid1a.Add(self.specieschoice,pos=(6,2))
    self.Bind(wx.EVT_CHOICE , self.choicefunc1a , self.specieschoice)

    ensemble_options = ["","NVT_MC","NVT_MIN","NPT_MC","GCMC","GEMC","NVT_MC_Fragment","NVT_MC_Ring_Fragment"]
    self.ensemblechoice = wx.Choice(self,-1,choices=ensemble_options,name="ensemblechoice")
    self.ensemblechoice.SetToolTipString("Select the simulation ensemble.  For further information, consult user_guide.pdf found in Cassandra/Documentation/")
    grid1a.Add(self.ensemblechoice,pos=(5,2))
    self.Bind(wx.EVT_CHOICE , self.choicefunc1a , self.ensemblechoice)

    # simulation directory selection
    self.simdirbtn = wx.Button(self,label="Select Directory",name="simdir")
    self.simdirbtn.SetToolTipString("Select the directory where the input file is to be placed.")
    grid1a.Add(self.simdirbtn,pos=(4,2),flag=wx.EXPAND)
    self.Bind(wx.EVT_BUTTON,self.btnfunc1a,self.simdirbtn)
    self.simdirdisp = wx.TextCtrl(self,value="",size=(200,-1),name="simdirdisp",style=wx.TE_READONLY)
    self.simdirdisp.SetToolTipString("Displays the selected directory.")
    grid1a.Add(self.simdirdisp,pos=(4,3))

  def textfunc1a(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "run_name":
      p1a.runname_store = text_data

  def choicefunc1a(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    choice_data = event.GetString()
    chosen_index = event.GetSelection()
    if widgetname == "specieschoice":
      if chosen_index == 0:
        choice_data = 1
      pub.sendMessage("species_msg", message = choice_data)
      p1a.nbrspecies_store = choice_data
    elif widgetname == "ensemblechoice":
      p1a.ensembletype_store = choice_data
      p1a.ensembleindex_store = chosen_index
      pub.sendMessage("ensemble_msg", message = choice_data)
      if choice_data == "GEMC":
        p1a.nbr_boxes = 2
      else: p1a.nbr_boxes = 1

  def btnfunc1a(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    dlg = wx.DirDialog(self,"Select Simulation Directory",style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      selected_path = dlg.GetPath()
      dir_path = selected_path
      self.dirname = os.path.split(dir_path)
      pub.sendMessage("directory_msg",message=dir_path)
      self.simdirdisp.SetValue(("/%s/" %self.dirname[1]))
    self.Layout()

class p1b(wx.Panel):
  angs = u'\u212B'.encode('utf-8') #angstroms unit
  mixingruletype_store = "LB"
  pairenergy_store = "TRUE"
  min_cutoffstore = "1.0"
  temperature_store = ""
  pressure_store = "0" #if zero, will be rejected on save and keyword not used
  chempot_store = ["","","","","",""]
  box1_shape = "CUBIC"
  box2_shape = "CUBIC"
  box1_length = ""
  box2_length = ""
  seed1_store = ""
  seed2_store = ""
  kappa_ins = ""
  kappa_rot = "0"
  kappa_dih = ""
  rcut_cbmc = ["",""]
  
  ensemble_received = ""
  nbr_spec_received = "1"
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid1b = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid1b)

    self.mixing_options = ["LB","geometric"]
    MixingRuleList = ["Lorentz-Berthelot","Geometric"]
    PairStorageList = ["TRUE","FALSE"]
    BoxShapeList = ["CUBIC"]

    pub.subscribe(self.ensemblefunc,"ensemble_msg")
    pub.subscribe(self.speciesfunc,"species_msg")

    self.templbl = wx.StaticText(self,label="Temperature (K): ")
    grid1b.Add(self.templbl, pos=(1,1))
    self.mixrulelbl = wx.StaticText(self,label="Mixing Rule: ")
    grid1b.Add(self.mixrulelbl,pos=(2,1))
    self.mincutlbl = wx.StaticText(self,label=("Minimum Cutoff (%s): " %(p1b.angs)))
    grid1b.Add(self.mincutlbl,pos=(3,1))
    self.pairstorelbl = wx.StaticText(self,label="Pair Storage:")
    grid1b.Add(self.pairstorelbl,pos=(4,1))
    self.seed1lbl = wx.StaticText(self,label="Seed 1:")
    grid1b.Add(self.seed1lbl,pos=(5,1))
    self.seed2lbl = wx.StaticText(self,label="Seed 2:")
    grid1b.Add(self.seed2lbl,pos=(6,1))
    self.pressurelbl = wx.StaticText(self,label="Pressure (bar): ")
    grid1b.Add(self.pressurelbl,pos=(7,1))
    self.chempotlbl = wx.StaticText(self,label="Chemical Potential (kJ/mol)")
    grid1b.Add(self.chempotlbl,pos=(1,6),flag=wx.ALIGN_CENTER)

    self.speclbl = []
    for i in range(6):
      self.speclbl.append(wx.StaticText(self,label = ("Species %d:" %(i+1))))
      grid1b.Add(self.speclbl[i],pos=(i+2,5))
    self.boxinfolbl = wx.StaticText(self,label="Box Information")
    grid1b.Add(self.boxinfolbl,pos=(9,1))
    self.boxshapelbl = wx.StaticText(self,label="Box Shape")
    grid1b.Add(self.boxshapelbl,pos=(9,2))
    self.boxedgelbl = wx.StaticText(self,label=("Box Edge Length (%s)" %(p1b.angs)))
    grid1b.Add(self.boxedgelbl,pos=(9,3),span=(1,2))
    self.box1lbl = wx.StaticText(self,label="Box 1: ")
    grid1b.Add(self.box1lbl,pos=(10,1))
    self.box2lbl = wx.StaticText(self,label="Box 2: ")
    grid1b.Add(self.box2lbl,pos=(11,1))

    self.textwidgets1b = []
    for i in range(13):
      self.textwidgets1b.append((wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name = ("textwidget %d" %i))))
      self.Bind(wx.EVT_TEXT, self.textfunc1b, self.textwidgets1b[i])
    grid1b.Add(self.textwidgets1b[0],pos=(1,2)) #index 0 is temperature widget (position determines widget prompt)
    grid1b.Add(self.textwidgets1b[1],pos=(3,2)) #index 1 is min cutoff
    grid1b.Add(self.textwidgets1b[2],pos=(5,2)) #index 2 is seed 1
    grid1b.Add(self.textwidgets1b[3],pos=(6,2)) #index 3 is seed 2
    grid1b.Add(self.textwidgets1b[4],pos=(7,2)) #index 4 is pressure
    grid1b.Add(self.textwidgets1b[5],pos=(10,3)) #index 5 is box 1 edge length
    grid1b.Add(self.textwidgets1b[6],pos=(11,3)) #index 6 is box 2 edge length
    
    for i in range(6):
      grid1b.Add(self.textwidgets1b[i+7],pos=(i+2,6),flag=wx.EXPAND) #species 1-6 chemical potential, widgets 7-12
      self.textwidgets1b[i+7].SetToolTipString("Valid entries are floating point numerics")
    self.choicewidgets1b = []
    self.choicewidgets1b.append(wx.Choice(self,-1,choices=MixingRuleList,name="mixingchoice"))
    grid1b.Add(self.choicewidgets1b[0],pos=(2,2),span=(1,2))
    self.choicewidgets1b.append(wx.Choice(self,-1,choices=PairStorageList,name="pairstoragechoice"))
    grid1b.Add(self.choicewidgets1b[1],pos=(4,2))
    self.choicewidgets1b.append(wx.Choice(self,-1,choices=BoxShapeList,name="box1shapechoice"))
    grid1b.Add(self.choicewidgets1b[2],pos=(10,2))
    self.choicewidgets1b.append(wx.Choice(self,-1,choices=BoxShapeList,name="box2shapechoice"))
    grid1b.Add(self.choicewidgets1b[3],pos=(11,2))
    for item in self.choicewidgets1b:
      self.Bind(wx.EVT_CHOICE, self.choicefunc1b, item)

    self.cbmclbl = wx.StaticText(self,label="CBMC Parameters")
    grid1b.Add(self.cbmclbl,pos=(13,1))
    
    self.kappainslbl = wx.StaticText(self,label="Trial Insertions:")
    grid1b.Add(self.kappainslbl,pos=(14,1))
    self.kapparotlbl = wx.StaticText(self,label="Rotational Bias:")
    grid1b.Add(self.kapparotlbl,pos=(15,1))
    self.kappadihlbl = wx.StaticText(self,label="Trial Orientations:")
    grid1b.Add(self.kappadihlbl,pos=(16,1))
    self.rcut_cbmc1 = wx.StaticText(self,label=("Cutoff (%s) Box 1:" %(p1b.angs)))
    self.rcut_cbmc2 = wx.StaticText(self,label=("Cutoff (%s) Box 2:" %(p1b.angs)))
    grid1b.Add(self.rcut_cbmc1,pos=(17,1))
    grid1b.Add(self.rcut_cbmc2,pos=(18,1))
    self.moretexts = []
    for i in range(5):
##############
############## delete the next three lines in later versions to allow for kappa_rot to be nonzero
      if i == 1:
        self.moretexts.append(wx.TextCtrl(self,value="",name=("CBMC %d" %i),style=wx.TE_READONLY))
        self.moretexts[i].SetValue("0")
      else:
        self.moretexts.append(wx.TextCtrl(self,value="",name=("CBMC %d" %i)))
      grid1b.Add(self.moretexts[i],pos=(14+i,2))
      self.Bind(wx.EVT_TEXT, self.textfunc1b, self.moretexts[i])
    # some preliminary hiding
    self.rcut_cbmc2.Hide(),self.moretexts[4].Hide()
    self.box2lbl.Hide(),self.textwidgets1b[6].Hide(),self.choicewidgets1b[3].Hide()

  def textfunc1b(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "textwidget 0":
      p1b.temperature_store = text_data
    elif widgetname == "textwidget 1":
      p1b.min_cutoffstore = text_data
    elif widgetname == "textwidget 2":
      p1b.seed1_store = text_data
    elif widgetname == "textwidget 3":
      p1b.seed2_store = text_data
    elif widgetname == "textwidget 4":
      p1b.pressure_store = text_data
    elif widgetname == "textwidget 5":
      p1b.box1_length = text_data
    elif widgetname == "textwidget 6":
      p1b.box2_length = text_data
    elif widgetname == "textwidget 7":
      p1b.chempot_store[0] = text_data
    elif widgetname == "textwidget 8":
      p1b.chempot_store[1] = text_data
    elif widgetname == "textwidget 9":
      p1b.chempot_store[2] = text_data
    elif widgetname == "textwidget 10":
      p1b.chempot_store[3] = text_data
    elif widgetname == "textwidget 11":
      p1b.chempot_store[4] = text_data
    elif widgetname == "textwidget 12":
      p1b.chempot_store[5] = text_data
    elif widgetname == "CBMC 0":
      p1b.kappa_ins = text_data
    elif widgetname == "CBMC 1":
      p1b.kappa_rot = text_data
    elif widgetname == "CBMC 2":
      p1b.kappa_dih = text_data
    elif widgetname == "CBMC 3":
      p1b.rcut_cbmc[0] = text_data
    elif widgetname == "CBMC 4":
      p1b.rcut_cbmc[1] = text_data

  def choicefunc1b(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    choice_data = event.GetString()
    chosen_index = event.GetSelection() 
    if widgetname == "mixingchoice":
      p1b.mixingruletype_store = self.mixing_options[chosen_index]
    elif widgetname == "pairstoragechoice":
      p1b.pairenergy_store = choice_data
    elif widgetname == "box1shapechoice":
      p1b.box1_shape = choice_data
    elif widgetname == "box2shapechoice":
      p1b.box2_shape = choice_data

  def ensemblefunc(self,message,arg2=None):
    msg = message
    p1b.ensemble_received = msg
    if msg == "GEMC":
      self.rcut_cbmc2.Show(),self.moretexts[4].Show()
      self.box2lbl.Show(),self.textwidgets1b[6].Show(),self.choicewidgets1b[3].Show()
      if int(p1b.nbr_spec_received)>1:
        self.pressurelbl.Show(),self.textwidgets1b[4].Show()
      self.chempotlbl.Hide()
      for i in range(6):
        self.speclbl[i].Hide()
        self.textwidgets1b[i+7].Hide()
    elif msg == "GCMC":
      for i in range(int(p1b.nbr_spec_received)):
        self.speclbl[i].Show()
        self.textwidgets1b[i+7].Show()
      self.chempotlbl.Show()
      self.rcut_cbmc2.Hide(),self.moretexts[4].Hide()
      self.box2lbl.Hide(),self.textwidgets1b[6].Hide(),self.choicewidgets1b[3].Hide()
      self.pressurelbl.Hide(),self.textwidgets1b[4].Hide()
    elif msg == "NPT_MC":
      self.pressurelbl.Show(),self.textwidgets1b[4].Show()
      for i in range(6):
        self.speclbl[i].Hide()
        self.textwidgets1b[i+7].Hide()
      self.chempotlbl.Hide()
      self.rcut_cbmc2.Hide(),self.moretexts[4].Hide()
      self.box2lbl.Hide(),self.textwidgets1b[6].Hide(),self.choicewidgets1b[3].Hide()
    else:
      self.pressurelbl.Hide(),self.textwidgets1b[4].Hide()
      for i in range(6):
        self.speclbl[i].Hide()
        self.textwidgets1b[i+7].Hide()
      self.chempotlbl.Hide()
      self.rcut_cbmc2.Hide(),self.moretexts[4].Hide()
      self.box2lbl.Hide(),self.textwidgets1b[6].Hide(),self.choicewidgets1b[3].Hide()
    self.Layout()

  def speciesfunc(self,message,arg2=None):
    msg = int(message)
    p1b.nbr_spec_received = message
    if p1b.ensemble_received == "GEMC" and msg > 1:
      self.pressurelbl.Show(),self.textwidgets1b[4].Show()
    elif p1b.ensemble_received == "GEMC" and msg <2:
      self.pressurelbl.Hide(),self.textwidgets1b[4].Hide()
    if p1b.ensemble_received == "GCMC":
      for i in range(msg):
        self.speclbl[i].Show()
        self.textwidgets1b[i+7].Show()
      self.chempotlbl.Show()
      for i in range(6-msg):
        self.speclbl[i+msg].Hide()
        self.textwidgets1b[i+msg+7].Hide()
    else:
      for i in range(msg):
        self.speclbl[i].Hide()
        self.textwidgets1b[i+7].Hide()
      self.chempotlbl.Hide()

# notebook holding p1a and p1b
class PanelOne(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    nested_notebook1 = wx.Notebook(self)
    nested_page1 = p1a(nested_notebook1)
    nested_page2 = p1b(nested_notebook1)
  
    nested_notebook1.AddPage(nested_page1,"Page 1")
    nested_notebook1.AddPage(nested_page2,"Page 2")

    nested_sizer = wx.BoxSizer()
    nested_sizer.Add(nested_notebook1,1,wx.EXPAND)
    self.SetSizer(nested_sizer)

###################################################################################
###################################################################################

# Intermolecular
class p2a(wx.Panel):
  # define the class variables
  angs = u'\u212B'.encode('utf-8')
  box1_vdw = ["","","","",""]
  box2_vdw = ["","","","",""]
  box1_charge = ["","","",""]
  box2_charge = ["","","",""]
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid2_a = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid2_a)

    pub.subscribe(self.ensemblefunc,"ensemble_msg")

    self.vdwlbl = wx.StaticText(self,label="van der Waals Style")
    self.chargelbl = wx.StaticText(self,label="Charge Style")
    self.box1lbl1 = wx.StaticText(self,label="Box 1:")
    self.box2lbl1 = wx.StaticText(self,label="Box 2:")
    self.box1lbl2 = wx.StaticText(self,label="Box 1:")
    self.box2lbl2 = wx.StaticText(self,label="Box 2:")

    vdwchoice1 = ["","Lennard Jones 12-6","None"]
    vdwchoice2 = ["","cut","cut_tail","cut_switch","cut_shift"]
    chargechoice1 = ["","Coulombic","NONE"]
    chargechoice2 = ["","Ewald","cut"]

    self.box1vdw_1 = wx.Choice(self,-1,choices=vdwchoice1,name="box1vdw 1")
    self.box1vdw_2 = wx.Choice(self,-1,choices=vdwchoice2,name="box1vdw 2")
    self.box2vdw_1 = wx.Choice(self,-1,choices=vdwchoice1,name="box2vdw 1")
    self.box2vdw_2 = wx.Choice(self,-1,choices=vdwchoice2,name="box2vdw 2")
    
    self.box1charge_1 = wx.Choice(self,-1,choices=chargechoice1,name="box1charge 1")
    self.box1charge_2 = wx.Choice(self,-1,choices=chargechoice2,name="box1charge 2")
    self.box2charge_1 = wx.Choice(self,-1,choices=chargechoice1,name="box2charge 1")
    self.box2charge_2 = wx.Choice(self,-1,choices=chargechoice2,name="box2charge 2")
    
    self.Bind(wx.EVT_CHOICE, self.choicefunc2a,self.box1vdw_1)
    self.Bind(wx.EVT_CHOICE, self.choicefunc2a,self.box1vdw_2)
    self.Bind(wx.EVT_CHOICE, self.choicefunc2a,self.box2vdw_1)
    self.Bind(wx.EVT_CHOICE, self.choicefunc2a,self.box2vdw_2)

    self.Bind(wx.EVT_CHOICE, self.choicefunc2a, self.box1charge_1)
    self.Bind(wx.EVT_CHOICE, self.choicefunc2a, self.box1charge_2)
    self.Bind(wx.EVT_CHOICE, self.choicefunc2a, self.box2charge_1)
    self.Bind(wx.EVT_CHOICE, self.choicefunc2a, self.box2charge_2)

    self.vdwcol1lbl = wx.StaticText(self,label="Functional Form")
    self.vdwcol2lbl = wx.StaticText(self,label="Tail Correction")
    grid2_a.Add(self.vdwcol1lbl,pos=(2,2),flag=wx.ALIGN_CENTER)
    grid2_a.Add(self.vdwcol2lbl,pos=(2,3),flag=wx.ALIGN_CENTER)

    self.vdw2col1lbl = wx.StaticText(self,label="Functional Form")
    self.vdw2col2lbl = wx.StaticText(self,label="Tail Correction")
    grid2_a.Add(self.vdw2col1lbl,pos=(5,2),flag=wx.ALIGN_CENTER)
    grid2_a.Add(self.vdw2col2lbl,pos=(5,3),flag=wx.ALIGN_CENTER)
    
    self.vdw1col3lbl = wx.StaticText(self,label=("Cutoff (%s)" %(p2a.angs)))
    self.vdw1col4lbl = wx.StaticText(self,label=("Spline Off (%s)" %(p2a.angs)))
    grid2_a.Add(self.vdw1col3lbl,pos=(2,4),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.vdw1col4lbl,pos=(2,5),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    
    self.vdw2col3lbl = wx.StaticText(self,label=("Cutoff (%s)" %(p2a.angs)))
    self.vdw2col4lbl = wx.StaticText(self,label=("Spline Off (%s)" %(p2a.angs)))
    grid2_a.Add(self.vdw2col3lbl,pos=(5,4),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.vdw2col4lbl,pos=(5,5),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    
    self.vdw1col5lbl = wx.StaticText(self,label="Logical (Optional)")
    grid2_a.Add(self.vdw1col5lbl,pos=(2,6),flag=wx.ALIGN_CENTER)
    self.vdw2col5lbl = wx.StaticText(self,label="Logical (Optional)")
    grid2_a.Add(self.vdw2col5lbl,pos=(5,6),flag=wx.ALIGN_CENTER)
    self.vdw1col5lbl.Hide(),self.vdw2col5lbl.Hide()
    
    options_col5 = ["","TRUE"]
    self.box1vdw_5 = wx.Choice(self,-1,choices=options_col5,name="box1vdw 5")
    self.box2vdw_5 = wx.Choice(self,-1,choices=options_col5,name="box2vdw 5")
    grid2_a.Add(self.box1vdw_5,pos=(3,6))
    grid2_a.Add(self.box2vdw_5,pos=(6,6))
    self.box1vdw_5.Hide(),self.box2vdw_5.Hide()
    self.Bind(wx.EVT_CHOICE,self.choicefunc2a,self.box1vdw_5)
    self.Bind(wx.EVT_CHOICE,self.choicefunc2a,self.box2vdw_5)
    
    # create the textcontrols - will need 8
    self.textwidgets = []
    for i in range(8):
      self.textwidgets.append((wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name=("textwidget %d" %i))))
      self.Bind(wx.EVT_TEXT, self.textfunc2_a, self.textwidgets[i])
      self.textwidgets[i].Hide() #hide until called for
    
    # textwidgets 0,1 = box1_vdw; 2,3 = box2_vdw; 4,5 = box1_charge; 6,7 = box2_charge

    # placement of widgets
    grid2_a.Add(self.vdwlbl,pos=(1,1),span=(1,2), flag = wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box1lbl1,pos=(3,1),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box2lbl1,pos=(6,1),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.chargelbl,pos=(8,1),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box1lbl2,pos=(10,1),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box2lbl2,pos=(13,1),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

    grid2_a.Add(self.box1vdw_1,pos=(3,2),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box1vdw_2,pos=(3,3),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box2vdw_1,pos=(6,2),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box2vdw_2,pos=(6,3),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
  
    grid2_a.Add(self.box1charge_1,pos=(10,2),flag=wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box1charge_2,pos=(10,3),flag=wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box2charge_1,pos=(13,2),flag=wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.box2charge_2,pos=(13,3),flag=wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN) #expand is for aesthetic purposes
    for i in range(8):
      if i < 2:
        grid2_a.Add(self.textwidgets[i],pos=(3,i+4),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
      elif i == 2 or i ==3:
        grid2_a.Add(self.textwidgets[i],pos=(6,i+2),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
      elif i == 4 or i == 5:
        grid2_a.Add(self.textwidgets[i],pos=(10,i),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
      elif i >5:
        grid2_a.Add(self.textwidgets[i],pos=(13,i-2),flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
      
    self.col2c1 = wx.StaticText(self,label="Functional Form")
    self.col2c2 = wx.StaticText(self,label="Functional Form")
    grid2_a.Add(self.col2c1,pos=(9,2),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.col2c2,pos=(12,2),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    
    self.col3c1 = wx.StaticText(self,label="Method")
    self.col3c2 = wx.StaticText(self,label="Method")
    grid2_a.Add(self.col3c1,pos=(9,3),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.col3c2,pos=(12,3),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    
    self.col4c1 = wx.StaticText(self,label=("Cutoff (%s)" %(p2a.angs)))
    self.col4c2 = wx.StaticText(self,label=("Cutoff (%s)" %(p2a.angs)))
    grid2_a.Add(self.col4c1,pos=(9,4),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.col4c2,pos=(12,4),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    
    self.col5c1 = wx.StaticText(self,label="Accuracy")
    self.col5c2 = wx.StaticText(self,label="Accuracy")
    grid2_a.Add(self.col5c1,pos=(9,5),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    grid2_a.Add(self.col5c2,pos=(12,5),flag=wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    
  # hide all widgets except initial choice functions and their labels
    self.emptyornone_select1 = [self.vdwcol2lbl,self.vdw1col3lbl,self.vdw1col4lbl,self.box1vdw_2,self.textwidgets[0],self.textwidgets[1],self.vdw1col5lbl,self.box1vdw_5]
    self.emptyornone_select2 = [self.vdw2col2lbl,self.vdw2col3lbl,self.vdw2col4lbl,self.vdw2col5lbl,self.box2vdw_2,self.textwidgets[2],self.textwidgets[3],self.box2vdw_5]
    ## hide all widgets except initial choice functions and their labels
    self.eorn_1c = [self.col3c1,self.col4c1,self.col5c1,self.textwidgets[4],self.textwidgets[5],self.box1charge_2]
    self.eorn_2c = [self.col3c2,self.col4c2,self.col5c2,self.textwidgets[6],self.textwidgets[7],self.box2charge_2]
    for item in self.emptyornone_select1:
      item.Hide()
    for item in self.emptyornone_select2:
      item.Hide()
    for item in self.eorn_1c:
      item.Hide()
    for item in self.eorn_2c:
      item.Hide()
      
  #define function to place values in to class variables and allow dynamic show/hide
  def textfunc2_a(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "textwidget 0":
      p2a.box1_vdw[2] = text_data
    elif widgetname == "textwidget 1":
      p2a.box1_vdw[3] = text_data
    elif widgetname == "textwidget 2":
      p2a.box2_vdw[2] = text_data
    elif widgetname == "textwidget 3":
      p2a.box2_vdw[3] = text_data
    elif widgetname == "textwidget 4":
      p2a.box1_charge[2] = text_data
    elif widgetname == "textwidget 5":
      p2a.box1_charge[3] = text_data
    elif widgetname == "textwidget 6":
      p2a.box2_charge[2] = text_data
    elif widgetname == "textwidget 7":
      p2a.box2_charge[3] = text_data
  
  def choicefunc2a(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    choice_data = event.GetString()
    chosen_index = event.GetSelection()
    if choice_data == "Lennard Jones 12-6":
      choice_data = "LJ"
    elif choice_data == "Coulombic":
      choice_data = "coul"
    for i in range(1,3):
      if widgetname == ("box1vdw %d" %i):
        p2a.box1_vdw[i-1] = str(choice_data)
      elif widgetname == ("box2vdw %d" %i):
        p2a.box2_vdw[i-1] = str(choice_data)
      elif widgetname == ("box1charge %d" %i):
        p2a.box1_charge[i-1] = str(choice_data)
      elif widgetname == ("box2charge %d" %i):
        p2a.box2_charge[i-1] = str(choice_data)
      elif widgetname == "box1vdw 5":
        p2a.box1_vdw[4] = choice_data
      elif widgetname == "box2vdw 5":
        p2a.box2_vdw[4] = choice_data

    #dynamic hiding, relabeling, and ensuring irrelevant class variables are empty strings when needed
    #box 1, empty string or "none" selected - clear all box 1 widgets, reset choices, clear class vars
    # box1vdw1 and chosen_index = (0 or 2) as one line statement had trouble working, separated seems to be OK
    if widgetname == "box1vdw 1" and chosen_index == 0:
      for item in self.emptyornone_select1:
        item.Hide()
      self.box1vdw_2.SetSelection(0),self.textwidgets[0].SetValue(""),self.textwidgets[1].SetValue(""),self.box1vdw_5.SetSelection(0)
      for i in range(1,5): p2a.box1_vdw[i] = ""
    elif widgetname == "box1vdw 1" and chosen_index == 2:
      for item in self.emptyornone_select1:
        item.Hide()
      self.box1vdw_2.SetSelection(0),self.textwidgets[0].SetValue(""),self.textwidgets[1].SetValue(""),self.box1vdw_5.SetSelection(0)
      for i in range(1,5): p2a.box1_vdw[i] = ""
    # box 2, empty string or "none" selected - clear all box 2 widgets, reset choices, clear class vars
    elif widgetname == "box2vdw 1" and chosen_index == 0:
      for item in self.emptyornone_select2:
        item.Hide()
      self.box2vdw_2.SetSelection(0),self.textwidgets[2].SetValue(""),self.textwidgets[3].SetValue(""),self.box2vdw_5.SetSelection(0)
      for i in range(1,5): p2a.box2_vdw[i] = ""
    elif widgetname == "box2vdw 1" and chosen_index == 2:
      for item in self.emptyornone_select2:
        item.Hide()
      self.box2vdw_2.SetSelection(0),self.textwidgets[2].SetValue(""),self.textwidgets[3].SetValue(""),self.box2vdw_5.SetSelection(0)
      for i in range(1,5): p2a.box2_vdw[i] = ""
    # box 1, LJ selected from choice 1 - show column 2 for box 1
    elif widgetname == "box1vdw 1" and chosen_index == 1:
      self.vdwcol2lbl.Show(), self.box1vdw_2.Show()
    # box 2, LJ selected from choice 1- show column 2 for box 2
    elif widgetname == "box2vdw 1" and chosen_index == 1:
      self.vdw2col2lbl.Show(), self.box2vdw_2.Show()
    # box 1, select empty string from choice 2 - hide box 1 columns 3 4 5, clear widgets, clear class variables
    elif widgetname == "box1vdw 2" and chosen_index == 0:
      self.textwidgets[0].SetValue(""),self.textwidgets[1].SetValue(""),self.box1vdw_5.SetSelection(0)
      p2a.box1_vdw[2] = ""
      p2a.box1_vdw[3] = ""
      p2a.box1_vdw[4] = ""
      self.vdw1col3lbl.Hide(),self.vdw1col4lbl.Hide(),self.textwidgets[0].Hide(),self.textwidgets[1].Hide(),self.vdw1col5lbl.Hide(),self.box1vdw_5.Hide()
    # box 1, select cut from choice 2 - hide box 1 column 4 & 5, clear textwidget[3,4], show column 3, set label to Cutoff (A) column 3
    elif widgetname == "box1vdw 2" and chosen_index == 1:
      self.textwidgets[0].Show(),self.vdw1col3lbl.Show(),self.textwidgets[1].Hide(),self.vdw1col4lbl.Hide(),self.vdw1col5lbl.Hide(),self.box1vdw_5.Hide()
      self.textwidgets[1].SetValue(""),self.box1vdw_5.SetSelection(0)
      p2a.box1_vdw[3] = ""
      p2a.box1_vdw[4] = ""
      self.vdw1col3lbl.SetLabel(("Cutoff (%s)" %(p2a.angs)))
    # box 1, select cut_tail from choice 2 - show box 1 col 3 & 5, hide col 4, clear textwidget[3], clear class var at index 3,
    elif widgetname == "box1vdw 2" and chosen_index == 2:
      self.textwidgets[0].Show(),self.vdw1col3lbl.Show(),self.textwidgets[1].Hide(),self.vdw1col4lbl.Hide(),self.vdw1col5lbl.Show(),self.box1vdw_5.Show()
      p2a.box1_vdw[3] = ""
      self.textwidgets[1].SetValue(""),self.vdw1col3lbl.SetLabel(("Cutoff (%s)" %(p2a.angs)))
    # box 1, select cut_switch from choice 2 - show box 1 col 3 & 4, hide col 5, clear class var index 4, set label col 4
    elif widgetname == "box1vdw 2" and chosen_index == 3:
      self.textwidgets[0].Show(),self.vdw1col3lbl.Show(),self.textwidgets[1].Show(),self.vdw1col4lbl.Show(),self.vdw1col5lbl.Hide(),self.box1vdw_5.Hide()
      self.vdw1col3lbl.SetLabel(("Spline On (%s)" %(p2a.angs))),self.vdw1col4lbl.SetLabel(("Spline Off (%s)" %(p2a.angs)))
      p2a.box1_vdw[4] = ""
      self.box1vdw_5.SetSelection(0)
    # box 1, select cut shift from choice 2- show box 1 col 3,hide col 4 & 5, clear class var index 4, set label col 4
    elif widgetname == "box1vdw 2" and chosen_index == 4:
      self.textwidgets[0].Show(),self.vdw1col3lbl.Show(),self.textwidgets[1].Hide(),self.vdw1col4lbl.Hide(),self.vdw1col5lbl.Hide(),self.box1vdw_5.Hide()
      self.vdw1col3lbl.SetLabel(("Cutoff (%s)" %(p2a.angs)))
      p2a.box1_vdw[3] = ""
      p2a.box1_vdw[4] = ""
      self.textwidgets[1].SetValue(""),self.box1vdw_5.SetSelection(0)
    ##### begin bulk of box 2 vdw widgets hiding/showing features
    # box 2, select empty string from choice 2 - hide box 2 column 3 4 5, clear widgets, clear class variables
    elif widgetname == "box2vdw 2" and chosen_index == 0:
      self.textwidgets[2].Hide(),self.textwidgets[3].Hide(),self.box2vdw_5.Hide(),self.vdw2col3lbl.Hide(),self.vdw2col4lbl.Hide(),self.vdw2col5lbl.Hide()
      self.textwidgets[2].SetValue(""),self.textwidgets[3].SetValue(""),self.box2vdw_5.SetSelection(0)
      p2a.box2_vdw[2] = ""
      p2a.box2_vdw[3] = ""
      p2a.box2_vdw[4] = ""
    # box 2, select cut from choice 2 - hide col 4 & 5, clear widget 4 & 5, clear class vars at index  3 & 4, relabel col 3
    # ^^^^^as above statement says, indices between storage variables and widgets do not align. caution.
    elif widgetname == "box2vdw 2" and chosen_index == 1:
      self.textwidgets[2].Show(),self.textwidgets[3].Hide(),self.box2vdw_5.Hide(),self.vdw2col3lbl.Show(),self.vdw2col4lbl.Hide(),self.vdw2col5lbl.Hide()
      self.vdw2col3lbl.SetLabel(("Cutoff (%s)" %(p2a.angs)))
      p2a.box2_vdw[3] = ""
      p2a.box2_vdw[4] = ""
      self.textwidgets[3].SetValue(""),self.box2vdw_5.SetSelection(0)
    # box 2, select cut_tail from choice 2 - show box 2 col 3 & 5, hide col 4, clear class var index 3, set label col 3
    elif widgetname == "box2vdw 2" and chosen_index == 2:
      self.textwidgets[2].Show(),self.vdw2col3lbl.Show(),self.textwidgets[3].Hide(),self.vdw2col4lbl.Hide(),self.vdw2col5lbl.Show(),self.box2vdw_5.Show()
      p2a.box1_vdw[3] = ""
      self.textwidgets[3].SetValue(""),self.vdw2col3lbl.SetLabel(("Cutoff (%s)" %(p2a.angs)))
    # box 2, select cut switch from choice 2 - show box 2 col 3 & 4, hide col 5, clear class var index 4, set labels
    elif widgetname == "box2vdw 2" and chosen_index == 3:
      self.textwidgets[2].Show(),self.vdw2col3lbl.Show(),self.textwidgets[3].Show(),self.vdw2col4lbl.Show(),self.vdw2col5lbl.Hide(),self.box2vdw_5.Hide()
      self.vdw2col3lbl.SetLabel(("Spline On (%s)" %(p2a.angs))),self.vdw2col4lbl.SetLabel(("Spline Off (%s)" %(p2a.angs)))
      p2a.box1_vdw[4] = ""
      self.box2vdw_5.SetSelection(0)
    # box 2, select cut shift from choice 2 - show box 2 col 3, hide 4 & 5, clear class var index 3 & 4 and widgets, set label
    elif widgetname == "box2vdw 2" and chosen_index == 4:
      self.textwidgets[2].Show(),self.vdw2col3lbl.Show(),self.textwidgets[3].Hide(),self.vdw2col4lbl.Hide(),self.vdw2col5lbl.Hide(),self.box2vdw_5.Hide()
      self.vdw2col3lbl.SetLabel(("Cutoff (%s)" %(p2a.angs))),
      p2a.box2_vdw[3] = ""
      p2a.box2_vdw[4] = ""
      self.textwidgets[3].SetValue(""),self.box2vdw_5.SetSelection(0)
    
    ############ charge_style dynamics
    # box 1 and 2 empty string or NONE selections (4 if statements)
    elif widgetname == "box1charge 1" and chosen_index == 0:
      for item in self.eorn_1c:
        item.Hide()
      for i in range(3): p2a.box1_charge[i+1] = ""
      self.textwidgets[4].SetValue(""),self.textwidgets[5].SetValue(""),self.box1charge_2.SetSelection(0)
    elif widgetname == "box2charge 1" and chosen_index == 0:
      for item in self.eorn_2c:
        item.Hide()
      for i in range(3): p2a.box2_charge[i+1] = ""
      self.textwidgets[6].SetValue(""),self.textwidgets[7].SetValue(""),self.box2charge_2.SetSelection(0)
    elif widgetname == "box1charge 1" and chosen_index == 2:
      for item in self.eorn_1c:
        item.Hide()
      for i in range(3): p2a.box1_charge[i+1] = ""
      self.textwidgets[4].SetValue(""),self.textwidgets[5].SetValue(""),self.box1charge_2.SetSelection(0)
    elif widgetname == "box2charge 1" and chosen_index == 2:
      for item in self.eorn_2c:
        item.Hide()
      for i in range(3): p2a.box2_charge[i+1] = ""
      self.textwidgets[6].SetValue(""),self.textwidgets[7].SetValue(""),self.box2charge_2.SetSelection(0)
    # box 1 coulomb selected
    elif widgetname == "box1charge 1" and chosen_index == 1:
      self.col3c1.Show(),self.box1charge_2.Show()
    elif widgetname == "box2charge 1" and chosen_index == 1:
      self.col3c2.Show(),self.box2charge_2.Show()
    # box 1 select empty string from choice 2 - hide textctrls, hide labels 4 and 5, clear index 2 & 3
    elif widgetname == "box1charge 2" and chosen_index == 0:
      self.col4c1.Hide(),self.col5c1.Hide(),self.textwidgets[5].Hide(),self.textwidgets[4].Hide()
      self.textwidgets[4].SetValue(""),self.textwidgets[5].SetValue("")
      p2a.box1_charge[2] = ""
      p2a.box1_charge[3] = ""

    elif widgetname == "box2charge 2" and chosen_index == 0:
      self.col4c2.Hide(),self.col5c2.Hide(),self.textwidgets[6].Hide(),self.textwidgets[7].Hide()
      self.textwidgets[6].SetValue(""),self.textwidgets[7].SetValue("")
      p2a.box2_charge[2] = ""
      p2a.box2_charge[3] = ""
    elif widgetname == "box1charge 2" and chosen_index == 1:
      self.col4c1.Show(),self.col5c1.Show(),self.textwidgets[4].Show(),self.textwidgets[5].Show()
    elif widgetname == "box2charge 2" and chosen_index == 1:
      self.col4c2.Show(),self.col5c2.Show(),self.textwidgets[6].Show(),self.textwidgets[7].Show()
    elif widgetname == "box1charge 2" and chosen_index == 2:
      self.col4c1.Show(),self.textwidgets[4].Show(),self.col5c1.Hide(),self.textwidgets[5].Hide()
      self.textwidgets[5].SetValue("")
      p2a.box1_charge[3] = ""
    elif widgetname == "box2charge 2" and chosen_index == 2:
      self.col4c2.Show(),self.textwidgets[6].Show(),self.col5c2.Hide(),self.textwidgets[7].Hide()
      self.textwidgets[7].SetValue("")
      p2a.box2_charge[3] = ""
    self.Layout()

  def ensemblefunc(self,message,arg2=None):
    msg = message
    if msg == "GEMC":
      self.box2lbl1.Show(),self.box2lbl2.Show(),self.vdw2col1lbl.Show(),self.col2c2.Show(),self.box2vdw_1.Show(),self.box2charge_1.Show()
    else:
      self.box2lbl1.Hide(),self.box2lbl2.Hide(),self.vdw2col1lbl.Hide(),self.col2c2.Hide(),self.box2vdw_1.Hide(),self.box2charge_1.Hide()

# begin Intramolecular panel
class p2b(wx.Panel):
  # define the class variables
  sp1vdw = ["","","",""]
  sp1c =   ["","","",""]
  sp2vdw = ["","","",""]
  sp2c =   ["","","",""]
  sp3vdw = ["","","",""]
  sp3c =   ["","","",""]
  sp4vdw = ["","","",""]
  sp4c =   ["","","",""]
  sp5vdw = ["","","",""]
  sp5c =   ["","","",""]
  sp6vdw = ["","","",""]
  sp6c =   ["","","",""]
  nbr_spec_received = "1" # default value
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent) 
    grid2_b = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid2_b)
    pub.subscribe(self.speciesfunc,"species_msg")
    self.infolbl1 = wx.StaticText(self,label="Select a scaling style:")
    self.infolbl2 = wx.StaticText(self,label="or enter custom values below.")
    
    grid2_b.Add(self.infolbl1,pos=(1,1),span=(1,2))
    grid2_b.Add(self.infolbl2,pos=(2,1),span=(1,2))
    self.amberchk = wx.CheckBox(self,-1,"AMBER",name="amber_param")
    self.charmmchk = wx.CheckBox(self,-1,"CHARMM",name="charmm_param")
    self.Bind(wx.EVT_CHECKBOX,self.checkfunc2b,self.amberchk)
    self.Bind(wx.EVT_CHECKBOX,self.checkfunc2b,self.charmmchk)
    grid2_b.Add(self.amberchk,pos=(1,3))
    grid2_b.Add(self.charmmchk,pos=(1,4))
    
    self.widgets2_b = []
    ii = 12
    jj = 4
    for i in range(ii):
      for j in range(jj):
        self.widgets2_b.append((wx.TextCtrl(self,value="",name=("intrascale %d %d" %(i,j)))))
        grid2_b.Add(self.widgets2_b[jj*i+j],pos=(i+5,j+3))
        self.Bind(wx.EVT_TEXT, self.textfunc2b, self.widgets2_b[4*i+j])

    # create the required labels list
    self.static12lbl = wx.StaticText(self,label="1-2 Scaling")
    self.static13lbl = wx.StaticText(self,label="1-3 Scaling")
    self.static14lbl = wx.StaticText(self,label="1-4 Scaling")
    self.static1nlbl = wx.StaticText(self,label="1-N Scaling")
    self.static1lbls = []
    self.static2lbls = []
    self.static1widgets =[]
    self.static2widgets = []
    kk = 6
    # below, create the static labels vdW and Coulombic iteratively, followed by Species 1,2,..,6
    for i in range(ii):
      if i%2 == 0:
        self.static1lbls.append("van der Waals")
      elif i%2 !=0:
        self.static1lbls.append("Coulombic") 
      self.static1widgets.append(wx.StaticText(self,label=self.static1lbls[i]))
      grid2_b.Add(self.static1widgets[i],pos=(i+5,2))
    for k in range(kk):
      self.static2lbls.append(("Species %d" %(k+1)))
      self.static2widgets.append(wx.StaticText(self,label=self.static2lbls[k]))  
      grid2_b.Add(self.static2widgets[k],pos=(2*k+5,1))
     
    grid2_b.Add(self.static12lbl,pos=(4,3),flag=wx.ALIGN_CENTER)
    grid2_b.Add(self.static13lbl,pos=(4,4),flag=wx.ALIGN_CENTER)
    grid2_b.Add(self.static14lbl,pos=(4,5),flag=wx.ALIGN_CENTER)
    grid2_b.Add(self.static1nlbl,pos=(4,6),flag=wx.ALIGN_CENTER)

  def speciesfunc(self,message,arg2=None):
    p2b.nbr_spec_received = message
    msg = int(message)
    kk = 2*msg
    ii = 12
    jj = 4
    for i in range(ii):
      for j in range(jj):
        if i<kk:
          self.widgets2_b[jj*i+j].Show()
          self.static1widgets[i].Show()
        else:
          self.widgets2_b[jj*i+j].Hide()
          self.static1widgets[i].Hide()
    for i in range(6):
      if i<msg:
        self.static2widgets[i].Show()
      else: self.static2widgets[i].Hide()

  def textfunc2b(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    for j in range(4):
      if widgetname == ("intrascale 0 %d" %j):
        p2b.sp1vdw[j] = text_data
      elif widgetname == ("intrascale 1 %d" %j):
        p2b.sp1c[j] = text_data
      elif widgetname == ("intrascale 2 %d" %j):
        p2b.sp2vdw[j] = text_data
      elif widgetname == ("intrascale 3 %d" %j):
        p2b.sp2c[j] = text_data
      elif widgetname == ("intrascale 4 %d" %j):
        p2b.sp3vdw[j] = text_data
      elif widgetname == ("intrascale 5 %d" %j):
        p2b.sp3c[j] = text_data
      elif widgetname == ("intrascale 6 %d" %j):
        p2b.sp4vdw[j] = text_data
      elif widgetname ==("intrascale 7 %d" %j):
        p2b.sp4c[j] = text_data
      elif widgetname == ("intrascale 8 %d" %j):
        p2b.sp5vdw[j] = text_data
      elif widgetname == ("intrascale 9 %d" %j):
        p2b.sp5c[j] = text_data
      elif widgetname == ("intrascale 10 %d" %j):
        p2b.sp6vdw[j] = text_data
      elif widgetname == ("intrascale 11 %d" %j):
        p2b.sp6c[j] = text_data

  def checkfunc2b(self,event):
    event_object=event.GetEventObject()
    widgetname = event_object.GetName()
    chk_data = event.IsChecked()
    n = int(p2b.nbr_spec_received)
    set_blank = ["","","",""]
    lj_settings = ["","","",""]
    elec_settings = ["","","",""]
    vdw_list = [p2b.sp1vdw , p2b.sp2vdw , p2b.sp3vdw , p2b.sp4vdw , p2b.sp5vdw , p2b.sp6vdw]
    coul_list = [p2b.sp1c , p2b.sp2c , p2b.sp3c , p2b.sp4c , p2b.sp5c , p2b.sp6c]
    if chk_data == True:
      if widgetname == "amber_param":
        self.charmmchk.SetValue(False)
        lj_settings = ["0.0","0.0","0.5","1.0"]
        elec_settings = ["0.0","0.0","0.8333","1.0"]
        for i in range(n):
          vdw_list[i] = lj_settings
          coul_list[i] = elec_settings
        for i in range(6-n):
          vdw_list[i+n] = set_blank
          coul_list[i+n] = set_blank
        for i in range(2*n):
          for j in range(4):
            if i%2 == 0:
              self.widgets2_b[4*i+j].SetValue(lj_settings[j])
            if i%2 != 0:
              self.widgets2_b[4*i+j].SetValue(elec_settings[j])
      elif widgetname == "charmm_param":
        self.amberchk.SetValue(False)
        lj_settings = ["0.0","0.0","0.0","1.0"]
        elec_settings = ["0.0","0.0","0.0","1.0"]
        for i in range(n):
          vdw_list[i] = lj_settings
          coul_list[i] = elec_settings
        for i in range(6-n):
          vdw_list[i+n] = set_blank
          coul_list[i+n] = set_blank
        for i in range(2*n):
          for j in range(4):
            if i%2 == 0:
              self.widgets2_b[4*i+j].SetValue(lj_settings[j])
            if i%2 != 0:
              self.widgets2_b[4*i+j].SetValue(elec_settings[j])
    if chk_data == False:
      for i in range(6):
        vdw_list[i] = set_blank
        coul_list[i] = set_blank
      for i in range(12):
        for j in range(4):
          self.widgets2_b[4*i+j].SetValue("")
    self.Layout()


# PanelTwo is the parent panel for the Inter- and Intra-molecular Parameters tabs
class PanelTwo(wx.Panel):
  # initialize the panel
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

    # Define the notebook and add the inter- and intra-molecular panels
    nested_notebook2 = wx.Notebook(self)
    nested_page1 = p2a(nested_notebook2)
    nested_page2 = p2b(nested_notebook2)

    nested_notebook2.AddPage(nested_page1,"Intermolecular")
    nested_notebook2.AddPage(nested_page2,"Intramolecular")

    nested_sizer = wx.BoxSizer()
    nested_sizer.Add(nested_notebook2,1,wx.EXPAND)
    self.SetSizer(nested_sizer)

###################################################################################
###################################################################################

# begin Translation Probabilities notebook tab

class p3a(wx.Panel):
  # class variables for data storage
  prob_trans = "0"
  ptb1s = ["","","","","",""]
  ptb2s = ["","","","","",""]

  ensemble_received = ""
  nbr_spec_received = "1"
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3_a = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3_a)
    
    pub.subscribe(self.ensemblefunc,"ensemble_msg")
    pub.subscribe(self.speciesfunc,"species_msg")

    self.widgetnames = [] #list of names for text control widgets
    self.widgets3_a = [] #list of text control widgets
    ii = 2
    jj = 6
    for i in range(ii):
      for j in range(jj):
        self.widgetnames.append("transprob %d %d" %(i,j))
        self.widgets3_a.append((wx.TextCtrl(self,value="",name=(self.widgetnames[jj*i+j]))))
        grid3_a.Add(self.widgets3_a[jj*i+j],pos=(j+5,i+3),flag=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_TEXT, self.textfunc3a, self.widgets3_a[jj*i+j])

    self.splbl_list = []
    for j in range(jj):
      self.splbl_list.append((wx.StaticText(self,label=("Species %d:" %(j+1)))))
      grid3_a.Add(self.splbl_list[j],pos=(j+5,2),flag=wx.ALIGN_RIGHT)

    self.noticelbl = wx.StaticText(self,label="Please note that the sum of the move probabilities for all move types used must sum to 1.")
    grid3_a.Add(self.noticelbl,pos=(1,2),span=(1,6))
    
    self.box1lbl = wx.StaticText(self,label="Box 1")
    self.box2lbl = wx.StaticText(self,label="Box 2")
    grid3_a.Add(self.box1lbl,pos=(4,3),flag=wx.ALIGN_CENTER)
    grid3_a.Add(self.box2lbl,pos=(4,4),flag=wx.ALIGN_CENTER)

    self.transproblbl = wx.StaticText(self,label="Move Probability:")
    grid3_a.Add(self.transproblbl,pos=(0,2),flag=wx.ALIGN_RIGHT)

    self.transprobquery = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name="transprobquery")
    grid3_a.Add(self.transprobquery,pos=(0,3))
    self.Bind(wx.EVT_TEXT, self.textfunc3a,self.transprobquery)

    string3a = "Enter the maximum displacement (Angstroms) allowed for each species in each box below."
    self.instructionlbl = wx.StaticText(self,label=string3a)
    grid3_a.Add(self.instructionlbl,pos=(3,2),span=(1,6))

  def textfunc3a(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    for k in range(6):
      if widgetname == ("transprob 0 %d" %k):
        p3a.ptb1s[k] = text_data
      elif widgetname ==("transprob 1 %d" %k):
        p3a.ptb2s[k] = text_data
    if widgetname == "transprobquery":
      p3a.prob_trans = text_data

  def ensemblefunc(self,message,arg2=None):
    msg = message
    ii = 2
    jj = int(p3a.nbr_spec_received)
    p3a.ensemble_received = message
    if msg == "GEMC":
      for j in range(jj):
        self.splbl_list[j].Show()
        self.widgets3_a[j].Show()
        self.widgets3_a[j+6].Show()
      for j in range(6-jj):
        self.splbl_list[j+jj].Hide()
        self.widgets3_a[j+jj].Hide()
        self.widgets3_a[j+6+jj].Hide()
      self.box2lbl.Show()
    else:
      for j in range(6):
        self.widgets3_a[j+6].Hide()
      self.box2lbl.Hide()
      for j in range(jj):
        self.splbl_list[j].Show()
        self.widgets3_a[j].Show()
      for j in range(6-jj):
        self.splbl_list[j+jj].Hide()
        self.widgets3_a[j+jj].Hide()
    self.Layout()

  def speciesfunc(self,message,arg2=None):
    msg = int(message)
    p3a.nbr_spec_received = message
    for i in range(msg):
      self.splbl_list[i].Show()
      self.widgets3_a[i].Show()
    for i in range(6-msg):
      self.splbl_list[i+msg].Hide()
      self.widgets3_a[i+msg].Hide()
    
    if p3a.ensemble_received !="GEMC":
      for i in range(6):
        self.widgets3_a[6+i].Hide()
      self.box2lbl.Hide(),
    else:
      for i in range(msg):
        self.widgets3_a[6+i].Show()
      for i in range(6-msg):
        self.widgets3_a[i+msg].Hide()
        self.widgets3_a[i+msg+6].Hide()
      self.box2lbl.Show()
    self.Layout()


# begin Rotation Probabilities notebook tab
class p3b(wx.Panel):
  # class variables for data storage
  prob_rot = "0"
  prb1s = ["","","","","",""]
  prb2s = ["","","","","",""]

  ensemble_received = ""
  nbr_spec_received = "1"
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    pub.subscribe(self.ensemblefunc,"ensemble_msg")
    pub.subscribe(self.speciesfunc,"species_msg")
    grid3_b = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3_b)
    self.widgetnames = [] #list of names for text control widgets
    self.widgets3_b = [] #list of text control widgets
    ii = 2
    jj = 6
    for i in range(ii):
      for j in range(jj):
        self.widgetnames.append("rotprob %d %d" %(i,j))
        self.widgets3_b.append((wx.TextCtrl(self,value="",name = (self.widgetnames[jj*i+j]))))
        grid3_b.Add(self.widgets3_b[jj*i+j],pos=(j+5,i+3),flag=wx.ALIGN_LEFT)
        self.Bind(wx.EVT_TEXT, self.textfunc3b, self.widgets3_b[jj*i+j])

    self.splbl_list = []
    for j in range(jj):
      self.splbl_list.append((wx.StaticText(self,label=("Species %d:" %(j+1)))))
      grid3_b.Add(self.splbl_list[j],pos=(j+5,2),flag=wx.ALIGN_RIGHT)

    self.noticelbl = wx.StaticText(self,label="Please note that the sum of the move probabilities for all move types used must sum to 1.")
    grid3_b.Add(self.noticelbl,pos=(1,2),span=(1,6))
    
    self.box1lbl = wx.StaticText(self,label="Box 1")
    self.box2lbl = wx.StaticText(self,label="Box 2")
    grid3_b.Add(self.box1lbl,pos=(4,3),flag=wx.ALIGN_CENTER)
    grid3_b.Add(self.box2lbl,pos=(4,4),flag=wx.ALIGN_CENTER)

    self.rotproblbl = wx.StaticText(self,label="Move Probability:")
    grid3_b.Add(self.rotproblbl,pos=(0,2),flag=wx.ALIGN_RIGHT)

    self.rotprobquery = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name="rotprobquery")
    grid3_b.Add(self.rotprobquery,pos=(0,3))
    self.Bind(wx.EVT_TEXT, self.textfunc3b,self.rotprobquery)

    string3b = "Enter the maximum rotational width in degrees for each species in each box below."
    self.instructionlbl = wx.StaticText(self,label=string3b)
    grid3_b.Add(self.instructionlbl,pos=(3,2),span=(1,6))

  def textfunc3b(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    for k in range(6):
      if widgetname == ("rotprob 0 %d" %k):
        p3b.prb1s[k] = text_data
      elif widgetname ==("rotprob 1 %d" %k):
        p3b.prb2s[k] = text_data
    if widgetname == "rotprobquery":
      p3b.prob_rot = text_data

  def ensemblefunc(self,message,arg2=None):
    msg = message
    ii = 2
    jj = int(p3b.nbr_spec_received)
    p3b.ensemble_received = message
    if msg == "GEMC":
      for j in range(jj):
        self.splbl_list[j].Show()
        self.widgets3_b[j].Show()
        self.widgets3_b[j+6].Show()
      for j in range(6-jj):
        self.splbl_list[j+jj].Hide()
        self.widgets3_b[j+jj].Hide()
        self.widgets3_b[j+6+jj].Hide()
      self.box2lbl.Show()
    else:
      for j in range(6):
        self.widgets3_b[j+6].Hide()
      self.box2lbl.Hide()
      for j in range(jj):
        self.splbl_list[j].Show()
        self.widgets3_b[j].Show()
      for j in range(6-jj):
        self.splbl_list[j+jj].Hide()
        self.widgets3_b[j+jj].Hide()
    self.Layout()

  def speciesfunc(self,message,arg2=None):
    msg = int(message)
    p3b.nbr_spec_received = message
    for i in range(msg):
      self.splbl_list[i].Show()
      self.widgets3_b[i].Show()
    for i in range(6-msg):
      self.splbl_list[i+msg].Hide()
      self.widgets3_b[i+msg].Hide()
    
    if p3b.ensemble_received !="GEMC":
      for i in range(6):
        self.widgets3_b[6+i].Hide()
      self.box2lbl.Hide(),
    else:
      for i in range(msg):
        self.widgets3_b[6+i].Show()
      for i in range(6-msg):
        self.widgets3_b[i+msg].Hide()
        self.widgets3_b[i+msg+6].Hide()
      self.box2lbl.Show()
    self.Layout()

# begin Regrowth Probabilities notebook tab
class p3c(wx.Panel):

  prob_regrowth = "0"
  pregrow_s = ["","","","","",""]
  nbr_spec_received = "1"
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    pub.subscribe(self.speciesfunc,"species_msg")
    grid3_c = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3_c)
    
    self.regrowproblbl = wx.StaticText(self,label="Move Probability:")
    grid3_c.Add(self.regrowproblbl,pos=(0,2),flag=wx.ALIGN_RIGHT)
    self.regrowprobquery = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name="regrowprobquery")
    grid3_c.Add(self.regrowprobquery,pos=(0,3))
    self.Bind(wx.EVT_TEXT, self.textfunc3c,self.regrowprobquery)

    string3c = "Enter the relative probability of regrowth for each species below."
    self.instructionlbl = wx.StaticText(self,label=string3c)
    grid3_c.Add(self.instructionlbl,pos=(3,2),span=(1,6))
    
    string3c_2 = "Note that the relative probabilities below must sum to 1."
    self.string3c_2lbl = wx.StaticText(self,label=string3c_2)
    grid3_c.Add(self.string3c_2lbl,pos=(4,2),span=(1,6))
    
    self.noticelbl = wx.StaticText(self,label="Please note that the sum of the move probabilities for all move types used must sum to 1.")
    grid3_c.Add(self.noticelbl,pos=(1,2),span=(1,6))
    
    self.staticlbl = []
    self.textwidgets3_c = []
    for i in range(6):
      self.textwidgets3_c.append((wx.TextCtrl(self,value="",name=("regprob %d" %(i)))))
      grid3_c.Add(self.textwidgets3_c[i],pos=(i+5,3))
      self.Bind(wx.EVT_TEXT, self.textfunc3c,self.textwidgets3_c[i])
      self.staticlbl.append((wx.StaticText(self,label=("Species %d:" %(i+1)))))
      grid3_c.Add(self.staticlbl[i],pos=(i+5,2),flag=wx.ALIGN_RIGHT)

  def textfunc3c(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    for i in range(6):
      if widgetname == ("regprob %d" %i):
        if not text_data:
          text_data = ""
        p3c.pregrow_s[i] = text_data
    if widgetname == "regrowprobquery":
      p3c.prob_regrowth = text_data

  def speciesfunc(self,message,arg2=None):
    msg = int(message)
    for i in range(msg):
      self.textwidgets3_c[i].Show()
      self.staticlbl[i].Show()
    for i in range(6-msg):
      self.textwidgets3_c[i+msg].Hide()
      self.staticlbl[i+msg].Hide()
    self.Layout()

# begin Volume Probabilities notebook tab
class p3d(wx.Panel):

  prob_vol = "0"
  pvol_box = ["",""]

  ensemble_received = ""
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3_d = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3_d)

    pub.subscribe(self.ensemblefunc,"ensemble_msg")
    
    self.volproblbl = wx.StaticText(self,label="Move Probability:")
    grid3_d.Add(self.volproblbl, pos=(0,2),flag=wx.ALIGN_RIGHT)
    self.volprobquery = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name = "volprobquery")
    grid3_d.Add(self.volprobquery,pos=(0,3))
    self.Bind(wx.EVT_TEXT, self.textfunc3d, self.volprobquery)

    string3d_1 = "Enter the maximum volume displacements in Angstroms^3 for the simulation box(es) below."
    self.instructionlbl = wx.StaticText(self,label=string3d_1)
    grid3_d.Add(self.instructionlbl,pos=(3,2),span=(1,6))
    
    string3d_2 = "This flag is required for NPT-MC, GEMC-NPT, and GEMC-NVT simulations, and should not be used for other simulation types."
    self.instruction2lbl = wx.StaticText(self,label=string3d_2)
    grid3_d.Add(self.instruction2lbl,pos=(4,2),span=(1,8))

    self.noticelbl = wx.StaticText(self,label="Please note that the sum of the move probabilities for all move types used must sum to 1.")
    grid3_d.Add(self.noticelbl,pos=(1,2),span=(1,6))
    self.noticelbl.Hide()
    
    self.box1lbl = wx.StaticText(self,label="Box 1: ")
    self.box2lbl = wx.StaticText(self,label="Box 2: ")
    self.box1query = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name="box1disp")
    self.box2query = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name="box2disp")
    grid3_d.Add(self.box1lbl,pos=(5,2),flag=wx.ALIGN_RIGHT)
    grid3_d.Add(self.box2lbl,pos=(6,2),flag=wx.ALIGN_RIGHT)
    grid3_d.Add(self.box1query,pos=(5,3))
    grid3_d.Add(self.box2query,pos=(6,3))

    self.notneedlbl = wx.StaticText(self,label="This flag is allowed only for NPT-MC, GEMC-NPT, and GEMC-NVT, and is not permitted for the currently selected ensemble.")
    grid3_d.Add(self.notneedlbl,pos=(6,4),span=(1,6))
    self.Bind(wx.EVT_TEXT, self.textfunc3d, self.box1query)
    self.Bind(wx.EVT_TEXT, self.textfunc3d, self.box2query)
 
  def textfunc3d(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "box1disp":
      p3d.pvol_box[0] = text_data
    elif widgetname =="box2disp":
      p3d.pvol_box[1] = text_data
    elif widgetname == "volprobquery":
      p3d.prob_vol = text_data
      
  def ensemblefunc(self,message,arg2=None):
    msg = message
    itemlist = [self.volproblbl,self.volprobquery,self.instructionlbl,self.instruction2lbl,self.box1lbl,self.box1query,self.box2lbl,self.box2query]
    if msg == "GEMC":
      self.noticelbl.Show()
      for i in range(8):
        itemlist[i].Show()
      self.notneedlbl.Hide()
    elif msg == "NPT_MC":
      self.noticelbl.Show()
      for i in range(6):
        itemlist[i].Show()
      for i in range(2):
        itemlist[i+6].Hide()
      self.notneedlbl.Hide()
    else:
      self.noticelbl.Hide()
      for item in itemlist:
        item.Hide()
      self.notneedlbl.Show()
    self.Layout()

# begin Insertion (and Deletion by extension) Probabilities notebook tab
class p3e(wx.Panel):
  prob_insert = "0"
  pinsert_spec = ["","","","","",""]
  
  nbr_spec_received = "1"
  ensemble_received = ""
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3_e = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3_e)

    pub.subscribe(self.ensemblefunc,"ensemble_msg")
    pub.subscribe(self.speciesfunc,"species_msg")

    infostring1 = "Select the insertion method for each relevant species in the simulation from the box(es) below."
    self.infostring1lbl = wx.StaticText(self,label=infostring1)
    grid3_e.Add(self.infostring1lbl,pos=(3,2),span=(1,6))

    self.insproblbl = wx.StaticText(self,label="Move Probability:")
    self.insprobquery = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name = "insprobquery")
    grid3_e.Add(self.insproblbl, pos=(0,2),flag=wx.ALIGN_RIGHT)
    grid3_e.Add(self.insprobquery,pos=(0,3))
    self.Bind(wx.EVT_TEXT, self.textfunc3e, self.insprobquery)

    self.noticelbl = wx.StaticText(self,label="Please note that the sum of the move probabilities for all move types used must sum to 1.")
    grid3_e.Add(self.noticelbl,pos=(1,2),span=(1,6))
    self.notice2lbl = wx.StaticText(self,label="Additionally, insertion moves define an equal probability of deletion, and so its probability should be counted twice when summing to 1.")
    grid3_e.Add(self.notice2lbl,pos=(2,2),span=(1,6))
    
    self.notneedlbl = wx.StaticText(self,label="This flag is allowed only for GCMC simulations, and is not permitted for the currently selected ensemble.")
    grid3_e.Add(self.notneedlbl,pos=(6,4),span=(1,6))
    self.notneedlbl.Hide()
    
    ins_options = ["","reservoir","none"]
    self.textwidgets3_e = []
    self.staticlbl = []
    for i in range(6):
      self.textwidgets3_e.append((wx.Choice(self,choices=ins_options,name=("insprob %d" %(i)))))
      grid3_e.Add(self.textwidgets3_e[i],pos=(i+5,3))
      self.Bind(wx.EVT_CHOICE, self.choicefunc3e,self.textwidgets3_e[i])
      self.staticlbl.append((wx.StaticText(self,label=("Species %d:" %(i+1)))))
      grid3_e.Add(self.staticlbl[i],pos=(i+5,2),flag=wx.ALIGN_RIGHT)

  def textfunc3e(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "insprobquery":
      p3e.prob_insert = text_data

  def choicefunc3e(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    choice_data = event.GetString()
    for i in range(6):
      if widgetname == ("insprob %d" %i):
        p3e.pinsert_spec[i] = choice_data

  def ensemblefunc(self,message,arg2=None):
    p3e.ensemble_received = message
    ii = int(p3e.nbr_spec_received)
    if message == "GCMC":
      self.noticelbl.Show()
      self.notice2lbl.Show()
      for i in range(ii):
        self.textwidgets3_e[i].Show()
        self.staticlbl[i].Show()
      for i in range(6-ii):
        self.textwidgets3_e[i+ii].Hide()
        self.staticlbl[i+ii].Hide()
      self.insproblbl.Show(),self.insprobquery.Show(),self.infostring1lbl.Show(),self.notneedlbl.Hide()
    else:
      self.noticelbl.Hide(), self.notice2lbl.Hide()
      for i in range(6):
        self.textwidgets3_e[i].Hide(),self.staticlbl[i].Hide()
      self.insproblbl.Hide(),self.insprobquery.Hide(),self.infostring1lbl.Hide(),self.notneedlbl.Show()
    self.Layout()

  def speciesfunc(self,message,arg2=None):
    p3e.nbr_spec_received = message
    ii = int(message)
    if p3e.ensemble_received == "GCMC":
      self.noticelbl.Show(),self.notice2lbl.Show()
      for i in range(ii):
        self.textwidgets3_e[i].Show()
        self.staticlbl[i].Show()
      for i in range(6-ii):
        self.textwidgets3_e[i+ii].Hide()
        self.staticlbl[i+ii].Hide()
      self.insproblbl.Show(),self.insprobquery.Show(),self.infostring1lbl.Show(),self.notneedlbl.Hide()
    else:
      self.noticelbl.Hide(),self.notice2lbl.Hide()
      for i in range(6):
        self.textwidgets3_e[i].Hide(),self.staticlbl[i].Hide()
      self.insproblbl.Hide(),self.insprobquery.Hide(),self.infostring1lbl.Hide(),self.notneedlbl.Show()
    self.Layout()

# begin Swap Probabilities notebook tab
class p3f(wx.Panel):
  prob_swap = "0"
  pswap_spec = ["","","","","",""]

  ensemble_received = ""
  nbr_spec_received = "1"
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

    pub.subscribe(self.ensemblefunc,"ensemble_msg")
    pub.subscribe(self.speciesfunc,"species_msg")

    grid3_f = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3_f)

    infostring1 = "Select the swap method for each relevant species in the simulation from the box(es) below."
    self.infostring1lbl = wx.StaticText(self,label=infostring1)
    grid3_f.Add(self.infostring1lbl,pos=(3,2),span=(1,6))

    self.swapproblbl = wx.StaticText(self,label="Move Probability:")
    self.swapprobquery = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name = "swapprobquery")
    grid3_f.Add(self.swapproblbl, pos=(0,2),flag=wx.ALIGN_RIGHT)
    grid3_f.Add(self.swapprobquery,pos=(0,3))
    self.Bind(wx.EVT_TEXT, self.textfunc3f, self.swapprobquery)
    
    self.noticelbl = wx.StaticText(self,label="Please note that the sum of the move probabilities for all move types used must sum to 1.")
    grid3_f.Add(self.noticelbl,pos=(1,2),span=(1,6))
    
    self.notneedlbl = wx.StaticText(self,label="This flag is allowed only for GEMC simulations, and is not permitted for the currently selected ensemble.")
    grid3_f.Add(self.notneedlbl,pos=(6,4),span=(1,6))
    self.notneedlbl.Hide()
    
    swap_options = ["","reservoir","none"]
    self.choicewidgets3_f = []
    self.staticlbl = []
    for i in range(6):
      self.choicewidgets3_f.append((wx.Choice(self,choices=swap_options,name=("swapprob %d" %(i)))))
      grid3_f.Add(self.choicewidgets3_f[i],pos=(i+5,3))
      self.Bind(wx.EVT_CHOICE, self.choicefunc3f,self.choicewidgets3_f[i])
      self.staticlbl.append((wx.StaticText(self,label=("Species %d:" %(i+1)))))
      grid3_f.Add(self.staticlbl[i],pos=(i+5,2),flag=wx.ALIGN_RIGHT)

  def textfunc3f(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "swapprobquery":
      p3f.prob_swap = text_data

  def choicefunc3f(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    choice_data = event.GetString()
    for i in range(6):
      if widgetname == ("swapprob %d" %i):
        p3f.pswap_spec[i] = choice_data

  def speciesfunc(self,message,arg2=None):
    p3f.nbr_spec_received = message
    msg = int(message)
    if p3f.ensemble_received == "GEMC":
      self.notneedlbl.Hide()
      self.noticelbl.Show()
      for i in range(msg):
        self.staticlbl[i].Show(),self.choicewidgets3_f[i].Show()
      self.infostring1lbl.Show(),self.swapproblbl.Show(),self.swapprobquery.Show()
      for i in range(6-msg):
        self.staticlbl[i+msg].Hide(),self.choicewidgets3_f[i+msg].Hide()
    else:
      self.noticelbl.Hide()
      for item in self.choicewidgets3_f:
        item.Hide()
      for item in self.staticlbl:
        item.Hide()
      self.infostring1lbl.Hide(),self.swapproblbl.Hide(),self.swapprobquery.Hide()
      self.notneedlbl.Show()
    self.Layout()

  def ensemblefunc(self,message,arg2=None):
    p3f.ensemble_received = message
    ii = int(p3f.nbr_spec_received)
    if message == "GEMC":
      self.notneedlbl.Hide()
      self.noticelbl.Show()
      for i in range(ii):
        self.staticlbl[i].Show(),self.choicewidgets3_f[i].Show()
      self.infostring1lbl.Show(),self.swapproblbl.Show(),self.swapprobquery.Show()
      for i in range(6-ii):
        self.staticlbl[i+ii].Hide(),self.choicewidgets3_f[i+ii].Hide()
    else:
      self.noticelbl.Hide()
      for item in self.choicewidgets3_f:
        item.Hide()
      for item in self.staticlbl:
        item.Hide()
      self.infostring1lbl.Hide(),self.swapproblbl.Hide(),self.swapprobquery.Hide()
      self.notneedlbl.Show()
    self.Layout()

# begin Ring Probabilities notebook tab
class p3g(wx.Panel):
  prob_flip = "0"
  max_angle = ""
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

    grid3_g = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3_g)

    self.fliplbl = wx.StaticText(self,label="Move Probability:")
    self.flipquery = wx.TextCtrl(self,value="",name="flipquery")
    grid3_g.Add(self.fliplbl,pos=(0,2),flag=wx.ALIGN_RIGHT)
    grid3_g.Add(self.flipquery,pos=(0,3))

    infostring1 = "Enter the maximum angular displacement for the move in degrees below."
    self.infostringlbl = wx.StaticText(self,label=infostring1)
    grid3_g.Add(self.infostringlbl,pos=(3,2),span=(1,6))

    self.textwidget = wx.TextCtrl(self,value="",name="maxdisp_angle")
    grid3_g.Add(self.textwidget,pos=(5,4))
    self.Bind(wx.EVT_TEXT,self.textfunc3g,self.textwidget)

    self.textlbl = wx.StaticText(self,label="Max Displacement:")
    grid3_g.Add(self.textlbl,pos=(5,2),span=(1,2),flag=wx.ALIGN_RIGHT)

  def textfunc3g(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "flipquery":
      p3g.prob_flip = text_data
    if widgetname == "maxdisp_angle":
      p3g.max_angle = text_data

# PanelThree is the parent panel for the probabilities tabs
class PanelThree(wx.Panel):
  # initialize the panel
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

    # Define the notebook
    nest = wx.Notebook(self)
    nest1 = p3a(nest)
    nest2 = p3b(nest)
    nest3 = p3c(nest)
    nest4 = p3d(nest)
    nest5 = p3e(nest)
    nest6 = p3f(nest)
    nest7 = p3g(nest)

    nest.AddPage(nest1,"Translation")
    nest.AddPage(nest2,"Rotation")
    nest.AddPage(nest3,"Regrowth")
    nest.AddPage(nest4,"Volume")
    nest.AddPage(nest5,"Insertion")
    nest.AddPage(nest6,"Swap")
    nest.AddPage(nest7,"Ring Flip")

    nested_sizer = wx.BoxSizer()
    nested_sizer.Add(nest,1,wx.EXPAND)
    self.SetSizer(nested_sizer)

###################################################################################
###################################################################################

#molecule files
class p4a(wx.Panel):
  mcf_s1 = ["",""]
  mcf_s2 = ["",""]
  mcf_s3 = ["",""]
  mcf_s4 = ["",""]
  mcf_s5 = ["",""]
  mcf_s6 = ["",""]
  
  simulation_dirname = ""
  nbr_spec_received = "1"
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    self.grid4a = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(self.grid4a)
    pub.subscribe(self.speciesfunc,"species_msg")
    pub.subscribe(self.directoryfunc,"directory_msg")
    
    self.dirname= ""
    self.mcflbl = wx.StaticText(self,label="Select the MCF files below.  Enter the maximum number of anticipated molecules for each species - this \nnumber will be used for memory allocation purposes.")
    self.grid4a.Add(self.mcflbl,pos=(1,1),span=(2,5))
    #begin in row 3
    self.editlbl = wx.StaticText(self,label="                     Selection                     ") #space to allow large filenames to be seen
    self.grid4a.Add(self.editlbl,pos=(3,3),flag=wx.ALIGN_CENTER,span=(1,1))
    
    self.selectlbl = wx.StaticText(self,label="Select")
    self.grid4a.Add(self.selectlbl,pos=(3,2),flag=wx.ALIGN_CENTER)
  
    self.nbr_moleculeslbl = wx.StaticText(self,label="# Molecules")
    self.grid4a.Add(self.nbr_moleculeslbl,pos=(3,4),flag=wx.ALIGN_CENTER,span=(1,1))
  
    self.labels4a= []
    self.buttons4a = []
    self.textwidgets4a_disp = []
    self.textwidgets4a_mol = []

    alertstring1 = "Please select a simulation directory on the page Basic Information > Page 1.\nThis is required to continue."
    self.alertstring1lbl = wx.StaticText(self,label=alertstring1)
    self.grid4a.Add(self.alertstring1lbl,pos=(6,5),span=(2,4))
    
    for i in range(6):
      self.labels4a.append(wx.StaticText(self,label=("Species %d:" %(i+1))))
      self.grid4a.Add(self.labels4a[i],pos=(i+4,1))
      self.buttons4a.append(wx.Button(self,label="Select MCF File",name=("buttonspec %d" %(i+1))))
      self.grid4a.Add(self.buttons4a[i],pos=(i+4,2))
      self.textwidgets4a_disp.append(wx.TextCtrl(self,value="",name=("view %d" %(i+1)),style=wx.TE_READONLY))
      self.grid4a.Add(self.textwidgets4a_disp[i],pos=(i+4,3),flag=wx.EXPAND,span=(1,1))
      self.textwidgets4a_mol.append(wx.TextCtrl(self,value="",name=("nbr %d" %(i+1))))
      self.grid4a.Add(self.textwidgets4a_mol[i],pos=(i+4,4),flag=wx.EXPAND)
      self.Bind(wx.EVT_BUTTON,self.selectMCF,self.buttons4a[i])
      self.Bind(wx.EVT_TEXT,self.textfunc4a,self.textwidgets4a_mol[i])
    for i in range(1,6):
      self.labels4a[i].Hide(),self.buttons4a[i].Hide(),self.textwidgets4a_disp[i].Hide(),self.textwidgets4a_mol[i].Hide()

    self.labels4a[0].Hide(),self.buttons4a[0].Hide(),self.textwidgets4a_disp[0].Hide(),self.textwidgets4a_mol[0].Hide(),self.editlbl.Hide(),self.selectlbl.Hide(),self.nbr_moleculeslbl.Hide(),self.mcflbl.Hide()

  def directoryfunc(self,message,arg2=None):
    p4a.simulation_dirname = message
    ii = int(p4a.nbr_spec_received)
    for i in range(ii):
      self.labels4a[i].Show,self.buttons4a[i].Show(),self.textwidgets4a_disp[i].Show(),self.textwidgets4a_mol[i].Show()
    self.editlbl.Show(),self.selectlbl.Show(),self.nbr_moleculeslbl.Show(),self.mcflbl.Show(),self.alertstring1lbl.Hide()
    for i in range(ii,6):
        self.labels4a[i].Hide(),self.buttons4a[i].Hide(),self.textwidgets4a_disp[i].Hide(),self.textwidgets4a_mol[i].Hide()
    self.Layout()
  
  def speciesfunc(self,message,arg2=None):
    p4a.nbr_spec_received = message
    if not p4a.simulation_dirname:
      self.alertstring1lbl.Show()
      self.labels4a[0].Hide(),self.buttons4a[0].Hide(),self.textwidgets4a_disp[0].Hide(),self.textwidgets4a_mol[0].Hide(),self.editlbl.Hide(),self.selectlbl.Hide(),self.nbr_moleculeslbl.Hide(),self.mcflbl.Hide()
    else:
      msg = int(message)
      for i in range(msg):
        self.labels4a[i].Show(),self.buttons4a[i].Show(),self.textwidgets4a_disp[i].Show(),self.textwidgets4a_mol[i].Show()
      self.editlbl.Show(),self.selectlbl.Show(),self.nbr_moleculeslbl.Show(),self.mcflbl.Show(),self.alertstring1lbl.Hide()
      for i in range(msg,6):
        self.labels4a[i].Hide(),self.buttons4a[i].Hide(),self.textwidgets4a_disp[i].Hide(),self.textwidgets4a_mol[i].Hide()
    self.Layout()

  def selectMCF(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    dlg = wx.FileDialog(self,"Select MCF File",self.dirname,"","*.mcf",wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
      self.filepathdata = dlg.GetPath()
      self.current_file_name = dlg.GetFilename()
      f = file(self.filepathdata)
      nfrags_ind = ""
      nfrags_data = ""
      for ind,line in enumerate(f):
        if '#' in line and 'Fragment_Info' in line:
          nfrags_ind = ind+1
        if nfrags_ind:
          if ind == nfrags_ind:
            line_data = line.split()
            nfrags_data = int(line_data[0])
      if not nfrags_data:
        print "Number of fragments not found.  Please enter fragment information manually upon creation of input file."
        nfrags_data = 0
      file_data = os.path.relpath(self.filepathdata, p4a.simulation_dirname)
      for i in range(6):
        if widgetname == ("buttonspec %d" %(i+1)):
          self.textwidgets4a_disp[i].SetValue(file_data)
      if widgetname == "buttonspec 1":
        p4a.mcf_s1[0] = str(file_data)
        pub.sendMessage("nfrags_spec1", message = nfrags_data)
      elif widgetname == "buttonspec 2":
        p4a.mcf_s2[0] = str(file_data)
        pub.sendMessage("nfrags_spec2",message = nfrags_data)
      elif widgetname == "buttonspec 3":
        p4a.mcf_s3[0] = str(file_data)
        pub.sendMessage("nfrags_spec3",message = nfrags_data)
      elif widgetname == "buttonspec 4":
        p4a.mcf_s4[0] = str(file_data)
        pub.sendMessage("nfrags_spec4",message=nfrags_data)
      elif widgetname == "buttonspec 5":
        p4a.mcf_s5[0] = str(file_data)
        pub.sendMessage("nfrags_spec5",message=nfrags_data)
      elif widgetname == "buttonspec 6":
        p4a.mcf_s6[0] = str(file_data)
        pub.sendMessage("nfrags_spec6",message=nfrags_data)

  def textfunc4a(self,event):
    event_object=event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "nbr 1":
      p4a.mcf_s1[1] = text_data
    elif widgetname == "nbr 2":
      p4a.mcf_s2[1] = text_data
    elif widgetname == "nbr 3":
      p4a.mcf_s3[1] = text_data
    elif widgetname == "nbr 4":
      p4a.mcf_s4[1] = text_data
    elif widgetname == "nbr 5":
      p4a.mcf_s5[1] = text_data
    elif widgetname == "nbr 6":
      p4a.mcf_s6[1] = text_data


#Fragment Files 1
class p4b(wx.Panel):
  s1_nfrags = ""
  s2_nfrags = ""
  s3_nfrags = ""
  s1_frags = ["","","","","","","","","","","","","","",""]
  s2_frags = ["","","","","","","","","","","","","","",""]
  s3_frags = ["","","","","","","","","","","","","","",""]
  
  frag_index = []
  for i in range(90):
    frag_index.append(int(i+1))
  
  simulation_dirname = ""
  nbr_spec_received = "1"
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid4b = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid4b)

    pub.subscribe(self.directoryfunc,"directory_msg")
    pub.subscribe(self.s1_nfragsfunc,"nfrags_spec1")
    pub.subscribe(self.s2_nfragsfunc,"nfrags_spec2")
    pub.subscribe(self.s3_nfragsfunc,"nfrags_spec3")
    
    infostr ="Please complete the Molecule Files tab to continue."
    self.infolbl = wx.StaticText(self,label=infostr)
    grid4b.Add(self.infolbl,pos=(1,1),span=(1,6))
    infostr2 = "Select the fragment files for each species below.  Please note that species with greater\nthan 15 fragments are not currently supported by the GUI and must be entered in to the input file manually."
    self.info2lbl = wx.StaticText(self,label=infostr2)
    grid4b.Add(self.info2lbl,pos=(3,1),span=(1,6))
    self.info2lbl.Hide()

    self.specieslbl2 = []
    self.buttonselect1 = []
    self.buttonselect2 = []
    self.buttonselect3 = []
    self.viewfrags1 = []
    self.viewfrags2 = []
    self.viewfrags3 = []
    
    for i in range(3):
      self.specieslbl2.append(wx.StaticText(self,label=("Species %d" %(i+1))))
      grid4b.Add(self.specieslbl2[i],pos=(5,(2*i)+1),span=(1,2),flag=wx.ALIGN_CENTER)
      self.specieslbl2[i].Hide()
    for i in range(15):
      self.buttonselect1.append(wx.Button(self,label="Select Fragment File",name=("sp1f %d" %(i+1))))
      self.buttonselect2.append(wx.Button(self,label="Select Fragment File",name=("sp2f %d" %(i+1))))
      self.buttonselect3.append(wx.Button(self,label="Select Fragment File",name=("sp3f %d" %(i+1))))
      grid4b.Add(self.buttonselect1[i],pos=(6+i,1))
      grid4b.Add(self.buttonselect2[i],pos=(6+i,3))
      grid4b.Add(self.buttonselect3[i],pos=(6+i,5))
      self.Bind(wx.EVT_BUTTON,self.buttonfunc4b,self.buttonselect1[i])
      self.Bind(wx.EVT_BUTTON,self.buttonfunc4b,self.buttonselect2[i])
      self.Bind(wx.EVT_BUTTON,self.buttonfunc4b,self.buttonselect3[i])
      self.viewfrags1.append(wx.TextCtrl(self,value="",name=("v1f %d" %(i+1)),style=wx.TE_READONLY))
      self.viewfrags2.append(wx.TextCtrl(self,value="",name=("v2f %d" %(i+1)),style=wx.TE_READONLY))
      self.viewfrags3.append(wx.TextCtrl(self,value="",name=("v3f %d" %(i+1)),style=wx.TE_READONLY))
      grid4b.Add(self.viewfrags1[i],pos=(6+i,2))
      grid4b.Add(self.viewfrags2[i],pos=(6+i,4))
      grid4b.Add(self.viewfrags3[i],pos=(6+i,6))
      self.buttonselect1[i].Hide(),self.buttonselect2[i].Hide(),self.buttonselect3[i].Hide()
      self.viewfrags1[i].Hide(),self.viewfrags2[i].Hide(),self.viewfrags3[i].Hide()

  def buttonfunc4b(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    dlg = wx.FileDialog(self,"Select Fragment Data",p4b.simulation_dirname,"","*.dat",wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
      self.dirname = dlg.GetDirectory()
      self.filepathdata = dlg.GetPath()
      file_data = os.path.relpath(self.filepathdata,p4b.simulation_dirname)
      for i in range(15):
        if widgetname == ("sp1f %d" %(i+1)):
          p4b.s1_frags[i] = str(file_data)
          self.viewfrags1[i].SetValue(str(file_data))
        elif widgetname == ("sp2f %d" %(i+1)):
          p4b.s2_frags[i] = str(file_data)
          self.viewfrags2[i].SetValue(str(file_data))
        elif widgetname== ("sp3f %d" %(i+1)):
          p4b.s3_frags[i] = str(file_data)
          self.viewfrags3[i].SetValue(str(file_data))
    self.Layout()

  def directoryfunc(self,message,arg2=None):
    p4b.simulation_dirname = message

  def s1_nfragsfunc(self,message,arg2=None):
    p4b.s1_nfrags = message
    self.infolbl.Hide(), self.info2lbl.Show()
    self.specieslbl2[0].Show()
    for i in range(15):
      if i < int(message):
        self.buttonselect1[i].Show(),self.viewfrags1[i].Show()
      if i>=int(message):
        self.buttonselect1[i].Hide(),self.viewfrags1[i].Hide()
    self.Layout()

  def s2_nfragsfunc(self,message,arg2=None):
    p4b.s2_nfrags = message
    if int(message) != 0:
      self.specieslbl2[1].Show()
    for i in range(15):
      if i<int(message):
        self.buttonselect2[i].Show(),self.viewfrags2[i].Show()
      if i>=int(message):
        self.buttonselect2[i].Hide(),self.viewfrags2[i].Hide()
    self.Layout()

  def s3_nfragsfunc(self,message,arg2=None):
    p4b.s3_nfrags = message
    if int(message)!=0:
      self.specieslbl2[2].Show()
    for i in range(15):
      if i<int(message):
        self.buttonselect3[i].Show(),self.viewfrags3[i].Show()
      if i>=int(message):
        self.buttonselect3[i].Hide(),self.viewfrags3[i].Hide()
    self.Layout()


#Fragment Files 2
class p4c(wx.Panel):
  s4_nfrags = ""
  s5_nfrags = ""
  s6_nfrags = ""
  s4_frags = ["","","","","","","","","","","","","","",""]
  s5_frags = ["","","","","","","","","","","","","","",""]
  s6_frags = ["","","","","","","","","","","","","","",""]
  simulation_dirname = ""
  nbr_spec_received = "1"
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid4c = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid4c)

    pub.subscribe(self.directoryfunc,"directory_msg")
    pub.subscribe(self.s4_nfragsfunc,"nfrags_spec4")
    pub.subscribe(self.s5_nfragsfunc,"nfrags_spec5")
    pub.subscribe(self.s6_nfragsfunc,"nfrags_spec6")
    
    infostr = "Please ignore this tab if there are less than 4 species in simulation.\nElse, provide fragment library information for species 4-6 below.  The Molecule Files tab must be completed first."
    self.infolbl = wx.StaticText(self,label=infostr)
    grid4c.Add(self.infolbl,pos=(1,1),span=(2,5))
    self.info2lbl = wx.StaticText(self,label="Species with more than 15 fragments are not currently supported by the GUI and must be entered manually to the input file.")
    grid4c.Add(self.info2lbl,pos=(3,1),span=(1,7))
    self.specieslbl2 = []
    self.buttonselect4 = []
    self.buttonselect5 = []
    self.buttonselect6 = []
    self.viewfrags4 = []
    self.viewfrags5 = []
    self.viewfrags6 = []
  
    for i in range(3):
      self.specieslbl2.append(wx.StaticText(self,label=("Species %d" %(i+4))))
      grid4c.Add(self.specieslbl2[i],pos=(5,(2*i)+1),span=(1,2),flag=wx.ALIGN_CENTER)
      self.specieslbl2[i].Hide()
    for i in range(15):
      self.buttonselect4.append(wx.Button(self,label="Select Fragment File",name=("sp4f %d" %(i+1))))
      self.buttonselect5.append(wx.Button(self,label="Select Fragment File",name=("sp5f %d" %(i+1))))
      self.buttonselect6.append(wx.Button(self,label="Select Fragment File",name=("sp6f %d" %(i+1))))
      grid4c.Add(self.buttonselect4[i],pos=(6+i,1))
      grid4c.Add(self.buttonselect5[i],pos=(6+i,3))
      grid4c.Add(self.buttonselect6[i],pos=(6+i,5))
      self.Bind(wx.EVT_BUTTON,self.buttonfunc4c,self.buttonselect4[i])
      self.Bind(wx.EVT_BUTTON,self.buttonfunc4c,self.buttonselect5[i])
      self.Bind(wx.EVT_BUTTON,self.buttonfunc4c,self.buttonselect6[i])
      self.viewfrags4.append(wx.TextCtrl(self,value="",name=("v4f %d" %(i+1)),style=wx.TE_READONLY))
      self.viewfrags5.append(wx.TextCtrl(self,value="",name=("v5f %d" %(i+1)),style=wx.TE_READONLY))
      self.viewfrags6.append(wx.TextCtrl(self,value="",name=("v6f %d" %(i+1)),style=wx.TE_READONLY))
      grid4c.Add(self.viewfrags4[i],pos=(6+i,2))
      grid4c.Add(self.viewfrags5[i],pos=(6+i,4))
      grid4c.Add(self.viewfrags6[i],pos=(6+i,6))
      self.buttonselect4[i].Hide(),self.buttonselect5[i].Hide(),self.buttonselect6[i].Hide()
      self.viewfrags4[i].Hide(),self.viewfrags5[i].Hide(),self.viewfrags6[i].Hide()

  def buttonfunc4c(self,event):
    event_object = event.GetEventObject()
    widgetname = event_object.GetName()
    dlg = wx.FileDialog(self,"Select Fragment Data",p4c.simulation_dirname,"","*.dat",wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
      self.dirname = dlg.GetDirectory()
      self.filepathdata = dlg.GetPath()
      file_data = os.path.relpath(self.filepathdata,p4c.simulation_dirname)
      for i in range(15):
        if widgetname == ("sp4f %d" %(i+1)):
          p4c.s4_frags[i] = str(file_data)
          self.viewfrags4[i].SetValue(str(file_data))
        elif widgetname == ("sp5f %d" %(i+1)):
          p4c.s5_frags[i] = str(file_data)
          self.viewfrags5[i].SetValue(str(file_data))
        elif widgetname== ("sp6f %d" %(i+1)):
          p4c.s6_frags[i] = str(file_data)
          self.viewfrags6[i].SetValue(str(file_data))
    self.Layout()

  def directoryfunc(self,message,arg2=None):
    p4c.simulation_dirname = message

  def s4_nfragsfunc(self,message,arg2=None):
    p4c.s4_nfrags = message
    self.infolbl.Hide(), self.info2lbl.Show()
    self.specieslbl2[0].Show()
    for i in range(15):
      if i < int(message):
        self.buttonselect4[i].Show(),self.viewfrags4[i].Show()
      if i>=int(message):
        self.buttonselect4[i].Hide(),self.viewfrags4[i].Hide()
    self.Layout()

  def s5_nfragsfunc(self,message,arg2=None):
    p4c.s5_nfrags = message
    if int(message) != 0:
      self.specieslbl2[1].Show()
    for i in range(15):
      if i<int(message):
        self.buttonselect5[i].Show(),self.viewfrags5[i].Show()
      if i>=int(message):
        self.buttonselect5[i].Hide(),self.viewfrags5[i].Hide()
    self.Layout()

  def s6_nfragsfunc(self,message,arg2=None):
    p4c.s6_nfrags = message
    if int(message)!=0:
      self.specieslbl2[2].Show()
    for i in range(15):
      if i<int(message):
        self.buttonselect6[i].Show(),self.viewfrags6[i].Show()
      if i>=int(message):
        self.buttonselect6[i].Hide(),self.viewfrags6[i].Hide()
    self.Layout()

#Input File - start type, run type
class p4d(wx.Panel):

  nbr_species_received = "1"
  ensemble_received = ""
  simulation_dirname = ""

  start_type = ""
  mkconfig_dat = ["","","","","","","","","","","",""]
  box1readold = ""
  box2readold = ""
  chkpoint = ""
  run_type = ["","",""]
  
  def __init__(self,parent):
    self.dirname = ""
    wx.Panel.__init__(self,parent)
    grid4d = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid4d)
    
    pub.subscribe(self.speciesfunc,"species_msg")
    pub.subscribe(self.ensemblefunc,"ensemble_msg")
    pub.subscribe(self.directoryfunc,"directory_msg")
    
    descript_string = "Provide information on the Start Type and Run Type for the simulation below."
    self.descript_stringlbl = wx.StaticText(self,label=descript_string)
    grid4d.Add(self.descript_stringlbl,pos=(1,1),span=(1,4))

    self.required_alert = wx.StaticText(self,label="Please select a simulation directory on the page Basic Information > Page 1. This is required to continue.")
    grid4d.Add(self.required_alert,pos=(2,1),span=(1,4))

    self.starttypelbl = wx.StaticText(self,label="Start Type:")
    grid4d.Add(self.starttypelbl,pos=(3,1))
    self.starttypelbl.Hide()
    
    choices_start = ["","make_config","checkpoint","read_old"]
    self.start_choice = wx.Choice(self,-1,choices=choices_start,name="start_choice")
    grid4d.Add(self.start_choice,pos=(3,2))
    self.Bind(wx.EVT_CHOICE,self.choicefunc4d,self.start_choice)
    self.start_choice.Hide()
    
    self.tw_start = []
    self.tw_speclbl = []
    for i in range(12):
      self.tw_start.append(wx.TextCtrl(self,value="",name=("twstart %d" %(i))))
      if i<6:
        grid4d.Add(self.tw_start[i],pos=(4+i,5))
      else: grid4d.Add(self.tw_start[i],pos=(i-2,6))
      self.tw_start[i].Hide()
      self.Bind(wx.EVT_TEXT,self.textfunc4d,self.tw_start[i])
    for i in range(6):
      self.tw_speclbl.append(wx.StaticText(self,label=("# Molecules Species %d:" %(i+1))))
      grid4d.Add(self.tw_speclbl[i],pos=(4+i,4),flag=wx.ALIGN_RIGHT)
      self.tw_speclbl[i].Hide()
    self.box1lbl = wx.StaticText(self,label="Box 1")
    self.box2lbl = wx.StaticText(self,label="Box 2")
    grid4d.Add(self.box1lbl,pos=(3,5),flag=wx.ALIGN_CENTER)
    grid4d.Add(self.box2lbl,pos=(3,6),flag=wx.ALIGN_CENTER)
    self.box1lbl.Hide(),self.box2lbl.Hide()
  
    self.chk_button = wx.Button(self,label="Select Checkpoint File",name="chk_btn")
    self.chk_txt = wx.TextCtrl(self,value="",name="chk_txt",style=wx.TE_READONLY)
    grid4d.Add(self.chk_button,pos=(4,1),span=(1,2),flag=wx.EXPAND)
    grid4d.Add(self.chk_txt,pos=(4,3))
    self.chk_button.Hide(),self.chk_txt.Hide()
    self.Bind(wx.EVT_BUTTON,self.buttonfunc4d,self.chk_button)
  
    self.readoldb1 = wx.Button(self,label="Select Box 1 Read Old",name="b1_ro")
    self.readoldb2 = wx.Button(self,label="Select Box 2 Read Old",name="b2_ro")
    grid4d.Add(self.readoldb1,pos=(5,1),span=(1,2),flag=wx.EXPAND)
    grid4d.Add(self.readoldb2,pos=(6,1),span=(1,2),flag=wx.EXPAND)
    self.chkrob1 = wx.TextCtrl(self,value="",name="rob1_txt",style=wx.TE_READONLY)
    self.chkrob2 = wx.TextCtrl(self,value="",name="rob2_txt",style=wx.TE_READONLY)
    grid4d.Add(self.chkrob1,pos=(5,3))
    grid4d.Add(self.chkrob2,pos=(6,3))
    self.Bind(wx.EVT_BUTTON,self.buttonfunc4d,self.readoldb1)
    self.Bind(wx.EVT_BUTTON,self.buttonfunc4d,self.readoldb2)
    self.readoldb1.Hide(),self.readoldb2.Hide(),self.chkrob1.Hide(),self.chkrob2.Hide()
  
    #### Run Type
    runtype_choice = ["","Equilibration","Production"]
    self.rt1lbl = wx.StaticText(self,label="Run Type:")
    grid4d.Add(self.rt1lbl,pos=(11,1))
    self.rt2lbl = wx.StaticText(self,label="Frequency at which Acceptance Ratios for moves output to Log File (MC Steps):")
    grid4d.Add(self.rt2lbl,pos=(12,1),span=(1,3))
    self.rt3lbl = wx.StaticText(self,label="Moves Before Volume Diplacement Output to Log File (MC Steps):")
    grid4d.Add(self.rt3lbl,pos=(13,1),span=(1,3))
    
    self.rt1_choose = wx.Choice(self,-1,choices=runtype_choice,name="runtype_choice")
    self.rt2_text = wx.TextCtrl(self,value="",name="rt_mdw")
    self.rt3_text = wx.TextCtrl(self,value="",name="rt_mvd")
    grid4d.Add(self.rt1_choose,pos=(11,2))
    grid4d.Add(self.rt2_text,pos=(12,4))
    grid4d.Add(self.rt3_text,pos=(13,4))
    self.Bind(wx.EVT_CHOICE,self.choicefunc4d,self.rt1_choose)
    self.Bind(wx.EVT_TEXT,self.textfunc4d,self.rt2_text)
    self.Bind(wx.EVT_TEXT,self.textfunc4d,self.rt3_text)

  def choicefunc4d(self,event):
    event_object=event.GetEventObject()
    widgetname = event_object.GetName()
    choice_data = event.GetString()
    choice_index = event.GetSelection()
    if widgetname == "start_choice":
      if choice_data:
        p4d.start_type = choice_data
        if choice_data == "make_config":
          self.box1lbl.Show()
          self.chk_button.Hide(),self.chk_txt.Hide()
          for i in range(int(p4d.nbr_species_received)):
            self.tw_start[i].Show()
            self.tw_speclbl[i].Show()
          if p4d.ensemble_received == "GEMC":
            self.box2lbl.Show()
            for i in range(int(p4d.nbr_species_received)):
              self.tw_start[i+6].Show()
          else: self.box2lbl.Hide()
          self.readoldb1.Hide(),self.readoldb2.Hide(),self.chkrob1.Hide(),self.chkrob2.Hide()
        elif choice_data == "checkpoint":
          self.chk_button.Show(),self.chk_txt.Show()
          for item in self.tw_start:
            item.Hide()
          for item in self.tw_speclbl:
            item.Hide()
          self.box1lbl.Hide(),self.box2lbl.Hide()
          self.readoldb1.Hide(),self.readoldb2.Hide(),self.chkrob1.Hide(),self.chkrob2.Hide()
        elif choice_data == "read_old":
          self.chk_button.Hide(),self.chk_txt.Hide()
          for item in self.tw_start:
            item.Hide()
          for item in self.tw_speclbl:
            item.Hide()
            self.box1lbl.Hide(),self.box2lbl.Hide()
          self.readoldb1.Show(),self.chkrob1.Show()
          if p4d.ensemble_received == "GEMC":
            self.readoldb2.Show(),self.chkrob2.Show()
      else:
        self.box1lbl.Hide(),self.box2lbl.Hide()
        self.readoldb1.Hide(),self.readoldb2.Hide(),self.chkrob1.Hide(),self.chkrob2.Hide()
        self.chk_button.Hide(),self.chk_txt.Hide()
        for item in self.tw_start:
          item.Hide()
        for item in self.tw_speclbl:
          item.Hide()
    elif widgetname == "runtype_choice":
      p4d.run_type[0] = choice_data
    self.Layout()
      
  def buttonfunc4d(self,event):
    event_object=event.GetEventObject()
    widgetname = event_object.GetName()
    if widgetname == "chk_btn":
      dlg = wx.FileDialog(self,"Select Checkpoint",p4d.simulation_dirname,"","*.chk",wx.OPEN)
      if dlg.ShowModal() == wx.ID_OK:
        self.filepathdata = dlg.GetPath()
        file_data = os.path.relpath(self.filepathdata,p4d.simulation_dirname)
        self.chk_txt.SetValue(str(file_data))
        p4d.chkpoint = str(file_data)
    elif widgetname == "b1_ro":
      dlg = wx.FileDialog(self,"Select Box 1 Read Old File",p4d.simulation_dirname,"","*.xyz",wx.OPEN)
      if dlg.ShowModal() == wx.ID_OK:
        self.filepathdata = dlg.GetPath()
        file_data = os.path.relpath(self.filepathdata,p4d.simulation_dirname)
        self.chkrob1.SetValue(str(file_data))
        p4d.box1readold = str(file_data)
    elif widgetname == "b2_ro":
      dlg = wx.FileDialog(self,"Select Box 2 Read Old File",p4d.simulation_dirname,"","*.xyz",wx.OPEN)
      if dlg.ShowModal() == wx.ID_OK:
        self.filepathdata = dlg.GetPath()
        file_data = os.path.relpath(self.filepathdata,p4d.simulation_dirname)
        self.chkrob2.SetValue(str(file_data))
        p4d.box2readold = str(file_data)
    self.Layout()

  def directoryfunc(self,message,arg2=None):
    self.starttypelbl.Show()
    self.start_choice.Show()
    self.required_alert.Hide()
    p4d.simulation_dirname = message
    self.Layout()

  def ensemblefunc(self,message,arg2=None):
    p4d.ensemble_received = message
    self.Layout()
  
  def speciesfunc(self,message,arg2=None):
    p4d.nbr_species_received = message
    self.Layout()

  def textfunc4d(self,event):
    event_object=event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    if widgetname == "rt_mdw":
      p4d.run_type[1] = text_data
    elif widgetname == "rt_mvd":
      p4d.run_type[2] = text_data
    for i in range(12):
      if widgetname == ("twstart %d" %i):
        p4d.mkconfig_dat[i] = text_data
    self.Layout()

# Output File - frequency, average, property output
class p4e(wx.Panel):
  freq_info = ["","","",""]
  average_dat = "1"
  freq_keywords = ["","",""]
  prop_b1 = ["","","","","","","","",""]
  prop_b2 = ["","","","","","","","",""]
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid4e = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid4e)
    
    pub.subscribe(self.ensemblefunc,"ensemble_msg")
  # frequency info requires one choice ctrl for first option followed by changing lbls and text controls
  # create required widgets and bind them
    self.freqtypelbl = wx.StaticText(self,label="Frequency Type:")
    grid4e.Add(self.freqtypelbl,pos=(1,1))
    freqtypechoices = ["","Timed","none"]
    self.freqtypechoice = wx.Choice(self,-1,choices=freqtypechoices,name ="freqtype_choice")
    grid4e.Add(self.freqtypechoice,pos=(1,2))
    self.Bind(wx.EVT_CHOICE,self.choicefunc4e,self.freqtypechoice)
    
    self.freqtxt = []
    for i in range(3):
      self.freqtxt.append(wx.TextCtrl(self,value="",name=("freqtxt %d" %i)))
      grid4e.Add(self.freqtxt[i],pos=(i+2,3))
      self.Bind(wx.EVT_TEXT,self.textfunc4e,self.freqtxt[i])
      self.freqtxt[i].Hide()
    self.freqlbl1 = wx.StaticText(self,label="Thermodynamic Output Frequency (minutes):")
    self.freqlbl2 = wx.StaticText(self,label="Coordinate Output Frequency (minutes):")
    self.freqlbl3 = wx.StaticText(self,label="Simulation Run Time (minutes):")
    grid4e.Add(self.freqlbl1,pos=(2,1),span=(1,2))
    grid4e.Add(self.freqlbl2,pos=(3,1),span=(1,2))
    grid4e.Add(self.freqlbl3,pos=(4,1),span=(1,2))
    self.freqlbl1.Hide(),self.freqlbl2.Hide(),self.freqlbl3.Hide()
    
  # average is set to 1
    self.avglbl = wx.StaticText(self,label="Average Information:")
    self.avgtxt = wx.TextCtrl(self,value="1",name="avg_text",style=wx.TE_READONLY)
    grid4e.Add(self.avglbl,pos=(6,1))
    grid4e.Add(self.avgtxt,pos=(6,2))
    
  # property output is a list of checkboxes for box 1, box2,
    self.proplbl = wx.StaticText(self,label="Desired Property Outputs")
    grid4e.Add(self.proplbl,pos=(8,1),span=(1,5),flag=wx.ALIGN_CENTER)
    self.propb1lbl = wx.StaticText(self,label="Box 1")
    self.propb2lbl = wx.StaticText(self,label="Box 2")
    grid4e.Add(self.propb1lbl,pos=(9,2))
    grid4e.Add(self.propb2lbl,pos=(9,4))

    self.prop_options = ["Energy_Total","Energy_LJ","Energy_Elec","Energy_Intra","Enthalpy","Pressure","Volume","Nmols","Density"]
    self.propb1_chks = []
    self.propb2_chks = []
    for i in range(9):
      self.propb1_chks.append(wx.CheckBox(self,-1,self.prop_options[i],name=("b1 %d" %i)))
      self.propb2_chks.append(wx.CheckBox(self,-1,self.prop_options[i],name=("b2 %d" %i)))
      grid4e.Add(self.propb1_chks[i],pos=(10+i,2))
      grid4e.Add(self.propb2_chks[i],pos=(10+i,4))
      self.Bind(wx.EVT_CHECKBOX,self.checkfunc4e,self.propb1_chks[i])
      self.Bind(wx.EVT_CHECKBOX,self.checkfunc4e,self.propb2_chks[i])
      for item in self.propb2_chks: #show only if GEMC selected
        item.Hide()
      self.propb2lbl.Hide()
      
  def ensemblefunc(self,message,arg2=None):
    msg = message
    if msg == "GEMC":
      for item in self.propb2_chks:
        item.Show()
      self.propb2lbl.Show()
    else:
      for item in self.propb2_chks:
        item.Hide()
      self.propb2lbl.Hide()
    self.Layout()
    
  def checkfunc4e(self,event):
    event_object=event.GetEventObject()
    widgetname = event_object.GetName()
    chk_data = event.IsChecked()
    for i in range(9):
      if widgetname == ("b1 %d" %i) and chk_data == True:
        p4e.prop_b1[i] = self.prop_options[i]
      elif widgetname == ("b1 %d" %i) and chk_data == False:
        p4e.prop_b1[i] = ""
      elif widgetname == ("b2 %d" %i) and chk_data == True:
        p4e.prop_b2[i] = self.prop_options[i]
      elif widgetname == ("b2 %d" %i) and chk_data == False:
        p4e.prop_b2[i] = ""

  def choicefunc4e(self,event):
    event_object=event.GetEventObject()
    widgetname = event_object.GetName()
    choice_data = event.GetString()
    choice_index = event.GetSelection()
    if widgetname == "freqtype_choice":
      if choice_index == 0:
        self.freqlbl1.Hide(),self.freqlbl2.Hide(),self.freqlbl3.Hide()
        for i in range(3):
          self.freqtxt[i].Hide()
          self.freqtxt[i].SetValue("")
        for i in range(4):
          p4e.freq_info[i] = ""
        p4e.freq_keywords[0]=""
        p4e.freq_keywords[1] = ""
        p4e.freq_keywords[2] = ""
      elif choice_index == 1:
        self.freqlbl1.Show(),self.freqlbl2.Show(),self.freqlbl3.Show()
        for i in range(3):
          self.freqtxt[i].Show()
        self.freqlbl1.SetLabel("Thermodynamic Output Frequency (minutes):")
        self.freqlbl2.SetLabel("Coordinate Output Frequency (minutes):")
        self.freqlbl3.SetLabel("Simulation Run Time (minutes):")
        p4e.freq_info[0] = choice_data
        p4e.freq_keywords[0] = "thermofreq"
        p4e.freq_keywords[1] = "coordfreq"
        p4e.freq_keywords[2] = "Stop"
      elif choice_index == 2:
        self.freqlbl1.Show(),self.freqlbl2.Show(),self.freqlbl3.Show()
        self.freqlbl1.SetLabel("Thermodynamic Output Frequency (MC Steps):")
        self.freqlbl2.SetLabel("Coordinate Output Frequency (MC Steps):")
        self.freqlbl3.SetLabel("Simulation Run Time (MC Steps):")
        for i in range(3):
          self.freqtxt[i].Show()
        p4e.freq_info[0] = choice_data
        p4e.freq_keywords[0] = "Nthermofreq"
        p4e.freq_keywords[1] = "Ncoordfreq"
        p4e.freq_keywords[2] = "MCsteps"

    self.Layout()

  def textfunc4e(self,event):
    event_object=event.GetEventObject()
    widgetname = event_object.GetName()
    text_data = event.GetString()
    for i in range(3):
      if widgetname == ("freqtxt %d" %i):
        p4e.freq_info[i+1] = text_data

# PanelFour is the parent panel for the file-handling tabs
class PanelFour(wx.Panel):
  # Declare the class variables

  # initialize the panel
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

  # define the notebook
    nest = wx.Notebook(self)
    nest1 = p4a(nest)
    nest2 = p4b(nest)
    nest3 = p4c(nest)
    nest4 = p4d(nest)
    nest5 = p4e(nest)

    nest.AddPage(nest1,"Molecule Files")
    nest.AddPage(nest2,"Fragment Files 1")
    nest.AddPage(nest3,"Fragment Files 2")
    nest.AddPage(nest4,"Input File")
    nest.AddPage(nest5,"Output File")
    nested_sizer = wx.BoxSizer()
    nested_sizer.Add(nest,1,wx.EXPAND)
    self.SetSizer(nested_sizer)

###################################################################################
###################################################################################


class MainFrame(wx.Frame):
  def __init__(self):
    wx.Frame.__init__(self,None,title="Cassandra Input File Editor v1.1")
    self.SetInitialSize((900,620))

    # Create the panel hosting the notebook
    panel_n = wx.Panel(self)
    notebook_m = wx.Notebook(panel_n)

    # Add PanelOne:PanelFour to notebook as pages
    page1 = PanelOne(notebook_m)
    page2 = PanelTwo(notebook_m)
    page3 = PanelThree(notebook_m)
    page4 = PanelFour(notebook_m)

    notebook_m.AddPage(page1, "Basic Information")
    notebook_m.AddPage(page2, "Interaction Parameters")
    notebook_m.AddPage(page3, "Probabilities Information")
    notebook_m.AddPage(page4, "File Handling")

    # Create a unique sizer for the notebook
    NBSizer = wx.BoxSizer()
    NBSizer.Add(notebook_m,1,wx.EXPAND)
    panel_n.SetSizer(NBSizer)

    # Create the menu to create the save file
    menubar = wx.MenuBar()
    file_menu = wx.Menu()
    save_select = file_menu.Append(wx.ID_ANY,"Save As...")
    self.Bind(wx.EVT_MENU,self.savefunc,save_select)
    menubar.Append(file_menu,"&File")
    self.SetMenuBar(menubar)

  def savefunc(self,event):
    s1_nf = int(p4b.s1_nfrags)
    s2_nf, s3_nf, s4_nf, s5_nf, s6_nf = 0,0,0,0,0
    if int(p1a.nbrspecies_store)>1:
      s2_nf = int(p4b.s2_nfrags)
    if int(p1a.nbrspecies_store)>2:
      s3_nf = int(p4b.s3_nfrags)
    if int(p1a.nbrspecies_store)>3:
      s4_nf = int(p4c.s4_nfrags)
    if int(p1a.nbrspecies_store)>4:
      s5_nf = int(p4c.s5_nfrags)
    if int(p1a.nbrspecies_store)>5:
      s6_nf = int(p4c.s6_nfrags)
    string1 = str(p1a.runname_store)+".inp"
    newsection = "!------------------"
    self.dirname = ''
    dialog = wx.FileDialog(self,"Save File As",p1a.simulation_dirname,string1,"*.*",wx.SAVE)
    if dialog.ShowModal() == wx.ID_OK:
      path = dialog.GetPath()
      f = file(path,'w')
      f.write("! This is the input file for %s, a %s-species %s simulation at %s K \n\n\n" %(p1a.runname_store,p1a.nbrspecies_store,p1a.ensembletype_store,p1b.temperature_store))
      f.write("# Run_Name \n%s \n%s \n\n" %(p1a.runname_store,newsection))
      f.write("# Sim_type \n%s \n%s \n\n" %(p1a.ensembletype_store,newsection))
      f.write("# Nbr_Species \n%s \n%s \n\n" %(p1a.nbrspecies_store,newsection))
      f.write("# Rcutoff_Low \n%s \n%s\n" %(p1b.min_cutoffstore,newsection))
      f.write("\n# Temperature_Info\n")
      if p1a.nbr_boxes == 1:
        f.write("%s\n%s\n\n" %(p1b.temperature_store,newsection))
      elif p1a.nbr_boxes == 2:
        f.write("%s %s\n%s\n\n" %(p1b.temperature_store,p1b.temperature_store,newsection))
      if p1a.ensembletype_store == "NPT_MC":
        f.write("# Pressure_Info\n%s\n%s \n\n" %(p1b.pressure_store,newsection))
      elif p1a.ensembletype_store == "GEMC" and int(p1a.nbrspecies_store)>1:
        f.write("# Pressure_Info\n%s %s\n%s\n\n" %(p1b.pressure_store,p1b.pressure_store,newsection))
      f.write("# Seed_Info\n%s %s\n%s\n\n" %(p1b.seed1_store,p1b.seed2_store,newsection))
      f.write("# Mixing_Rule \n%s \n%s \n\n" %(p1b.mixingruletype_store,newsection))
      f.write("# Charge_Style \n")
      f.write("%s %s %s %s \n" %(p2a.box1_charge[0],p2a.box1_charge[1],p2a.box1_charge[2],p2a.box1_charge[3]))
      if p1a.nbr_boxes >1:
        f.write("%s %s %s %s \n" %(p2a.box2_charge[0],p2a.box2_charge[1],p2a.box2_charge[2],p2a.box2_charge[3]))
      f.write("%s \n\n# VDW_Style\n" %(newsection))
      f.write("%s %s %s %s %s\n" %(p2a.box1_vdw[0],p2a.box1_vdw[1],p2a.box1_vdw[2],p2a.box1_vdw[3],p2a.box1_vdw[4]))
      if p1a.nbr_boxes >1:
        f.write("%s %s %s %s %s\n" %(p2a.box2_vdw[0],p2a.box2_vdw[1],p2a.box2_vdw[2],p2a.box2_vdw[3],p2a.box2_vdw[4]))
      f.write("%s \n\n" %(newsection))
      f.write("# Intra_Scaling \n%s %s %s %s\n%s %s %s %s\n" %(p2b.sp1vdw[0],p2b.sp1vdw[1],p2b.sp1vdw[2],p2b.sp1vdw[3],p2b.sp1c[0],p2b.sp1c[1],p2b.sp1c[2],p2b.sp1c[3]))
      if int(p1a.nbrspecies_store) > 1:
        f.write("%s %s %s %s \n%s %s %s %s\n" %(p2b.sp2vdw[0],p2b.sp2vdw[1],p2b.sp2vdw[2],p2b.sp2vdw[3],p2b.sp2c[0],p2b.sp2c[1],p2b.sp2c[2],p2b.sp2c[3]))
      if int(p1a.nbrspecies_store) > 2:
        f.write("%s %s %s %s \n%s %s %s %s\n" %(p2b.sp3vdw[0],p2b.sp3vdw[1],p2b.sp3vdw[2],p2b.sp3vdw[3],p2b.sp3c[0],p2b.sp3c[1],p2b.sp3c[2],p2b.sp3c[3]))
      if int(p1a.nbrspecies_store) > 3:
        f.write("%s %s %s %s \n%s %s %s %s\n" %(p2b.sp4vdw[0],p2b.sp4vdw[1],p2b.sp4vdw[2],p2b.sp4vdw[3],p2b.sp4c[0],p2b.sp4c[1],p2b.sp4c[2],p2b.sp4c[3]))
      if int(p1a.nbrspecies_store) > 4:
        f.write("%s %s %s %s \n%s %s %s %s\n" %(p2b.sp5vdw[0],p2b.sp5vdw[1],p2b.sp5vdw[2],p2b.sp5vdw[3],p2b.sp5c[0],p2b.sp5c[1],p2b.sp5c[2],p2b.sp5c[3]))
      if int(p1a.nbrspecies_store) > 5:
        f.write("%s %s %s %s \n%s %s %s %s\n" %(p2b.sp6vdw[0],p2b.sp6vdw[1],p2b.sp6vdw[2],p2b.sp6vdw[3],p2b.sp6c[0],p2b.sp6c[1],p2b.sp6c[2],p2b.sp6c[3]))
      f.write("%s \n\n" %(newsection))
      f.write("# Box_Info \n%s\n%s\n%s %s %s\n" %(p1a.nbr_boxes,p1b.box1_shape, p1b.box1_length,p1b.box1_length, p1b.box1_length))
      if p1a.nbr_boxes > 1:
        f.write("\n%s\n%s %s %s\n" %(p1b.box2_shape,p1b.box2_length,p1b.box2_length,p1b.box2_length))
      f.write("%s \n\n" %(newsection))
      f.write("# Move_Probability_Info\n")
      if p3a.prob_trans !="0":
        f.write("\n# Prob_Translation\n%s\n%s %s %s %s %s %s\n" %(p3a.prob_trans,p3a.ptb1s[0],p3a.ptb1s[1],p3a.ptb1s[2],p3a.ptb1s[3],p3a.ptb1s[4],p3a.ptb1s[5]))
        if p1a.nbr_boxes>1:
          f.write("%s %s %s %s %s %s\n" %(p3a.ptb2s[0],p3a.ptb2s[1],p3a.ptb2s[2],p3a.ptb2s[3],p3a.ptb2s[4],p3a.ptb2s[5]))
      if p3b.prob_rot !="0":
        f.write("\n# Prob_Rotation\n%s\n%s %s %s %s %s %s\n" %(p3b.prob_rot,p3b.prb1s[0],p3b.prb1s[1],p3b.prb1s[2],p3b.prb1s[3],p3b.prb1s[4],p3b.prb1s[5]))
        if p1a.nbr_boxes>1:
          f.write("%s %s %s %s %s %s\n" %(p3b.prb2s[0],p3b.prb2s[1],p3b.prb2s[2],p3b.prb2s[3],p3b.prb2s[4],p3b.prb2s[5]))
      if p3c.prob_regrowth !="0":
        f.write("\n# Prob_Regrowth\n%s\n%s %s %s %s %s %s\n" %(p3c.prob_regrowth,p3c.pregrow_s[0],p3c.pregrow_s[1],p3c.pregrow_s[2],p3c.pregrow_s[3],p3c.pregrow_s[4],p3c.pregrow_s[5]))
      if p3d.prob_vol !="0":
        f.write("\n# Prob_Volume\n%s\n%s\n%s\n" %(p3d.prob_vol,p3d.pvol_box[0],p3d.pvol_box[1]))
      if p3e.prob_insert !="0":
        f.write("\n# Prob_Insertion\n%s\n%s\n%s\n"%(p3e.prob_insert,"insertion method",p3e.pinsert_spec[0]))
        if int(p1a.nbrspecies_store)>1:
          f.write("%s\n%s\n" %("insertion method",p3e.pinsert_spec[1]))
        if int(p1a.nbrspecies_store)>2:
          f.write("%s\n%s\n" %("insertion method",p3e.pinsert_spec[2]))
        if int(p1a.nbrspecies_store)>3:
          f.write("%s\n%s\n" %("insertion method",p3e.pinsert_spec[3]))
        if int(p1a.nbrspecies_store)>4:
          f.write("%s\n%s\n" %("insertion method",p3e.pinsert_spec[4]))
        if int(p1a.nbrspecies_store)>5:
          f.write("%s\n%s\n" %("insertion method",p3e.pinsert_spec[5]))
      if p3e.prob_insert !="0":
        f.write("\n# Prob_Deletion\n%s\n" %(p3e.prob_insert))
      if p3f.prob_swap !="0":
        f.write("\n# Prob_Swap\n%s\n%s\n%s" %(p3f.prob_swap,"insertion method",p3f.pswap_spec[0]))
        if int(p1a.nbrspecies_store)>1:
          f.write("\n%s\n%s\n" %("insertion method",p3f.pswap_spec[1]))
        if int(p1a.nbrspecies_store)>2:
          f.write("%s\n%s\n" %("insertion method",p3f.pswap_spec[2]))
        if int(p1a.nbrspecies_store)>3:
          f.write("%s\n%s\n" %("insertion method",p3f.pswap_spec[3]))
        if int(p1a.nbrspecies_store)>4:
          f.write("%s\n%s\n" %("insertion method",p3f.pswap_spec[4]))
        if int(p1a.nbrspecies_store)>5:
          f.write("%s\n%s\n" %("insertion method",p3f.pswap_spec[5]))
      if p3g.prob_flip !="0":
        f.write("\n\n# Prob_Ring\n%s %s" %(p3g.prob_flip,p3g.max_angle))
      f.write("\n# Done_Probability_Info\n%s\n" %(newsection))
      f.write("\n# Molecule_Files\n%s %s\n" %(p4a.mcf_s1[0],p4a.mcf_s1[1]))
      if int(p1a.nbrspecies_store)>1:
        f.write("%s %s\n" %(p4a.mcf_s2[0],p4a.mcf_s2[1]))
      if int(p1a.nbrspecies_store)>2:
        f.write("%s %s\n" %(p4a.mcf_s3[0],p4a.mcf_s3[1]))
      if int(p1a.nbrspecies_store)>3:
        f.write("%s %s\n" %(p4a.mcf_s4[0],p4a.mcf_s4[1]))
      if int(p1a.nbrspecies_store)>4:
        f.write("%s %s\n" %(p4a.mcf_s5[0],p4a.mcf_s5[1]))
      if int(p1a.nbrspecies_store)>5:
        f.write("%s %s\n" %(p4a.mcf_s6[0],p4a.mcf_s6[1]))
      f.write("%s\n" %(newsection))
      f.write("\n# Fragment_Files\n")
      for i in range(int(p4b.s1_nfrags)):
        f.write("%s %s\n" %(p4b.s1_frags[i],p4b.frag_index[i]))
      if int(p1a.nbrspecies_store)>1:
        for i in range(int(p4b.s2_nfrags)):
          f.write("%s %s\n" %(p4b.s2_frags[i],p4b.frag_index[i+s1_nf]))
      if int(p1a.nbrspecies_store)>2:
        for i in range(int(p4b.s3_nfrags)):
          f.write("%s %s\n" %(p4b.s3_frags[i],p4b.frag_index[i+s1_nf+s2_nf]))
      if int(p1a.nbrspecies_store)>3:
        for i in range(int(p4c.s4_nfrags)):
          f.write("%s %s\n" %(p4c.s4_frags[i],p4b.frag_index[i+s1_nf+s2_nf+s3_nf]))
      if int(p1a.nbrspecies_store)>4:
        for i in range(int(p4c.s5_nfrags)):
          f.write("%s %s\n" %(p4c.s5_frags[i],p4b.frag_index[i+s1_nf+s2_nf+s3_nf+s4_nf]))
      if int(p1a.nbrspecies_store)>5:
        for i in range(int(p4c.s6_nfrags)):
          f.write("%s %s\n" %(p4c.s6_frags[i],p4b.frag_index[i+s1_nf+s2_nf+s3_nf+s4_nf+s5_nf]))
      f.write("%s\n\n" %(newsection))
      f.write("# Pair_Energy\n%s\n%s\n\n" %(p1b.pairenergy_store,newsection))
      if p1a.ensembletype_store == "GCMC":
        f.write("# Chemical_Potential_Info\n%s %s %s %s %s %s\n\n%s\n\n" %(p1b.chempot_store[0],p1b.chempot_store[1],p1b.chempot_store[2],p1b.chempot_store[3],p1b.chempot_store[4],p1b.chempot_store[5],newsection))
      f.write("# Start_Type\n")
      if p4d.start_type == "make_config":
        f.write("make_config\n")
        if p1a.nbr_boxes == 1:
          for i in range(int(p1a.nbrspecies_store)):
            f.write("%s\n" %(p4d.mkconfig_dat[i]))
          f.write("%s\n\n" %(newsection))
        elif p1a.nbr_boxes == 2:
          for i in range(int(p1a.nbrspecies_store)):
            f.write("%s %s\n" %(p4d.mkconfig_dat[i],p4d.mkconfig_dat[i+6]))
          f.write("%s\n\n" %(newsection))
      elif p4d.start_type == "checkpoint":
        f.write("checkpoint\n%s\n%s\n\n" %(p4d.chkpoint,newsection))
      elif p4d.start_type == "read_old":
        f.write("read_old\n")
        if p1a.nbr_boxes == 1:
          f.write("%s\n%s\n\n" %(p4d.box1readold,newsection))
        elif p1a.nbr_boxes == 2:
          f.write("%s\n%s\n%s\n\n" %(p4d.box1readold,p4d.box2readold,newsection))
      f.write("# Run_Type\n%s %s %s\n%s\n\n" %(p4d.run_type[0],p4d.run_type[1],p4d.run_type[2],newsection))
      f.write("# Frequency_Info\n%s %s\n" %("freq_type",p4e.freq_info[0]))
      for i in range(3):
        f.write("%s %s\n" %(p4e.freq_keywords[i],p4e.freq_info[i+1]))
      f.write("# Done_Frequency_Info\n%s\n\n" %(newsection))
      f.write("# Average_Info\n%s\n%s\n\n" %(p4e.average_dat,newsection))
      self.b1propslist = []
      self.b2propslist = []
      for item in p4e.prop_b1:
        if item:
          self.b1propslist.append(item)
      for item in p4e.prop_b2:
        if item:
          self.b2propslist.append(item)
      if self.b1propslist:
        f.write("# Property_Info 1\n")
        for item in self.b1propslist:
          f.write("%s\n" %(item))
        f.write("%s\n\n" %(newsection))
      if int(p1a.nbr_boxes)>1:
        if self.b2propslist:
          f.write("# Property_Info 2\n")
          for item in self.b2propslist:
            f.write("%s\n" %(item))
          f.write("%s\n\n" %(newsection))
      f.write("# CBMC_Info\n%s %s\n%s %s\n%s %s\n%s %s %s\n%s\n" %("kappa_ins",p1b.kappa_ins,"kappa_rot",p1b.kappa_rot,"kappa_dih",p1b.kappa_dih,"rcut_cbmc",p1b.rcut_cbmc[0],p1b.rcut_cbmc[1],newsection))
      f.write("\n\nEND")

      f.close
if __name__ =="__main__":
  app = wx.App()
  MainFrame().Show()
  app.MainLoop()