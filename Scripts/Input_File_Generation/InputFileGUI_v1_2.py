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
#   Revision history: 
#   23 July 2015 - (BK) rewrote
#
#********************************************************************************

import re, wx, os, sys
from wx.lib.pubsub import setupkwargs
from wx.lib.pubsub import pub

# dictionary where we will put the user data
myDict = {}
obj_counter = 0
# angstroms unit
angstrom = u'\u212B'.encode('utf-8')

# max number of species allowed by the GUI
maxnbrspecies = 6
# max number of boxes allowed by the GUI
maxnbrboxes = 2

class p1a(wx.Panel):
  # create the widgets to be displayed on this panel, and bind them to 
  # relevant methods
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid1a = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid1a)

    self.dirname = ""
    self.staticlbls = []
  
    # place the instructional label
    alertstring = "Please provide the following information before continuing."
    self.alertlbl = wx.StaticText(self,label=alertstring)
    self.alertlbl.SetFont(wx.Font(14,wx.DEFAULT, wx.NORMAL, wx.BOLD))
    grid1a.Add(self.alertlbl,pos=(1,2),span=(1,4))
    
    # create and place the label widgets
    labelstr1 = ["Run Name: ", "Simulation Directory: ", "Ensemble: ",\
                 "Number of Species: "]
    for i in range(len(labelstr1)):
      self.staticlbls.append(wx.StaticText(self,label=labelstr1[i]))
      grid1a.Add(self.staticlbls[i],pos=(i+3,1))
    
    # asks for the run name of the simulation
    self.runnamequery = wx.TextCtrl(self,value="",name="runName")
    self.Bind(wx.EVT_TEXT,self.datafunc,self.runnamequery)

    # choice widget for the number of species in simulation
    nbrspec_options = ["%s" %(i) for i in range(maxnbrspecies+1)]
    nbrspec_options[0] = ""
    self.specieschoice=wx.Choice(self,-1,choices=nbrspec_options,\
                                 name="nbrSpecies",style=wx.ALIGN_CENTER)
    self.Bind(wx.EVT_CHOICE, self.datafunc, self.specieschoice)

    # create the choice widget for the ensemble
    ensemble_options = ["","NVT_MC","NVT_MIN","NPT_MC","GCMC","GEMC","GEMC_NPT"]
    self.ensemblechoice = wx.Choice(self,-1,choices=ensemble_options,name="ensemble")
    self.Bind(wx.EVT_CHOICE,self.datafunc,self.ensemblechoice)

    # button prompting for the directory of the simulation
    self.simdirbtn = wx.Button(self,label="Select Directory",name="simDir")
    self.Bind(wx.EVT_BUTTON,self.datafunc,self.simdirbtn)

    # displays the current directory as indicated by button selection
    self.simdirdisp= wx.TextCtrl(self,value="",size=(200,-1),\
                                 name="simdirdisp",style=wx.TE_READONLY)
    grid1a.Add(self.simdirdisp,pos=(4,3))

    # a list of the interactive widgets on this panel
    widgetslist = [self.runnamequery,self.simdirbtn,self.ensemblechoice,\
                   self.specieschoice]
    flagList = [(wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN\
                           | wx.ALIGN_CENTER)]*len(widgetslist)
   
    # place the widgets
    for index, item in enumerate(widgetslist):
      grid1a.Add(item,pos=(index+3,2),flag=flagList[index])

    # create the button for saving the files
    self.saveBtn = wx.Button(self,label="Create Input File",name="save")
    self.Bind(wx.EVT_BUTTON,self.saveFunction,self.saveBtn)
    grid1a.Add(self.saveBtn,pos=(12,2))

    # create the label describing the save button
    descString1 = "Have you completed all prompts on all pages?" 
    descString2 = "If so, click here:"
    self.saveBtnDesc = wx.StaticText(self,label=descString1)
    grid1a.Add(self.saveBtnDesc,pos=(11,1),span=(1,3))
    self.saveBtnDesc2 = wx.StaticText(self,label=descString2)
    grid1a.Add(self.saveBtnDesc2,pos=(12,1),span=(1,1))
    

  # places data as dictionary key value pairs; key is the name of 
  # the interactive widget associated with the event
  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    if (obj == "simDir"):
      val = str(self.button(obj))
      self.dirname = os.path.split(val)
      if val:
        self.simdirdisp.SetValue(("/%s/" %self.dirname[1]))
    data(obj,val)

    if obj == 'ensemble':
      if val[0:4] == 'GEMC':
        data('nbrBoxes','2')
      else:
        data('nbrBoxes','1')

  def button(self,obj):
    dlg = wx.DirDialog(self,"Select Simulation Directory", \
                       style = wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      val = dlg.GetPath()
      return val      
    else:
      val = ''
      return val

  def saveFunction(self,event):
   
   # loop through the dictionary and check what values may or may not be missing
   # allow an option to continue if missing values are found
   # so that the user may use the GUI to make a template

   ####### loop here, before creating string1 below

   # but for now just do a try/except KeyError. 
    try:
      string1 = str(myDict['runName'])+".inp"
      newsection = "!--------------------------------------"
      self.dirname = ''
      dialog = wx.FileDialog(self,"Save File As ...",myDict['simDir'],string1,"*.*",\
                           style = wx.SAVE)
  
      if dialog.ShowModal() == wx.ID_OK:
   
        path = dialog.GetPath()
        f = file(path,'w')
        f.write("! This is the input file for %s,\n! a %s-species %s simulation at %s K\n" \
             %(myDict['runName'],myDict['nbrSpecies'],\
               myDict['ensemble'],myDict['temperature']))
        f.write("! Created with Cassandra Editor v1.1\n\n")
        f.write("# Run_Name \n%s \n%s \n\n" %(myDict['runName'],newsection))
        f.write("# Sim_Type \n%s \n%s \n\n" %(myDict['ensemble'],newsection))
        f.write("# Nbr_Species \n%s \n%s \n\n" %(myDict['nbrSpecies'],newsection))
        f.write("# Rcutoff_Low \n%s \n%s\n\n" %(myDict['minCutoff'],newsection))
        f.write("# Temperature_Info\n")
        if myDict['ensemble'][0:4] == 'GEMC':
          f.write("%s %s\n%s\n\n" %(myDict['temperature'],myDict['temperature'],newsection))
        else:
          f.write("%s\n%s\n\n" %(myDict['temperature'],newsection))
        if 'pressure' in myDict:
          if myDict['ensemble'] == "GEMC" and int(myDict['nbrSpecies']) > 1:
            f.write('# Pressure_Info\n%s %s\n%s\n\n' %(myDict['pressure'],\
                  myDict['pressure'],newsection))
          elif myDict['ensemble'] == "GEMC_NPT":
            f.write('# Pressure_Info\n%s %s\n%s\n\n' %(myDict['pressure'],\
                  myDict['pressure'],newsection))
          else:
            f.write('# Pressure_Info\n%s\n%s\n\n' %(myDict['pressure'],newsection))
        f.write('# Seed_Info\n%s %s\n%s\n\n' %(myDict['seed1'],myDict['seed2'],newsection))
        f.write('# Mixing_Rule\n%s\n%s\n\n' %(myDict['mixingChoice'],newsection))
        f.write('# Charge_Style\n')
        # we don't know if all of these values are here; try to write them anyways
        tryToWrite(f,'box 1 charge functional','box 1 charge method','box 1 elec cutoff',\
                     'box 1 ewald accuracy')
        if myDict['ensemble'][0:4] == 'GEMC':
          tryToWrite(f,'box 2 charge functional','box 2 charge method','box 2 elec cutoff',\
                       'box 2 ewald accuracy')
        f.write("%s\n\n# VDW_Style\n" %newsection)
        tryToWrite(f,'box 1 vdw functional','box 1 vdw tail','box 1 vdw spline on',\
                     'box 1 vdw spline off','box 1 logical')
        if myDict['ensemble'][0:4] == 'GEMC':
          tryToWrite(f,'box 2 vdw functional','box 2 vdw tail','box 2 vdw spline on',\
                       'box 2 vdw spline off','box 2 logical')
        f.write('%s\n\n' %newsection)
        f.write('# Intra_Scaling \n')
        for i in range(int(myDict['nbrSpecies'])):
          f.write('%s %s %s %s\n%s %s %s %s\n' %(myDict['s%s vdw 1-2' %(i+1)], \
                  myDict['s%s vdw 1-3' %(i+1)],myDict['s%s vdw 1-4' %(i+1)],\
                  myDict['s%s vdw 1-N' %(i+1)],myDict['s%s coul 1-2' %(i+1)],\
                  myDict['s%s coul 1-3' %(i+1)],myDict['s%s coul 1-4' %(i+1)],\
                  myDict['s%s coul 1-N' %(i+1)]))
        f.write('%s\n\n' %newsection)

        if myDict['box1Shape'] == 'CUBIC':
          f.write('# Box_Info\n%s\n%s\n%s %s %s\n' %(myDict['nbrBoxes'],myDict['box1Shape'],\
                  myDict['box1Length'],myDict['box1Length'],myDict['box1Length']))
        elif myDict['box1Shape'] == 'NON-CUBIC':
          f.write('# Box_Info\n%s\n' %(myDict['nbrBoxes']))
          # check whether the off diagonal components are zero or not.
          # if they are all zero, then we have an orthorhombic non-cubic
          # else, use 'cell matrix' keyword
          checkForThese = ['xy 1','xz 1','yx 1','yz 1','zx 1','zy 1']
          thisKey = 'ORTHOGONAL'
          for index, item in enumerate(checkForThese):
            if item in myDict:
              if float(myDict[item]) != 0.0:
                thisKey = 'CELL_MATRIX'
                break
          if thisKey == 'ORTHOGONAL':
            f.write('%s\n%s %s %s\n' %(thisKey,myDict['xx 1'],myDict['yy 1'],myDict['zz 1']))
          elif thisKey == 'CELL_MATRIX':
            f.write('%s\n%s %s %s\n%s %s %s\n%s %s %s\n' %(thisKey, myDict['xx 1'],myDict['xy 1'],\
                    myDict['xz 1'],myDict['yx 1'],myDict['yy 1'],myDict['yz 1'],myDict['zx 1'],\
                    myDict['zy 1'],myDict['zz 1']))

        if int(myDict['nbrBoxes']) > 1:
          f.write('\n%s\n%s %s %s \n' %(myDict['box2Shape'],myDict['box2Length'],\
                                        myDict['box2Length'],myDict['box2Length']))
        f.write('%s\n\n' %newsection)
        f.write('# Move_Probability_Info\n')
        if 'prob translation' in myDict:
          f.write('\n# Prob_Translation\n%s\n' %myDict['prob translation'])
          tryToWrite(f,'prob trans s1 b1', 'prob trans s2 b1', \
                       'prob trans s3 b1','prob trans s4 b1',\
                       'prob trans s5 b1', 'prob trans s6 b1')
          if int(myDict['nbrBoxes'])>1:
            tryToWrite(f,'prob trans s1 b2','prob trans s2 b2',\
                         'prob trans s3 b2','prob trans s4 b2',\
                         'prob trans s5 b2','prob trans s6 b2')
        if 'prob rotation' in myDict:
          f.write('\n# Prob_Rotation\n%s\n' %(myDict['prob rotation']))
          tryToWrite(f,'prob rot s1 b1','prob rot s2 b1',\
                       'prob rot s3 b1','prob rot s4 b1',\
                       'prob rot s5 b1','prob rot s6 b1')
          if int(myDict['nbrBoxes'])>1:
            tryToWrite(f,'prob rot s1 b2','prob rot s2 b2',\
                         'prob rot s3 b2','prob rot s4 b2',\
                         'prob rot s5 b2','prob rot s6 b2')
        if 'prob regrowth' in myDict:
          f.write('\n# Prob_Regrowth\n%s\n' %(myDict['prob regrowth']))
          tryToWrite(f,'prob regrowth s1','prob regrowth s2',\
                       'prob regrowth s3','prob regrowth s4',\
                      'prob regrowth s5','prob regrowth s6')
        if 'prob vol' in myDict:
          f.write('\n# Prob_Volume\n%s\n' %myDict['prob vol'])
          f.write('%s\n' %(myDict['prob vol b1']))
          if int(myDict['nbrBoxes'])>1:
            f.write('%s\n' %(myDict['prob vol b2']))
        if 'prob insertion' in myDict:
          f.write('\n# Prob_Insertion\n%s\n' %(myDict['prob insertion']))
          f.write('%s' %(myDict['insertion method s%s' %(1)]))
          for i in range(1,int(myDict['nbrSpecies'])):
            f.write(' %s' %(myDict['insertion method s%s' %(i+1)]))
          f.write('\n')
        if 'prob insertion' in myDict:
          f.write('\n# Prob_Deletion\n%s\n' %(myDict['prob insertion']))
        if 'prob swap' in myDict:
          f.write('\n# Prob_Swap\n%s\n' %(myDict['prob swap']))
          f.write('%s' %(myDict['insertion method s%s' %(1)]))
          for i in range(1,int(myDict['nbrSpecies'])):
            f.write(' %s' %(myDict['insertion method s%s' %(i+1)]))
          f.write('\n')
        f.write('\n# Done_Probability_Info\n%s\n\n' %(newsection))
        f.write('# Molecule_Files\n')
        for i in range(int(myDict['nbrSpecies'])):
          f.write('%s %s\n' %(myDict['MCF s%s' %(i+1)],\
                  myDict['max nmols s%s' %(i+1)]))
        f.write('%s\n\n' %(newsection))
        if myDict['fragFilesCheck'] == 'Yes':
          f.write('# Fragment_Files\n')
          frag_count = 0
          for i in range(int(myDict['nbrSpecies'])):
            for index, item in enumerate(myDict['frag files s%s' %(i+1)]):
              frag_count = frag_count + 1
              f.write('%s %s\n' %(str(item),str(frag_count)))
          f.write('%s\n\n' %(newsection)) 

        f.write('# Pair_Energy\n%s\n' %(myDict['pairStorage']))
        f.write('%s\n\n' %(newsection))
        if myDict['ensemble'] == 'GCMC':
          f.write('# Chemical_Potential_Info\n')
          tryToWrite(f,'chempot 1','chempot 2','chempot 3',\
                       'chempot 4','chempot 5','chempot 6')
          f.write('%s\n\n' %(newsection))
        f.write('# Start_Type\n')
        if myDict['start type'] == 'make_config':
          for j in range(int(myDict['nbrBoxes'])):
            f.write('make_config')
            for i in range(int(myDict['nbrSpecies'])):
              f.write(' %s' %(myDict['makeConfig s%s b%s' %(i+1,j+1)]))
            f.write('\n')
          f.write('%s\n\n' %(newsection))
        elif myDict['start type'] == 'checkpoint':
          f.write('checkpoint\n%s\n%s\n\n' %(myDict['checkpoint file'],newsection))
        elif myDict['start type'] == 'read_config':
          for j in range(int(myDict['nbrBoxes'])):
            f.write('read_config')
            for i in range(int(myDict['nbrSpecies'])):
              f.write(' %s' %(myDict['readConfig s%s b%s' %(i+1,j+1)]))
            f.write(' %s\n' % (myDict['read old box %s' %(j+1)]))
          f.write('%s\n\n' %(newsection))
        f.write('# Run_Type\n')
        tryToWrite(f,'run type','acceptance ratio output','volume displacement output')
        f.write('%s\n\n' %(newsection)) 

        f.write('# Simulation_Length_Info\n%s %s\n' %('Units',myDict['freqTypeChoice']))
        if myDict['freqTypeChoice'] == 'Minutes':
          f.write('%s %s\n%s %s\n%s %s\n' %('Prop_Freq',myDict['thermoFreq'],\
                  'Coord_Freq',myDict['coordFreq'],'Run',myDict['simRunTime']))
          f.write('%s\n\n' %(newsection))
        elif myDict['freqTypeChoice'] == 'Steps':
          f.write('%s %s\n%s %s\n%s %s\n' %('Prop_Freq',myDict['thermoFreq'],\
                  'Coord_Freq',myDict['coordFreq'],'Run',myDict['simRunTime']))
          f.write('%s\n\n' %(newsection))

        f.write('# Average_Info\n%s\n%s\n\n' %(myDict['avgInfo'],newsection))
        f.write('# Property_Info 1\n')
        for index, item in enumerate(myDict['propList b1']):
          f.write('%s\n' %(item))
        f.write('%s\n\n' %(newsection))
        if int(myDict['nbrBoxes']) > 1:
          f.write('# Property_Info 2\n')
          for index, item in enumerate(myDict['propList b2']):
            f.write('%s\n' %(item))
          f.write('%s\n\n' %(newsection))
        f.write('# CBMC_Info\n%s %s\n%s %s\n%s %s\n' %('kappa_ins',myDict['kappains'],\
                'kappa_rot',myDict['kapparot'],'kappa_dih',myDict['kappadih']))
        if int(myDict['nbrBoxes']) == 1:
          f.write('%s %s\n' %('rcut_cbmc',myDict['rcutCBMC 1']))
        elif int(myDict['nbrBoxes']) == 2:
          f.write('%s %s %s\n' %('rcut_cbmc',myDict['rcutCBMC 1'],myDict['rcutCBMC 2']))
        f.write('%s\n\n' %(newsection))
        f.write('\n\nEND')
        f.close

    # this will identify what key threw an exception.
    # in the future, possibly link each key to a more user-friendly string,
    # to aid in identification of exactly what information was missing.
    except KeyError as err:
      cause = err.args[0]
      print 'Missing required info: %s \nInput file not written.' %(str(cause)) +\
            '  Please check all panels for completeness.'

      return
  # check the dictionary for all keys to be used in the save file function
  # feature to be added
  def checkDictionary(self,keepGoing):
    pass

class p1b(wx.Panel):
 
  def __init__(self,parent):
    pub.subscribe(self.showAndHide, "refresh")

    # explicity tell it the coordinates for every widget in
    # a list of tuples; same with spans, etc
    wx.Panel.__init__(self,parent)
    self.grid1b = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(self.grid1b)

    mixingOptions = ["","Lorentz-Berthelot", "Geometric"]
    pairNrgOptions = ["","TRUE", "FALSE"]

    # 
    self.boxShapeOptions = ["","CUBIC"]
    choicesList = [(pairNrgOptions), (self.boxShapeOptions), \
                   (self.boxShapeOptions),(mixingOptions)]

    lblList1 = ["Temperature (K): ",         "Mixing Rule: ",       
                "Minimum Cutoff (%s): " %(angstrom),    "Pair Storage: ",
                "Seed 1: ",                  "Seed 2: ",          
                "Pressure (bar): ",          "Box Information",
                "Box Shape",                 "Box Edge Length (%s)" %(angstrom),  
                "Box 1: ",                    "Box 2: ",
                "CBMC Parameters",           "Trial Insertions: ",    
                "Rotational Bias: ",         "Trial Orientations: ",
                "Cutoff (%s) Box 1: " %(angstrom),  "Cutoff (%s) Box 2: " %(angstrom),  
                "Chemical Potential (kJ/mol)"]
 
    lblListCoords = [(1,1),  (2,1),   (3,1),    (4,1),
                     (5,1),  (6,1),   (7,1),    (9,1),
                     (9,2),  (9,3),   (10,1),   (11,1),
                     (13,1), (14,1),  (15,1),   (16,1),
                     (17,1), (18,1),  (1,5)]   

    self.staticlbls1 = []
    spanList = [(1,1)]*(len(lblListCoords)-1) + [(1,2)]
    staticFlagList = [(wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.ALIGN_LEFT)]*8 + \
                     [(wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)]*2 + \
                     [(wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.ALIGN_LEFT)]*8 + \
                     [(wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)]

    for index, item in enumerate(lblList1):
      self.staticlbls1.append(wx.StaticText(self,label=str(item), name = item))
      self.grid1b.Add(self.staticlbls1[index],pos=lblListCoords[index],\
                      span=spanList[index],flag=staticFlagList[index])

    self.speclbl = []
    for i in range(maxnbrspecies):
      self.speclbl.append(wx.StaticText(self,label=("Species %d:" %(i+1))))
      self.grid1b.Add(self.speclbl[i],pos=(i+2,5),flag= \
                      (wx.ALIGN_RIGHT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN))
     
    chempotnames = []
    self.chempotwidgets = []
    for i in range(maxnbrspecies):
      chempotnames.append("chempot %s" %(i+1))
      self.chempotwidgets.append(wx.TextCtrl(self,value="",name=chempotnames[i]))
      self.grid1b.Add(self.chempotwidgets[i],pos=(i+2,6))
      self.Bind(wx.EVT_TEXT,self.datafunc,self.chempotwidgets[i])

    self.widgets1b = []
    widgetnames = ["temperature",        "minCutoff",   \
                   "seed1",              "seed2",       \
                   "pressure",           "box1Length",  \
                   "box2Length",         "kappadih",    \
                   "kappains",           "kapparot",    \
                   "rcutCBMC 1",         "rcutCBMC 2",  \
                   "pairStorage",        "box1Shape",   \
                   "box2Shape",          "mixingChoice" ]
   
    coords1 = [(1,2),   (3,2),   (5,2),   (6,2),   (7,2), 
               (10,3),  (11,3),  (16,2),  (14,2),  (15,2), 
               (17,2), (18,2),   (4,2),   (10,2),   (11,2),
               (2,2)]

    span1 = [(1,1)]*15 + [(1,3)]
    widgetFlags = [(wx.RESERVE_SPACE_EVEN_IF_HIDDEN)]*len(widgetnames)

    # add 'style' attribute
    styleList = [(wx.TE_LEFT)]*9 + [(wx.TE_READONLY)] + [(wx.TE_LEFT)]*2
    for index, item in enumerate(widgetnames):
      # create the widgets; last ones are the choice widgets
      if index < 12:
        self.widgets1b.append(wx.TextCtrl(self,value="",\
                        name=item,style=styleList[index]))
        self.grid1b.Add(self.widgets1b[index],\
                        pos=coords1[index],span=span1[index],\
                        flag = widgetFlags[index])
        self.Bind(wx.EVT_TEXT, self.datafunc, self.widgets1b[index])
      else:
        self.widgets1b.append(wx.Choice(self,-1,\
                    choices=choicesList[index-12],name=item))
        self.grid1b.Add(self.widgets1b[index],\
                        pos=coords1[index],span=span1[index],\
                        flag = widgetFlags[index])
        self.Bind(wx.EVT_CHOICE, self.datafunc, self.widgets1b[index])


    # hard code kappa rot as 0 for now.  
    # This will be changed in future versions.
    self.widgets1b[9].SetValue('0')

    # create h-matrix buttons
    hMatrixNames = ["h matrix 1","h matrix 2"]
    hMatrixpos = [(10,4),(11,4)]
    self.hMatrixBtns = []
    hMatrixFlag = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    for index, item in enumerate(hMatrixNames):
      self.hMatrixBtns.append(wx.Button(self,label="Create H Matrix",\
        name = item))
      self.grid1b.Add(self.hMatrixBtns[index],pos=hMatrixpos[index],\
                      flag=hMatrixFlag)
      self.Bind(wx.EVT_BUTTON, self.hMatrixFunc, self.hMatrixBtns[index])

  def datafunc(self,event):  
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())

    # adjust values for cassandra if needed
    if val == 'Lorentz-Berthelot':
      val = 'LB'
    if val == 'Gemoetric':
      val = 'geometric'

    if obj == 'box1Shape' or obj == 'box2Shape':
      # then we need to display the appropriate prompts
      # pass the widget name and the value of the string
      self.dispBoxPrompts(obj,val)

    data(obj,val)
  
  def dispBoxPrompts(self,obj,val):
    if obj == 'box1Shape':
      if val == 'CUBIC':
        # hide box 1 h matrix
        self.hMatrixBtns[0].Hide()
        # show the text widget and label
        self.staticlbls1[9].Show()
        self.widgets1b[5].Show()
      elif val == 'NON-CUBIC':
        self.widgets1b[5].Hide()
        self.hMatrixBtns[0].Show()
        if not 'box2Shape' in myDict:
          self.staticlbls1[9].Hide()
        elif myDict['box2Shape'] != 'CUBIC':
          self.staticlbls1[9].Hide()
        else:
          self.staticlbls1[9].Show()
      elif not val:
        self.hMatrixBtns[0].Hide()
        self.widgets1b[5].Hide()
        if not 'box2Shape' in myDict:
          self.staticlbls1[9].Hide()
        elif myDict['box2Shape'] != 'CUBIC':
          self.staticlbls1[9].Hide()
        else:
          self.staticlbls1[9].Show()

    elif obj == 'box2Shape':
      if val == 'CUBIC':
        self.widgets1b[6].Show()
        self.hMatrixBtns[1].Hide()
        self.staticlbls1[9].Show()
      elif val == 'NON-CUBIC':
        self.widgets1b[6].Hide()
        self.hMatrixBtns[1].Show()
        if 'box1Shape' in myDict:
          if myDict['box1Shape'] != 'CUBIC':
            self.staticlbls1[9].Hide()
          else:
            self.staticlbls1[9].Show()
        else:
          self.staticlbls1[9].Hide()
      elif not val:
        self.widgets1b[6].Hide()
        self.hMatrixBtns[1].Hide()
        if 'box1Shape' in myDict:
          if myDict['box1Shape'] != 'CUBIC':
            self.staticlbls1[9].Hide()
          else:
            self.staticlbls1[9].Show()
        else:
          self.staticlbls1[9].Hide()
    self.Layout()
        

    # now, hide the box edge length (angstrom) labels, textwidgets,
    # and the hMatrix buttons.
    #for index, item in enumerate(self.hMatrixBtns):
    #  item.Hide()
    #self.staticlbls1[9].Hide()
    #self.widgets1b[5].Hide()
    #self.widgets1b[6].Hide()



  # we can trigger this functionality using the pub module within wxPython. 
  # For some reason, simply calling a class method from the global 'data' function
  # without using the publisher module did not produce 
  # the desired show/hide functionality. Could be due to incorrect event propagation?
  # in any event, the publisher module is used across all panels and is prevalent
  # in all or nearly all classes of the GUI.
  def showAndHide(self, message, arg2 = None):
    # return without doing anything if missing a required value. 
    if not 'ensemble' in myDict.keys():
      return
    if not 'nbrSpecies' in myDict.keys():
      return

    # make a list of widgets only required by GEMC; in this case, widgets that allude
    # to properties belonging to a second box
    showOnlyIfGEMC = [self.widgets1b[6], self.widgets1b[11], \
                      self.widgets1b[14], self.staticlbls1[11], self.staticlbls1[17]]

    # these widgets pertain to chemical potentials
    showOnlyIfGCMC = [self.staticlbls1[18]] + [self.speclbl[i] \
                     for i in range(maxnbrspecies)] + \
                     [self.chempotwidgets[i] for i in range(maxnbrspecies)]
    
    # pressure widgets
    pressureWidgets = [self.widgets1b[4], self.staticlbls1[6]]

    # GEMC behavior
    if myDict['ensemble'][0:4] != "GEMC":
      for item in showOnlyIfGEMC:
        item.Hide()
    else:
      for item in showOnlyIfGEMC:
        item.Show()
       # Gibbs phase rule - if there are 2 or more species in solution, then we will need
       # a pressure.  If there is one species, we do not require a pressure.
      if int(myDict['nbrSpecies']) >= 2 or myDict['ensemble'] == 'GEMC_NPT':
        for item in pressureWidgets:
          item.Show()
      else:
        for item in pressureWidgets:
          item.Hide()

    # GCMC behavior
    if myDict['ensemble'] != "GCMC":
      for item in showOnlyIfGCMC:
        item.Hide()
    else:
      self.staticlbls1[18].Show()
      for index, item in enumerate(self.speclbl):
        if index in range(int(myDict['nbrSpecies'])):
          item.Show()
        else: item.Hide()
      for index, item in enumerate(self.chempotwidgets):
        if index in range(int(myDict['nbrSpecies'])):
          item.Show() 
        else: item.Hide()
      for item in pressureWidgets:
        item.Hide()
      self.widgets1b[4].SetValue('')

    # NPT_MC behavior
    if myDict['ensemble'] == "NPT_MC":
      for item in pressureWidgets:
        item.Show()

    # NVT_MC behavior; use this for all NVT
    if myDict['ensemble'][0:3] == "NVT":
      for item in pressureWidgets:
        item.Hide()

    # now, hide the box edge length (angstrom) labels, textwidgets,
    # and the hMatrix buttons.
    for index, item in enumerate(self.hMatrixBtns):
      item.Hide()
    self.staticlbls1[9].Hide()
    self.widgets1b[5].Hide()
    self.widgets1b[6].Hide()

    # we will also use this function to alter the list contained as choices
    # for box shape
    if myDict['ensemble'] == "GCMC" or myDict['ensemble'][0:3] == "NVT":
      boxShapeOptions = ["","CUBIC","NON-CUBIC"]
      # box shape choice widgets are index 13, 14 in self.widgets1b for boxes
      # 1 & 2, respectively
      # clear the current choices
      self.widgets1b[13].Clear()
      self.widgets1b[14].Clear()
      self.widgets1b[13].AppendItems(boxShapeOptions)
      self.widgets1b[14].AppendItems(boxShapeOptions)
    else:
      # clear the box shape options and set as ["","CUBIC"]
      boxShapeOptions = ["","CUBIC"]
      self.widgets1b[13].Clear()
      self.widgets1b[14].Clear()
      self.widgets1b[13].AppendItems(boxShapeOptions)
      self.widgets1b[14].AppendItems(boxShapeOptions)

      
  
    self.Layout()

  def hMatrixFunc(self,event):
    # create H object on self to allow hMatrixFrame to inherit the data
    # self.Hobj()
    self.Hobj = str(event.GetEventObject().GetName())
    # the box number is the last character of the object name that was clicked
    self.Hobj = str(self.Hobj[-1])
    # initialize the frame and pass self to the child frame
    thisWindow = hMatrixFrame(self)
    # at this point we assume the button is only available if the box shape option
    # for non-orthorhombic is a valid option
    thisWindow.Show(True)


# define the hMatrixFrame that will be standalone and called from the function
# hMatrixFunc(self,event) above
class hMatrixFrame(wx.Frame):
  def __init__(self,parent):

    # we can inherit the self.() objects from the calling method
    # except now they will have the parent.() syntax

    thisIndex = parent.Hobj 
    wx.Frame.__init__(self,parent, title = "H Matrix")
    panel = wx.Panel(self,-1)
    hGrid = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(hGrid)

    # hopefully this is ~correct size for all operating systems (?)
    self.SetInitialSize((400,300))
  
    lblList = ['Enter your H-matrix vectors below.','x','y','z','x','y','z']
    lblCoords = [(1,1)] + [(2,i+2) for i in range(3)] +\
                [(i+3,1) for i in range(3)]
    lblSpan = [(1,4)] + [(1,1)]*(len(lblList)-1)
    lblFlags = [(wx.ALIGN_CENTER)]*len(lblList)
    self.hMatrixLbls = []

    for index, item in enumerate(lblList):
      self.hMatrixLbls.append(wx.StaticText(self,label=item))
      hGrid.Add(self.hMatrixLbls[index],pos=lblCoords[index],\
                span = lblSpan[index],flag=lblFlags[index])
    
    # these will be text widgets forming the h-matrix
    #    xx      yx       zx
    #    xy      yy       zy
    #    xz      yz       zz
    
    widgetNames = ['xx %s' %(thisIndex),'xy %s' %(thisIndex),'xz %s' %(thisIndex),\
                   'yx %s' %(thisIndex),'yy %s' %(thisIndex),'yz %s' %(thisIndex),\
                   'zx %s' %(thisIndex),'zy %s' %(thisIndex),'zz %s' %(thisIndex)]
    widgetCoords = [(3,i+2) for i in range(3)] +\
                   [(4,i+2) for i in range(3)] +\
                   [(5,i+2) for i in range(3)]
    self.hMatrixWidgets = []

    for index, item in enumerate(widgetNames):
      if str(item) in myDict:
        val = myDict[str(item)]
      else:
        val = ''
      self.hMatrixWidgets.append(wx.TextCtrl(self,value=val,name=item))
      hGrid.Add(self.hMatrixWidgets[index],pos=widgetCoords[index])
      self.Bind(wx.EVT_TEXT,self.hDatafunc,self.hMatrixWidgets[index])


    # add a 'done' button that can close the window
    self.doneButton = wx.Button(self,label="Done")
    hGrid.Add(self.doneButton,pos=(7,4))
    self.Bind(wx.EVT_BUTTON, self.doneButtonFunc, self.doneButton)

  def hDatafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    data(obj,val)


  def doneButtonFunc(self,event):
    self.Close()


#########################
class PanelOne(wx.Panel):

  def __init__(self,parent):

    # this places p1a and p1b as children of PanelOne, 
    # allowing for the tabs within the "Basic Information" tab
    wx.Panel.__init__(self,parent)
    nested_notebook = wx.Notebook(self)
    self.p1a = p1a(nested_notebook)
    self.p1b = p1b(nested_notebook)

    nested_notebook.AddPage(self.p1a,"Page 1")
    nested_notebook.AddPage(self.p1b,"Page 2")

    nested_sizer = wx.BoxSizer()
    nested_sizer.Add(nested_notebook,1,wx.EXPAND)
    self.SetSizer(nested_sizer)

    # if changing the page, check to make sure all info is put in to dictionary
    self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)

  def OnPageChanging(self,event):
    oldPage = event.GetOldSelection()
    newPage = event.GetSelection()
    
    # if the page we are leaving is 'Page 1' of basic information, 
    # check that our dictionary has the required keywords
    if oldPage == 0:
      required_keys = ["runName", "nbrSpecies", "ensemble", "simDir"]
      count_num_keys = 0
      # check if these keys exist in the dictionary; if they are in the
      # dictionary, then the data we require is there
      for item in required_keys:
        if item in myDict:
          count_num_keys = count_num_keys + 1
      # if these keys are not in the dictionary, veto the page change request
      if count_num_keys < len(required_keys):
        event.Veto()
    
##########################################
# intermolecular panel
class p2a(wx.Panel):

  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid2a = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid2a)

    pub.subscribe(self.showAndHide, "refresh")

    # create the labels, coordinates, and spans of the label widgets
    lblList = ["van der Waals Style",      "Charge Style",       \
               "Box 1: ",                  "Box 2: ",            \
               "Box 1: ",                  "Box 2: ",            \
               "Functional Form" ,    "Tail Correction",         \
               "Spline On (%s)" %(angstrom), "Spline Off (%s)" %(angstrom),    \
               "Logical (Optional)",       "Functional Form",    \
               "Tail Correction",       "Spline On (%s)" %(angstrom),          \
               "Spline Off (%s)" %(angstrom),  "Logical (Optional)",           \
               "Functional Form",          "Method",             \
               "Cutoff (%s)" %(angstrom),    "Accuracy",         \
               "Functional Form",          "Method",             \
               "Cutoff (%s)" %(angstrom),  "Accuracy"]

    lblListCoords = [(1,1),   (8,1),    (3,1),    (6,1),    (10,1),
                     (13,1),  (2,2),    (2,3),    (2,4),    (2,5),
                     (2,6),   (5,2),    (5,3),    (5,4),    (5,5),
                     (5,6),   (9,2),    (9,3),    (9,4),    (9,5),
                     (12,2),  (12,3),   (12,4),   (12,5)]

    spanList = [(1,2)] + [(1,1)]*(len(lblListCoords)-1)
    staticFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

    self.staticLbls2a = []
  
    # place the labels on the panel
    for index, item in enumerate(lblList):
      self.staticLbls2a.append(wx.StaticText(self,label=item))
      grid2a.Add(self.staticLbls2a[index],\
                 pos=lblListCoords[index],span=spanList[index])

 

    # list of options for the choice widgets
    vdwFunctionalOptions = ["", "Lennard Jones 12-6", "None"]
    vdwTailCorrectionOptions = ["","cut","cut_tail","cut_switch","cut_shift"]
    chargeFunctionalOptions = ["","Coulombic","NONE"]
    chargeMethodOptions = ["","Ewald","cut"]
    logicalOptions = ["", "TRUE"]

    choicesList = [(vdwFunctionalOptions), (vdwTailCorrectionOptions)]*2 + \
                  [(chargeFunctionalOptions), (chargeMethodOptions)]*2 + \
                  [(logicalOptions)]*2
    choiceCoords = [(3,2),   (3,3),  (6,2), (6,3),
                    (10,2),  (10,3), (13,2), (13,3),
                    (3,6),   (6,6)]
    choiceNames = ["box 1 vdw functional", "box 1 vdw tail", \
                   "box 2 vdw functional", "box 2 vdw tail",
                   "box 1 charge functional", "box 1 charge method", \
                   "box 2 charge functional", "box 2 charge method", \
                   "box 1 logical",        "box 2 logical"]
    choiceFlags = (wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.ALIGN_CENTER | wx.EXPAND)



    self.choicewidgets2a = []
    # create the choice widgets, place to grid, and bind to function
    for index, item in enumerate(choiceNames):
      self.choicewidgets2a.append(wx.Choice(self,-1,\
                choices=choicesList[index],name=item))
      grid2a.Add(self.choicewidgets2a[index],\
                 pos=choiceCoords[index],flag=choiceFlags)
      self.Bind(wx.EVT_CHOICE, self.datafunc, self.choicewidgets2a[index])


    # note that 'cutoff' and 'spline on' will both be stored in 'box 1 spline on'
    self.textNames = ["box 1 vdw spline on", "box 1 vdw spline off", \
                 "box 2 vdw spline on", "box 2 vdw spline off", \
                 "box 1 elec cutoff",   "box 1 ewald accuracy", \
                 "box 2 elec cutoff",   "box 2 ewald accuracy"]
    textCoords = [(3,4), (3,5), (6,4), (6,5),
                  (10,4), (10,5), (13,4), (13,5)]
    textFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.EXPAND)



    
    self.textwidgets2a = []
    for index, item in enumerate(self.textNames):
      self.textwidgets2a.append(wx.TextCtrl(self,value="",name=item))
      grid2a.Add(self.textwidgets2a[index],pos=textCoords[index],flag=textFlags)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets2a[index])


    # make a list of the indices of labels relevant to 
    # box 2 - only display these if the selected ensemble is GEMC
    # use enumerate on the relevant lists to verify these are the correct indices
    self.onlyShowGEMC_staticlbls = [3, 5, 11, 12, 13, 14, 15, 20, 21, 22, 23]
    self.onlyShowGEMC_choicewidgets = [2, 3, 6, 7, 9]
    self.onlyShowGEMC_textwidgets = [2,3,6,7]

    # a list of items to initially hide
    self.initially_hide_labels = [7, 8, 9, 10, 17, 18, 19]
    self.initially_hide_choices = [1, 5, 8]
    self.initially_hide_textwidgets = [0, 1, 4, 5]

    for item in self.onlyShowGEMC_staticlbls:
      self.staticLbls2a[item].Hide()
    for item in self.onlyShowGEMC_choicewidgets:
      self.choicewidgets2a[item].Hide()
    for item in self.onlyShowGEMC_textwidgets:
      self.textwidgets2a[item].Hide() 
    for item in self.initially_hide_labels:
      self.staticLbls2a[item].Hide()
    for item in self.initially_hide_choices:
      self.choicewidgets2a[item].Hide()
    for item in self.initially_hide_textwidgets:
      self.textwidgets2a[item].Hide()

  def datafunc(self,event): 
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    if val == "Lennard Jones 12-6":
      val = "LJ"
    elif val == "Coulombic":
      val = "coul"
    data(obj,val)

    selection = event.GetSelection()
    # inform the function what object the user changed
    self.refreshThisPanel(obj,selection)
    # show and hide relevant widgets
    
  # dynamically updates the widgets displayed on the panel according to user's selections
  def refreshThisPanel(self,obj,selection):
    # if obj is text ctrl, we can return without doing anything
    # if the object is a choice control, proceed
    if obj in self.textNames:
      return
    # the logicals have no effect on what is displayed
    if obj == "box 1 logical" or obj == "box 2 logical":
      return

    labels_vdw = []
    labels_vdw.append([])
    labels_vdw.append([])
    # we need to make a list with: tail correction, spline on, spline off, logical
    # for both box 1 and box 2.

    for i in range(4):
      labels_vdw[0].append((i+7))
    for i in range(4):
      labels_vdw[1].append((i+12))



    labels_coul = []
    labels_coul.append([])
    labels_coul.append([])
    labels_coul[0] = [17,18,19]
    labels_coul[1] = [21,22,23]

    # we start with the functional form
    choice_vdw = []
    choice_vdw.append([])
    choice_vdw.append([])
    choice_coul = []
    choice_coul.append([])
    choice_coul.append([])
    text_vdw = []
    text_vdw.append([])
    text_vdw.append([])
    text_coul = []
    text_coul.append([])
    text_coul.append([])

    # these values are just the indices of the widgets in their appropriate object
    # e.g., choice_vdw[0] lists indices of box 1 choice widgets in self.choicewidgets2a
    # and choice_vdw[1] lists indices of box 2 choice widgets in self.choicewidgets2a,
    # specifically those in the 'vdw style' section (as opposed to charge style)

    # they can be checked by using enumerate() on the self.choicewidgets2a list
    choice_vdw[0] = [1,8]
    choice_vdw[1] = [3,9]
    choice_coul[0] = [5]
    choice_coul[1] = [7]

    text_vdw[0] = [0,1]
    text_vdw[1] = [2,3]
    text_coul[0] = [4,5]
    text_coul[1] = [6,7]

    if "vdw functional" in obj:
      box_id = int(obj[4]) - 1
      # selected empty value, or none, for vdw functional form
      if selection == 0 or selection == 2:
        for item in labels_vdw[box_id]:
          self.staticLbls2a[item].Hide()
        for item in choice_vdw[box_id]:
          self.choicewidgets2a[item].Hide()
          self.choicewidgets2a[item].SetSelection(0)
        for item in text_vdw[box_id]:
          self.textwidgets2a[item].Hide()
          self.textwidgets2a[item].SetValue('')
      else:
        # selection is "lennard Jones 12-6" for the functional form, so display the next option
        self.choicewidgets2a[choice_vdw[box_id][0]].Show()
        self.staticLbls2a[labels_vdw[box_id][0]].Show()

    # now, we look at the user's selection for vdw tail correction
    elif "vdw tail" in obj:
      box_id = int(obj[4]) - 1
      if selection == 0:
        for item in labels_vdw[box_id][1:]:
          self.staticLbls2a[item].Hide()
        for item in choice_vdw[box_id][1:]:
          self.choicewidgets2a[item].Hide()
          self.choicewidgets2a[item].SetSelection(0)
        for item in text_vdw[box_id]:
          self.textwidgets2a[item].Hide()
          self.textwidgets2a[item].SetValue('')
      # the user selected 'cut' OR 'cut_shift' as their tail correction
      elif selection == 1 or selection == 4:
        # note that cut and cut_shift have the same requirements according to user_guide.pdf
        # here, we change the label to Cutoff (angs); hide the spline off widget & label,
        # and hide the relevant logical and its label
        changeThisLabel = labels_vdw[box_id][1]
        self.staticLbls2a[changeThisLabel].SetLabel('Cutoff (%s)' %angstrom)
        self.staticLbls2a[changeThisLabel].Show()
        for item in labels_vdw[box_id][2:]:
          self.staticLbls2a[item].Hide()
        self.choicewidgets2a[choice_vdw[box_id][1]].Hide()
        self.choicewidgets2a[choice_vdw[box_id][1]].SetSelection(0)
        self.textwidgets2a[text_vdw[box_id][1]].Hide()
        self.textwidgets2a[text_vdw[box_id][1]].SetValue('')
        self.textwidgets2a[text_vdw[box_id][0]].Show()
      elif selection == 2:
        # here, we change the label to Cutoff (angs); hide the spline off widget & label,
        # and show the relevent logical and its label
        self.staticLbls2a[labels_vdw[box_id][1]].SetLabel('Cutoff (%s)' %angstrom)
        self.staticLbls2a[labels_vdw[box_id][1]].Show()
        self.staticLbls2a[labels_vdw[box_id][2]].Hide()
        self.staticLbls2a[labels_vdw[box_id][3]].Show()
        
        # show & hide the widgets now; always clear value of widgets that are being hidden
        self.textwidgets2a[text_vdw[box_id][1]].Hide()
        self.textwidgets2a[text_vdw[box_id][1]].SetValue('')
        self.choicewidgets2a[choice_vdw[box_id][1]].Show()
        self.textwidgets2a[text_vdw[box_id][0]].Show()

      elif selection == 3:
        # hide the logical and its label, and clear it; change 
        # the statics labels for the text widgets to
        # 'spline on' and 'spline off'; show the text widgets
        self.staticLbls2a[labels_vdw[box_id][1]].SetLabel('Spline On (%s)' %angstrom)
        self.staticLbls2a[labels_vdw[box_id][2]].SetLabel('Spline Off (%s)' %angstrom)
        self.staticLbls2a[labels_vdw[box_id][1]].Show()
        self.staticLbls2a[labels_vdw[box_id][2]].Show()
        
        # logical - hide, set to empty, hide label
        self.choicewidgets2a[choice_vdw[box_id][1]].Hide()
        self.choicewidgets2a[choice_vdw[box_id][1]].SetSelection(0)
        self.staticLbls2a[labels_vdw[box_id][3]].Hide()
 
        # show the text widgets
        for item in text_vdw[box_id]:
          self.textwidgets2a[item].Show() 
    
    elif "charge functional" in obj:
      # obtain the box id from the object name (passed in as obj)
      box_id = int(obj[4]) - 1
      if selection == 0 or selection == 2:
        # empty, or none; therefore, hide all the other box widgets and erase their values
        for item in labels_coul[box_id]:
          self.staticLbls2a[item].Hide()
      
        self.choicewidgets2a[choice_coul[box_id][0]].Hide()
        self.choicewidgets2a[choice_coul[box_id][0]].SetSelection(0)
        for item in text_coul[box_id]:
          self.textwidgets2a[item].Hide()
          self.textwidgets2a[item].SetValue('')
      elif selection == 1:
        # display the method widget and label for this box_id
        self.staticLbls2a[labels_coul[box_id][0]].Show()
        self.choicewidgets2a[choice_coul[box_id][0]].Show()
    elif "charge method" in obj:
      box_id = int(obj[4]) - 1
      if selection == 0:
        # hide the two text widgets and their labels for this box (and clear the widgets)
        for item in labels_coul[box_id][1:]:
          self.staticLbls2a[item].Hide()
        for item in text_coul[box_id]:
          self.textwidgets2a[item].Hide()
          self.textwidgets2a[item].SetValue('')
      elif selection == 1:
        # show the two text widgets and their labels for this box
        for item in labels_coul[box_id][1:]:
          self.staticLbls2a[item].Show()
        for item in text_coul[box_id]:
          self.textwidgets2a[item].Show()
      elif selection == 2:
        # show the 'cutoff' label and textwidget, but hide the 'accuracy' widget
        self.staticLbls2a[labels_coul[box_id][1]].Show()
        self.textwidgets2a[text_coul[box_id][0]].Show()
        self.staticLbls2a[labels_coul[box_id][2]].Hide()
        self.textwidgets2a[text_coul[box_id][1]].Hide()
        self.textwidgets2a[text_coul[box_id][1]].SetValue('')

    self.Layout()

  # our function called by global data() through pub.sendMessage at an appropriate time
  # serves to initialize the panel once parameters 'ensemble' and 'nbrSpecies' have been
  # determined. 
  def showAndHide(self,message,arg2=None):
    if 'ensemble' not in myDict.keys():
      return
    if 'nbrSpecies' not in myDict.keys():
      return

    if myDict['ensemble'][0:4] != "GEMC":
      for item in self.onlyShowGEMC_staticlbls:
        self.staticLbls2a[item].Hide()
      for item in self.onlyShowGEMC_choicewidgets:
        self.choicewidgets2a[item].Hide()
        self.choicewidgets2a[item].SetSelection(0) # re-initialize the selection
      for item in self.onlyShowGEMC_textwidgets:
        self.textwidgets2a[item].Hide()
        self.textwidgets2a[item].SetValue('') # erase the value
    else:
      # show the initial choice widgets with their labels
      #self.staticLbls2a[
      GEMC_init_display_lbls = [self.staticLbls2a[3], self.staticLbls2a[5],  \
                           self.staticLbls2a[20]]
      GEMC_init_display_choice =  [self.choicewidgets2a[2], self.choicewidgets2a[6]]
      for item in GEMC_init_display_lbls:
        item.Show()
      for item in GEMC_init_display_choice:
        item.Show()
    self.Layout()
 


class p2b(wx.Panel):
  
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid2b = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid2b)

    pub.subscribe(self.showAndHide, "refresh")

    self.lblList = ["Select a scaling style: ", "or enter custom values below.", \
                    "1-2 Scaling", "1-3 Scaling", \
                    "1-4 Scaling",   "1-N Scaling"] + \
                    ["van der Waals", "Coulombic"]*(maxnbrspecies) \
                    + ["Species %d" %(i+1) for i in range(maxnbrspecies)]
    flagStaticLbls = wx.ALIGN_CENTER
    lblCoords = [(1,1),    (2,1),     (4,3),     (4,4),
                 (4,5),    (4,6),     
                 (5,2),    (6,2),    (7,2),      (8,2),     
                 (9,2),    (10,2),   (11,2),     (12,2),   
                 (13,2),   (14,2),   (15,2),     (16,2), 
                 (5,1),    (7,1),    (9,1),      (11,1),
                 (13,1),   (15,1)] 
    spanListLbls = [(1,2)]*2 + [(1,1)]*4 + [(1,1)]*12 + [(1,1)]*maxnbrspecies

    self.staticLbls2b = []
    for index, item in enumerate(self.lblList):
      self.staticLbls2b.append(wx.StaticText(self,label=item))
      grid2b.Add(self.staticLbls2b[index],pos=lblCoords[index], \
                 flag=flagStaticLbls,span=spanListLbls[index])

    # create the names of the textwidgets
    strings_to_add = [" 1-2", " 1-3", " 1-4", " 1-N"]
    type_of_interaction = ["vdw","coul"]
    self.widgetnames2b = []
    for i in range(maxnbrspecies):
      for k in range(len(type_of_interaction)):
        for j in range(len(strings_to_add)):
          (self.widgetnames2b.append(("s%s" %(i+1)) + \
          (" %s" %(type_of_interaction[k]) + strings_to_add[j])))
 
    # create the list of coordinate tuples
    textCoords = []
    for i in range(2*maxnbrspecies):
      for j in range(len(strings_to_add)):
        textCoords.append((i+5,j+3))
    
    flagTextwidgets = (wx.ALIGN_CENTER | wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

    self.textwidgets2b = []
    for index, item in enumerate(self.widgetnames2b):
      self.textwidgets2b.append(wx.TextCtrl(self,value="",name=item))
      grid2b.Add(self.textwidgets2b[index],pos=textCoords[index],flag=flagTextwidgets)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets2b[index])

    # this list will be useful for changing data via checkboxes
    self.textwidgetsBySpecies = []
    for i in range(maxnbrspecies):
      self.textwidgetsBySpecies.append([])
    for item in self.textwidgets2b:
      thisName = item.GetName()
      indexThisWidget = int(thisName[1])
      self.textwidgetsBySpecies[indexThisWidget-1].append(item)

    # AMBER and CHARMM checkboxes below
    self.amberchk = wx.CheckBox(self,-1,"AMBER",name = "amber_param")
    self.charmmchk = wx.CheckBox(self,-1,"CHARMM",name="charmm_param")
    self.Bind(wx.EVT_CHECKBOX, self.datafunc, self.amberchk)
    self.Bind(wx.EVT_CHECKBOX, self.datafunc, self.charmmchk)
    grid2b.Add(self.amberchk,pos=(1,3))
    grid2b.Add(self.charmmchk,pos=(1,4))


  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())      
    val = str(event.GetString())
    if obj == "amber_param" or obj == "charmm_param":

      val = event.IsChecked()
      self.checkboxSelected(obj,val)

      return
    data(obj,val)

  def checkboxSelected(self,obj,val):
    if obj == "amber_param" and val == True:
      self.charmmchk.SetValue(False)
    elif obj == "charmm_param" and val == True:
      self.amberchk.SetValue(False)

    charmm_lj = ["0.0", "0.0", "0.0", "1.0"]
    charmm_elec = ["0.0", "0.0", "0.0", "1.0"]
    charmm_lj_elec = charmm_lj + charmm_elec #concatenate

    amber_lj = ["0.0","0.0","0.5","1.0"]
    amber_elec = ["0.0","0.0","0.8333","1.0"]
    amber_lj_elec = amber_lj + amber_elec 

    if (self.amberchk.GetValue()):
      for i in range(maxnbrspecies):
        for index, item in enumerate(self.textwidgetsBySpecies[i]):
          nameToPass = item.GetName()
          if i < int(myDict['nbrSpecies']):
            valToPass = amber_lj_elec[index]
          else:
            valToPass = ''
          item.SetValue(valToPass)
         # note that item.SetValue() causes its own event
         # this means that the value will be sent to the relevant widget for display,
         # and an event for that widget will be sent to the bound datafunc() method above
         # and therefore the value is sent to the global dictionary
    elif (self.charmmchk.GetValue()): 
      for i in range(maxnbrspecies):
        for index, item in enumerate(self.textwidgetsBySpecies[i]):
          nameToPass = item.GetName()
          if i < int(myDict['nbrSpecies']):
            valToPass = charmm_lj_elec[index]
          else:
            valToPass = ''
          item.SetValue(valToPass)
        
    else:
      for i in range(maxnbrspecies):
        for index, item in enumerate(self.textwidgetsBySpecies[i]):
          nameToPass = item.GetName()
          valToPass = ''
          item.SetValue(valToPass) 

  def showAndHide(self,message,arg2 = None):
    if 'nbrSpecies' not in myDict.keys():
      return

    # manually inspect label list to determine relevant indices;
    # expand using enumerate function
    # to easily determine label indices
    specLbls = self.staticLbls2b[18:]
    # 2*i since the vdw labels are alternating with the coulomb labels; 
    # first label located at index 6
    vdwLbls = [self.staticLbls2b[2*i+6] for i in range(maxnbrspecies)]
    # same reasoning
    coulLbls = [self.staticLbls2b[2*i+7] for i in range(maxnbrspecies)]

    # python counts the index starting at zero, therefore use a less than sign
    # the species labels
    for index, item in enumerate(specLbls):
      if index < int(myDict['nbrSpecies']):
        item.Show()
      else: item.Hide()
      
    # vdw labels
    for index, item in enumerate(vdwLbls):
      if index < int(myDict['nbrSpecies']):
        item.Show()
      else: item.Hide()
  
    # coulombic labels
    for index, item in enumerate(coulLbls):
      if index < int(myDict['nbrSpecies']):
        item.Show()
      else: item.Hide()

    # text widgets
    # the textwidgets contain the actual species number rather than python's index;
    # therefore, use a less than or equal to sign for sorting
    for index, item in enumerate(self.textwidgets2b):
      name = str(item.GetName())
      # note that the form of the widget's names are 's2 vdw 1-4', 's4 coul 1-2' etc.
      # splice the widget name to get what species it belongs to
      species = name[1] 
      if int(species) <= (int(myDict['nbrSpecies'])):
        item.Show()
      else: item.Hide()

    self.Layout()


class PanelTwo(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

    nested_notebook2 = wx.Notebook(self)
    self.p2a = p2a(nested_notebook2)
    self.p2b = p2b(nested_notebook2)

    nested_notebook2.AddPage(self.p2a,"Intermolecular")
    nested_notebook2.AddPage(self.p2b,"Intramolecular")

    nested_sizer = wx.BoxSizer()
    nested_sizer.Add(nested_notebook2,1,wx.EXPAND)
    self.SetSizer(nested_sizer)


####################################################################################

# Translation
class p3a(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3a = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3a)

    pub.subscribe(self.showAndHide, "refresh")

    # the static widgets
    lblList = ["Move Probability: ", "Box 1", "Box 2"] + \
              ["Species %s: " %(i+1) for i in range(maxnbrspecies)]
    flagList = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    lblCoords = [(0,2), (4,3), (4,4)] + [(i+5,2) for i in range(maxnbrspecies)]
    self.staticwidgets3a = []
    for index,item in enumerate(lblList):
      self.staticwidgets3a.append(wx.StaticText(self,label=item))
      grid3a.Add(self.staticwidgets3a[index],pos=lblCoords[index],flag=flagList)

    # probability of translation for species i in box j
    textwidgetNames = ["prob trans s%s b1" %(i+1) for i in range(maxnbrspecies)] + \
                      ["prob trans s%s b2" %(i+1) for \
                       i in range(maxnbrspecies)] + ["prob translation"]
    textwidgetCoords = [(i+5,3) for i in range(maxnbrspecies)] + [(i+5,4) \
                        for i in range(maxnbrspecies)] + [(0,3)]
    textwidgetFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.textwidgets3a = []
    for index,item in enumerate(textwidgetNames):
      self.textwidgets3a.append(wx.TextCtrl(self,value="",name=item))
      grid3a.Add(self.textwidgets3a[index],pos=textwidgetCoords[index],flag=textwidgetFlags)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets3a[index])

    # instructions and alert labels
    string3a1 = "Enter the maximum displacement " +\
                "(%s) allowed for each species in each box below." %(angstrom)
    self.instructionlbl = wx.StaticText(self,label=string3a1)
    grid3a.Add(self.instructionlbl,pos=(3,2),span=(1,6))

    string3a2 = "Please note that the sum of the move probabilities " +\
                "across all move types must sum to 1."
    self.noticelbl = wx.StaticText(self,label=string3a2)
    grid3a.Add(self.noticelbl,pos=(1,2),span=(1,6))


  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    data(obj,val)

  def showAndHide(self, message, arg2 = None):
    # return without doing anything if missing a required value. 
    if not 'ensemble' in myDict.keys():
      return
    if not 'nbrSpecies' in myDict.keys():
      return

    for i in range(maxnbrspecies):
      if i < int(myDict['nbrSpecies']):
        self.staticwidgets3a[i+3].Show()
        self.textwidgets3a[i].Show()
      else:
        self.staticwidgets3a[i+3].Hide()
        self.textwidgets3a[i].Hide()
        self.textwidgets3a[i].SetValue('')

    if myDict['ensemble'][0:4] == "GEMC":
      for i in range(maxnbrspecies):
        if i < int(myDict['nbrSpecies']):
          self.textwidgets3a[i+6].Show()
        else:
          self.textwidgets3a[i+6].Hide()
          self.textwidgets3a[i+6].SetValue('')
      self.staticwidgets3a[2].Show()
    else:
      for i in range(maxnbrspecies):
        self.textwidgets3a[i+6].Hide()
        self.textwidgets3a[i+6].SetValue('')
      self.staticwidgets3a[2].Hide()

    self.Layout()

# Rotation
class p3b(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3b = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3b)

    pub.subscribe(self.showAndHide, "refresh")

    # the static widgets
    lblList = ["Move Probability: ", "Box 1", "Box 2"] +\
              ["Species %s: " %(i+1) for i in range(maxnbrspecies)]
    flagList = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    lblCoords = [(0,2), (4,3), (4,4)] + [(i+5,2) for i in range(maxnbrspecies)]
    self.staticwidgets3b = []
    for index,item in enumerate(lblList):
      self.staticwidgets3b.append(wx.StaticText(self,label=item))
      grid3b.Add(self.staticwidgets3b[index],pos=lblCoords[index],flag=flagList)

    # probability of rotation for species i in box j
    textwidgetNames = ["prob rot s%s b1" %(i+1) for i in range(maxnbrspecies)] +\
                      ["prob rot s%s b2" %(i+1) for \
                       i in range(maxnbrspecies)] + ["prob rotation"]
    textwidgetCoords = [(i+5,3) for i in range(maxnbrspecies)] + \
                       [(i+5,4) for i in range(maxnbrspecies)] + [(0,3)]
    textwidgetFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.textwidgets3b = []
    for index,item in enumerate(textwidgetNames):
      self.textwidgets3b.append(wx.TextCtrl(self,value="",name=item))
      grid3b.Add(self.textwidgets3b[index],pos=textwidgetCoords[index],flag=textwidgetFlags)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets3b[index])

    # instructions and alert labels
    string3b1 = "Enter the maximum rotational width in " +\
                "degrees for each species in each box below."
    self.instructionlbl = wx.StaticText(self,label=string3b1)
    grid3b.Add(self.instructionlbl,pos=(3,2),span=(1,6))

    string3b2 = "Please note that the sum of the move probabilities " +\
                "across all move types must sum to 1."
    self.noticelbl = wx.StaticText(self,label=string3b2)
    grid3b.Add(self.noticelbl,pos=(1,2),span=(1,6))

  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    data(obj,val)

  def showAndHide(self, message, arg2 = None):
    # return without doing anything if missing a required value. 
    if not 'ensemble' in myDict.keys():
      return
    if not 'nbrSpecies' in myDict.keys():
      return

    for i in range(maxnbrspecies):
      if i < int(myDict['nbrSpecies']):
        self.staticwidgets3b[i+3].Show()
        self.textwidgets3b[i].Show()
      else:
        self.staticwidgets3b[i+3].Hide()
        self.textwidgets3b[i].Hide()
        self.textwidgets3b[i].SetValue('')

    if myDict['ensemble'][0:4] == "GEMC":
      for i in range(maxnbrspecies):
        if i < int(myDict['nbrSpecies']):
          self.textwidgets3b[i+6].Show()
        else:
          self.textwidgets3b[i+6].Hide()
          self.textwidgets3b[i+6].SetValue('')
      self.staticwidgets3b[2].Show()
    else:
      for i in range(maxnbrspecies):
        self.textwidgets3b[i+6].Hide()
        self.textwidgets3b[i+6].SetValue('')
      self.staticwidgets3b[2].Hide()

    self.Layout()

# Regrowth
class p3c(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3c = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3c)

    pub.subscribe(self.showAndHide, "refresh")

    # static labels
    lblList = ["Move Probability: "] + ["Species %s: " %(i+1) for i in range(maxnbrspecies)]
    flagList = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    lblCoords = [(0,2)] + [(i+5,2) for i in range(maxnbrspecies)]
    self.staticwidgets3c = []
    for index,item in enumerate(lblList):
      self.staticwidgets3c.append(wx.StaticText(self,label=item))
      grid3c.Add(self.staticwidgets3c[index],pos=lblCoords[index],flag=flagList)

    # text widgets
    textwidgetNames = ["prob regrowth s%s" %(i+1) for i in range(maxnbrspecies)] +\
                      ["prob regrowth"]
    textwidgetCoords = [(i+5,3) for i in range(maxnbrspecies)] + [(0,3)]
    textwidgetFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.textwidgets3c = []
    for index,item in enumerate(textwidgetNames):
      self.textwidgets3c.append(wx.TextCtrl(self,value="",name = item))
      grid3c.Add(self.textwidgets3c[index],\
                 pos=textwidgetCoords[index],flag=textwidgetFlags)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets3c[index])

    # instructions and alert labels
    string3c = ["Enter the relative probability of regrowth for each species below"]  + \
      ["Note that the relative probability below must sum to 1."] + \
      ["Please note that the sum of the move probabilities across all move types must sum to 1."]
    alertCoords = [(3,2), (4,2), (1,2)]
    spanList = [(1,6)]*3
    stringFlagList = wx.RESERVE_SPACE_EVEN_IF_HIDDEN
    for index,item in enumerate(string3c):
      self.staticwidgets3c.append(wx.StaticText(self,label=item))
      grid3c.Add(self.staticwidgets3c[-1],pos=alertCoords[index],\
                 span=spanList[index],flag=stringFlagList)


  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    data(obj,val)

  def showAndHide(self, message, arg2 = None):
    # return without doing anything if missing a required value. 
    if not 'ensemble' in myDict.keys():
      return
    if not 'nbrSpecies' in myDict.keys():
      return

    for i in range(maxnbrspecies):
      if i < int(myDict['nbrSpecies']):
        self.staticwidgets3c[i+1].Show()
        self.textwidgets3c[i].Show()
      else:
        self.staticwidgets3c[i+1].Hide()
        self.textwidgets3c[i].Hide()
        self.textwidgets3c[i].SetValue('')

    # when making any changes to which widgets are being shown/hidden, it is
    # important to call the Layout() method
    self.Layout()

# Volume
class p3d(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3d = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3d)

    pub.subscribe(self.showAndHide, "refresh")

    # static labels
    lblList = ["Move Probability: ", "Box 1: ", "Box 2: "]
    flagList = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    lblCoords = [(0,2), (5,2), (6,2)]
    self.staticwidgets3d = []
    for index,item in enumerate(lblList):
      self.staticwidgets3d.append(wx.StaticText(self,label=item))
      grid3d.Add(self.staticwidgets3d[index],pos=lblCoords[index],flag=flagList)

    # text widgets
    textwidgetNames = ["prob vol b1", "prob vol b2", "prob vol"]
    textwidgetCoords = [(5,3),(6,3),(0,3)]
    textwidgetFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.textwidgets3d = []
    for index,item in enumerate(textwidgetNames):
      self.textwidgets3d.append(wx.TextCtrl(self,value="",name=item))
      grid3d.Add(self.textwidgets3d[index],\
                 pos=textwidgetCoords[index],flag=textwidgetFlags)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets3d[index])


    self.string3d = ["Enter the maximum volume displacements in" +\
                     " %s^3 for the simulation box(es) below." %angstrom, \
                "This flag is required for NPT-MC, GEMC-NPT, and " +\
                "GEMC-NVT simulations, and may not be used for other simulation types.", \
                "Please note that the sum of the move probabilities " +\
                "across all move types must sum to 1."]
    stringCoords = [(3,2), (4,2), (1,2)]
    stringSpan = [(1,6), (1,8), (1,6)]
    stringFlags = (wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.stringwidgets3d = []
    for index, item in enumerate(self.string3d):
      self.stringwidgets3d.append(wx.StaticText(self,label=item))
      grid3d.Add(self.stringwidgets3d[index],pos=stringCoords[index],\
                 span = stringSpan[index], flag = stringFlags)

  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    data(obj,val)

  def showAndHide(self, message, arg2 = None):
    # return without doing anything if missing a required value. 
    if not 'ensemble' in myDict.keys():
      return
    if not 'nbrSpecies' in myDict.keys():
      return

    # here, Volume is a valid move only for NPT, GEMC simulations. Hide all widgets for others
    valid_ensembles = ["NPT_MC", "GEMC", "GEMC_NPT"]
    if myDict['ensemble'] not in valid_ensembles:
      for item in self.staticwidgets3d:
        item.Hide()
      for item in self.textwidgets3d:
        item.SetValue('')
        item.Hide()
      self.stringwidgets3d[0].Hide()
      self.stringwidgets3d[2].Hide()
    else:
      for item in self.staticwidgets3d:
        item.Show()
      for item in self.textwidgets3d:
        item.Show()
      for item in self.stringwidgets3d:
        item.Show()
      if myDict['ensemble'][0:4] != "GEMC":
        self.textwidgets3d[1].Hide()
        self.textwidgets3d[1].SetValue('')
        self.staticwidgets3d[2].Hide()
      else:
        self.textwidgets3d[1].Show()
        self.staticwidgets3d[2].Show()

    self.Layout()

# Insertion (& implicitly Deletion) Probabilities
class p3e(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3e = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3e)

    pub.subscribe(self.showAndHide, "refresh")

    lblList = ["Move Probability: "] + ["Species %s: " %(i+1) for i in range(maxnbrspecies)]
    flagList = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    lblCoords = [(0,2)] + [(i+5,2) for i in range(maxnbrspecies)]
    self.staticlbls3e = []
    for index, item in enumerate(lblList):
      self.staticlbls3e.append(wx.StaticText(self,label=item))
      grid3e.Add(self.staticlbls3e[index],pos=lblCoords[index],flag=flagList)

    # only one textwidget on this panel - asking for the probability of insertion move
    self.insprobquery = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name = "prob insertion")
    grid3e.Add(self.insprobquery,pos=(0,3),flag=flagList)
    self.Bind(wx.EVT_TEXT, self.datafunc, self.insprobquery)

    # need to ask for the insertion method for 
    # each species in simulation; use choice widgets to do so
    choiceOptions = ["", "cbmc","none"]
    self.choicewidgets3e = []
    choiceNames = ["insertion method s%s" %(i+1) for i in range(maxnbrspecies)]
    choiceCoords = [(i+5,3) for i in range(maxnbrspecies)]
    choiceFlags = (wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.ALIGN_CENTER)
    for index,item in enumerate(choiceNames):
      self.choicewidgets3e.append(wx.Choice(self,choices=choiceOptions,name=item))
      grid3e.Add(self.choicewidgets3e[index],pos=choiceCoords[index],flag=choiceFlags)
      self.Bind(wx.EVT_CHOICE, self.datafunc, self.choicewidgets3e[index])

    # some strings providing information about the panel
    string3e = ["Please note that the sum of the move " +\
                "probabilities across all move types must sum to 1.", \
      "Additionally, insertion moves define an equal probability of deletion," +\
      " and so this probability", \
      "should be counted twice when summing to 1.", \
      "This flag is allowed only for GCMC simulations."]
    stringCoords = [(1,2), (2,2), (3,2), (8,4)]
    stringSpan = [(1,6), (1,6), (1,6), (1,6)]
    stringFlags = (wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.ALIGN_LEFT)
    self.stringwidgets3e = []
    for index, item in enumerate(string3e):
      self.stringwidgets3e.append(wx.StaticText(self,label=item))
      grid3e.Add(self.stringwidgets3e[index],pos=stringCoords[index],\
                 span=stringSpan[index],flag=stringFlags)

  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    data(obj,val)

  def showAndHide(self, message, arg2 = None):
    # return without doing anything if missing a required value. 
    if not 'ensemble' in myDict.keys():
      return
    if not 'nbrSpecies' in myDict.keys():
      return

    if myDict['ensemble'] != "GCMC":
      for item in self.staticlbls3e:
        item.Hide()
      for item in self.choicewidgets3e:
        item.Hide()
        item.SetSelection(0)
      self.stringwidgets3e[0].Hide()
      self.stringwidgets3e[1].Hide()
      self.stringwidgets3e[2].Hide()
      self.insprobquery.Hide()
      self.insprobquery.SetValue('')
      self.stringwidgets3e[3].Show()

    else:
      for i in range(maxnbrspecies):
        if i < int(myDict['nbrSpecies']):
          self.staticlbls3e[i+1].Show()
          self.choicewidgets3e[i].Show()
        else:
          self.staticlbls3e[i+1].Hide()
          self.choicewidgets3e[i].Hide()
          self.choicewidgets3e[i].SetSelection(0)
      self.staticlbls3e[0].Show()
      self.stringwidgets3e[0].Show()
      self.stringwidgets3e[1].Show()
      self.stringwidgets3e[2].Show()
      self.insprobquery.Show()
      self.stringwidgets3e[3].Hide()

    self.Layout()


# Swap probabilities
class p3f(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid3f = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid3f)

    pub.subscribe(self.showAndHide, "refresh")

    lblList = ["Move Probability: "] + ["Species %s: " %(i+1) for i in range(maxnbrspecies)]
    flagList = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    lblCoords = [(0,2)] + [(i+5,2) for i in range(maxnbrspecies)]
    self.staticlbls3f = []
    for index, item in enumerate(lblList):
      self.staticlbls3f.append(wx.StaticText(self,label=item))
      grid3f.Add(self.staticlbls3f[index], pos=lblCoords[index],flag = flagList)

    # only one textwidget on this panel - asking for the probability of swap move
    self.swapprobquery = wx.TextCtrl(self,value="",style=wx.ALIGN_LEFT,name="prob swap")
    grid3f.Add(self.swapprobquery,pos=(0,3),flag=flagList)
    self.Bind(wx.EVT_TEXT, self.datafunc, self.swapprobquery)

    # ask for swap method for each species in simulation; use choice widgets to do so
    choiceOptions = ["", "reservoir", "none"]
    self.choicewidgets3f = []
    choiceNames = ["swap method s%s" %(i+1) for i in range(maxnbrspecies)]
    choiceCoords = [(i+5,3) for i in range(maxnbrspecies)]
    choiceFlags = (wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.ALIGN_CENTER)
    for index, item in enumerate(choiceNames):
      self.choicewidgets3f.append(wx.Choice(self,choices=choiceOptions,name=item))
      grid3f.Add(self.choicewidgets3f[index],pos=choiceCoords[index],flag=choiceFlags)
      self.Bind(wx.EVT_CHOICE, self.datafunc, self.choicewidgets3f[index])

    # strings to inform user what this panel is for
    string3f = ["Please note that the sum of the move probabilities " +\
                "across all move types must sum to 1.", \
                "This flag is allowed only for GEMC simulations.", \
                "Select the swap method for each relevant species in the simulation below."]
    stringCoords = [(1,2), (8,4), (3,2)]
    stringSpan = [(1,6), (1,6), (1,6)]
    stringFlags = (wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.stringwidgets3f = []
    for index, item in enumerate(string3f):
      self.stringwidgets3f.append(wx.StaticText(self,label=item))
      grid3f.Add(self.stringwidgets3f[index], pos=stringCoords[index], \
                 span=stringSpan[index], flag=stringFlags)

  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    data(obj,val)

  def showAndHide(self, message, arg2 = None):
    # return without doing anything if missing a required value. 
    if not 'ensemble' in myDict.keys():
      return
    if not 'nbrSpecies' in myDict.keys():
      return

    if myDict['ensemble'][0:4] != "GEMC":
      for item in self.staticlbls3f:
        item.Hide()
      for item in self.choicewidgets3f:
        item.Hide()
      for item in self.stringwidgets3f:
        item.Hide()
      self.swapprobquery.Hide()
      self.swapprobquery.SetValue('')
      self.stringwidgets3f[1].Show()
    else:
      for i in range(maxnbrspecies):
        if i < int(myDict['nbrSpecies']):
          self.staticlbls3f[i+1].Show()
          self.choicewidgets3f[i].Show()
        else:
          self.staticlbls3f[i+1].Hide()
          self.choicewidgets3f[i].Hide()
          self.choicewidgets3f[i].SetSelection(0)
      for item in self.stringwidgets3f:
        item.Show()
      self.stringwidgets3f[1].Hide()
      self.swapprobquery.Show()
      self.staticlbls3f[0].Show()

    self.Layout()


class PanelThree(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

    nest = wx.Notebook(self)
    # self. notation makes the panels accessible at the global level
    self.p3a = p3a(nest)
    self.p3b = p3b(nest)
    self.p3c = p3c(nest)
    self.p3d = p3d(nest)
    self.p3e = p3e(nest)
    self.p3f = p3f(nest)

    nest.AddPage(self.p3a,"Translation")
    nest.AddPage(self.p3b,"Rotation")
    nest.AddPage(self.p3c,"Regrowth")
    nest.AddPage(self.p3d,"Volume")
    nest.AddPage(self.p3e,"Insertion")
    nest.AddPage(self.p3f,"Swap")


    nested_sizer = wx.BoxSizer()
    nested_sizer.Add(nest,1,wx.EXPAND)
    self.SetSizer(nested_sizer)


# molecule files
class p4a(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid4a = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid4a)

    pub.subscribe(self.showAndHide, "refresh")

    # labels for the interactive widgets
    lblList = ["Select","                     Selection                     ", \
               "# Molecules",] + \
              ["Species %s: " %(i+1) for i in range(maxnbrspecies)]
    flagList = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    lblCoords = [(3,2), (3,3), (3,4)] + [(i+4,1) for i in range(maxnbrspecies)]
    self.staticlbls4a = []
    for index, item in enumerate(lblList):
      self.staticlbls4a.append(wx.StaticText(self,label=item))
      grid4a.Add(self.staticlbls4a[index],pos=lblCoords[index],flag=flagList)

    # buttons to select MCF file for appropriate species
    self.buttonNames = ["MCF s%s" %(i+1) for i in range(maxnbrspecies)]
    buttonLabel = "Select MCF File"
    buttonCoords = [(i+4,2) for i in range(maxnbrspecies)]
    buttonFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.EXPAND)
    self.buttonwidgets4a = []
    for index, item in enumerate(self.buttonNames):
      self.buttonwidgets4a.append(wx.Button(self,label=buttonLabel,name = item))
      grid4a.Add(self.buttonwidgets4a[index],pos=buttonCoords[index],flag=buttonFlags)
      self.Bind(wx.EVT_BUTTON, self.datafunc, self.buttonwidgets4a[index])

    self.textNames = ["max nmols s%s" %(i+1) for i in range(maxnbrspecies)]
    textCoords = [(i+4,4) for i in range(maxnbrspecies)]
    textFlags = (wx.EXPAND | wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.textwidgets4a = []
    for index, item in enumerate(self.textNames):
      self.textwidgets4a.append(wx.TextCtrl(self,value="",name=item))
      grid4a.Add(self.textwidgets4a[index],pos=textCoords[index],flag=textFlags)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets4a[index])

    # some read only text widgets to display the user's selected mcf
    dispCoords = [(i+4,3) for i in range(maxnbrspecies)]
    dispFlags = (wx.EXPAND | wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.dispwidgets4a = []
    for index, item in enumerate(dispCoords):
      self.dispwidgets4a.append(wx.TextCtrl(self,value="",style = wx.TE_READONLY))
      grid4a.Add(self.dispwidgets4a[index],pos=item,flag=dispFlags)

    # instructional label
    string4a = ["Select the MCF files below.  Enter the maximum number" +\
                " of anticipated molecules for each species.", \
                "This number will be used for memory allocation purposes only."]
    stringCoords = [(1,1), (2,1)]
    stringSpan = [(1,5), (1,5)]
    stringFlags = (wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.stringwidgets4a = []
    for index, item in enumerate(string4a):
      self.stringwidgets4a.append(wx.StaticText(self,label=item))
      grid4a.Add(self.stringwidgets4a[index],pos=stringCoords[index],\
                 span=stringSpan[index],flag=stringFlags)

  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())

    # check if this is an MCF button
    if "MCF" in obj:
      val = str(self.button(obj))
      dispThisName = os.path.split(val)
      if val:
        # obj[-1] since the species number is the last character of the object name
        # int(obj[-1])-1 since python starts with 0 - therefore, it is matched up with self.dispwidgets4a
        self.dispwidgets4a[int(obj[-1])-1].SetValue(("/%s/" %dispThisName[1]))
    data(obj,val) 

  def button(self,obj):
    dlg = wx.FileDialog(self,"Select MCF File", myDict['simDir'],"","*.mcf",wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
      val = dlg.GetPath()
      f = file(val)
      nfrags_ind = ""
      nfrags_data = ""
      for ind, line in enumerate(f):
        if '#' in line and 'Fragment_Info' in line:
          nfrags_ind = ind+1
        if nfrags_ind:
          if ind == nfrags_ind:
            line_data = line.split()
            nfrags_data = int(line_data[0])
      if not nfrags_data:
        print "Number of fragments could not be identified for this species."
        nfrags_data = 0
      file_data = os.path.relpath(val,str(myDict['simDir']))
      val = file_data

      # pass the expected number of fragments in to the data function:
      # first, retrieve the species number from the original object
      thisSpeciesNum = int(obj[-1])
      nameOfFragData = "nfrags expected s%s" %(thisSpeciesNum)
      # store the data we obtained for the expected number of fragments
      data(nameOfFragData,nfrags_data)

      # now, return the relative path to the mcf file that we came here for.
      return val      
    else:
      val = ''
      return val

      val = str(self.button(obj))
      self.dirname = os.path.split(val)
      if val:
        self.simdirdisp.SetValue(("/%s/" %self.dirname[1]))
    data(obj,val) 

  def showAndHide(self, message, arg2 = None):
    # return without doing anything if missing a required value. 
    if not 'ensemble' in myDict.keys():
      return
    if not 'nbrSpecies' in myDict.keys():
      return

    for i in range(maxnbrspecies):
      if i < int(myDict['nbrSpecies']):
        self.textwidgets4a[i].Show()
        self.staticlbls4a[i+3].Show()
        self.buttonwidgets4a[i].Show()
        self.dispwidgets4a[i].Show()
      else:
        self.textwidgets4a[i].Hide()
        self.textwidgets4a[i].SetValue('')
        self.staticlbls4a[i+3].Hide()
        self.buttonwidgets4a[i].Hide() 
        self.dispwidgets4a[i].Hide()
        self.dispwidgets4a[i].SetValue('')

        # buttons do not have SetValue attribute; call data(objName, '') to erase from dict
        thisName = self.buttonwidgets4a[i].GetName()
        data(thisName,'')

    self.Layout()

# Fragment Files
class p4b(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid4b = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid4b)

    pub.subscribe(self.showAndHide, "refresh")

    lblList = ["Are your fragment files prepared?"] + \
              ["Select the fragment files for each species below."] +\
              ["After creating your input file, consult the user guide" +\
               " regarding creation of your fragment files."]
   
    lblCoords = [(1,1)] + [(3,1)] + [(2,1)]
    lblSpan = [(1,2)] + [(1,3)] + [(1,6)]
    lblFlags = [(wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)]*3 

    self.staticlbls4b = []
    for index,item in enumerate(lblList):
      self.staticlbls4b.append(wx.StaticText(self,label=item))
      grid4b.Add(self.staticlbls4b[index],pos=lblCoords[index],\
                 span=lblSpan[index],flag=lblFlags[index])

    optionsPrepared = ["","Yes","No"]
    self.areFragFilesPrepared = wx.Choice(self,-1,name='fragFilesCheck',\
                              choices = optionsPrepared,style = wx.ALIGN_CENTER)
    self.Bind(wx.EVT_CHOICE, self.datafunc, self.areFragFilesPrepared)
    grid4b.Add(self.areFragFilesPrepared,pos=(1,3))
    
    self.speciesOptions = [""] + ["%s" %(i+1) for i in range(maxnbrspecies)]
    self.chooseFragSpecies = wx.Choice(self,-1,name='fragChoice',\
                             choices = self.speciesOptions,style=wx.ALIGN_CENTER)
    self.Bind(wx.EVT_CHOICE, self.datafunc,self.chooseFragSpecies)
    grid4b.Add(self.chooseFragSpecies,pos=(6,2))
    
    speciesString = "Select a Species:"
    self.speciesLabel=wx.StaticText(self,label=speciesString)
    grid4b.Add(self.speciesLabel,pos=(6,1))

    self.fragFileBtn = wx.Button(self,label="Select Fragment File",name='fragBtn')
    grid4b.Add(self.fragFileBtn,pos=(6,3))
    self.Bind(wx.EVT_BUTTON,self.fragSelect,self.fragFileBtn)
    
    # multi line textctrl to display the current fragment list of the 
    # selected species
    self.dispVals = wx.TextCtrl(self,value='',\
                    style=(wx.TE_MULTILINE | wx.TE_READONLY),size=(300,300))
    grid4b.Add(self.dispVals,pos=(7,3),span=(4,4))
    self.Bind(wx.EVT_TEXT, self.te_multiFunc, self.dispVals)

    self.thisFragIndex = ''

    # initially, we hide all widgets except for the question label and yes/no choice.
    hideThese = [self.staticlbls4b[1],self.staticlbls4b[2],self.chooseFragSpecies,\
                 self.speciesLabel, self.fragFileBtn,  self.dispVals]
    for index, item in enumerate(hideThese):
      item.Hide()

  # should we use this? keep for
  def te_multiFunc(self,event):
    pass

  def showAndHide(self,message,arg2=None):
    if 'nbrSpecies' not in myDict:
      return
    # in this class, we don't use this to show/hide anything
    # instead, use to clear and append relevant strings to the choice control
    # put the show/hide behavior inside of refreshThisPanel()
    # and initially hide all except the top question and choice widget
    stringsToUse = [""] + ["%s" %(i+1) for i in range(int(myDict['nbrSpecies']))]
    self.chooseFragSpecies.Clear()
    self.chooseFragSpecies.AppendItems(stringsToUse)
    self.Layout()

  def fragSelect(self,event):
    dlg = wx.FileDialog(self,"Select Fragment File", myDict['simDir'],"","*.dat",wx.OPEN)
    val = ''
    if dlg.ShowModal() == wx.ID_OK:
      thisPath = str(dlg.GetPath())
      val = str(os.path.relpath(thisPath, myDict['simDir']))
    if self.thisFragIndex:
      thisSpecies = str(self.thisFragIndex)
      thisName = 'frag files s%s' %(thisSpecies)
  # we want to store our files for each species as a list
  # because that is how it is written in the save function for ease of use
  # the enumerate function was used because parsing through the dictionary list
  # using 'for item in []' caused the strings to split in to individual letters (bad)
    if thisName in myDict:
      thisSpeciesFragList = []
      for index, item in enumerate(myDict[thisName]):   
        thisSpeciesFragList.append(str(item))
      if val:
        thisSpeciesFragList.append(val)
        data(thisName,thisSpeciesFragList)
      self.dispVals.SetValue(str(thisSpeciesFragList))
    else:
      thisSpeciesFragList = []
      if val:
        thisSpeciesFragList.append(val)
        data(thisName,thisSpeciesFragList)
      self.dispVals.SetValue(str(thisSpeciesFragList))

  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    if obj == 'fragFilesCheck':
      self.refreshThisPanel(obj,val)
    elif obj == 'fragChoice':
      # refresh the panel - we need to hide the button if not val
      self.refreshThisPanel(obj,val)
      # create an object to store the current species
      # allows for communication with the text control
    data(obj,val)

  # this function receives the object reporting and the value
  # from here, we decide what to show and hide
  def refreshThisPanel(self,obj,val):
    hideThese = [self.staticlbls4b[1],self.staticlbls4b[2],self.chooseFragSpecies,\
                 self.speciesLabel, self.fragFileBtn,  self.dispVals]
    if obj == 'fragFilesCheck':
      if not val:
        # if val is empty, hide everything. don't bother erasing data though. 
        for index, item in enumerate(hideThese):
          item.Hide()
      elif val == "No":
        for index, item in enumerate(hideThese):
          item.Hide()
        self.staticlbls4b[2].Show()
      elif val == "Yes":
        self.staticlbls4b[2].Hide()
        self.staticlbls4b[1].Show()
        self.chooseFragSpecies.Show()
        self.speciesLabel.Show()
    elif obj == 'fragChoice': 
      if not val:
        self.fragFileBtn.Hide()
        self.dispVals.Hide()
      else:
        self.thisFragIndex = val
        thisSpecies = str(self.thisFragIndex)
        thisName = 'frag files s%s' %(thisSpecies)

        self.fragFileBtn.Show()
        self.dispVals.Show()
        if thisName in myDict:
          self.dispVals.SetValue(str(myDict[thisName]))
        else:
          self.dispVals.SetValue('')
    self.Layout() 


# Input File
class p4c(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid4c = wx.GridBagSizer(hgap=5, vgap=5)
    self.SetSizer(grid4c)

    pub.subscribe(self.showAndHide, "refresh")

    # labels - static text widgets
    lblList = ["Provide information regarding the Start Type and Run Type for the simulation below.", \
               "Start Type: ", "Box 1", "Box 2"] + \
              ["# Molecules Species %s: " %(i+1) for i in range(maxnbrspecies)] +\
              ["Run Type: ", "Frequency at which Acceptance" +\
               " Ratios are output to Log File (MC Steps): ", \
               "Frequency at which Volume Acceptance Ratios are output to Log File (MC Steps): "]

    lblCoords = [(1,1), (3,1), (3,5), (3,6),] + \
                [(i+4,4) for i in range(maxnbrspecies)] + [(10,1), (11,1), (12,1)]

    lblSpan = [(1,6), (1,1), (1,1), (1,1)] + \
              [(1,1)]*maxnbrspecies + [(1,1), (1,4), (1,4)]

    lblFlags = [(wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)] + \
               [(wx.ALIGN_RIGHT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)] + \
               [(wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)]*2 + \
               [wx.ALIGN_RIGHT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN]*(len(lblList)-4)
    self.staticlbls4c = []
    for index, item in enumerate(lblList):
      self.staticlbls4c.append(wx.StaticText(self,label=item))
      grid4c.Add(self.staticlbls4c[index], pos=lblCoords[index], \
                 span=lblSpan[index], flag = lblFlags[index])
      
    # text widgets
    self.textwidgetnames = ["makeConfig s%s b%s" %(i+1,j+1) \
                            for j in range(maxnbrboxes) for i in range(maxnbrspecies)] + \
                           ["acceptance ratio output", "volume displacement output"]
    textwidgetCoords = [(i+4,5+j) for j in range(maxnbrboxes) \
                         for i in range(maxnbrspecies)] + [(11,5), (12,5)]
    textwidgetFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.textwidgets4c = []
    for index, item in enumerate(self.textwidgetnames):
      self.textwidgets4c.append(wx.TextCtrl(self,value="",name=item))
      grid4c.Add(self.textwidgets4c[index],pos=textwidgetCoords[index],flag=textwidgetFlags)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets4c[index])
      
    # choice widgets
    runTypeOptions = ["", "Equilibration", "Production"]
    startTypeOptions = ["", "make_config", "checkpoint", "read_config"]

    choiceOptionsList = [startTypeOptions, runTypeOptions]
    self.choicewidgetnames = ["start type","run type"]
    choiceCoords = [(3,2), (10,2)]
    choiceFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    self.choicewidgets4c = []
    for index, item in enumerate(self.choicewidgetnames):
      self.choicewidgets4c.append(wx.Choice(self,choices=choiceOptionsList[index],name=item))
      grid4c.Add(self.choicewidgets4c[index],pos=choiceCoords[index], flag = choiceFlags)
      self.Bind(wx.EVT_CHOICE, self.datafunc, self.choicewidgets4c[index])

    # button widgets
    self.buttonwidgetnames = ["checkpoint file", "read old box 1", "read old box 2"]
    buttonwidgetLabels = ["Select Checkpoint File",\
                          "Select Box 1 Read Config", "Select Box 2 Read Config"]
    buttonwidgetCoords = [(4,1), (5,1), (6,1)]
    buttonwidgetSpan = (1,2)
    buttonwidgetFlags = (wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.EXPAND)
    self.buttonwidgets4c = []
    for index, item in enumerate(self.buttonwidgetnames):
      self.buttonwidgets4c.append(wx.Button(self,label=buttonwidgetLabels[index],name = item))
      grid4c.Add(self.buttonwidgets4c[index],pos=buttonwidgetCoords[index],\
                 span=buttonwidgetSpan,flag = buttonwidgetFlags)
      self.Bind(wx.EVT_BUTTON, self.datafunc, self.buttonwidgets4c[index])

    # text controls to display selection from buttons
    self.textdispnames = ["checkpoint disp", "read old b1 disp", "read old b2 disp"]
    textdispCoords = [(4,3), (5,3), (6,3)]
    textdispFlags = (wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN | wx.EXPAND)
    self.textdisp4c = []
    for index, item in enumerate(self.textdispnames):
      self.textdisp4c.append(wx.TextCtrl(self,value="",name=item,style=wx.TE_READONLY))
      grid4c.Add(self.textdisp4c[index],pos=textdispCoords[index],flag = textdispFlags)

    # some widgets will initially be hidden.  These include:
    # make config widgets; read old widgets; checkpoint widgets
    self.lbls_init_hide = [2,3] + [(i+4) for i in range(maxnbrspecies)]
    self.text_init_hide = [i for i in range(maxnbrspecies*2)]
    self.button_init_hide = [0,1,2]
    self.disp_init_hide = [0,1,2]

    for item in self.lbls_init_hide:
      self.staticlbls4c[item].Hide()
    for item in self.text_init_hide:
      self.textwidgets4c[item].Hide()
      self.textwidgets4c[item].SetValue('')
    for item in self.button_init_hide:
      self.buttonwidgets4c[item].Hide()
      thisWidgetName = self.buttonwidgets4c[item].GetName()
      data(thisWidgetName,'')
    for item in self.disp_init_hide:
      self.textdisp4c[item].Hide()
      self.textdisp4c[item].SetValue('')

    # also, hide the widgets for run type other than the choice widget and label
    self.staticlbls4c[11].Hide()
    self.staticlbls4c[12].Hide()
    self.textwidgets4c[12].Hide()
    self.textwidgets4c[13].Hide()


  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    if obj == "start type":
      self.startTypeChange(val)
    elif obj == 'checkpoint file':
      ext = 'chk'
      self.buttonClick(obj,val,ext)
      return
    elif 'read old' in obj:
      ext = 'xyz'
      self.buttonClick(obj,val,ext)
      return
    elif obj == "run type":
      self.runTypeChange(val)
    data(obj,val)


  def buttonClick(self,obj,val,ext):
    dlg = wx.FileDialog(self,'Select File',myDict['simDir'],'','*.%s' %(ext),wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
      thisFilePath = dlg.GetPath()
      thisFilePath = os.path.relpath(thisFilePath,myDict['simDir'])
      if obj == 'checkpoint file':
        self.textdisp4c[0].SetValue(thisFilePath)
      elif obj == 'read old box 1':
        self.textdisp4c[1].SetValue(thisFilePath)
      elif obj == 'read old box 2':
        self.textdisp4c[2].SetValue(thisFilePath)
      val = thisFilePath
      data(obj,val)
      return obj, val

  def runTypeChange(self,val):
    # some show and hide features with a label change
    if val == "Equilibration":
      # then change the labels around. check what ensemble we're in
      equilMov = "Frequency at which Maximum Displacement Width Updated (MCsteps):"
      equilVol = "Frequency at which Volume Displacement Width Updated (MCsteps):"
      if myDict['ensemble'] == "GCMC" or myDict['ensemble'][0:3] == "NVT":
        # hide the volume lbl and widget regardless
        self.staticlbls4c[12].Hide()
        self.textwidgets4c[13].Hide()
        
        self.staticlbls4c[11].SetLabel(equilMov)
        self.staticlbls4c[11].Show()
        self.textwidgets4c[12].Show()
      else:
        self.staticlbls4c[12].SetLabel(equilVol)
        self.staticlbls4c[12].Show()
        self.textwidgets4c[13].Show()
    
        self.staticlbls4c[11].SetLabel(equilMov)
        self.staticlbls4c[11].Show()  
        self.textwidgets4c[12].Show()

    elif val == "Production":
      prodMov = "Frequency at which Acceptance Ratios are Output to Log File (MCsteps):"
      prodVol = "Frequency at which Volume Acceptance Ratios are Output to Log File (MCsteps):"
      if myDict['ensemble'] == "GCMC" or myDict['ensemble'][0:3] == "NVT":
        self.staticlbls4c[12].Hide()
        self.textwidgets4c[13].Hide()
    
        self.staticlbls4c[11].SetLabel(prodMov)
        self.staticlbls4c[11].Show()
        self.textwidgets4c[12].Show()
      else:
        self.staticlbls4c[12].SetLabel(prodVol)
        self.staticlbls4c[12].Show()
        self.textwidgets4c[13].Show()
  
        self.staticlbls4c[11].SetLabel(prodMov)
        self.staticlbls4c[11].Show()
        self.textwidgets4c[12].Show()
    else:
      # empty string, hide everything
      self.staticlbls4c[11].Hide()
      self.staticlbls4c[12].Hide()
      self.textwidgets4c[12].Hide()
      self.textwidgets4c[13].Hide()

    self.Layout()


  def startTypeChange(self,val):
    if not val:
      for item in self.lbls_init_hide:
        self.staticlbls4c[item].Hide()
      for item in self.text_init_hide:
        self.textwidgets4c[item].Hide()
        self.textwidgets4c[item].SetValue('')
      for item in self.button_init_hide:
        self.buttonwidgets4c[item].Hide()
        thisWidgetName = self.buttonwidgets4c[item].GetName()
        data(thisWidgetName,'')
      for item in self.disp_init_hide:
        self.textdisp4c[item].Hide()
        self.textdisp4c[item].SetValue('')
    elif val == 'make_config' or val == 'read_config':
      # hide the widgets not pertaining to make config
      for item in self.disp_init_hide:
        self.textdisp4c[item].Hide()
        self.textdisp4c[item].SetValue('')
      for item in self.button_init_hide:
        self.buttonwidgets4c[item].Hide()
        thisWidgetName = self.buttonwidgets4c[item].GetName()
        data(thisWidgetName,'')
      self.staticlbls4c[2].Show()
      for i in range(maxnbrspecies):
        if i < int(myDict['nbrSpecies']):
          self.staticlbls4c[i+4].Show()
          self.textwidgets4c[i].Show()
        else:
          self.staticlbls4c[i+4].Hide()
          self.textwidgets4c[i].Hide()
          self.textwidgets4c[i].SetValue('')
      if myDict['ensemble'][0:4] == "GEMC":
        self.staticlbls4c[3].Show()
        for i in range(maxnbrspecies):
          if i < int(myDict['nbrSpecies']):
            self.textwidgets4c[i+6].Show()
          else:
            self.textwidgets4c[i+6].Hide()
            self.textwidgets4c[i+6].SetValue('')
      else:
        for i in range(maxnbrspecies):
          self.textwidgets4c[i+6].Hide()
          self.textwidgets4c[i+6].SetValue('')
        self.staticlbls4c[3].Hide()
      if val == 'read_config':
        self.buttonwidgets4c[0].Hide(),self.textdisp4c[0].Hide()
        self.buttonwidgets4c[1].Show(), self.textdisp4c[1].Show()
        if myDict['ensemble'][0:4] == "GEMC":
          self.buttonwidgets4c[2].Show(), self.textdisp4c[2].Show()
        else:
          self.buttonwidgets4c[2].Hide(), self.textdisp4c[2].Hide()
          thisName = self.buttonwidgets4c[2].GetName()
          data(thisName,'')
    elif val == 'checkpoint':
      # hide and clear the make config widgets
      for item in self.lbls_init_hide:
        self.staticlbls4c[item].Hide()
      for item in self.text_init_hide:
        self.textwidgets4c[item].Hide()
        self.textwidgets4c[item].SetValue('')
      # hide and clear read old widgets
      self.buttonwidgets4c[1].Hide(), self.textdisp4c[1].Hide()
      self.buttonwidgets4c[2].Hide(), self.textdisp4c[2].Hide()
      thisName = self.buttonwidgets4c[1].GetName()
      data(thisName,'')
      thisName = self.buttonwidgets4c[2].GetName()
      data(thisName,'')

      # show the checkpoint widgets
      self.buttonwidgets4c[0].Show(),self.textdisp4c[0].Show()
      self.buttonwidgets4c[1].Hide()
      self.buttonwidgets4c[1].Hide()
      self.textdisp4c[1].Hide(), self.textdisp4c[2].Hide()
      thisName = self.buttonwidgets4c[1].GetName()
      data(thisName,'')
      thisName = self.buttonwidgets4c[2].GetName()
      data(thisName,'')
      # hide and clear the make config widgets
      # hide and clear checkpoint widgets
      # show the read old widget(s) - only box 2 if ensemble is GEMC

    # then done. except need to program the buttons too.
    self.Layout()

  def showAndHide(self,message,arg2=None):
    # make an event in the nbr of species or ensemble reset the start type
    try:
      thisSelect = self.choicewidgets4c[0].GetString()
      self.choicewidgets4c[0].SetString(thisSelect)
      self.startTypeChange(str(thisSelect))
      self.Layout()
    except TypeError:
      pass


# Output File
class p4d(wx.Panel):
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)
    grid4d = wx.GridBagSizer(hgap=5,vgap=5)
    self.SetSizer(grid4d)

    pub.subscribe(self.showAndHide, "refresh")

    self.lblList = ["Units: ", "MC Steps", "MC Steps",\
                    "MC Steps", "Average Information: ", \
                    "Desired Property Outputs", "Box 1", "Box 2",\
                    "Thermodynamic Output Frequency",
                    "Coordinate Output Frequency", "Simulation Run Length"]
    lblCoords = [(2,1), (3,4), (4,4), (5,4), (7,1), (9,2), (10,2), (10,4), (3,1), (4,1), (5,1)]
    lblFlags = [(wx.ALIGN_RIGHT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)] + \
               [(wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)]*3 + \
               [(wx.ALIGN_RIGHT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)] + \
               [(wx.ALIGN_CENTER | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)]*3 + \
               [(wx.ALIGN_RIGHT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)]*3
    lblSpan = [(1,1)]*5 + [(1,3)] + [(1,1)]*2 + [(1,2)]*3
    self.staticlbls4d = []
    for index, item in enumerate(self.lblList):
      self.staticlbls4d.append(wx.StaticText(self,label=item))
      grid4d.Add(self.staticlbls4d[index],pos=lblCoords[index],\
                 span=lblSpan[index], flag=lblFlags[index])

    choiceOptions = ["", "Minutes", "Steps"]
    self.freqtypechoice = wx.Choice(self,-1,choices=choiceOptions,name="freqTypeChoice")
    grid4d.Add(self.freqtypechoice,pos=(2,2))
    self.Bind(wx.EVT_CHOICE, self.datafunc,self.freqtypechoice)

    self.textwidgets4d = []
    textCoords = [(3,3), (4,3), (5,3), (7,2)]
    textFlags = (wx.EXPAND | wx.RESERVE_SPACE_EVEN_IF_HIDDEN)
    textNames = ["thermoFreq", "coordFreq", "simRunTime", "avgInfo"]
    styleList = [(wx.DEFAULT)]*3 + [(wx.TE_READONLY)]
    for index, item in enumerate(textNames):
      self.textwidgets4d.append(wx.TextCtrl(self,value="",name=item,\
                                style = styleList[index]))
      grid4d.Add(self.textwidgets4d[index],pos=textCoords[index],flag=textFlags)
      self.Bind(wx.EVT_TEXT, self.datafunc, self.textwidgets4d[index])

    self.simLengthLbl = wx.StaticText(self,label="Simulation Length Information")
    grid4d.Add(self.simLengthLbl,pos=(1,1),span=(1,4),flag=wx.ALIGN_CENTER)
    

    #######
    # hard code avg info to be '1' - to be changed in later versions.
    #######
    self.textwidgets4d[3].SetValue('1')


    self.prop_options = ["Energy_Total", "Energy_LJ", "Energy_Elec", "Energy_Intra", \
                         "Enthalpy", "Pressure", "Volume", "Nmols", "Density"]
    self.propb1_chks = []
    self.propb2_chks = []
    for i in range(len(self.prop_options)):
      self.propb1_chks.append(wx.CheckBox(self,-1,self.prop_options[i],name=('b1 %s' \
                                          %(self.prop_options[i]))))
      self.propb2_chks.append(wx.CheckBox(self,-1,self.prop_options[i],name=('b2 %s' \
                                          %(self.prop_options[i]))))
      grid4d.Add(self.propb1_chks[i],pos=(11+i,2),flag = \
                 (wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN))
      grid4d.Add(self.propb2_chks[i],pos=(11+i,4),flag= \
                 (wx.ALIGN_LEFT | wx.RESERVE_SPACE_EVEN_IF_HIDDEN))
      self.Bind(wx.EVT_CHECKBOX,self.checkfunc,self.propb1_chks[i])
      self.Bind(wx.EVT_CHECKBOX,self.checkfunc,self.propb2_chks[i])

    # the values we will be passing to our dictionary
    self.propb1 = []
    self.propb2 = []
    # initially, hide the MC steps labels, three textwidgets, and 
    # three labels preceding MC steps labels
    self.labels_to_hide = [1,2,3,8,9,10]
    self.text_to_hide = [0,1,2]
    self.step_lbls = [1,2,3]
    for item in self.labels_to_hide:
      self.staticlbls4d[item].Hide()
    for item in self.text_to_hide:
      self.textwidgets4d[item].Hide()


  def checkfunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = event.IsChecked()
    thisBox = obj[1] 
    valToPass = []
    objToPass = "propList b%s" %(thisBox)
    if objToPass == "propList b1":
      for index, item in enumerate(self.propb1_chks):
        if item.IsChecked():
          thisProp = str(item.GetName()[3:])
          valToPass.append(thisProp) 
        else:
          pass
    elif objToPass == "propList b2":
      for index, item in enumerate(self.propb2_chks):
        if item.IsChecked():
          thisProp = str(item.GetName()[3:])
          valToPass.append(thisProp)
        else:
          pass
    data(objToPass,valToPass)

  def datafunc(self,event):
    obj = str(event.GetEventObject().GetName())
    val = str(event.GetString())
    if obj == "freqTypeChoice":
      self.showFreqTypes(val)
    data(obj,val)


  def showFreqTypes(self,val):
    if not val:
      for  item in self.text_to_hide:
        self.textwidgets4d[item].Hide()
        self.textwidgets4d[item].SetValue('')
      for item in self.labels_to_hide:
        self.staticlbls4d[item].Hide()
    elif val == 'Steps':
      for item in self.text_to_hide:
        self.textwidgets4d[item].Show()
      for item in self.labels_to_hide:
        self.staticlbls4d[item].Show()
      for item in self.step_lbls:
        self.staticlbls4d[item].SetLabel('MC Steps')
    elif val == 'Minutes':
      for item in self.text_to_hide:
        self.textwidgets4d[item].Show()
      for item in self.labels_to_hide:
        self.staticlbls4d[item].Show()
      for item in self.step_lbls:
        self.staticlbls4d[item].SetLabel('minutes')

    self.Layout()

  def showAndHide(self,message,arg2=None):
    if not 'ensemble' in myDict:
      return
    if myDict['ensemble'][0:4] == "GEMC":
      for index, item in enumerate(self.propb2_chks):
        item.Show()
      self.staticlbls4d[7].Show()
    else:
      for index, item in enumerate(self.propb2_chks):
        item.Hide()
      self.staticlbls4d[7].Hide()
    self.Layout()


class PanelFour(wx.Panel):
  # Declare the class variables

  # initialize the panel
  def __init__(self,parent):
    wx.Panel.__init__(self,parent)

  # define the notebook
    nest = wx.Notebook(self)
    self.nest1 = p4a(nest)
    self.nest2 = p4b(nest)
    self.nest3 = p4c(nest)
    self.nest4 = p4d(nest)

    nest.AddPage(self.nest1,"Molecule Files")
    nest.AddPage(self.nest2,"Fragment Files")
    nest.AddPage(self.nest3,"Input File")
    nest.AddPage(self.nest4,"Output File")
    nested_sizer = wx.BoxSizer()
    nested_sizer.Add(nest,1,wx.EXPAND)
    self.SetSizer(nested_sizer)


  
################################################################################################
#
# The main frame, of which all other panels are children.
# Controls the size of the window, and the primary notebook housing the panels.
# Also contains the 'File > Save As..' feature.
#
################################################################################################
class MainFrame(wx.Frame):
  def __init__(self):
    wx.Frame.__init__(self,None,title="Cassandra Input File Editor v1.1")
    self.SetInitialSize((900,620))

    # Create the panel hosting the notebook
    panel_n = wx.Panel(self)
    notebook_m = wx.Notebook(panel_n)

    # Add PanelOne:PanelFour to notebook as pages
    self.P1 = PanelOne(notebook_m)
    self.P2 = PanelTwo(notebook_m)
    self.P3 = PanelThree(notebook_m)
    self.P4 = PanelFour(notebook_m)

    notebook_m.AddPage(self.P1, "Basic Information")
    notebook_m.AddPage(self.P2, "Interaction Parameters")
    notebook_m.AddPage(self.P3, "Probabilities Information")
    notebook_m.AddPage(self.P4, "File Handling")

    # Create a unique sizer for the notebook
    NBSizer = wx.BoxSizer()
    NBSizer.Add(notebook_m,1,wx.EXPAND)
    panel_n.SetSizer(NBSizer)


    # if changing the page, check to make sure all info is put in to dictionary
    self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)
   # uncomment above self.Bind statement to enable page changing restrictions
   # between the top level notebook pages

  def OnPageChanging(self,event):
    oldPage = event.GetOldSelection()
    newPage = event.GetSelection()
    
    # if the page we are leaving is 'Page 1' of basic information,
    # check that our dictionary has the required keywords
    if oldPage == 0:
      required_keys = ["runName", "nbrSpecies", "ensemble", "simDir"]
      count_num_keys = 0
      # check if these keys exist in the dictionary; if they are in 
      # the dictionary, then the data we require is there
      for item in required_keys:
        if item in myDict:
          count_num_keys = count_num_keys + 1
      # if all required keys are not present, the user is 
      # not allowed to change pages
      if count_num_keys < len(required_keys):
        event.Veto()
    

def tryToWrite(thisFile,*args):
  # writes an arbitrary number of arguments that may or may not have values to a single line
  # useful for, e.g. our probability of translation
  # we don't know ahead of time how many species are in the system
  # and therefore some of the keywords will not be in the dictionary.
  # this function prevents raising an error and allows the file to be written
  listToDel = []
  for item in args:
    if item not in myDict.keys():
      myDict[item] = ''
      listToDel.append(item)
  strToWrite = ''
  for index, item in enumerate(args):
    if index == 0:
      strToWrite = str(myDict[item])
    else:
      strToWrite = strToWrite + "  " +  str(myDict[item])
  thisFile.write('%s\n' %(strToWrite))
  # take it back out of our dictionary - it was added as an empty value
  for item in listToDel:
    data(item,'')
  return

# function used to populate myDict as needed
# if we pass an object (usually the widget's name) and the
# value received is empty, the key is removed from the dictionary
def data(obj,val,listStyle=False):
  myDict[obj] = val
  if not val:
    del myDict[obj]
  
   
  # if the object being changed is the ensemble or number
  # of species, call the show and hide protocols
  # these only need to be called after the counter is
  # greater than or equal to 2, as prior to this the user will
  # be stuck on the initial panel. 
  if obj == "ensemble" or obj == "nbrSpecies":
    global obj_counter
    obj_counter += 1
    if obj_counter >= 2:
      try:
        pub.sendMessage("refresh",message=None)
      except: 
        pub.sendMessage("refresh")



if __name__ =="__main__":
  app = wx.App()
  MainFrame().Show()
  app.MainLoop()
