#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gs
import glob
xlim = [0,5] #This is the y limits of your subplots
ylim = [-10,10] #The x limits of your subplots
colorkeyfile = "/home/adam/processing/output/color.key" #To add color to your lines
fontsize = 16
 
def Symmetries(fstring):
  f = open(fstring,'r')
  x = np.zeros(0)
  for i in f:
    if "high-symmetry" in i:
      x = np.append(x,float(i.split()[-1]))
  f.close()
  return x
 
def dosplot(directory,subplot, fermi): 
  filelist = glob.glob(directory + "/*wfc*") #This is all the files in the directory
  keys = {}
  for i in filelist: #This loops over all the files in the filelist and generates the dictionary with all the numpy arrays
    l = i.find('(') + 1
    r = i.find(')',l)
    atom = i[l:r] # Simply the atom name
    if atom not in keys:
      t = np.loadtxt(i) 
      keys[atom] = t[:,[0,1]] # We only want the energy and ldos
    else:
      t = np.loadtxt(i)
      t = t[:,[0,1]]
      keys[atom][:,1] = np.maximum(keys[atom][:,1],t[:,1]) # Since we already had a numpy array here, we take the maximum of each point
  for i in keys:
    keys[i] = np.fliplr(keys[i]) # We flip the array for vertical plotting of the pdos
  Compiled = []
  filekey = open(colorkeyfile,'r')
  colorkey = {'null' : 12345}
  for i in filekey:
    colorkey[i.split()[0]] = i.split()[1] #reading in the colorkey
  filekey.close()
  for i in keys:
    subplot.plot(keys[i][:,0], keys[i][:,1], colorkey[i]) #Plots the curve
    subplot.fill_betweenx(keys[i][:,1],keys[i][:,0], where=keys[i][:,1] < fermi, color=colorkey[i],alpha=0.25) # This fills in the plot
  subplot.plot([0,1000],[fermi,fermi],color='red') #Shows the Fermi Level
  subplot.set_xlabel('DOS [arb]',fontsize=fontsize) #Labels the x-axis
def bndplot(datafile,fermi,symmetryfile,subplot,label):
  z = np.loadtxt(datafile) #This loads the bandx.dat.gnu file
  x = np.unique(z[:,0]) #This is all the unique x-points
  bands = []
  bndl = len(z[z[:,0]==x[1]]) #This gives the number of bands in the calculation
  Fermi = float(fermi)
  axis = [min(x),max(x),Fermi - 4, Fermi + 4]
  for i in range(0,bndl):
    bands.append(np.zeros([len(x),2])) #This is where we store the bands
  for i in range(0,len(x)):
    sel = z[z[:,0] == x[i]]  #Here is the energies for a given x
    test = []
    for j in range(0,bndl): #This separates it out into a single band
      bands[j][i][0] = x[i]
      bands[j][i][1] = np.multiply(sel[j][1],13.605698066)
  for i in bands: #Here we plots the bands
    subplot.plot(i[:,0],i[:,1],color="black")
  temp = Symmetries(symmetryfile)
  for j in temp: #This is the high symmetry lines
    x1 = [j,j]
    x2 = [axis[2],axis[3]]
    subplot.plot(x1,x2,'--',lw=0.55,color='black',alpha=0.75)
  subplot.plot([min(x),max(x)],[Fermi,Fermi],color='red',)
  subplot.set_xticklabels([])
  subplot.set_ylim([axis[2],axis[3]])
  subplot.set_xlim([axis[0],axis[1]])
  subplot.set_xlabel('K-Path',fontsize=fontsize)
  subplot.set_ylabel('Energy [eV]',fontsize=fontsize)
  subplot.tick_params(axis='both',which='major',labelsize=14)
 
 
gs1 = gs.GridSpec(1,2,width_ratios=[3,1]) # Creating the subplots
gs1.update(wspace=0,hspace=0.0) # I want them right next to eachother
BND = plt.subplot(gs1[0,0]) #My band diagram
DOS = plt.subplot(gs1[0,1]) #My DOS plot
bndplot('bands.dat.gnu',0.0,'bands.out',BND,'test') #This is my bndplot with *.gnu coming from the output of bands.x, 2.1082 is the fermi level, bandx.out is the symmetries file, and 'test' is the title
ylim = BND.get_ylim()
dosplot('.',DOS,0.0) #This is my DOS plot with 2.1082 being the fermi level
DOS.set_ylim(BND.get_ylim()) #Placing the same y-limits
DOS.set_xlim(xlim) #This was defined in the header it was [0,10]
DOS.axes.get_yaxis().set_ticklabels([]) #No Tick Labels
DOS.axes.get_xaxis().set_ticklabels([]) #No Tick Lavels
plt.suptitle('Test',fontsize=24) #A Title for the Plot
plt.show()
