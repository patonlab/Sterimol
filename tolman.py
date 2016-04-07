#!/usr/bin/python

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# kelvin.jackson@chem.ox.ac.uk or robert.paton@chem.ox.ac.uk

#######################################################################
#                          tolman.py                                  #
#                                                                     #
#    Outline of the program                                           #
#######################################################################
#######  Written by:  Kelvin Jackson and Robert Paton #################
#######  Last modified:  Apr 07, 2016 #################################
#######################################################################

#Script Settings
minmax=1
verbose=0

#Python Libraries 
import subprocess, sys, os
from numpy import *
from scipy import *
from math import *
import numpy as np

# For reading Gaussian formatted input/output files
from ccParse import *

#Useful Arrays
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

metals = ["Li","Be","Na","Mg","Al","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv"]

heteroatoms=["Cl","P"]

def bondiRadius(massno):
   #Bondi van der Waals radii for all atoms from: Bondi, A. J. Phys. Chem. 1964, 68, 441-452, except hydrogen, which is taken from Rowland, R. S.; Taylor, R. J. Phys. Chem. 1996, 100, 7384-7391. Radii that are not available in either of these publications have RvdW = 2.00 Angstrom
   bondi = [0.0,1.09, 1.40, 1.82,2.00,2.00,1.70,1.55,1.52,1.47,1.54,2.27,1.73,2.00,2.10,1.80,1.80,1.75,1.88,2.75,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.40,1.39,1.87,2.00,1.85,1.90,1.85,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.72,1.58,1.93,2.17,2.00,2.06,1.98,2.16,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.72,1.66,1.55,1.96,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.86]
   if massno<len(bondi): radius = bondi[massno]
   else: radius = 2.0
   return radius

def vdwr(atom):
   atomno=periodictable.index(atom)
   return bondiRadius(atomno)

def getfragment(atom,molcart):
   bondlist=[atom]
   for a in range(len(molcart)):
      if distcalc(molcart[atom],molcart[a])<1.6 and a not in bondlist:bondlist.append(a)

      for b in range(len(bondlist)):
         for c in range(len(molcart)):
            if distcalc(molcart[bondlist[b]],molcart[c])<1.6 and c not in bondlist:bondlist.append(c)
   return bondlist

def dotproduct(v1, v2):
        return sum((a*b) for a, b in zip(v1, v2))

def length(v):
        return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
        return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def calcopposite(atom1,atom2,angle,molcart):
   h=distcalc(molcart[atom1],molcart[atom2])
   d=h*math.sin(angle)
   return d

def getcoords(atom,molcart):
   coords=[]
   for i in range(3): coords.append(molcart[atom][i])
   return coords

def calcvector(atom1,atom2):
   x=atom1[0]-atom2[0]; y=atom1[1]-atom2[1]; z=atom1[2]-atom2[2]
   return x, y, z

def avpoints(atomnos,molcart):
   xcoords=[]; ycoords=[]; zcoords=[]
   for a in atomnos:
      xcoords.append(molcart[a][0]); ycoords.append(molcart[a][1]); zcoords.append(molcart[a][2])
   syslength=len(xcoords)
   x=0;y=0;z=0
   for i in range(syslength):
      x=x+xcoords[i]
      y=y+ycoords[i]
      z=z+zcoords[i]
   x=x/syslength; y=y/syslength; z=z/syslength
   return round(x,8),round(y,8),round(z,8)

def distcalc(atom1,atom2):
   x=atom1[0]-atom2[0]; y=atom1[1]-atom2[1]; z=atom1[2]-atom2[2]
   dist = (x**2+y**2+z**2)**0.5
   return dist

def dihedral(atoma,atomb,atomc,atomd):
   x1=atoma[0]; y1=atoma[1]; z1=atoma[2]
   x2=atomb[0]; y2=atomb[1]; z2=atomb[2]
   x3=atomc[0]; y3=atomc[1]; z3=atomc[2]
   x4=atomd[0]; y4=atomd[1]; z4=atomd[2]
   ax= (y2-y1)*(z2-z3)-(z2-z1)*(y2-y3)
   ay= (z2-z1)*(x2-x3)-(x2-x1)*(z2-z3)
   az= (x2-x1)*(y2-y3)-(y2-y1)*(x2-x3)
   bx= (y3-y2)*(z3-z4)-(z3-z2)*(y3-y4)
   by= (z3-z2)*(x3-x4)-(x3-x2)*(z3-z4)
   bz= (x3-x2)*(y3-y4)-(y3-y2)*(x3-x4)
   nbx= (y2-y3)*(z4-z3)-(z2-z3)*(y4-y3)
   nby= (z2-z3)*(x4-x3)-(x2-x3)*(z4-z3)
   nbz= (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3)
   torsion=180.0/math.pi*math.acos((ax*bx+ay*by+az*bz)/(math.sqrt(ax*ax+ay*ay+az*az)*math.sqrt(bx*bx+by*by+bz*bz)))
   sign=180.0/math.pi*math.acos((nbx*(x2-x1)+nby*(y2-y1)+nbz*(z2-z1))/(math.sqrt(nbx*nbx+nby*nby+nbz*nbz)*math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))))
   if sign<90.0: torsion=torsion*-1.0
   return torsion

def concheck(conpar,val):
   cons=[]
   for a in range(len(conpar)):
      for b in range(len(conpar[a])):
         if val ==conpar[a][0]:
            for c in range(len(conpar[a])-1): cons.append(conpar[a][c+1])
            return cons

class calcTolman:
   def __init__(self, file):
   
      ## Use ccParse to get the Cartesian coordinates from Gaussian input/output files or pdb files
      if len(file.split(".com"))>1 or len(file.split(".gjf"))>1: fileData = getinData(file)
      if len(file.split(".out"))>1 or len(file.split(".log"))>1: fileData = getoutData(file)
      if len(file.split(".pdb"))>1: fileData = getpdbData(file)


      # Find the atom numbers of any metal centers
      metalatoms=[]   
      for i in range(len(fileData.ATOMTYPES)):
         if fileData.ATOMTYPES[i] in metals:metalatoms.append(i)
      if verbose==1:print "Detected metal atoms", metalatoms

      # Looks for C-C bonds
      ivals=[]; jvals=[]
      for i in range(len(fileData.ATOMTYPES)):
         for j in range(len(fileData.ATOMTYPES)):
            if distcalc(fileData.CARTESIANS[i],fileData.CARTESIANS[j]) <1.511 and fileData.ATOMTYPES[j] == "C" and fileData.ATOMTYPES[i] == "C": #Change criteria for analysis of heterocycles
               ivals.append(i); jvals.append(j)

      conpar=[]
      for a in range(len(ivals)):
         rar=[]
         rar.append(ivals[a])
                 
         for b in range(len(ivals)):
            if ivals[a]==ivals[b]:rar.append(jvals[b])
         if rar not in conpar:conpar.append(rar)

      allrings=[]
      for a in range(len(conpar)):
         z = conpar[a][0]
         for b in concheck(conpar,z):
            y = b
            for c in concheck(conpar,y):
               x = c
               for d in concheck(conpar,x):
                  w = d
                  for e in concheck(conpar,w):
                     v = e
                     rar = []
                     rar.extend([z,y,x,w,v])
                     if z in concheck(conpar,v) and sorted(rar) not in allrings and len(set(rar))==5:allrings.append(sorted(rar))
                     for f in concheck(conpar,v):
                        u = f
                        tar = []
                        tar.extend([z,y,x,w,v,u])
                        if z in concheck(conpar,u) and sorted(tar) not in allrings and len(set(tar))==6:allrings.append(sorted(tar))
      #Define centroid vector
      mcdists=[]; mcdist=9999
      for a in allrings:
         #Start point for rings of interest
         if len(a)==5:
            
            cent=avpoints(a,fileData.CARTESIANS)
            m=fileData.CARTESIANS[metalatoms[0]]
            tempmcdist=mcdist
            mcdist=distcalc(m,cent)

            for b in metalatoms:
               m=fileData.CARTESIANS[b]
               if mcdist>=distcalc(m,cent):moi=b;mcdist=distcalc(m,cent);metal=b
            mcdists.append([mcdist,moi])
            frag=getfragment(a[0],fileData.CARTESIANS)
            vect1=calcvector(getcoords(metal,fileData.CARTESIANS),cent)

            if tempmcdist==mcdist:break #Stops if dealing with identical ring system as before
            
            self.tolman = []; opplist=[]; alpha=[]; beta=[]; ringang=[]; self.nontolman=[]
            theta=[] # Candidate Tolman angle substituent
            omega=[] # standardised atom "dihedral" orientation
            eta=[] # Candidate non-Tolman cone angle substituent
            for i in frag:
               vect2=calcvector(getcoords(metal,fileData.CARTESIANS),getcoords(i,fileData.CARTESIANS))
               oppdist=calcopposite(metal,i,angle(vect1,vect2),fileData.CARTESIANS)
               opplist.append(oppdist+vdwr(fileData.ATOMTYPES[i]))
               alpha.append(angle(vect1,vect2))
               adj=distcalc(getcoords(i,fileData.CARTESIANS),getcoords(metal,fileData.CARTESIANS))
               beta.append(math.asin(vdwr(fileData.ATOMTYPES[i])/adj))
               theta.append(alpha[-1]+beta[-1])
               eta.append(alpha[-1])
               omega.append(dihedral([10,10,10],getcoords(metal,fileData.CARTESIANS),cent,getcoords(i,fileData.CARTESIANS)))
               if i in a:ringang.append(dihedral([10,10,10],getcoords(metal,fileData.CARTESIANS),cent,getcoords(i,fileData.CARTESIANS)))

            interval=180/len(a)
            for k in ringang:
               tlist=[];tang=[];tfrag=[];ntang=[]

               for h in range(len(frag)):
                  if k-interval<omega[h]<k+interval:tlist.append(frag[h])
                  if k>(180-interval) and k-interval<omega[h]+360<k+interval:tlist.append(frag[h])
                  if k<-(180-interval) and k-interval<omega[h]-360<k+interval:tlist.append(frag[h])
               
               for t in range(len(frag)):
                  if frag[t] in tlist: tang.append(theta[t]);tfrag.append(frag[t]);ntang.append(eta[t])
               self.tolman.append(math.degrees(max(tang)))
               self.nontolman.append(math.degrees(max(ntang)))
                  
            self.dist_to_centroid = distcalc(cent,getcoords(metal,fileData.CARTESIANS))


if __name__ == "__main__":
   if verbose==0 and minmax==0:print "FileName","MCentDist","TolmanCA","NonTolmanCA",
   if verbose==0 and minmax==1:print "FileName","MCentDist","TolmanCA","MinTolmanCA","MaxTolmanCA","NonTolmanCA",
   filetypes=["out","log","com","gjf","pdb"]
   files = []
   if len(sys.argv) > 1:
      for i in range(1,len(sys.argv)):
         if len(sys.argv[i].split(".")) > 1:
            if sys.argv[i].split(".")[1] in filetypes:
            #if sys.argv[i].split(".")[1] == "out" or sys.argv[i].split(".")[1] == "log" or sys.argv[i].split(".")[1] == "com" or sys.argv[i].split(".")[1] == "gjf" or sys.argv[i].split(".")[1] == "pdb":
               files.append(sys.argv[i])
   else: print "\nWrong number of arguments used. Correct format: tolman.py file(s)\n"; sys.exit()
   print ""
   for file in files:
      file_Params = calcTolman(file)
      if verbose==1:print file_Params.tolman;print file_Params.nontolman
      max_tolman_angle = round(max(file_Params.tolman),3)
      min_tolman_angle = round(min(file_Params.tolman),3)
      av_tolman_angle = round(2 * sum(file_Params.tolman) / float(len(file_Params.tolman)),3)
      mc_dist = round(file_Params.dist_to_centroid,3)
      av_nontolman=round(2 * sum(file_Params.nontolman) / float(len(file_Params.nontolman)),3)

      if verbose==1:print file, "--- Tolman angle", av_tolman_angle, "degrees; M-Centroid_Distance", mc_dist, "Angstrom; Min Tolman angle",min_tolman_angle, "degrees; Max Tolman angle", max_tolman_angle, "degrees; Non-Tolman angle",av_nontolman, ;print "\n"
      if verbose==0 and minmax==0:print file, mc_dist, av_tolman_angle, av_nontolman
      if verbose==0 and minmax==1:print file, mc_dist, av_tolman_angle, min_tolman_angle,max_tolman_angle,av_nontolman
