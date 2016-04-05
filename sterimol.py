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

###############################################################
#                         sterimol.py                         #
#                                                             #
###############################################################
#######  Written by:  Kelvin Jackson ##########################
#######  Last modified:  Mar 20, 2016 #########################
###############################################################

#Python Libraries 
import subprocess, sys, os
from numpy import *
from scipy import *
from math import *
import numpy as np
#from vpython import *

# For reading Gaussian formatted input/output files
from ccParse import *

# Dependent on parameter file
from pars import *

#Avoid number error warnings
import warnings
warnings.filterwarnings("ignore")

#Useful Arrays
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

def bondiRadius(massno):
   #Bondi van der Waals radii for all atoms from: Bondi, A. J. Phys. Chem. 1964, 68, 441-452, except hydrogen, which is taken from Rowland, R. S.; Taylor, R. J. Phys. Chem. 1996, 100, 7384-7391. Radii that are not available in either of these publications have RvdW = 2.00 Angstrom
   bondi = [0.0,1.09, 1.40, 1.82,2.00,2.00,1.70,1.55,1.52,1.47,1.54,2.27,1.73,2.00,2.10,1.80,1.80,1.75,1.88,2.75,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.40,1.39,1.87,2.00,1.85,1.90,1.85,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.72,1.58,1.93,2.17,2.00,2.06,1.98,2.16,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.72,1.66,1.55,1.96,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.86]
   if massno<len(bondi): radius = bondi[massno]
   else: radius = 2.0
   return radius

# Verloop's original Sterimol parameters use CPK atomic VdW radii based on atom-type definitions
sterimol_atomtypes = ["C", "C2", "C3", "C4", "C5/N5", "C6/N6", "C7", "C8", "H", "N", "C66", "N4", "O", "O2", "P", "S", "S1", "F", "C1", "S4", "B1", "I"]

# CPK VdW radii in pm
cpk_radii = [150,160,160,150,170,170,170,150,100,150,170,145,135,135,140,170,100,135,180,140,195,215]

# Generate Sterimol atom type from connectivity data
def generate_atom_types(atomtype, cn):
   st_types = []
   for i in range(0,len(atomtype)):
      atom = atomtype[i]
      if atom == "H": st_types.append("H")
      elif atom == "P": st_types.append("P")
      elif atom == "F": st_types.append("F")
      elif atom == "Cl": st_types.append("C1")
      elif atom == "Br": st_types.append("B1")
      elif atom == "I": st_types.append("I")
      elif atom == "O": #Sterimol distinguishes between "normal", and double-bonded O atoms
         if cn[i] < 1.5: st_types.append("O2")
         if cn[i] > 1.5: st_types.append("O")
      elif atom == "S": #Sterimol distinguishes between "normal", tetrahedral, and octohedral S atoms
         if cn[i] < 2.5: st_types.append("S")
         if 5.5 > cn[i] > 2.5: st_types.append("S4")
         if cn[i] > 5.5: st_types.append("S1")
      elif atom == "N": #Sterimol distinguishes between tetrahedral and planar (amide) N atoms
         if cn[i] > 2.5: st_types.append("N")
         if cn[i] < 2.5: st_types.append("C6/N6")
      elif atom == "C": #Sterimol distinguishes between myriad types of C atoms ...
         if cn[i] < 2.5: st_types.append("C3")
         if 3.5 > cn[i] > 2.5: # need to differentiate between sp2 carbon and aromatic carbon ...
            st_types.append("C6/N6") # assumes aromatic rather than sp2
         if cn[i] > 3.5: st_types.append("C")
   return st_types

# Calculation of atomic coordination numbers (taken from Grimme's DFTD3 definitions)
def ncoord(natom, rcov, atomtype, coords):
   max_elem = 94
   k1 = 16.0
   k2 = 4.0/3.0
   cn =[]
   for i in range(0,natom):
      xn = 0.0
      for iat in range(0,natom):
         if iat != i:
            dx = coords[iat][0] - coords[i][0]
            dy = coords[iat][1] - coords[i][1]
            dz = coords[iat][2] - coords[i][2]
            r2 = dx*dx+dy*dy+dz*dz
            r = math.pow(r2,0.5)
            r = r
            for k in range(0,max_elem):
               if atomtype[i].find(elements[k])>-1:Zi=k
               if atomtype[iat].find(elements[k])>-1:Ziat=k
            
            rco = rcov[Zi]+rcov[Ziat]
            rco = rco*k2
            rr=rco/r
            damp=1.0/(1.0+math.exp(-k1*(rr-1.0)))
            xn=xn+damp
      cn.append(xn)
   return cn

def calcdist(atoma,atomb,coords):
   x1=coords[atoma][0]
   y1=coords[atoma][1]
   z1=coords[atoma][2]
   x2=coords[atomb][0]
   y2=coords[atomb][1]
   z2=coords[atomb][2]
   ba = [x1-x2, y1-y2, z1-z2]
   dist = math.sqrt(ba[0]*ba[0]+ba[1]*ba[1]+ba[2]*ba[2])
   return dist

#3D rotation about normalized vector perpendicular to both vector arguments; rotation parallelizes vect1 with vect2 - vect 3 is rotated relative
def rotrel(vect1,vect2,vect3):
   ax=np.cross(vect1,vect2)###Find axis of rotation
   ang=math.acos((np.dot(vect1,vect2))/(np.linalg.norm(vect1)*np.linalg.norm(vect2)))###Find angle between vectors
   norm=1/(np.linalg.norm(ax))
   axnorm=np.dot(ax,norm)
   ux=axnorm[0]
   uy=axnorm[1]
   uz=axnorm[2]
   a=math.cos(ang)+((ux*ux)*(1-math.cos(ang)))
   b=(ux*uy*(1-math.cos(ang)))-(uz*math.sin(ang))
   c=(ux*uz*(1-math.cos(ang)))+(uy*math.sin(ang))
   d=(uy*ux*(1-math.cos(ang)))+(uz*math.sin(ang))
   e=(math.cos(ang))+(uy*uy*(1-math.cos(ang)))
   f=(uy*uz*(1-math.cos(ang)))-(ux*math.sin(ang))
   g=(uz*ux*(1-math.cos(ang)))-(uy*math.sin(ang))
   h=(uz*uy*(1-math.cos(ang)))+(ux*math.sin(ang))
   i=math.cos(ang)+(uz*uz*(1-math.cos(ang)))
   bigmat=([[a,b,c],[d,e,f,],[g,h,i]])
   vect=np.dot(bigmat,vect3)
   #vecttest=np.cross(vect,vect2) ###Since vect has been made parallel with vect2, this should produce a vector of zero magnitude
   return vect

def calcopposite(atom1,atom2,angle,molcart):
   h=calcdist(atom1,atom2,molcart)
   d=h*math.sin(angle)
   return d

def calcadj(atom1,atom2,angle,molcart):
   h=calcdist(atom1,atom2,molcart)
   d=h*math.cos(angle)
   return d

def getcoords(atom,molcart):
   coords=[]
   for i in range(3): coords.append(molcart[atom][i])
   return coords

def calcvector(atom1,atom2):
   x=atom1[0]-atom2[0]
   y=atom1[1]-atom2[1]
   z=atom1[2]-atom2[2]
   return x, y, z

def dist(atom1,atom2,molcart):#atom 2 is defined by cart
        x=molcart[atom1][0]-atom2[0]
        y=molcart[atom1][1]-atom2[1]
        z=molcart[atom1][2]-atom2[2]
        dist = (x**2+y**2+z**2)**0.5
        return dist

def dprod(v1, v2): return sum((a*b) for a, b in zip(v1, v2))

def length(v): return math.sqrt(dprod(v, v))

def angle(v1, v2):
   val = dprod(v1, v2) / length(v1) / length(v2)
   if val > 0.999999: val = 1.0
   if val < -0.999999: val = -1.0
   return math.acos(val)

def twod_dist(a,b,c):#Returns distance from a to line vector bc
   vect1=np.subtract(a,b)
   vect2=np.subtract(b,c)
   ang=angle(vect1,vect2)
   return math.sin(ang)*np.linalg.norm(vect1)

def twod_vect(a,b,c):#Returns normal vector intersect to line
   vect1=np.subtract(a,b)
   vect2=np.subtract(b,c)
   ang=angle(vect1,vect2)
   nvect2=vect2/np.linalg.norm(vect2)
   return ((math.cos(ang)*np.linalg.norm(vect1))*nvect2)+b

def twod_rot(vect,theta):#counterclockwise rotation of vector
   a=math.cos(theta)
   b=math.sin(theta)
   mat=[[a,-b],[b,a]]
   vect=np.dot(mat,vect)
   return vect
def linearcheck(carts):
   ans=0;xgrad=[];ygrad=[]
   for row in carts:xgrad.append(round(np.gradient(row)[0],4));ygrad.append(round(np.gradient(row)[1],4))
   if min(xgrad)==max(xgrad) and min(ygrad)==max(ygrad):ans=1
   return ans
# main class calculates the L, B1 and B5 parameters
class calcSterimol:
   def __init__(self, file, radii, atomA, atomB):
      
      if len(file.split(".com"))>1 or len(file.split(".gjf"))>1: fileData = getinData(file)
      if len(file.split(".out"))>1 or len(file.split(".log"))>1: fileData = getoutData(file)
      
      # initialize the array of atomic vdw radii
      molcart = fileData.CARTESIANS; atomtype = fileData.ATOMTYPES; natoms = len(molcart); vdw_radii = []
      
      if radii == "cpk":
         print "\n   STERIMOL: using original CPK Van der Waals parameters"
         atomic_co_no = ncoord(natoms, rcov, atomtype, molcart)
         sterimol_types = generate_atom_types(atomtype, atomic_co_no)
         #print sterimol_types
         for i in range(0,natoms):
            for j in range(0,len(sterimol_atomtypes)):
               if sterimol_types[i] == sterimol_atomtypes[j]: vdw_radii.append(cpk_radii[j]/100.00)

      if radii == "bondi":
         print "\n   STERIMOL: using Bondi Van der Waals parameters"
         for i in range(0,natoms): vdw_radii.append(bondiRadius(periodictable.index(fileData.ATOMTYPES[i])))
      
      # Define vector along the L-axis connecting base atom and the next attached atom
      # subtract one since the array starts from zero not one
      atomA = atomA - 1; atomB = atomB - 1
      base_atom = molcart[atomA]
      next_atom = molcart[atomB]
      vect1=calcvector(getcoords(atomA,molcart),next_atom)
      print "   Atoms", atomA, "and", atomB, "define the L-axis and direction", vect1

      print "\n", "   Atom ".ljust(9), "  Xco/A".rjust(9), "  Yco/A".rjust(9), "  Zco/A".rjust(9), " VdW/pm".rjust(9)
      print "   ##############################################"
      # Remove the base atom from the list of atoms to be considered for sterics (after printing all)
      atomlist = list(xrange(0,natoms))
      for atom in atomlist:
         if radii == "cpk": print "  ", sterimol_types[atom].ljust(6),
         if radii == "bondi": print "  ", atomtype[atom].ljust(6),
         for coord in molcart[atom]:
            if coord < 0.0: print "   %.3f".rjust(6) % coord,
            else: print "    %.3f".rjust(6) % coord,
         print "    %.1f" % round(vdw_radii[atom]*100)
      atomlist.remove(atomA)

      adjlist=[]; opplist=[]; theta=[]
      for i in atomlist:
         vect2=calcvector(getcoords(atomA,molcart),getcoords(i,molcart))
         oppdist=calcopposite(atomA,i,angle(vect1,vect2),molcart)
         opplist.append(oppdist+vdw_radii[i])
         adjdist=calcadj(atomA,i,angle(vect1,vect2),molcart)
         #minadjlist.append(adjdist-vdw_radii[i])
         adjlist.append(adjdist+vdw_radii[i])
         #if i == atomB: minval = adjdist-vdw_radii[i]
      #print adjlist
      #print minadjlist
      #print minval
      B5=max(opplist)
      #self.lval=max(adjlist)-minval
      # A bit weird, but seems like original sterimol adds on the difference between the bond length and vdw radius of atom B. For a C-H bond this is 1.50 - 1.10 = 0.40 Angstrom)
      self.lval=max(adjlist)+0.40

      ###Useful - do not delete!
      #print "   B5 atom", atomlist[opplist.index(max(opplist))]+1, "distance", max(opplist)
      #print "   Highest atom", atomlist[adjlist.index(max(adjlist))]+1,"distance", max(adjlist),"\n   Lowest atom", atomlist[minadjlist.index(min(minadjlist))]+1,"distance", min(minadjlist)
  
      zcarts=[]#zeroed carts
      for i in atomlist: zcarts.append(np.subtract(molcart[i],molcart[atomA]))
      zvect=[0,0,1]
      zcent=np.subtract(next_atom,molcart[atomA])
      for cart in range(len(zcarts)):
         zcoord= rotrel(zcent,zvect,zcarts[cart])
         zcarts[cart]=zcoord
      twodcarts=[]
      for row in zcarts: twodcarts.append([row[0],row[1]])
      fragrad=[]#radii of fragment atoms
      for t in atomlist: fragrad.append(vdw_radii[t])
      singledist=[]
      for t in range(len(fragrad)):
         d=np.linalg.norm(twodcarts[t])#;print d
         d=d+fragrad[t]
         singledist.append(d)
      self.newB5=max(singledist) #This is the same as the 3D calculated value from above

      center=[0,0]
      vlist=[]#list of distances from the origin to the tangential vectors
      alist=[]#list of atoms between which the tangential vectors pass through no other atoms
      iav=[]#interatomic vectors

      for x in range(len(twodcarts)):
         for y in range(len(twodcarts)):
            if x!=y:
               try:nvect= (twod_vect(center,twodcarts[x],twodcarts[y]))#origin normal vector to connecting atomic centers vector
               except ValueError:nvect=[0,0]
               iav=np.subtract(twodcarts[x],twodcarts[y])#interatomic vector
               iad=np.linalg.norm(iav)#interatomic distance
               try:theta=math.asin((fragrad[y]-fragrad[x])/iad)#calculates angle by which to rotate vdw radii before adding
               except ValueError: theta=np.pi/2
               try:unvect=nvect/np.linalg.norm(nvect)
               except RuntimeWarning:pass#unvect=[0,0]
               xradv=twod_rot(unvect*fragrad[x],theta)
               yradv=twod_rot(unvect*fragrad[y],theta)
               nvect= (twod_vect(center,twodcarts[x]+xradv,twodcarts[y]+yradv))#origin normal vector to connecting atomic surfaces tangential vector
               newx=twodcarts[x]+xradv
               newy=twodcarts[y]+yradv
               if np.cross(nvect,xradv)<0.000000001 and theta!=np.pi/2:
                  satpoint=[]#Satisfied points not within range of tangential vector
                  for z in range(len(twodcarts)):
                     pvdist=twod_dist(twodcarts[z],newx,newy)
                     if z!=x and z!=y and pvdist>fragrad[z]:satpoint.append(pvdist)
                  if len(satpoint)==len(atomlist)-2:vlist.append(np.linalg.norm(nvect));alist.append([x,y]);#print x,y
      #print vlist, len(atomlist), min(vlist),alist[vlist.index(min(vlist))]
      if linearcheck(twodcarts)==1:self.B1 = max(fragrad)
      elif len(vlist) > 0: self.B1=min(vlist)
      else: self.B1 = max(fragrad)

if __name__ == "__main__":
   files = []
   radii = "cpk"; atom1 = 1; atom2 = 2
   if len(sys.argv) > 1:
      for i in range(1,len(sys.argv)):
         if sys.argv[i] == "-radii": radii = (sys.argv[i+1])
         elif sys.argv[i] == "-a1": atom1 = int(sys.argv[i+1])
         elif sys.argv[i] == "-a2": atom2 = int(sys.argv[i+1])
         else:
            if len(sys.argv[i].split(".")) > 1:
               if sys.argv[i].split(".")[1] == "out" or sys.argv[i].split(".")[1] == "log" or sys.argv[i].split(".")[1] == "com" or sys.argv[i].split(".")[1] == "gjf":
                  files.append(sys.argv[i])
   else: print "\nWrong number of arguments used. Correct format: sterimol.py (-radii cpk/bondi) (-a1 atomid) (-a2 atomid) file(s)\n"; sys.exit()

   for file in files:
      file_Params = calcSterimol(file, radii, atom1, atom2)
      lval = file_Params.lval; B1 = file_Params.B1; B5 = file_Params.newB5
      print "\n   Structure               L1        B1        B5 (Ang)"
      print "  ", file.ljust(16), "  %.2f".rjust(9) % lval, "  %.2f".rjust(9) % B1, "  %.2f".rjust(9) % B5
      print ""
