#!/usr/bin/python

###############################################################
#                         ccParse.py                          #
#                Reads compchem job file(s)                   #
###############################################################

#Python Libraries 
import subprocess, sys, os

## Check for integer when parsing ##
def is_number(s):
   try: int(s); return True
   except ValueError: return False

#Some useful arrays for chemists
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr", "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl", "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo","Bq"]

def elementID(massno):
   if massno < len(periodictable): return periodictable[massno]
   else: return "XX"


#Read Cartesian coordinate data from a PDB file
class getpdbData:
   def __init__(self, file):
      if not os.path.exists(file):
         print ("\nFATAL ERROR: Input file [ %s ] does not exist"%file); sys.exit()

      def getATOMS(self, inlines):
         self.ATOMTYPES = []
         self.CARTESIANS = []
         for i in range(0,len(inlines)):
            if inlines[i].find("ATOM") > -1:
               self.ATOMTYPES.append(int(inlines[i].split()[1]))
               self.CARTESIANS.append(float(inlines[i].split()[2:4]))
            if inlines[i].find("HETATM")>-1:
               self.ATOMTYPES.append(inlines[i].split()[-1])
               self.CARTESIANS.append([float(inlines[i].split()[-4]), float(inlines[i].split()[-3]), float(inlines[i].split()[-2])])

      infile = open(file,"r")
      inlines = infile.readlines()
      getATOMS(self, inlines)


#Read Cartesian data from a Gaussian formatted input file (*.com or *.gjf)
class getinData:
   def __init__(self, file):
      if not os.path.exists(file):
         print ("\nFATAL ERROR: Input file [ %s ] does not exist"%file); sys.exit()

      def getATOMTYPES(self, inlines):
         self.ATOMTYPES = []
         for i in range(0,len(inlines)):
            if inlines[i].find("#") > -1:
               if len(inlines[i+1].split()) == 0: start = i+5
               if len(inlines[i+2].split()) == 0: start = i+6
               break
         for i in range(start,len(inlines)):
            if len(inlines[i].split()) ==0: break
            else: self.ATOMTYPES.append(inlines[i].split()[0])

      def getCARTESIANS(self, inlines, natoms):
         self.CARTESIANS = []
         for i in range(0,len(inlines)):
            if inlines[i].find("#") > -1:
               start = i+5
               break

         for i in range(start,len(inlines)):
            if len(inlines[i].split()) == 0: break
            elif len(inlines[i].split()) == 4: self.CARTESIANS.append([float(inlines[i].split()[1]), float(inlines[i].split()[2]), float(inlines[i].split()[3])])

      def getMETHOD(self, inlines):
         self.FUNCTIONAL = "none"
         # looks for a selected group of methods (some of my favourites...)
         for i in range(0,len(inlines)):
            if inlines[i].find("#") > -1:
               if inlines[i].find("B3LYP") > -1 or inlines[i].find("b3lyp") > -1 : self.FUNCTIONAL = "B3LYP"
               if inlines[i].find("TPSSTPSS") > -1 or inlines[i].find("tpsstpss") > -1 : self.FUNCTIONAL = "TPSSTPSS"

      def getBONDINDEX(self,inlines,natoms):
         conn=[]
         connectivity = 0

         for j in range(0,len(inlines)):
            if " 1 " in inlines[j]:
               startconn = j
               connectivity  = 1
               break

         if connectivity == 1:
            for j in range(startconn,len(inlines)): conn.append(inlines[j])

         self.BONDINDEX=[]
         for j in range(0,natoms):
            self.BONDINDEX.append([0])
            for k in range(0,natoms): self.BONDINDEX[j].append(0)
         
         for j in range(0,natoms):
            if connectivity == 1:
               for bonded in conn[j].split():
                  if is_number(bonded) ==True:
                     if int(bonded)-1!=j:
                        self.BONDINDEX[j][int(bonded)-1]=1
                        self.BONDINDEX[int(bonded)-1][j]=1

      infile = open(file,"r")
      inlines = infile.readlines()
      getATOMTYPES(self, inlines)
      self.NATOMS=len(self.ATOMTYPES)
      getMETHOD(self,inlines)
      getBONDINDEX(self,inlines,self.NATOMS)
      getCARTESIANS(self, inlines, self.NATOMS)


#Read Cartesian data from a Gaussian formatted output file (*.log or *.out)
class getoutData:
   def __init__(self, file):
      if not os.path.exists(file):
         print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file); sys.exit()

      def getATOMTYPES(self, outlines):
         self.ATOMTYPES = []
         self.CARTESIANS = []
         self.ATOMICTYPES = []

         anharmonic_geom=0
         for i in range(0,len(outlines)):
            if outlines[i].find("Input orientation") > -1: standor = i
            if outlines[i].find("Standard orientation") > -1: standor = i
            if outlines[i].find("Vib.Av.Geom.") > -1: standor = i; anharmonic_geom=1
            if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1: self.NATOMS = i-standor-6

         try: standor
         except NameError: pass
         else:
            for i in range (standor+5,standor+5+self.NATOMS):
               self.ATOMTYPES.append(elementID(int(outlines[i].split()[1])))
               #self.ATOMICTYPES.append(int(outlines[i].split()[2]))

               if anharmonic_geom==0:
                  if len(outlines[i].split()) > 5: self.CARTESIANS.append([float(outlines[i].split()[3]),float(outlines[i].split()[4]),float(outlines[i].split()[5])])
                  else: self.CARTESIANS.append([float(outlines[i].split()[2]),float(outlines[i].split()[3]),float(outlines[i].split()[4])])
               if anharmonic_geom==1:self.CARTESIANS.append([float(outlines[i].split()[2]),float(outlines[i].split()[3]),float(outlines[i].split()[4])])
      
      def getMETHOD(self, outlines):
         self.FUNCTIONAL = "none"
         # looks for a selected group of methods (some of my favourites...)
         for i in range(0,len(outlines)):
            if outlines[i].find("B3LYP)") > -1: self.FUNCTIONAL = "B3LYP"
            if outlines[i].find("CAM-B3LYP)") > -1: self.FUNCTIONAL = "CAM-B3LYP"
            if outlines[i].find("B-P86)") > -1: self.FUNCTIONAL = "BP86"
            if outlines[i].find("B2PLYP)") > -1: self.FUNCTIONAL = "B2PLYP"
            if outlines[i].find("M06)") > -1: self.FUNCTIONAL = "M06"
            if outlines[i].find("M062X)") > -1: self.FUNCTIONAL = "M06-2X"
            if outlines[i].find("M06L)") > -1: self.FUNCTIONAL = "M06L"
            if outlines[i].find("B97D)") > -1: self.FUNCTIONAL = "B97D"
            if outlines[i].find("TPSS-TPSS)") > -1: self.FUNCTIONAL = "TPSSTPSS"

      if os.path.exists(file):outfile = open(file,"r")
      outlines = outfile.readlines()
      getMETHOD(self, outlines)
      getATOMTYPES(self, outlines)