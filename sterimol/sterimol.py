#!/usr/bin/python
from __future__ import print_function, absolute_import

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome:
# kelvin.jackson@chem.ox.ac.uk or robert.paton@chem.ox.ac.uk

###############################################################
#                         sterimol.py                         #
#                                                             #
###############################################################
####  Written by: Dr Kelvin Jackson and Prof Robert Paton #####
####  Last modified:  Nov 09, 2016 ############################
###############################################################


from .sterimoltools import *

def main():
   files = []
   jobtype=2
   radii = "cpk"; atom1 = 1; atom2 = 2
   if len(sys.argv) > 1:
      for i in range(1,len(sys.argv)):
         if sys.argv[i] == "-radii": radii = (sys.argv[i+1])
         elif sys.argv[i] == "-a1": atom1 = int(sys.argv[i+1]);jobtype=1
         elif sys.argv[i] == "-a2": atom2 = int(sys.argv[i+1])
         else:
            if len(sys.argv[i].split(".")) > 1:
               if sys.argv[i].split(".")[1] == "out" or sys.argv[i].split(".")[1] == "log" or sys.argv[i].split(".")[1] == "com" or sys.argv[i].split(".")[1] == "gjf":
                  files.append(sys.argv[i])
   else: print("\nWrong number of arguments used. Correct format: sterimol.py (-radii cpk/bondi) (-a1 atomid) (-a2 atomid) file(s)\n"); sys.exit()
   if jobtype==1:
      for file in files:
         file_Params = calcSterimol(file, radii, atom1, atom2, True)
         lval = file_Params.lval; B1 = file_Params.B1; B5 = file_Params.newB5
         print("\n   STERIMOL: using", radii, "van der Waals parameters")
         print("\n","   Structure".ljust(25),"L".rjust(9),"B1".rjust(9),"B5".rjust(9))
         print("   "+file.ljust(22), "  %.2f".rjust(9) % lval, "  %.2f".rjust(9) % B1, "  %.2f".rjust(9) % B5)
         print("")
   if jobtype==2:
      print("\n   Sandwich Analysis\n   STERIMOL: using original CPK Van der Waals parameters\n")
      print("   "+"Structure".ljust(25), "Tolman_CA".rjust(9),"MC_dist".rjust(9),"L".rjust(9),"B1".rjust(9),"B5".rjust(9))
      for file in files:
         calcSandwich(file)
      print("\n")

if __name__ == "__main__":
    main()

