from sterimol import *
example_folder = "examples"
for radii in ["cpk", "bondi"]:
    print "\n   Structure               L1        B1        B5 (Ang): from", radii, "definitions"

    for structure in ["H", "Me", "Et", "iPr", "nBu", "CH2iPr", "cHex", "nPr", "Ad", "tBu", "CH2tBu", "CHEt2", "CHiPr2","CHPr2", "CEt3", "Ph","Bn", "4ClPh","4MePh","4MeOPh","35diMePh","1Nap"]:
        file = example_folder+"/"+structure+".com"
        atom1 = 1; atom2 = 2
        file_Params = calcSterimol(file, radii, atom1, atom2, False)
        lval = file_Params.lval; B1 = file_Params.B1; B5 = file_Params.newB5
        print "  ", file.ljust(16), "  %.2f".rjust(9) % lval, "  %.2f".rjust(9) % B1, "  %.2f".rjust(9) % B5
