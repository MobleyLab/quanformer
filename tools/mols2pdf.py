#!/usr/bin/env python
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

#############################################################################
# Converts molecules into a multi-page PDF document.
#############################################################################

import sys
from openeye import oechem
from openeye import oedepict


def main(argv=[__name__]):

    itf = oechem.OEInterface(InterfaceData)
    oedepict.OEConfigureReportOptions(itf)
    oedepict.OEConfigurePrepareDepictionOptions(itf)
    oedepict.OEConfigure2DMolDisplayOptions(itf)

    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    iname = itf.GetString("-in")
    ifs = oechem.oemolistream()
    ifs.SetConfTest(oechem.OEAbsoluteConfTest()) # VTL
    if not ifs.open(iname):
        oechem.OEThrow.Fatal("Cannot open input file!")

    oname = itf.GetString("-out")
    ext = oechem.OEGetFileExtension(oname)
    if ext != "pdf":
        oechem.OEThrow.Fatal("Output must be PDF format.")

    ofs = oechem.oeofstream()
    if not ofs.open(oname):
        oechem.OEThrow.Fatal("Cannot open output file!")

    if itf.HasString("-ringdict"):
        rdfname = itf.GetString("-ringdict")
        if not oechem.OEInit2DRingDictionary(rdfname):
            oechem.OEThrow.Warning("Cannot use user-defined ring dictionary!")

    ropts = oedepict.OEReportOptions()
    oedepict.OESetupReportOptions(ropts, itf)
    ropts.SetFooterHeight(25.0)
    report = oedepict.OEReport(ropts)

    popts = oedepict.OEPrepareDepictionOptions()
    oedepict.OESetupPrepareDepictionOptions(popts, itf)

    dopts = oedepict.OE2DMolDisplayOptions()
    oedepict.OESetup2DMolDisplayOptions(dopts, itf)
    dopts.SetDimensions(report.GetCellWidth(), report.GetCellHeight(), oedepict.OEScale_AutoScale)

    for mol in ifs.GetOEMols(): # VTL ignore confs; dont use GetOEGraphMols
        print(mol.GetTitle()) # VTL
        cell = report.NewCell()
        oedepict.OEPrepareDepiction(mol, popts)
        disp = oedepict.OE2DMolDisplay(mol, dopts)
        oedepict.OERenderMolecule(cell, disp)

    font = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Bold, 12,
                           oedepict.OEAlignment_Center, oechem.OEBlack)
    for pagenum, footer in enumerate(report.GetFooters()):
        text = "Page %d of %d" % (pagenum + 1, report.NumPages())
        oedepict.OEDrawTextToCenter(footer, text, font)

    oedepict.OEWriteReport(ofs, ext, report)

    return 0


#############################################################################
# INTERFACE
#############################################################################

InterfaceData = '''
!BRIEF [-in] <input> [-out] <output pdf> [-ringdict] <rd file>

!CATEGORY "input/output options :"

  !PARAMETER -in
    !ALIAS -i
    !TYPE string
    !REQUIRED true
    !KEYLESS 1
    !VISIBILITY simple
    !BRIEF Input filename
  !END

  !PARAMETER -out
    !ALIAS -o
    !TYPE string
    !REQUIRED true
    !KEYLESS 2
    !VISIBILITY simple
    !BRIEF Output filename
  !END

  !PARAMETER -ringdict
    !ALIAS -rd
    !TYPE string
    !REQUIRED false
    !VISIBILITY simple
    !BRIEF User-defined 2D ring dictionary
    !DETAIL
        2D ring dictionaries can be generated by the following OEChem examples:
        C++    - createringdict.cpp
        Python - createringdict.py
        Java   - CreateRingDict.java
        C#     - CreateRingDict.cs
    !END

!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
