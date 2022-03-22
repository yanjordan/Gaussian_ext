# Gaussian_ext
## Fortran codes

### gau_xtb_V2.f90 gau_xtb_V3.f90
    This is a Gaussian external to xTB code

### GauOPENMP.f90  
    This is a program to run Gaussian jobs on parralel on one node 
    Please define the Gauexe before: g16 or g09
    Usage: gjf1 gjf2 ... gjfn

### Gau_ANI.py
    B3LYP/6-31g(d) EmpiricalDispersion=GD3 --> external='Gau_ANI.py 2x’
    1x	ANI-1x(H C N O elements wB97X/6-31G(d))
	1ccx	ANI-1ccx(H C N O elements CCSD(T)*/CBS (CCSD(T))
	2x	ANI-2x(H C N O F S Cl elements wB97X/6-31G(d)) default

### Gaudftb.py
    This is a Gaussian external to DFTB+
    modify the parameter file path if needed
    
### gau2mpxyz.f90 & call_molpro
    This is a Gaussian external to xyz files (for Molpro in my case)

### GauCenters_list.f90
    split to two or more frag and computed by Gaussian. Fragments defined by a parameter files.
    see [QR_ONIOM_Mode](https://github.com/yanjordan/ONIOM_QR_mod) Manual.pdf for details

    Usage: External='GauCenters list parafile'
    parafile is a text file to define atom lists and other information.
    2 #number of fragments
    1-17,36 #atom list (including link atoms) of fragment 1
    18-35 #atom list (including link atoms) of fragment 2
    head1.txt tail11.txt #head and tail files for Gaussian computation of fragment 1
    head2.txt tail12.txt #head and tail files for Gaussian computation of fragment 2
    Here, the atoms list is dfferent to the label in gjf file of whole system. Thus, the tool to obatin the correct
    atomic label is implemented in GauCenters list .
    Usage: External='GauCenters list'
    Call the external program GauCenters list in Gaussian without extra parameter, the job will be terminated
    (error) with generating testtmp.gjf file, which are the high layer with link atoms. Therefore, users can use
    the Atom selection tools in GaussianView to obtain the atom list string.

### Gau_orca.f90
    This is a Gaussian external to ORCA, need a head file to define ORCA parameter 

### Gau_srcf.f90
    extract part of system using solvation model, since large system with solvation model requires memory
    The result seems not good

### Getscanxyz.sh
    get coordination of scan (xyz format)

in continue...

