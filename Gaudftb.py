#!/usr/bin/python3
# -*- coding: utf-8 -*-
#lzx 20190706

import os
import sys
import numpy as np
#mod parameter file path
SKF_dir = '/share/home/yanz/Apps/DFTB/3ob-3-1'
G09_ANGS_TO_BOHR = 0.52917721092
chemical_symbols = [
    # 0
    'X',
    # 1
    'H', 'He',
    # 2
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    # 3
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    # 4
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    # 5
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    # 6
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn',
    # 7
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc',
    'Lv', 'Ts', 'Og']

def parse_ifile(ifile):
    with open(ifile) as f:
        lines = f.readlines()
        tokens = lines[0].split()
        natoms = int(tokens[0])
        deriv = int(tokens[1])
        charge = int(tokens[2])
        spin = int(tokens[3])
        coords = np.zeros((natoms,3))
        atomlist = []
        atomtype = []
        atomcode = []
        for i,line in enumerate(lines[1:1+natoms]):
            tokens = line.split()
            a = chemical_symbols[int(tokens[0])]
            c = np.array([float(tokens[1]), float(tokens[2]), float(tokens[3])])*G09_ANGS_TO_BOHR
            for j in range(3): coords[i][j] = c[j]
            atomlist.append(a)
        atomtype, atomcode = list2type(atomlist)
        #('atomlist: ',atomlist)  ['C','H','H','H','H']
        #('atomtype: ',atomtype)  ['C','H']
        #('atomcode: ',atomcode)  [ 1 , 2 , 2 , 2 , 2 ]  
    return natoms,deriv,charge,spin,coords,atomlist,atomtype,atomcode     

def list2type(atomlist):
    atomtype=[]
    atomcode=[]
    for i in range(len(atomlist)):
        if atomlist[i] not in atomtype:
            atomtype.append(atomlist[i])
        atomcode.append(atomtype.index(atomlist[i])+1)
    return atomtype,atomcode

def write_hsd(natoms,deriv,charge,spin,coords,atomlist,atomtype,atomcode):
    MAM = {
    'H'  : 'H  = "s"' ,
    'C'  : 'C  = "p"' ,
    'N'  : 'N  = "p"' ,
    'O'  : 'O  = "p"' ,
    'F'  : 'F  = "p"' ,
    'Na' : 'Na = "p"' ,
    'Mg' : 'Mg = "p"' ,
    'P'  : 'P  = "d"' ,
    'S'  : 'S  = "d"' ,
    'Cl' : 'Cl = "d"' ,
    'K'  : 'K  = "p"' ,
    'Ca' : 'Ca = "p"' ,
    'Zn' : 'Zn = "d"' ,
    'Br' : 'Br = "d"' ,
    'I'  : 'I  = "d"'
    }
    HD = {    
    'H'  : 'H  = -0.1857' ,
    'C'  : 'C  = -0.1492' ,
    'N'  : 'N  = -0.1535' ,
    'O'  : 'O  = -0.1575' ,
    'F'  : 'F  = -0.1623' ,
    'Na' : 'Na  = -0.0454' ,
    'Mg' : 'Mg  = -0.02' ,
    'P'  : 'P  = -0.14' ,
    'S'  : 'S  = -0.11' ,
    'Cl' : 'Cl  = -0.0697' ,
    'K'  : 'K  = -0.0339' ,
    'Ca' : 'Ca  = -0.0340' ,
    'Zn' : 'Zn  = -0.03' ,
    'Br' : 'Br  = -0.0573' ,
    'I'  : 'I  = -0.0433'
    }
    if  not os.path.exists("geometry.gen"): os.mknod("geometry.gen")
    with open("geometry.gen",'w') as gg:
        gg.write(str(natoms))
        gg.write('  C\n')
        for i in atomtype: gg.write(i+'  ')
        gg.write('\n')
        for i in range(natoms):
            gg.write(" %3i  %-3s   %20.12f %20.12f %20.12f" % \
                    (i+1, atomcode[i], coords[i][0], coords[i][1], coords[i][2]))
            gg.write('\n')            
    if not os.path.exists('dftb_in.hsd'): os.mknod('dftb_in.hsd')
    with open('dftb_in.hsd','w') as dftbin:
        dftbin.write('Geometry = GenFormat {\n')
        dftbin.write('<<<"geometry.gen"\n')
        dftbin.write('}\n\n')
        if deriv == 2:
            dftbin.write('Driver = SecondDerivatives {\n')
            dftbin.write('  Atoms = 1:-1\n')
            dftbin.write('  Delta = 1.0E-004\n')
            dftbin.write('}\n\n')        
        dftbin.write('Hamiltonian = DFTB {\n')
        dftbin.write('  SCC = Yes\n')
        dftbin.write('  MaxSCCIterations = 500\n')
        if deriv == 2 : dftbin.write('  SCCTolerance = 1e-7\n')
        dftbin.write('  MaxAngularMomentum={\n') 
        for i in atomtype: 
            dftbin.write('  ')
            dftbin.write(MAM[i])
            dftbin.write('\n')
        dftbin.write('  }\n')
        dftbin.write('  Charge = '+str(charge)+'\n')
        dftbin.write('  Eigensolver = QR {}\n') # change to Solvre in 19.1
        dftbin.write('  SlaterKosterFiles = Type2FileNames {\n')
        dftbin.write('    Prefix = '+SKF_dir+'\n')
        dftbin.write('    Separator = "-"\n')
        dftbin.write('    Suffix = ".skf"\n')
        dftbin.write('  }\n')
        dftbin.write('  ThirdOrderFull = Yes\n')
        dftbin.write('  HubbardDerivs = {\n')
        for i in atomtype:
            dftbin.write('  ')
            dftbin.write(HD[i])
            dftbin.write('\n')
        dftbin.write('  }\n')
        dftbin.write('  HCorrection = H5 {}\n')
        dftbin.write('}\n\n')    
        if deriv >= 1:
            dftbin.write('Analysis = {\n')
            dftbin.write('  CalculateForces = Yes\n')
            dftbin.write('}\n')
        dftbin.write('Options = {\n')
        #dftbin.write('  WriteResultsTag = Yes\n') #maybe in next version
        dftbin.write('  WriteDetailedOut = Yes\n')
        dftbin.write('}\n')

def get_energy(): # get energy from detailed.out
    with open('detailed.out') as detail:
        energy = 0.0
        lines = detail.readlines()
        imark = 0
        for i,line in enumerate(lines):
            if 'Total energy:' in line: imark = i
        tokens = lines[imark].split()
        energy = float(tokens[2])
    return energy

def get_forces(natoms): # get forces from detailed.out
    with open('detailed.out') as detail:
        forces = np.zeros((natoms,3))
        lines = detail.readlines()
        imark = 0
        for i,line in enumerate(lines):
            if 'Total Forces' in line: imark = i+1
        for iatom in range(natoms):
            tokens = lines[imark+iatom].split()
            for jcoor in range(3):
                forces[iatom][jcoor] = float(tokens[jcoor])
    return forces

def get_hessian(natoms): # get hessian from hessian.out
    with open('hessian.out') as hessian:
        lines = hessian.readlines() 
        hess_list = []
        hess_matrix = np.zeros((3*natoms,3*natoms))
        for line in lines:
            if line: tokens = line.split()
            for token in tokens: hess_list.append(token)
        for i in range(len(hess_list)): hess_matrix[i//(3*natoms)][i%(3*natoms)]=float(hess_list[i])
    return hess_matrix

def write_ofile(natoms,Outputfile,deriv): 
    with open(Outputfile,'a') as ofile:
        energy = get_energy() # always write energy
        ofile.write("%20.12f"%energy)
        for i in range(3): ofile.write("%20.12f"%0.0)
        ofile.write('\n')   
        if deriv > 0: # if forces is needed
            forces = get_forces(natoms)
            for i in range(natoms):
                for j in range(3):
                    ofile.write("%20.12f"%((-1)*forces[i][j]))
                ofile.write('\n')            
            # write polarizability
            polarizability=np.zeros((2,3))
            for line in polarizability:
                for i in range(3):
                    ofile.write("%20.12f"%line[i])
                ofile.write('\n')
            polarizability=np.zeros((3*natoms,3))
            for line in polarizability:
                for i in range(3):
                    ofile.write("%20.12f"%line[i])
                ofile.write('\n')            
        if deriv > 1: # if hessian is needed
            total = 0
            num = 0
            hessian = get_hessian(natoms)
            for i in range(3*natoms):
                for j in range(3*natoms):
                    if j<=i:
                        ofile.write("%20.12f"%hessian[i][j])
                        num += 1
                        if (num%3)==0: ofile.write('\n')
                        total += 1

if __name__ == '__main__':
    InputFile = sys.argv[2]
    OutputFile = sys.argv[3]
    kw = parse_ifile(InputFile)
    natoms = kw[0]
    deriv = kw[1]
    write_hsd(*kw)
    os.system('dftb+ dftb_in.hsd > dftbout') # don't add & !!!
    write_ofile(natoms,OutputFile,deriv)
    os.system('rm -f autotest.tag band.out charges.bin geo_end.gen geo_end.xyz')
    os.system('rm -f geometry.gen hessian.out dftbout detailed.out dftb_in.hsd')
    os.system('rm -f dftb_pin.hsd') # you can keep dftb_pin.hsd
