#!/usr/bin/env python3
import sys
import os
import numpy as np
import torch
import torchani
import math
import time
import datetime

def ANI_run(eles, coordlist, method, derivatives):
    # ANI
    #device = torch.device('cpu')
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    if method == '1x':
        print('Method '+method+' was used')
        model = torchani.models.ANI1x(periodic_table_index=True).to(device).double()
    elif method == '1ccx':
        print('Method '+method+' was used')
        model = torchani.models.ANI1ccx(periodic_table_index=True).to(device).double()
    elif method == '2x':
        print('Method '+method+' was used')
        model = torchani.models.ANI2x(periodic_table_index=True).to(device).double()
    else: 
        print('Method '+method+' was unkown')
        sys.exit()
    coordinates = torch.from_numpy(coordlist).requires_grad_(True).unsqueeze(0)
    species = torch.from_numpy(eles).unsqueeze(0)
    masses = torchani.utils.get_atomic_masses(species)

    # Now let's compute energy and force:
    energy = model((species, coordinates)).energies
    if derivatives == 2: #first derivatives
        hessian = torchani.utils.hessian(coordinates, energies=energy)
        hess_mat = hessian.numpy()[0] * 0.52917721092 * 0.52917721092
    
    if derivatives == 1 or derivatives == 2: #first derivatives
        derivative = torch.autograd.grad(energy.sum(), coordinates)[0]
        force = derivative
        forcearray=force.squeeze().numpy() * 0.52917721092


    #freq, modes, fconstants, rmasses = torchani.utils.vibrational_analysis(masses, hessian, mode_type='MDU')
    #torch.set_printoptions(precision=12, sci_mode=False)
    if derivatives == 0:
        return energy.item()
    elif derivatives == 1:
        return energy.item(), forcearray
    elif derivatives == 2:
        return energy.item(), forcearray, hess_mat


time_start = time.time()

if len(sys.argv) == 7:
    filein = sys.argv[2]
    fileout = sys.argv[3]
    method='2x'
elif len(sys.argv) == 8:
    method= sys.argv[1]
    filein = sys.argv[3]
    fileout = sys.argv[4]
else:
    print('Gaussian external interface using ANI-1x, ANI-1ccx, ANI-2x')
    print('')
    print('ANI-1x(H C N O elements wB97X/6-31G(d))')
    print('ANI-1ccx(H C N O elements CCSD(T)*/CBS (CCSD(T))')
    print('ANI-2x(H C N O F S Cl elements wB97X/6-31G(d)) default')
    print('')
    print('7 inputs: method Gaussian external 6 input')
    print('Gau_External layer InputFile OutputFile MsgFile FChkFile MatElFile')
    print('More infomation see Gaussian website')
    sys.exit()


with open(filein, 'r') as f:
    [atoms, deriva, charge ,spin] = [int(s) for s in f.readline().split()]
    ele = np.zeros((atoms,), dtype = np.int)
    

    atomlist = np.zeros((atoms,3), dtype = np.float32)

    for i in range(atoms):
        
        axestr = f.readline().split()
        ele[i] = int(axestr[0])
        atomlist[i][0] = float(axestr[1])*0.52917724
        atomlist[i][1] = float(axestr[2])*0.52917724
        atomlist[i][2] = float(axestr[3])*0.52917724

    if deriva == 0:
        ene = ANI_run(ele, atomlist, method, deriva)
    elif deriva == 1:
        ene, force= ANI_run(ele, atomlist, method, deriva)
    elif deriva == 2:
        ene, force, hess = ANI_run(ele, atomlist, method, deriva)

time_end = time.time()
timstot=time_end - time_start
print(str(datetime.timedelta(seconds=timstot)))

fw = open(fileout,"w")
fw.write('%20.12E%20.12E%20.12E%20.12E\n' % (ene, 0.0, 0.0, 0.0))
if deriva == 1 or deriva == 2: #first derivatives
    for i in range(atoms):       
        fw.write('%20.12E%20.12E%20.12E\n' % (force[i][0], force[i][1],force[i][2]))
    if deriva == 2: #first derivatives
        fw.write('%20.12E%20.12E%20.12E\n' % (0,0,0))
        fw.write('%20.12E%20.12E%20.12E\n' % (0,0,0))
        for i in range(3*atoms):    
            fw.write('%20.12E%20.12E%20.12E\n' % (0,0,0))

        hess_lt = hess[np.tril_indices(3*atoms)]
        k = int(len(hess_lt)/3)
        for i in range(k):
            a = hess_lt[i*3]
            b = hess_lt[i*3+1]
            c = hess_lt[i*3+2]
            fw.write('%20.12E%20.12E%20.12E\n' % (hess_lt[i*3], hess_lt[i*3+1],hess_lt[i*3+2]))
fw.close()

    
