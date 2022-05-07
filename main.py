"""
Author: Xiao Wang (Marco liviz)
Edit time: 2022/3/8
Function: Synthesis potential function of Binary alloy
Source: https://www.ctcms.nist.gov/potentials/system/Ta/#Cu-Ta Zhou04_create_v2.f
"""
import os
import numpy as np
import math
import pandas as pd
import csv
import tkinter as tk
from tkinter.messagebox import *


def calc_electro_density(par,r,it):
    """
    calcualte electro density
    :param par:
    :param r:
    :param it:
    :return: f
    """
    f = par[it-1][1]*math.exp(-par[it-1][23]*(r/par[it-1][0]-1.0))
    f = f/(1.0+(r/par[it-1][0]-par[it-1][24])**20)
    return f
def calc_pair_potential(par,it1,it2,r):
    """
    calculate pair potential
    :param par:
    :param it1:
    :param it2:
    :param r:
    :return: psi
    """
    if it1 == it2:
        psi1 = par[it1 - 1][6] * math.exp(-par[it1 - 1][4] * (r / par[it1 - 1][0] - 1.0))
        psi1 = psi1 / (1.0 + (r / par[it1 - 1][0] - par[it1 - 1][8]) ** 20)
        psi2 = par[it1 - 1][7] * math.exp(-par[it1 - 1][5] * (r / par[it1 - 1][0] - 1.0))
        psi2 = psi2 / (1.0 + (r / par[it1 - 1][0] - par[it1 - 1][9]) ** 20)
        psi = psi1 - psi2
    else:
        psi1 = par[it1 - 1][6] * math.exp(-par[it1 - 1][4] * (r / par[it1 - 1][0] - 1.0))
        psi1 = psi1 / (1.0 + (r / par[it1 - 1][0] - par[it1 - 1][8]) ** 20)
        psi2 = par[it1 - 1][7] * math.exp(-par[it1 - 1][5] * (r / par[it1 - 1][0] - 1.0))
        psi2 = psi2 / (1.0 + (r / par[it1 - 1][0] - par[it1 - 1][9]) ** 20)
        psia = psi1 - psi2
        psi1 = par[it2 - 1][6] * math.exp(-par[it2 - 1][4] * (r / par[it2 - 1][0] - 1.0))
        psi1 = psi1 / (1.0 + (r / par[it2 - 1][0] - par[it2 - 1][8]) ** 20)
        psi2 = par[it2 - 1][7] * math.exp(-par[it2 - 1][5] * (r / par[it2 - 1][0] - 1.0))
        psi2 = psi2 / (1.0 + (r / par[it2 - 1][0] - par[it2 - 1][9]) ** 20)
        psib = psi1 - psi2
        f1 = calc_electro_density(par, r, it1)
        f2 = calc_electro_density(par, r, it2)
        psi = 0.5 * (f2/f1*psia+f1/f2*psib)
    return psi
def calc_embedding_energy(par,it,rho):
    """
    calculate embed energy
    :param par:
    :param it:
    :param rho:
    :return: emb
    """
    rhoin = np.zeros(2, dtype=float)
    rhoout = np.zeros(2, dtype=float)
    rhoin[it-1] = par[it-1][25]*par[it-1][2]
    rhoout[it-1] = par[it-1][26]*par[it-1][2]
    if rho == par[it-1][2]:
        Fm33 = par[it-1][17]
    else:
        Fm33 = par[it-1][22]
    if rho == 0:
        emb = 0.0
    elif rho < rhoin[it-1]:
        emb = par[it-1][10]+par[it-1][11]*(rho/rhoin[it-1]-1.0)+par[it-1][12]*(rho/rhoin[it-1]-1.0)**2+par[it-1][13]*(rho/rhoin[it-1]-1.0)**3
    elif rho < rhoout[it-1]:
        emb = par[it-1][14]+par[it-1][15]*(rho/(par[it-1][2]-1.0))+par[it-1][16]*(rho/(par[it-1][2]-1.0))**2+Fm33*(rho/(par[it-1][2]-1.0))**3
    else:
        emb = par[it-1][19]*(1.0-par[it-1][18]*math.log(rho/par[it-1][3]))*(rho/par[it-1][3])**par[it-1][18]
    return emb
def data_process(elementa, elementb):
    """
    process the input information & find data in EAM_code
    :param elementa:
    :param elementb:
    :return: par_elementa, par_elementb, ntypes
    """
    ntypes = 0  # No. of atom type
    par_elemta = []
    par_elemtb = []
    # data = pd.read_csv('../EAM_code', sep='\t', header=None)#pyinstaller use
    data = pd.read_csv('EAM_code', sep='\t', header=None)
    data = data[0].tolist()
    if data.index(elementa):
        par_elemta = data[data.index(elementa) + 1: data.index(elementa) + 28]
        par_elemta = list(map(float, par_elemta))
        ntypes = ntypes + 1
    else:
        os.system('exit')
    if data.index(elementb):
        par_elemtb = data[data.index(elementb) + 1: data.index(elementb) + 28]
        par_elemtb = list(map(float, par_elemtb))
        ntypes = ntypes + 1
    else:
        os.system('exit')

    return par_elemta, par_elemtb, ntypes
def calc_output(parameter):
    """
    calculate the output information
    :param parameter:
    :return: Fr, rhor, z2r, extra_par, blat
    """
    z2r = np.zeros((2000, 2, 2), dtype=float)
    Fr = np.zeros((2000, 2), dtype=float)
    rhor = np.zeros((2000, 2), dtype=float)
    blat_a = float(math.sqrt(2.0))*parameter[0][0]
    blat_b = float(math.sqrt(2.0))*parameter[1][0]
    nr = 2000
    nrho = 2000
    alatmax = max(blat_a, blat_b)
    rhoemax = max(parameter[0][2], parameter[1][2])
    rc = math.sqrt(10.0)/2.0*alatmax
    rst = 0.5
    dr = rc/(nr-1.0)
    fmax = -1.0
    for i1 in range(1, parameter[2]+1):
        for i2 in range(1, i1+1):
            if i1 == i2:
                for i in range(1, nr+1):
                    r = (i-1)*dr
                    if r < rst:
                        r = rst
                    fvalue = calc_electro_density(parameter, r, i1)
                    if fmax < fvalue:
                        fmax = fvalue
                    rhor[i-1][i1-1] = calc_electro_density(parameter, r, i1)
                    z2r[i-1][i1-1][i2-1] = r*calc_pair_potential(parameter, i1, i2, r)
            else:
                for i in range(1, nr+1):
                    r = (i-1)*dr
                    if r < rst:
                        r = rst
                    z2r[i-1][i1-1][i2-1] = r*calc_pair_potential(parameter, i1, i2, r)
                    z2r[i-1][i2-1][i1-1] = z2r[i-1][i1-1][i2-1]
    rhom = fmax
    if rhom < 2.0*rhoemax:
        rhom = 2.0*rhoemax
    if rhom < 100.0:
        rhom = 100.0
    drho = rhom/(nrho-1.0)
    for it in range(1, parameter[2]+1):
        for i in range(1, nrho+1):
            rhoF = (i-1)*drho
            if i == 1:
                rhoF = 0.0
            Fr[i-1][it-1] = calc_embedding_energy(parameter, it, rhoF)
    extra_par = np.array([nrho, drho, nr, dr, rc])
    blat = np.array([blat_a, blat_a])
    return Fr, rhor, z2r, extra_par, blat
def write_out(elementa, elementb):
    """
    define the output information format
    :param elementa, elementb:
    :return:
    """
    parameter = data_process(elementa, elementb)
    print(parameter)
    file_name = "./"+elementa+elementb+"_zhou04_wang22.eam.alloy"
    extra_info = (' DATE: 2018-03-30 CONTRIBUTOR: Xiaowang Zhou xzhou@sandia.gov and Lucas Hale lucas.hale@nist.gov '
              'CITATION:X. W. Zhou, R. A. Johnson, H. N. G. Wadley, Phys. Rev. B, 69, 144113(2004)',
              ' Generated from main.py Edit by Xiao Wang (Marco liviz) 2022/3/8', ' Fixes precision issues with older version')

    head_info = ()
    struc="fcc"
    outfile = "./"+"_zhou04.eam.alloy"
    extra_info = pd.DataFrame(extra_info)
    print(extra_info)
    f = calc_output(parameter)
    output_head0 = pd.DataFrame(np.array([2, elementa, elementb]).reshape(-1, 3))
    output_head1 = pd.DataFrame(f[3].reshape(-1, 5))
    output_head_elementa = pd.DataFrame(np.array([parameter[0][20], parameter[0][21], f[4][0], 'fcc']).reshape(-1, 4))
    output_head_elementb = pd.DataFrame(np.array([parameter[1][20], parameter[1][21], f[4][1], 'fcc']).reshape(-1, 4))
    output_part0 = pd.DataFrame(f[0].reshape(-1, 2).transpose().reshape(-1, 5))
    output_part0_0 = output_part0[0:400]
    output_part0_1 = output_part0[400:]
    output_part1 = pd.DataFrame(f[1].reshape(-1, 2).transpose().reshape(-1, 5))
    output_part1_0 = output_part1[0:400]
    output_part1_1 = output_part1[400:]
    output_part2 = pd.DataFrame(f[2].reshape(-1, 2).transpose().reshape(-1, 5))
    output_part2 = output_part2[0:1200]
    extra_info.to_csv(file_name, sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)
    output_head0.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)
    output_head1.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)
    output_head_elementa.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)
    output_part0_0.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)
    output_part1_0.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)
    output_head_elementb.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)
    output_part0_1.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)
    output_part1_1.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)
    output_part2.to_csv(file_name, sep='\t', header=None, index=False, mode='a', quoting=csv.QUOTE_NONE)


    return None
def my_GUI():
    """
    GUI interface for user
    :return:
    """
    def callback(a, b):
        if (a == b) | (a == 'elementa') | (b == 'elementb'):
            askretrycancel(title='retry/cancel dialog', message='please reselect, element is NULL or SAME!')
        else:
            write_out(v1.get(), v2.get())

    exe = tk.Tk()
    exe.title("potential syn.exe")
    frame1 = tk.Frame(exe)
    frame2 = tk.Frame(exe)
    frame3 = tk.Frame(exe)

    v1 = tk.StringVar()
    v2 = tk.StringVar()
    v1.set('elementa')
    v2.set('elementb')
    element_list = ['Cu', 'Ag', 'Au', 'Ni', 'Pd', 'Pt', 'Al', 'Pb', 'Fe', 'Mo', 'Ta', 'W', 'Mg', 'Co', 'Ti', 'Zr']
    tk.OptionMenu(frame1, v1, *element_list).pack()
    tk.OptionMenu(frame2, v2, *element_list).pack()
    tk.Button(frame3, text="create potential file！", width=20, command=lambda: callback(v1.get(), v2.get())).pack()

    frame1.pack()
    frame2.pack()
    frame3.pack()
    tk.mainloop()
    return None


if __name__ == '__main__':
    my_GUI()

