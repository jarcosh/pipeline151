# -*- coding: utf-8 -*-
"""
Created on Fri May 28 10:57:16 2021

@author: Laboratorio
"""
        
#note script is optimized for Windows, might try to make linux-compatible or linux version eventually
##to do list:
#Input species since there are small differences between b values depending on coil used; based on input search for gtab changes
#bruker2nifti text output leaves b=0 as ~21 compare with MyMapAnalyzer to see if substracting 21 to all other b values gives most similar results
#see if flag for processing can also be used to rename folder for more ease of finding relevant nifty file
#see if from dipy.data import get_fnames is better for iteration than what you got here
#magnetization transfer should be pretty easy as it's only load both as nifti, check method or even better just average values and higher mean = on
#and use NewNiftyFile = 100 * (1- (NiftyMTOn/NiftyMTOff)), save nifti
#but group often has more than 1 MT nifty file per raw data folder, if it can't be automated (if one that one, if several last) you need extra user input
#unresolved issue with curve fitting for T2 (and by extension, T1 and T2*) in separate script
#might be simple mistake due to still learning, check edx 6.0002 tutorial
#use https://stackoverflow.com/questions/9845292/a-tool-to-convert-matlab-code-to-python (last comit 2018 so might not work) on https://github.com/npnl/T2-Maps
#try to get adc to work with acquisition directions, not the sphere in default function 
#we need R^2 curve fit for each voxel, also to warn if it's very low etc.
#output parametric map with color scales using dipy and if not fsl functions; print default scale, minimum and maximum values and ask if wish to change
#maybe print to inform of progress inside each function, e.g. processing DTI of rat 48, processing DTI rat 49 etc.

from bruker2nifti.converter import Bruker2Nifti
import os
import shutil
import numpy as np
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
import dipy.reconst.dti as dti
import math
from dipy.reconst.dti import fractional_anisotropy, color_fa
#last one unused so far, must be compared with fsl script


print("RECUERDA: tus sujetos tienen que estar en la carpeta inmediamente inferior")
print("por ejemplo, las carpetas de número que nos da bruker deben estar en")
print("tudirectorio/20200825Rata 11, tudirectorio/20200825Rata 12, etc.")
print("no en una subcarpeta como tu directiorio/20200825/20200825_150640_J250820_Rat12_day13_1_1/Rata11")

def setfolder(question, default_no=True):
    choices = ' [y/N]: ' if default_no else ' [Y/n]: '
    default_answer = 'n' if default_no else 'y'
    '''
    this is a simple initial input for the user
    first question asks if the current python working directory is where
    the data is located, if not, it asks your to input data location as text
    important to note folders with raw data acquisition number 
    must be inmediate subfolders
    '''
    reply = str(input(question + choices)).lower().strip() or default_answer
    if reply[0] == 'y':
        os.chdir(os.getcwd())
    if reply[0] == 'n':
         userinput = input("¿Directorio de los datos?")
         os.chdir(userinput)
    
setfolder("¿Es el current working de python dónde están ya los datos? Y/n") 

seleccion = os.getcwd()
estudios = os.listdir(seleccion)
estudios.remove('rata.bval')
#just removed from object list to avoid issues later with shutil etc.
##remember to expand when you add mice bval and bvec!!!
estudios.remove('rata.bvec')
estudios.remove('pipeline1.51lab.py')
print(estudios)
seleccion = seleccion + '\\'

#bruker2nifti all raw data studies
for i in estudios:
    root_dir =  (seleccion)
    pfo_study_in = (seleccion + i)
    pfo_study_out = (seleccion + 'convertidos')
    '''
    takes for the given directory all folders as in
    if not out folder exists, it's created under the name 'convertidos''
    '''
    if not os.path.exists(seleccion + 'convertidos'):
        os.makedirs(seleccion + 'convertidos')
    bru = Bruker2Nifti(pfo_study_in, pfo_study_out, study_name= 'convertido' + i)
    '''
    by default everything unnecesary for our work is on false but might eventually be put as true
    verbose can also be changed to 1-2 by user preference, no issue
    '''
    bru.verbose = 0
    bru.correct_slope = True
    bru.get_acqp = False
    bru.get_method = True
    #^^DON'T TOUCH THIS!!11º ╚(•⌂•)╝
    bru.get_reco = False
    bru.nifti_version = 1
    bru.qform_code = 1
    bru.sform_code = 2
    bru.save_human_readable = True
    bru.save_b0_if_dwi = False
    bru.convert()

del root_dir

print('Obteniendo lista de directorios convertidos')

def GetterSupremo(dirName):
    listOfFile = os.listdir(dirName)
    allFiles = list()
    '''
    not complicated either, as shutil needs full directory to copy something we must first get al sub-sub directories 
    and filter from there
    '''
    for entry in listOfFile:
        fullPath = os.path.join(dirName, entry)
        if os.path.isdir(fullPath):
            allFiles = allFiles + GetterSupremo(fullPath)
        else:
            allFiles.append(fullPath)           
    return allFiles

root = seleccion + 'convertidos'

listaconvertidos = list(GetterSupremo(root))
#prob overly long but helps keep track of what is being created, filtered moved etc.
del root

methodfiles = []
niftifiles = []

mtfolders = []

for i in listaconvertidos:
   if 'method.txt' in i:
       methodfiles.append(i)
       
for i in listaconvertidos:
    if '.nii' in i:
        niftifiles.append(i)

transfer = methodfiles + niftifiles

temp = []

for string in transfer:
    bark = string.replace('convertidos', 'procesados')
    temp.append(bark)
del bark   
if not os.path.exists(seleccion + 'procesados'):
    os.makedirs(seleccion + 'procesados')    
 
dirpros = []

for string in temp:
    meow = string.replace('convertido', 'procesada')
    dirpros.append(meow)
del meow
del temp    

print('Creando carpeta de procesados y moviendo archivos clave')

ForShutil0 = []

for string in dirpros:
    quack = string.split('\\')
    ForShutil0.append(quack)
del quack

ForShutil1 = []

for i in dirpros:
    for num in range(7,8):
        oink = '\\'.join(i.split("\\")[0:num])
        ForShutil1.append(oink)
del oink

for i in ForShutil1:
    if not os.path.exists(i):
        os.makedirs(i)
#maybe rename since it gets used for more than shutil later       
ForShutil2 = []

for i in dirpros:
    for num in range(8,9):
        woof = '\\'.join(i.split("\\")[0:num])
        ForShutil2.append(woof)
del woof

for i in ForShutil2:
    if not os.path.exists(i):
        os.makedirs(i)

for i, j in zip(transfer, dirpros):
    shutil.copy(i, j)
#all folders are crated, now define functions for each type of processing

UniqueDirectories = list(set(ForShutil2))

#must use return/global if you want dti object for something

Unitchange = 1000_000


def ProcesarDTIs(dirName):
    os.chdir(dirName)
    carpetaactual = os.getcwd()
    getmethod = 'acquisition_method.txt'
    '''this just looks in the text file that bruker2nifti created and moved
    a string that is unique for dti methods
    '''
    with open(getmethod) as f:
        lines = f.readlines()
        isdti = lines[0]
    if isdti.startswith('Dti') == True:
        getnifty = carpetaactual.split("\\")[7] + '_subscan_0.nii.gz'
        dti_fname, dti_bval_fname, dti_bvec_fname = (getnifty, seleccion + 'rata.bval', seleccion + 'rata.bvec')
        data, affine = load_nifti(dti_fname)
        bvals, bvecs = read_bvals_bvecs(dti_bval_fname, dti_bvec_fname)
        '''
        note b values taken from text files that must be in same folder as subjects subfolders
        '''
        gtab = gradient_table(bvals, bvecs)
        tenmodel = dti.TensorModel(gtab) 
        tenfit = tenmodel.fit(data)
        #todo: check export of color_fa whether what we want or need other scales too
        FA = fractional_anisotropy(tenfit.evals)
        FA = FA / (math.sqrt(3))
        '''
        there appears to be more than one FA formula in the literature
        our team uses the one that can be found in, say, Alexander et. al 2007  10.1016/j.nurt.2007.05.011
        dipy uses one that can be found in Gavin P Winston 2012,  10.3978/j.issn.2223-4292.2012.12.05
        if you use the one in Nucifora et. al 2007 10.1148/radiol.2452060445 just change 3 for (4/3)
        '''
        save_nifti('FAdipy.nii.gz', FA.astype(np.float32), affine)
        MD = tenfit.md * Unitchange 
        '''seems dipy gives data as mm^2 so all diffusitivities must be multiplied to be converted to µm^2
        conversion not only for ease of reading. ImageJ will export too low values as 0'''
        save_nifti('MDdipy', MD.astype(np.float32), affine)
        AD = tenfit.ad * Unitchange 
        AD = AD * Unitchange 
        save_nifti('ADdipy', AD.astype(np.float32), affine)
        RD = tenfit.rd 
        RD = RD * Unitchange 
        save_nifti('RDdipy', MD.astype(np.float32), affine)

#no mask is used because resolution of data and number of slices isn't very computationally demanding
#ImageJ can easily remove ugly background anyway

print('procesando DTIs')


for i in UniqueDirectories:
    ProcesarDTIs(i)
    
