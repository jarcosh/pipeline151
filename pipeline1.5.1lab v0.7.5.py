# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 09:59:03 2021

@author: Laboratorio
"""

#r2, mirar y saltar labeles, no abrir varias veces pnj de colorbars

#Before reporting an error, if it's related with multiprocessing and exits with 1 during Tmap processing, check line 852-853

#This pipeline uses two modules associated with scientific publications
#If you use it, both of them would be thankful if you could cite
#For any diffusion related data, Garyfallidis E, Brett M, Amirbekian B, Rokem A, van der Walt S, Descoteaux M, Nimmo-Smith I and Dipy Contributors (2014). DIPY, a library for the analysis of diffusion MRI data. Frontiers in Neuroinformatics, vol.8, no.8.
#Link at https://www.frontiersin.org/articles/10.3389/fninf.2014.00008/full if you need it for correct citation formatting in your context
#For any T_ map data, "Multi-parametric quantitative in vivo spinal cord MRI with unified signal readout and image denoising". Grussu F, Battiston M, Veraart J, Schneider T, Cohen-Adad J, Shepherd TM, Alexander DC, Fieremans E, Novikov DS, Gandini Wheeler-Kingshott CAM; NeuroImage 2020, 217: 116884 (DOI: 10.1016/j.neuroimage.2020.116884).
#Link at https://linkinghub.elsevier.com/retrieve/pii/S1053811920303700 if you need it for correct formatting in your context
#It should also be noted that although full credit goes to Francesco Grussu, this pipeline uses slight modifications of his GetT2T2star.py script
#This includes a separate file, getT1TR, with slight modifications so it can get a exponential curve fitting on varying TR instead of decay
#Anyone interested to estimate T1 based in inversion recovery methos is welcome to use his script at https://github.com/fragrussu/MyRelax instead
#Version of script for Windows

#User modifications: line 850 if you don't wish to use avaliable CPUs -1 for multiprocessing
#Todo: should user be allowed to skip mask and parametric map previews and go 'straight to save'? (pending conversation with PI)
#Todo: if folder already has labeled subfolders, skip 

import os
import shutil
import numpy as np
import math
import multiprocessing
from bruker2nifti.converter import Bruker2Nifti
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
import dipy.reconst.dti as dti
from dipy.core.sphere import Sphere
from dipy.reconst.dti import apparent_diffusion_coef
from dipy.reconst.dti import fractional_anisotropy
from myrelax import getT2T2star
from myrelax import getT1TR
import matplotlib.cm as cm
import seaborn as sns
import matplotlib.pyplot as plt
from PIL import Image

import warnings
warnings.filterwarnings("ignore")
def EditableBrukerConverter(directory, str):
    '''
    This function is defined early on before use so an user can modify the options easily
    To understand them, you can check https://github.com/SebastianoF/bruker2nifti/wiki/Example:-use-bruker2nifti-in-a-python-(Ipython)-session
    That said, setting get_method or save_ human readable to False is extremely not recommended
    '''
    pfo_study_in = (directory + str)
    pfo_study_out = (directory + 'convertidos')
    if not os.path.exists(directory + 'convertidos'):
        os.makedirs(directory + 'convertidos')
    bru = Bruker2Nifti(pfo_study_in, pfo_study_out, study_name=('convertido' + str))
    bru.verbose = 0
    bru.correct_slope = True
    bru.get_acqp = False
    bru.get_method = True
    # ^^DON'T TOUCH THIS!!11º ╚(•⌂•)╝
    bru.get_reco = False
    bru.nifti_version = 1
    bru.qform_code = 1
    bru.sform_code = 2
    bru.save_human_readable = True
    bru.save_b0_if_dwi = False
    bru.convert()
    print('Actualizando lista de directorios convertidos')
    
    
print("RECUERDA: tus sujetos tienen que estar en la carpeta inmediamente inferior")
print("por ejemplo, las carpetas de número que nos da bruker deben estar en")
print("tudirectorio/20200825Rata 11, tudirectorio/20200825Rata 12, etc.")
print("no en una subcarpeta como tu directiorio/20200825/20200825_150640_J250820_Rat12_day13_1_1/Rata11")

def AskUser(question):
    '''
    a very simple function that requests a y or n input, correcting for all caps
    it will insist if given different input
    made simply because before any need for deeper user engagement user is typically asked
    '''
    reply1 = input(question)
    while reply1.lower() not in ('y', 'n'):
        print('Por favor introduce una de las dos opciones')
        reply1 = input(question)
    if reply1.lower() == 'y':
        reply1 = 'y'
        return reply1
    elif reply1.lower() == 'n':
        reply1 = 'n'
        return reply1

answer = AskUser("¿Es el current working directory de python dónde están ya los datos? Y/n")

if answer =='y':
    os.chdir(os.getcwd())
elif answer =='n':
    while True:
        answer2 = input("¿Directorio de los datos?")
        try:
            os.chdir(answer2)
            break
        except FileNotFoundError:
            print('Por favor, introduce un directorio válido')
#this allows user to change directory either by using python ide manually or pasting directory as text
os.chdir(answer2)
  
seleccion = os.getcwd()
seleccion = seleccion + '\\'
estudios = os.listdir(seleccion); estudios = [s for s in estudios if s not in ('procesados', 'convertidos', 'supplfiles')]
#we are left only with the raw data folders in this variable

if not os.path.exists(seleccion + 'convertidos'):
    os.makedirs(seleccion + 'convertidos')

convertdesraw = os.listdir(seleccion + 'convertidos')

#the following checks if convertidos folder is empty 
#if it is converts to nifti all raw data folders
#if it's not, checks if all raw data folders have a converted in that folder, if so notifies and skips
#if there are extra raw data folders not yet converted to nifty, asks user if they want those converted
if convertdesraw == []:
    print(estudios)
    for i in estudios:
        EditableBrukerConverter(seleccion, i)
elif len(estudios) == len(convertdesraw):
    for j in estudios:
            if ('convertido' + j) in convertdesraw:
                    print('Ya convertiste ' + j)
else:
    copyestudios = estudios[:]
    #need to add [:] so estudios itself is not mutated
    copyestudios = [s for s in copyestudios if('convertido' + s) not in convertdesraw]
    for j in copyestudios: 
        reply2 = AskUser('Tienes ' + j + 'y otros en el directorio principal pero no en convertidos, ¿convertirlos?')
        if reply2.lower() == 'y':
            for j in copyestudios:
                EditableBrukerConverter(seleccion, j)
            break
        elif reply2.lower() == 'n':
            break 
    

def GetterSupremo(dirName):
    listOfFile = os.listdir(dirName)
    allFiles = list()
    '''
    as shutil needs full directory to copy something we must first get al sub-sub directories 
    and filter from there
    '''
    for entry in listOfFile:
        fullPath = os.path.join(dirName, entry)
        if os.path.isdir(fullPath):
            allFiles = allFiles + GetterSupremo(fullPath)
        else:
            allFiles.append(fullPath)
    return allFiles


root_dir = seleccion + 'convertidos'

listaconvertidos = list(GetterSupremo(root_dir))

#all final files will be in a separate folder generated by script, titled 'procesados' by default
if not os.path.exists(seleccion + 'procesados'):
    os.makedirs(seleccion + 'procesados')
    
def CreaProcesados():
    '''this function is structured with three purposes:
    ease to edit by filtering other files generated by bruker2nifti into procesados,
    generating the full directories for shutil to copy the filtered list
    returning this filtered list as seen below'''
    methodfiles = []
    niftifiles = []
    for i in listaconvertidos:
        if 'method.txt' in i:
            methodfiles.append(i)
    for i in listaconvertidos:
        if '.nii' in i:
            niftifiles.append(i)
    transfer = methodfiles + niftifiles
    temp0 = []
    for string in transfer:
        x = string.replace('convertidos', 'procesados')
        temp0.append(x)
    dirpros = []
    for string in temp0:
        y = string.replace('convertido', 'procesada')
        dirpros.append(y)        
    print('Creando carpeta de procesados y moviendo archivos clave')
    for i in dirpros:
        if not os.path.exists('\\'.join(i.split("\\")[0:-2])):
             os.makedirs('\\'.join(i.split("\\")[0:-2]))
    for i in dirpros:
        if not os.path.exists('\\'.join(i.split("\\")[0:-1])):
            os.makedirs('\\'.join(i.split("\\")[0:-1]))
    for i, j in zip(transfer, dirpros):
        shutil.copy(i, j)
    finalist = []
    for i in dirpros:
        z = '\\'.join(i.split("\\")[0:-1])
        finalist.append(z)
    return finalist
             
FilteredList = CreaProcesados()

UniqueDirectories = list(set(FilteredList))
#this gets rid of any duplicates, hence name
print('Etiquetando las carpetas relevantes con su método')


def SlapDTI(dirname):
    '''
    all these functions have identical purpose, they parse the method files
    in order to rename the folder for ease of location later and feeding
    into further processing functions
    originally it was one single big function, but couldn't fix a bug
    that gave permission errors on a non-deterministic folder and didn't edit
    despite specific call for fclose()
    expected this won't happen on linux version
    '''
    
    try:
        os.chdir(dirname + '\\')
        carpeta = os.getcwd()
        getflag = 'acquisition_method.txt'
        with open(getflag, 'r') as f:
            check = f.readlines()
            f.close()
        if check[0] == 'DtiEpi':
            slap = '\\'.join(carpeta.split('\\')[0:-1]) + '\\DT' + carpeta.split("\\")[-1]
            os.chdir(seleccion)
            os.rename(dirname, slap)
            dirname.replace (dirname, slap)
    except FileNotFoundError:
        placeholder = 'skipping renamed folder'

for i in UniqueDirectories:
    i.replace('\\\\', '\\')
    SlapDTI(i)
    

def SlapT2(dirname):
    try:
        os.chdir(dirname)
        carpeta = os.getcwd()
        getflag = 'acquisition_method.txt'
        getmethod = carpeta.split("\\")[-1] + '_method.txt'
        with open(getmethod, 'r') as f:
            check2 = f.readlines()
            T2flag = []
            for i in check2:
                if i.startswith('EffectiveTE ='):
                    T2flag.append(i)
            f.close()
            with open(getflag, 'r') as f:
                check = f.readlines()
                f.close()
        if check[0] == 'MSME':
            if len(T2flag[0]) > 35:
                slap = '\\'.join(carpeta.split('\\')[0:-1]) + '\\T2' + carpeta.split("\\")[-1]
                os.chdir(seleccion)
                os.rename(dirname, slap)
    except FileNotFoundError:
        placeholder = 'skipping renamed folder'
        
for i in UniqueDirectories:
    i.replace('\\\\', '\\')
    SlapT2(i)
    
def SlapMT(dirname):
        try:
            os.chdir(dirname)
            carpeta = os.getcwd()
            getflag = 'acquisition_method.txt'
            getmethod = carpeta.split("\\")[-1] + '_method.txt'
            with open(getmethod, 'r') as f:
                check2 = f.readlines()
                MTflag = []
                for i in check2:
                    if i.startswith('MagTransOnO'):
                        MTflag.append(i)
                f.close()
                with open(getflag, 'r') as f:
                    check = f.readlines()
                    f.close()
            if check[0] == 'MSME':
                if (MTflag[0]).startswith('MagTransOnOff = On'):
                    slap = '\\'.join(carpeta.split('\\')[0:-1]) + '\\MT' + carpeta.split("\\")[-1]
                    os.chdir(seleccion)
                    os.rename(dirname, slap)
        except FileNotFoundError:
            placeholder = 'skipping renamed folder'
        
for i in UniqueDirectories:
    i.replace('\\\\', '\\')
    SlapMT(i)
    

def SlapM0(dirname):
    try:
        os.chdir(dirname)
        carpeta = os.getcwd()
        getflag = 'acquisition_method.txt'
        getmethod = carpeta.split("\\")[-1] + '_method.txt'
     
        with open(getmethod, 'r') as f:
            check2 = f.readlines()
            T2flag = []
            M0flag = []
            for i in check2:
                if i.startswith('EffectiveTE ='):
                    T2flag.append(i)
            for i in check2:
                if i.startswith('DigFilte'):
                    M0flag.append(i)
            f.close()
        with open(getflag, 'r') as f:
            check = f.readlines()
            f.close()
        if check[0] == 'MSME':
            if len(T2flag[0]) < 50:
                if (M0flag[0]).startswith('DigFilter = Digital_Medi'):
                    #unorthodox workaround to avoid some T1 maps to show up as M0, will try to find something else
                    slap = '\\'.join(carpeta.split('\\')[0:-1]) + '\\M0' + carpeta.split("\\")[-1]
                    os.chdir(seleccion)
                    os.rename(dirname, slap)
    except FileNotFoundError:
        placeholder = 'skipping renamed folder'
                    
for i in UniqueDirectories:
    i.replace('\\\\', '\\')
    SlapM0(i)
    

def SlapT1(dirname):
    try:
        os.chdir(dirname)
        carpeta = os.getcwd()
        getflag = 'acquisition_method.txt'
        with open(getflag, 'r') as f:
            check = f.readlines()
            f.close()
        if check[0] == 'RAREVTR':
            slap = '\\'.join(carpeta.split('\\')[0:-1]) + '\\T1' + carpeta.split("\\")[-1]
            os.chdir(seleccion)
            os.rename(dirname, slap)
    except FileNotFoundError:
        placeholder = 'skipping renamed folder'
        
for i in UniqueDirectories:
    i.replace('\\\\', '\\')
    SlapT1(i)
    

def SlapT2star(dirname):
    try:
        os.chdir(dirname)
        carpeta = os.getcwd()
        getflag = 'acquisition_method.txt'
        with open(getflag, 'r') as f:
            check = f.readlines()
            f.close()
        if check[0] == 'MGE':
            slap = '\\'.join(carpeta.split('\\')[0:-1]) + '\\TS' + carpeta.split("\\")[-1]
            os.chdir(seleccion)
            os.rename(dirname, slap)
    except FileNotFoundError:
        placeholder = 'skipping renamed folder'
        
for i in UniqueDirectories:
    i.replace('\\\\', '\\')
    SlapT2star(i)
    
Labels = GetterSupremo(seleccion + 'procesados')
Labelsfiltered =[]

for i in Labels:
    x = '\\'.join(i.split('\\')[0:-1])
    Labelsfiltered.append(x)

#once again we filter to have unique directories 
Tobedone = list(set(Labelsfiltered))

#easiest way to rotate seaborn output is to simply rotate the original array
def rotated(array_2d):
    '''rotates an array coming from a nifti file 90 degrees, last line keeps it in array form
    note later functions still need to rotate x axis but your data might not have to
    this is the easiest way to deal with matplotlib/seaborn not getting the orientation 
    you want by default
    '''
    list_of_tuples = zip(*array_2d[::-1])
    return [list(elem) for elem in list_of_tuples]

def CreateMasks (dirname):
    ### pending to warn if R^2 below input threshold
    '''will produce figures to allow user to preview if the mask on a treshold method is okay
    if not, requests a new treshold value
    pending to see other alternatives including dipy's, since this is certainly not comfortable for data with dozens of slices
    '''
    os.chdir(dirname)
    carpetaactual = os.getcwd()
    try:
        (nfile, affine) = load_nifti((carpetaactual + '\\' + (carpetaactual.split("\\")[-1])[2:] + '_subscan_0.nii.gz'))
    except FileNotFoundError: 
        (nfile, affine) = load_nifti((carpetaactual + '\\' + (carpetaactual.split("\\")[-1])[2:] + '.nii.gz'))
    nfile = nfile[:,:,:,0]   
    nsl = nfile.shape[2]
    copyformask = nfile.copy()
    nfilter = 0
    method = (carpetaactual.split("\\")[-1])[0:2]
    if method == 'T2': 
            nfilter = 40
    if method == 'TS':
            nfilter = 20
    if method == 'DT':
            nfilter = 1000
    if method == 'T1':
            nfilter = 30
    copyformask[copyformask <=nfilter] = 0
    copyformask[copyformask >nfilter] = 1
    slicedmask =[]        
    for j in range(0,nsl):
        submask = copyformask[:,:,j]
        slicedmask.append(submask)
    slicednif = []
    for z in range(0,nsl):
        sl = nfile[:,:,z]
        slicednif.append(sl)
    dims = (8.7, 5.27)
    for i in range(0, nsl):
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=dims)
        sns.heatmap((rotated(slicedmask[i])), cmap= cm.gray, cbar= False, xticklabels= False, yticklabels= False, ax=ax1).invert_xaxis()
        plt.title('Mostrando selección para slice ' + str(i + 1) + ', ¿es correcto?', loc = 'right', fontsize=20)
        sns.heatmap((rotated(slicednif[i])), cmap = cm.gray, cbar = False, xticklabels= False, yticklabels= False, ax=ax2).invert_xaxis()
        fig.subplots_adjust(wspace=0.01)
        plt.show()
        plt.clf()
    while True:
        print(carpetaactual.split("\\")[-1])
        approval = AskUser('Responde Y si el archivo mask estaría correcto, N si quieres cambiar el filtro ')
        if approval == 'n':
            print('Por favor introduce el nuevo valor de filtrado. Valor usado: ' + str(nfilter))
            while True:
                try:
                    nfilter = int(input('Recuerda, valores mayores si coge demasiado fondo blanco y menores si hay puntos negros en región de interés. Nuevo valor: ' ))
                    break
                except ValueError:
                    print('Por favor introduce una cifra')
                   
            copyformask = nfile.copy()
            copyformask[copyformask <= nfilter] = 0
            copyformask[copyformask > nfilter] = 1
            slicedmask =[]        
            for j in range(0,nsl):
                submask = copyformask[:,:,j]
                slicedmask.append(submask)
            for i in range(0, nsl):
                fig, (ax1, ax2) = plt.subplots(1,2, figsize=dims)
                sns.heatmap((rotated(slicedmask[i])), cmap= cm.gray, cbar= False, xticklabels= False, yticklabels= False, ax=ax1).invert_xaxis()
                plt.title('Mostrando selección para slice ' + str(i + 1) + ', ¿es correcto?', loc = 'right', fontsize=20)
                sns.heatmap((rotated(slicednif[i])), cmap = cm.gray, cbar = False, xticklabels= False, yticklabels= False, ax=ax2).invert_xaxis()
                fig.subplots_adjust(wspace=0.01)
                plt.show()
                plt.clf()
        elif approval == 'y':
            break
    save_nifti((method + 'mask'), copyformask.astype(np.float32), affine)

while True:
    tasks=int(input('¿Cuántos procesamientos quieres hacer? (el máximo son 5)'))
    if tasks <6:
        break
    
ToDo= []

while tasks > 0:
    '''asking user which processing out of the options they want to do, then fills lists with thsoe with a counter'''
    Process = input("Introduce procesamiento " + str(tasks) + " Opciones válidas: DTI, MT, T2, TS, T1 ")
    if Process in ['DTI', 'MT', 'T1', 'T2', 'TS']:
        ToDo.append(Process)
        tasks -= 1
    else: 
        print('No introduciste procesamiento válido. Para salir, introduce b')
        ##ahora mismo hay que darle dos veces 
        Process = input("Introduce procesamiento " + str(tasks) + " Opciones válidas: DTI, MT, T2, TS, T1 ")
        if Process == 'b':
            break

KN95 = []
FFP2 = []
FFP3 = []
P100 = []


for i in Tobedone:
    if 'T2' in i and 'T2' in ToDo:
        i.replace('\\\\', '\\')
        FFP2.append(i)
for i in Tobedone:
    if 'TS' in i and 'TS' in ToDo:
        i.replace('\\\\', '\\')
        FFP3.append(i)
for i in Tobedone:
    if 'DT' in i and 'DTI' in ToDo:
        i.replace('\\\\', '\\')
        P100.append(i)
for i in Tobedone:
    if 'T1' in i and 'T1' in ToDo:
        i.replace('\\\\', '\\')
        KN95.append(i)      

#only those in the previously givenlist of tasks get masks generated   
if KN95 != []:
    for i in KN95:
        i.replace('\\\\', '\\')
        CreateMasks(i)
        
if FFP2 != []:
    for i in FFP2:
        i.replace('\\\\', '\\')
        CreateMasks(i)
        
if FFP3 != []:
    for i in FFP3:
        i.replace('\\\\', '\\')
        CreateMasks(i)
        
if P100 != []:
    for i in P100:
        i.replace('\\\\', '\\') 
        CreateMasks(i)
        
Total = KN95 + FFP2 + FFP3 + P100


def GibNifty(dirname):
    '''will get list of target nifti files'''
    funky = []
    for i in Total:
        os.chdir(i)
        carpetaactual = os.getcwd()
        try: 
            me_nifti = carpetaactual + '\\' + (carpetaactual.split("\\")[-1])[2:] + '_subscan_0.nii.gz'
            funky.append(me_nifti)
        except FileNotFoundError: 
            me_nifti = carpetaactual + '\\' + (carpetaactual.split("\\")[-1])[2:] + 'nii.gz'
            funky.append(me_nifti)
    return (funky)

funky = GibNifty(Total)

isfirst = True
ugib1 = None
ugib2 = None
ugib3 = 'nipy_spectral'

def previewheatmap(array3d, titlestring, robust = isfirst, vmin = ugib1, vmax = ugib2, cmap = ugib3): 
    '''to avoid repeating code, is called in saveheatmap twice, previous global variables necessary
    for function to be defined, but will catch any redefinition within saveheatmap calls'''
    method = titlestring 
    counter = 1 
    sls = array3d.shape[2] 
    slicing = [] 
    for j in range(0,sls): 
        x = array3d[:,:,j] 
        slicing.append(x)  
    for i in slicing: 
        mlem = plt.subplots() 
        sns.heatmap(rotated(i), robust = robust, vmin = vmin, vmax = vmax, cmap = cmap).invert_xaxis() 
        plt.title(method + ' slice ' + str(counter)) 
        mlem[1].set(xticklabels=[]) 
        mlem[1].set(yticklabels=[]) 
        counter += 1 
        plt.show() 
        plt.clf()


           
def saveheatmap (array3d, titlestring):
    '''similar to previous mask function, uses robust and default colormap to preview figures, ask user if fine like this
    if not, vmin, vmax and colormap are requested, opening repertoire from a .png in supplfiles'''
    carpeta = os.getcwd()
    method = titlestring
    counter = 1
    sls = array3d.shape[2]
    slicing = []
    for j in range(0,sls):
        x = array3d[:,:,j]
        slicing.append(x) 
    isfirst = True
    ugib1 = None
    ugib2 = None
    ugib3 = 'nipy_spectral'
    previewheatmap(array3d, titlestring, robust = isfirst, vmin = ugib1, vmax = ugib2, cmap = ugib3)
    print(carpeta.split('\\')[-2][9:])
    beforesave = AskUser("¿Guardar como está?")
    if beforesave == 'y':
        counter = 0
        for i in slicing:
            mlem = plt.subplots()
            sns.heatmap(rotated(i), robust = True, cmap = 'nipy_spectral').invert_xaxis()
            plt.title(titlestring + ' slice ' + str(counter))
            mlem[1].set(xticklabels=[])
            mlem[1].set(yticklabels=[])
            counter += 1
            mlem[0].savefig(method + 'slice' + str(counter) + '.png')
            plt.clf()
    elif beforesave == 'n':
        while True: 
            isfirst=False
            ugib1 = float(input('Escribe valor mínimo'))
            ugib2 = float(input('Escriba valor máximo'))
            img = Image.open(seleccion + 'supplfiles' + '\\colorbars.png')
            img.show()
            ugib3 = input('selecciona escala de color (presa atención a mayúsculas/minúsculas)')
            if ugib3 == '':
                ugib3 = 'nipy_spectral'
            counter = 1
            previewheatmap(array3d, titlestring, robust = isfirst, vmin = ugib1, vmax = ugib2, cmap = ugib3)
            print(carpeta.split('\\')[-2][9:])
            uadj = AskUser ('¿Quieres guardar ahora?')
            if uadj == 'y':
                counter = 0
                for i in slicing:
                    mlem = plt.subplots()
                    sns.heatmap(rotated(i), robust = False, vmin = ugib1, vmax =ugib2, cmap = ugib3).invert_xaxis()
                    plt.title(titlestring + ' slice ' + str(counter))
                    mlem[1].set(xticklabels=[])
                    mlem[1].set(yticklabels=[])
                    counter += 1
                    mlem[0].savefig(method + 'slice' + str(counter) + '.png')
                    plt.clf()
                break

#encadena bien procesar mt con calcula mt
temp = []
def ProcesarMT(directorylist):
    '''due to specific arrangements in our lab, correct nifti MT file is not predictable but user should know
    which one it is if present during adquisition
    this generates an extra folder with short filename and asks user to move the correct MT on MT off files there'''
    for i in directorylist:
        x = ('\\').join(i.split('\\')[0:-1])
        temp.append(x)
        meow = list(set(temp))
        counter = len(meow)
    while counter > 0:      
        for i in meow:
            if not os.path.exists(i + '\\MTCarpeta'):
                os.makedirs(i + '\\MTCarpeta')
            print('Procesamiento pausado en ' + str(i)) 
            print('Por favor, mueve los archivos de las carpetas MT y M0 de tus sujetos a la carpeta creada')
            print('Renombra los archivos MT.nii y M0.nii, sin dejarte el nii')  
            qmt = AskUser('Pulsa y cuando esté listo o n para no hacer este procesamiento ' )
            if qmt == 'y':
                counter -=1
                print('')
                print('Pasando a siguiente carpeta')
                print('')
            elif qmt == 'n':
                break
        return meow
##type error continua        
listMT = ProcesarMT(UniqueDirectories)


def CalculaMTR (dirname):
    '''after user request, estimates MT ratio, generates nifti file, calls save heatmap funciton'''
    os.chdir(dirname + '\\MTCarpeta')
    (MTon, affineo) = load_nifti(os.getcwd() + '\\MT.nii.gz')
    (MToff, affinef) = load_nifti(os.getcwd() + '\\M0.nii.gz')
    #arreglaestoconmáscara
    MTon[MTon <=20] = 0
    MToff[MToff <=20] = 0
    MTfinal = 100 * (1 - (MTon/MToff))
    MTfinal[MTfinal <0] = 0
    save_nifti('MTR.nii.gz', MTfinal.astype(np.float32), affineo)
    saveheatmap(MTfinal, 'MT')
    
for i in listMT:
    CalculaMTR(i)
    #to do: add raise nonetype error to avoid issue if user skips with n above

        
def AskBTable(question):
     '''
     to allow for flexible adc b values and directions in acquisition and still save files correctly
     number of basal images and number of b values for each direction is asked
     this avoids issues like 4th adc dimension with each direction duplicated
     '''
     DTIinfo = []
     Bdir = input("¿Número de b valores para cada dirección?" )
     Bdir = int(Bdir)
     DTIinfo.append(Bdir)
     BI = input("¿Número de imágenes basales?" )
     BI = int(BI)
     DTIinfo.append(BI)
     return DTIinfo
 

for i in funky:
    if 'DT' in i:
        DTIinfo = AskBTable("Solicitando información al usuario para imágenes DTI")
        Bdir = DTIinfo[0]
        BI = DTIinfo[1]     
        break
        

def ProcesarDTIs(filelist):
    os.chdir('\\'.join(filelist.split('\\')[0:-1]))
    #pending input from user to choose bval bdir files from repertoire in supplfiles, in our case rat and mice have different ones due to different coils
    carpetaactual = os.getcwd()
    try:
        getnifty = (carpetaactual.split("\\")[-1])[2:] + '_subscan_0.nii.gz'
    except FileNotFoundError:
        getnifty = (carpetaactual.split("\\")[-1])[2:] + '.nii.gz'
    dti_fname, dti_bval_fname, dti_bvec_fname = (getnifty, seleccion + ('supplfiles' + '\\' + 'rata.bval'), seleccion + ('supplfiles' + '\\' + 'rata.bvec'))
    data, affine = load_nifti(dti_fname)
    bvals, bvecs = read_bvals_bvecs(dti_bval_fname, dti_bvec_fname)
    (maskhere, affine) = load_nifti(carpetaactual + '\\' + 'DTmask.nii')
    maskhere = np.stack((maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere, maskhere), -1)
    #quite ugly but fastest way found to make dimensions match and allow for masking based on basal image with minimal relevant voxel loss
    data = data * maskhere
    '''
    note b values taken from text files that must be in same supplfiles subfolders
    '''
    gtab = gradient_table(bvals, bvecs)
    tenmodel = dti.TensorModel(gtab, fit_method='NLLS')
    tenfit = tenmodel.fit(data)
    ##todo: check export of color_fa whether what we want or need other scales too
    FA = fractional_anisotropy(tenfit.evals)
    FA = FA / (math.sqrt(3))
    '''
    there appears to be more than one FA formula in the literature
    our team uses the one that can be found in, say, Alexander et. al 2007  10.1016/j.nurt.2007.05.011
    dipy uses one that can be found in Gavin P Winston 2012,  10.3978/j.issn.2223-4292.2012.12.05
    if you use the one in Nucifora et. al 2007 10.1148/radiol.2452060445 just change 3 for (4/3)
    '''
    os.makedirs(carpetaactual + '\\FA')
    save_nifti('FAdipy.nii', FA.astype(np.float32), affine)
    os.chdir(carpetaactual + '\\FA'); saveheatmap(FA, 'FA');os.chdir(carpetaactual)
    Unitchange = 1_000_000
    MD = tenfit.md * Unitchange 
    '''seems dipy gives data as mm^2 so all diffusitivities must be multiplied to be converted to µm^2
     conversion not only for ease of reading. ImageJ will export too low values as 0'''
    os.makedirs(carpetaactual + '\\MD')
    save_nifti('MDdipy', MD.astype(np.float32), affine)
    os.chdir(carpetaactual + '\\MD'); saveheatmap(MD, 'MD');os.chdir(carpetaactual)
    AD = tenfit.ad * Unitchange 
    save_nifti('ADdipy', AD.astype(np.float32), affine)
    os.makedirs(carpetaactual + '\\AD')
    os.chdir(carpetaactual + '\\AD'); saveheatmap(AD, 'AD');os.chdir(carpetaactual)
    os.makedirs(carpetaactual + '\\RD')
    RD = tenfit.rd 
    RD = RD * Unitchange 
    save_nifti('RDdipy', MD.astype(np.float32), affine)
    os.chdir(carpetaactual + '\\RD'); saveheatmap(RD, 'RD');os.chdir(carpetaactual)
    ''' now this is optimized for our method, using directions and b values during acquisition instead
    of sphere directions as is default in '''
    my_sphere = Sphere(xyz=gtab.bvecs[~gtab.b0s_mask])
    adc = apparent_diffusion_coef(tenfit.quadratic_form, my_sphere)
    '''we have a fourth dimension "dupplicated" by default
    luckily this is easy to fix, previous input was made to ensure it was
    solved despite acqusition methods changing'''
    b = (len(bvecs[:,:]) - BI)
    #4th dimension minus basal images
    funkynator = [*range(1,b,Bdir)]
    #skips based on previous input. in this case for 14, 1, 3, 5... in others avoids duplicates depending on input
    adc = adc[:,:,:,funkynator]
    adc = adc * Unitchange 
    save_nifti('ADCdipy', adc.astype(np.float32), affine)
    #pending (asking PI) how to organize .png output into folders
    printadc = []
    for d in range (0, int(b/Bdir)):
        direc = adc[:,:,:,d]
        printadc.append(direc)
    counter = 1
    for n in printadc:
        os.makedirs (carpetaactual + '\\ADCdir' + str(counter))
        os.chdir(carpetaactual + '\\ADCdir' + str(counter)); saveheatmap(n, ('ADC dir' + str(counter)));  os.chdir(carpetaactual)
        counter += 1
        

print('procesando DTIs')

##scripts originales T1-T2 tardan menos?

for i in funky:
    if 'DT' in i:
        ProcesarDTIs(i)

'''all the following are different forms by PI request to generate the text files with echo or repetition times
you can make them get it directly from method files: just use i.startswith('EffectiveTE ='): as in line ~250 and turn all after = into array
similar por repetition times
'''   
def AskTR(question):
     TRlist = []
     TotalT1 = int(input("¿Cuántos tiempos de repetición se usaron en la acquisición? "))
     for i in range(0,TotalT1):
         T1 = input("Introduce tiempo " + str(i+1) + " de repetición ")
         TRlist.append(T1)
     return TRlist
  
def AskTE(question):
     startte = int(input("¿Primer tiempo de eco en T2? "))
     finaltte = int(input("¿Último tiempo de eco en T2? "))
     interval = int(input("¿Cuántos tiempos de eco hay en T2? "))
     TElist = list(range(startte, (finaltte + startte), finaltte//interval))
     return TElist
 

def AskTEstars(question):
    #this one has to be done in a different way to avoid floating point errors
     startte = float(input("¿Primer tiempo de eco en T2 estrella? "))
     interval = float(input("¿Cuántos tiempos de eco hay en T2 estrella? "))
     jumpt = float(input("¿Separación entre tiempos de eco? "))
     counter = 1
     TEstarlist = [startte]
     while counter < interval :
         TEstarlist.append(startte + (counter * jumpt))
         counter += 1
     return TEstarlist 
    
for i in funky:
    if 'T2' in i:
        try:
            f = open(seleccion + 'supplfiles\\' + "TiemposEco.txt", "x")
        except FileExistsError:
            print ('Usando archivo de sujeto anterior para tiempos de eco T2')
            te_file = seleccion + 'supplfiles\\' + "TiemposEco.txt"
            break
        TElist = AskTE('Comprobando si hay archivo previo')
        TEs = (seleccion + 'supplfiles\\' + 'TiemposEco.txt')
        with open(TEs, 'w') as f:
            for i in TElist:
                content = (str(i) + ' ')
                f.write(content)
                te_file = TEs
            break
        
for i in funky:
    if 'TS' in i:
        try:
            f = open(seleccion + 'supplfiles\\' + "TiemposEcoStar.txt", "x")
        except FileExistsError:
            print ('Usando archivo de sujeto anterior para tiempos de eco T2 estrella')
            ts_file = seleccion + 'supplfiles\\' + "TiemposEcoStar.txt"
            break
        TSlist = AskTEstars('Comprobando si hay archivo previo')
        TSs = (seleccion + 'supplfiles\\' + 'TiemposEcoStar.txt')
        with open(TEs, 'w') as f:
            for i in TSlist:
                content = (str(i) + ' ')
                f.write(content)
                ts_file = TSs
      
for i in funky:
    if 'T1' in i:
        try:
            f = open(seleccion + 'supplfiles\\' + "TiemposRepeticion.txt", "x")
        except FileExistsError:
            print ('Usando archivo de sujeto anterior para tiempos de repetición T1')
            tr_file = seleccion + 'supplfiles\\' + "TiemposRepeticion.txt"
            break
        TRlist = AskTR('Comprobando si hay archivo previo')
        TRs = (seleccion + 'supplfiles\\' + 'TiemposRepeticion.txt')
        with open(TRs, 'w') as f:
            for i in TRlist:
                content = (str(i) + ' ')
                f.write(content)
                tr_file = TRs


if 'T1' or 'T2' or 'TS' in Tobedone:
    algo = 'nonlinear'
    ncpu = multiprocessing.cpu_count() - 1


''' IF YOU GET AN ERROR BECAUSE YOU RAN THE SCRIPT DIRECTLY INSTEAD OF PASTING TO CONSOLE, COPY AND PASTE THIS LAST SECTION AND PRESS ENTER '''
'''SI TE DA ERROR POR HABER ARRANCADO EL SCRIPT DIRECTAMENTE EN VEZ DE COPIARLO EN LA CONSOLA, COPIA Y PEGA ESTA ÚLTIMA SECCIÓN EN LA CONSOLA Y PULSA INTRO '''
#will later attempt to build dictionary to turn this into just one repetition of the code instead of three

for i in funky:
    if 'T2' in i:
        os.chdir((('\\').join(i.split("\\")[0:-1])) + '\\')
        getT2T2star.TxyFitME(i, te_file, (('\\').join(i.split("\\")[0:-1]) + '\\'), algo, ncpu, ('\\').join(i.split("\\")[0:-1]) + '\\' + 'T2mask.nii')
        if not os.path.exists('otrosarchivos'):
            os.makedirs('otrosarchivos')
        otherfiles = os.listdir((('\\').join(i.split("\\")[0:-1])))
        tomove = []
        for z in otherfiles:
            if z.endswith('ME.nii'):
                tomove.append(z)
                for n in tomove:
                    if n.endswith('xyME.nii'):
                        (T2map, affine) = load_nifti(n)
                        saveheatmap(T2map, 'T2')
                        tomove.remove(n)
        for g in tomove:
            shutil.move(((('\\').join(i.split("\\")[0:-1])) + '\\' + g), ((('\\').join(i.split("\\")[0:-1])) + '\\' + 'otrosarchivos\\' + g))
    elif 'TS' in i:
        os.chdir((('\\').join(i.split("\\")[0:-1])) + '\\')
        getT2T2star.TxyFitME(i, ts_file, (('\\').join(i.split("\\")[0:-1]) + '\\'), algo, ncpu, ('\\').join(i.split("\\")[0:-1]) + '\\' + 'TSmask.nii')
        if not os.path.exists('otrosarchivos'):
            os.makedirs('otrosarchivos')
        otherfiles = os.listdir((('\\').join(i.split("\\")[0:-1])))
        tomove = []
        for z in otherfiles:
            if z.endswith('ME.nii'):
                tomove.append(z)
                for n in tomove:
                    if n.endswith('xyME.nii'):
                        (T2starmap, affine) = load_nifti(n)
                        saveheatmap(T2starmap, 'T2star')
                        tomove.remove(n)
        for g in tomove:
            shutil.move(((('\\').join(i.split("\\")[0:-1])) + '\\' + g), ((('\\').join(i.split("\\")[0:-1])) + '\\' + 'otrosarchivos\\' + g))
    elif 'T1' in i:
        os.chdir((('\\').join(i.split("\\")[0:-1])) + '\\')
        getT1TR.TxyFitME(i, tr_file, (('\\').join(i.split("\\")[0:-1]) + '\\'), algo, ncpu, ('\\').join(i.split("\\")[0:-1]) + '\\' + 'T1mask.nii')
        if not os.path.exists('otrosarchivos'):
            os.makedirs('otrosarchivos')
        otherfiles = os.listdir((('\\').join(i.split("\\")[0:-1])))
        tomove = []
        for z in otherfiles:
            if z.endswith('ME.nii'):
                tomove.append(z)
                for n in tomove:
                    if n.endswith('xyME.nii'):
                        (T1map, affine) = load_nifti(n)
                        saveheatmap(T1map, 'T1')
                        tomove.remove(n)
        for g in tomove:
            shutil.move(((('\\').join(i.split("\\")[0:-1])) + '\\' + g), ((('\\').join(i.split("\\")[0:-1])) + '\\' + 'otrosarchivos\\' + g))