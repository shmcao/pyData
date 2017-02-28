# -*- coding: utf-8 -*-
"""
Data manipulation procedures,
suits for data in Chenlab.

TODO: Add more comments and reconstruct some codes.
TODO: add some more procedure.
v0.9 @Feb 21 2017, by shmcao.

"""
"""
Created on Sat Aug 27 21:54:31 2016
Modified on Mon Nov 14 00:18:02 2016
Modified on Tue Feb 21 10:22:58 2017
@author: shmcao
"""

import re
import os,shutil
import numpy as np
import scipy as sp
from scipy import signal
#=====================================array manipulation API=======================
def zipper(array,value):
    '''
    This method puts val in array.
    The type of array is numpy nd array.
    The type of value is numpy (n-1)d array/element.
    '''
    return np.concatenate((array,[value]))

def insert(array,val,index):
    '''
    USAGE: array=insert(array,val,...)
    
    insert val in array[index]
    the type of array is numpy nd array.
    index can be a number or a list
    '''
    if(len(array)==0):return np.array([val])
    if(index>len(array)-1 or index<-len(array)):
        print("[insert]Index error! index="+str(index)+" len(array)="+str(len(array)))        
        return array
    if(index<0):index+=len(array)+1
    return np.concatenate((array[:index],[val],array[index:]))
    
def delete(array,index):
    '''
    USAGE: array=delete(array,index)
    
    delete array[index] in array
    the type of array is numpy nd array.
    index can be a number or a list
    '''
    array=np.array(array)
    res=np.array([])
    print(array.shape)
    if(type(index)==int):index=[index]
    for i in range(len(array)):
        if(not i in index):res=zipper(res,array[i])
    return res
    
def find_nearest(array,value):
    '''
    Find the index of nearest element to given value in 1d array,
    considering the continuous condition.
    cost[i] = abs(array[i]-value)+abs(array[i+1]+array[i-1]-2*array[i])
    '''
    if(array is None or len(array) is 0):return None
    if(len(array) is 1):return array[0]
    curr=abs(array[0]-value)+abs(array[1]-array[0])
    ind=0
    for i in range(len(array)):
        #print(i,len(array))
        if(i < len(array)-1):
            diff=abs(array[i]-value)+abs(array[i+1]+array[i-1]-2*array[i])
        else:diff=abs(array[i]-value)+abs(array[i-1]-array[i])
        if(diff<curr):
            curr=diff
            ind=i
    
    return ind

    
def groupbydiff(array,index,threshold):
    '''
    Continuously group an n dimension array by array[index], 
    using threshold of adjacent difference.
    returns a list of arrays.
    '''
    i,j=0,0
    #print("grouper:"+str(array.shape)+str(len(array[index])))
    res=[np.array([])]
    arrayp=array.T
    #print(len(array[index]))
    for i in range(len(array[index]))[0:-1]:
        #print(i)
        res[j]=insert(res[j],arrayp[i],-1)
        #insert a row
        if(abs(array[index][i]-array[index][i+1])>threshold):
            j+=1
            res+=[np.array([])]
    res[j]=insert(res[j],arrayp[-1],-1)
    return [e.T for e in res]
    
def groupbydir(array,index):
    '''
    Continuously group an n dimension array by array[index], 
    using the direction of elements.
    returns a list of arrays.
    '''
    i,j=0,0
    #print("grouper:"+str(array.shape)+str(len(array[index])))
    res=[np.array([])]
    arrayp=array.T
    #print(len(array[index]))
    a=array[index][1]-array[index][0]
    for i in range(len(array[index]))[1:-1]:
        #print(i)
        b=array[index][i]-array[index][i-1]
        if(a*b>0): #same direction
            res[j]=insert(res[j],arrayp[i],-1)
        elif(a==0 and b!=0 or a!=0 and b==0 or a*b<0):
        #insert a row
            j+=1
            res+=[np.array([])]
        a=b
        
    res[j]=insert(res[j],arrayp[-1],-1)
    return [e.T for e in res]
#=====================file I/O & string manipulation API=======================
    
#STRONGLY RECOMMAND: np.loadtxt()&np.savetxt()
    # np.loadtxt(filename) gives a horizontal array.
    # Use .T  attribute to get vertical array.
    # np.savetxt(filename,data,delimiter='\t',newline=os.linesep)
    # where data is saved horizontally.
def read_v(fname):
    '''
    using read_v you will get a column-major ordered array.
    eg. res[0] is the first column of data.
    '''
    return np.loadtxt(fname).T
def read_h(fname):
    '''
    using read_h you will get a row-major ordered array.
    eg. res[0] is the first row of data.
    '''
    return np.loadtxt(fname)
def write_v(fname,array):
    '''
    for vertically saved data
    '''
    return np.savetxt(fname,array.T,delimiter='\t',newline=os.linesep)
def write_h(fname,array):
    '''
    for horizontally saved data
    '''
    return np.savetxt(fname,array,delimiter='\t',newline=os.linesep)
def rm_str_in_name(dataDir,strList):
    '''
    For all files in data_dir, 
    remove certain string from filename,
    to make the files neat.
    Fhe strings are in strList.
    '''
    if(dataDir==None or dataDir==""):return
    if(type(strList) is not list):strList=[strList]
    for fName in os.listdir(dataDir):
        newName=fName
        for s in strList:newName=newName.replace(s,'')
        print('renaming as:'+newName)
        os.rename(dataDir+fName,dataDir+newName)
    return dataDir
def sscan_raw(re_string,string,prefix='',postfix=''):
    '''
    sort of a sscanf to scan the first match of re_string in  string
    '''
    searcher=prefix+re_string+postfix
    p=re.compile(searcher)
    m=p.search(string)
    if(m is None):return None
    res=m.group(0)
    if(len(res) is 0):return None
    
    start=len(prefix)
    end=len(res)-len(postfix)
    return res[start:end]
def sscan_double(string,prefix='',postfix=''):
    '''
    sort of a sscanf("%f") or sscanf("%lf")
    '''
    return sscan(float,string,prefix,postfix)
def sscan_int(string,prefix='',postfix=''):
    '''
    sort of a sscanf("%d")
    '''
    return sscan(int,string,prefix,postfix)
def sscan(v_type,string,prefix='',postfix=''):
    '''
    this simulating sscanf
    v_type: scan type, like float, int, hex, oct, str.
    string: the string to scan.
    pre/post fix: string, default to be ''.
    
    return type depends on v_type.
    '''
    if(v_type is float):re_string='[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
    elif(v_type is int):re_string='[-+]?\d+'
    elif(v_type is hex):re_string='[-+]?(0[xX])?[\dA-Fa-f]+'
    elif(v_type is oct):re_string='[-+]?[0-7]+'
    elif(v_type is str):re_string='\S+'
    else:
        print('[sscan]Unrecognized type')
        return None
    searcher=prefix+re_string+postfix
    p=re.compile(searcher)
    m=p.search(string)
    if(m is None):return None
    res=m.group(0)
    if(len(res) is 0):return None
    
    start=len(prefix)
    end=len(res)-len(postfix)
    if(v_type is float):return float(res[start:end])
    elif(v_type is int):return int(res[start:end])
    elif(v_type is hex):return hex(int(res[start:end]))
    elif(v_type is oct):return oct(int(res[start:end]))
    elif(v_type is str):return str(res[start:end])
    else:return None
#=====================================signal process API=======================
def runningMean(x, N):
    '''
    x is the array to be processed
    N is the size of meanning window
    '''
    # from stack overflow
    # using the fact that the running mean of an array is actually the convolution of the array and [1/N]*N
    return np.convolve(x, np.ones((N,))/N)[(N-1):]
    

    
#=====================================data process API=========================
# all API in this section use data_dir as first variable.
    
def medfilt_data(data_dir,index,kernel_size,fname_filter=''):
    if(data_dir==None or data_dir==""):return
    if(type(index)==int):index=[index]
    if(type(kernel_size)==int):kernel_size=[kernel_size]*len(index)
    
    new_dir=data_dir+"-mfilt"
    if not os.path.exists(new_dir):os.makedirs(new_dir)
        
    for f in os.listdir(data_dir):
        fOld=os.path.join(data_dir,f)
        fNew=os.path.join(new_dir,f)
        if(os.path.isdir(fNew)):continue
        if(f.find(fname_filter)<0):
            shutil.copy(fOld,fNew)
            continue
        data=read_v(fOld)
        print("Median filting "+ f +":"+str(data.T.shape))
        #vertical?
        for i in range(len(index)):
            if(index[i]>=0 and index[i]<len(data)):
                #print(index[i])
                data[index[i]]=runningMean(data[index[i]],kernel_size[i])
        fNew=fNew[:-4]+"-medfilted.txt"
        write_v(fNew,data)
    return new_dir

def MR(data_dir,H_col=0,Vxx_col=1,Vxy_col=2,to_calc='%',
       bench=0.0,fname_filter='',curr_in_nA=20,carrier_factor=1):
    '''
    compute the magnetoresistance
    'd' for delta resistance
    '%' for percent resistance
    'r' for just ratio
    and exactly as the order in to_calc

    new feature:
    calculate Hall effiecients and carrier density (in cm^(-2))
    H in Oe
    carrier factor can be either 1 or -1,
    meaning electron or hole.
    '''
    if(data_dir==None or data_dir==""):return
    #now we can only average on B(0) or Vg(1)
    z=0 #zero field index
    if(type(Vxx_col)!=list):Vxx_col=[Vxx_col]
    new_dir=data_dir+"-MR"
    if not os.path.exists(new_dir):os.makedirs(new_dir)
    to_calc=to_calc.lower()
    if(len(to_calc)<1):return(data_dir)
    
    for f in os.listdir(data_dir):
        fOld=os.path.join(data_dir,f)
        fNew=os.path.join(new_dir,f)
        if(os.path.isdir(fNew)):continue
        if(f.find(fname_filter)<0):
            shutil.copy(fOld,fNew)
            continue
        data=read_v(fOld)
        print("Finding MR of "+ f +":"+str(data.T.shape))

        z=find_nearest(data[H_col],bench)
        for i in range(len(Vxx_col)):
            if(Vxx_col[i]<0):Vxx_col[i]+=len(data)
        for r in Vxx_col:
            for operator in to_calc:
                if(operator is 'd'):data=zipper(data,data[r]-data[r][z])
                elif(operator is '%'):data=zipper(data,100*(data[r]-data[r][z])/data[r][z])
                elif(operator is 'r'):data=zipper(data,data[r]/data[r][z])
                elif(operator is not 'n'):
                    print("Wield variable in MR program:"+str(to_calc))
        
        Hall_fit_res=np.polyfit(data[H_col],data[Vxy_col],1)
        slope=Hall_fit_res[0]
        carrier_density=carrier_factor*curr_in_nA/slope/1.6*100
        print('k=',slope,', n=',carrier_density)
        if('n' in to_calc): fNew=os.path.join(new_dir,'n='+str(carrier_density)+'_'+f)
        write_v(fNew,data)     

    return new_dir
    
def reverse_data_col(data_dir,col=0,fname_filter=''):
    # reverse the order of some column
    if(data_dir==None or data_dir==""):return
    if(type(col)!=list):col=[col]
    new_dir=data_dir+"-rev"
    if not os.path.exists(new_dir):os.makedirs(new_dir)
    
    for f in os.listdir(data_dir):
        fOld=os.path.join(data_dir,f)
        fNew=os.path.join(new_dir,f)
        if(os.path.isdir(fOld)):continue
        if(f.find(fname_filter)<0):
            shutil.copy(fOld,fNew)
            continue
        data=read_v(fOld)
        print("Reversing " + f + ":" + str(data.T.shape))

        for r in col:
            data=zipper(data,data[r][::-1])
        write_v(fNew,data)
    return new_dir
    
def symm_data_col(data_dir,col=0,fname_filter='',operator='+',factor=1):
    # symmetrize some columns of data
    if(data_dir==None or data_dir==""):return
    if(type(col) is not list):col=[col]
    if(type(operator) is not list):operator=[operator]*len(col)
    new_dir=data_dir+"-symm"
    if not os.path.exists(new_dir):os.makedirs(new_dir)
    
    for f in os.listdir(data_dir):
        fOld=os.path.join(data_dir,f)
        fNew=os.path.join(new_dir,f)
        if(os.path.isdir(fOld)):continue
        if(f.find(fname_filter)<0):
            shutil.copy(fOld,fNew)
            continue
        data=read_v(fOld)
        print("Symmetrizing " + f + ":" + str(data.T.shape))

        for i in range(len(col)):
            c=col[i]
            if(operator[i] is '+'):data=zipper(data,factor*(data[c][::-1]+data[c])/2)
            elif(operator[i] is '-'):data=zipper(data,factor*(data[c]-data[c][::-1])/2)
        write_v(fNew,data)
    return new_dir
    
def sigma(data_dir,Vxx_col=1,Vxy_col=3,to_calc='xy',
            aspect_ratio=1,curr_in_nA=20,fname_filter='',n=1):
    # calculate sigmaxy in units of n*e2/h
    if(data_dir==None or data_dir==""):return
    new_dir=data_dir+"-sigma"
    if not os.path.exists(new_dir):os.makedirs(new_dir)
    to_calc=to_calc.lower()
    
    for f in os.listdir(data_dir):
        fOld=os.path.join(data_dir,f)
        fNew=os.path.join(new_dir,f)
        if(os.path.isdir(fOld)):continue
        if(f.find(fname_filter)<0):
            shutil.copy(fOld,fNew)
            continue
        
        data=read_v(fOld)
        print("Calculating sigma" + to_calc + " of " + f + ":" + str(data.T.shape))
        
        C=curr_in_nA*10**(-9)*25812.8/n
        A=aspect_ratio
        if('x' in to_calc):
            data=zipper(data,C*data[Vxx_col]/((data[Vxx_col]/A)**2+data[Vxy_col]**2))
        if('y' in to_calc):
            data=zipper(data,-C*data[Vxy_col]/((data[Vxx_col]/A)**2+data[Vxy_col]**2))
        write_v(fNew,data)
    return new_dir
    
def split_data(data_dir,key_item=0,prefix="",name_col=1,zero_col=0,threshold=0.1):
    if(data_dir==None or data_dir=="" or key_item<0):return data_dir
    new_dir=data_dir+"-split"
    if not os.path.exists(new_dir):os.makedirs(new_dir)
    for f in os.listdir(data_dir):
        fOld=os.path.join(data_dir,f)
        if(os.path.isdir(fOld)):continue
    
        data=read_v(fOld)
        print("Spliting "+ f +":"+str(data.shape))
        #for array in [arrays[i:i+41] for i in range(0,len(arrays),41)]:
        for block in groupbydiff(data,key_item,threshold):
            #print(data,len(block))
            if(len(block) is 0):continue
            #print(str(array[key_item][0]))
            z=find_nearest(block[zero_col],0)
            fNew=os.path.join(new_dir,prefix+str(block[name_col][z])+"-"+f)
            write_v(fNew,block)
    return new_dir

def select_data(data_dir,key_item,new_dir,fname_filter='',prefix=''):
    if(data_dir is None or len(data_dir) is 0):return data_dir
    if not os.path.exists(new_dir):os.makedirs(new_dir)
    if(type(key_item) is not list):key_item=[key_item]
    
    for f in os.listdir(data_dir):
        if(f.find(fname_filter)<0):continue
        fOld=os.path.join(data_dir,f)
        fNew=os.path.join(new_dir,prefix+f)
        data=read_v(fOld)
        print("Selecting "+ f +":"+str(data.shape)+str(key_item))
        new=np.array([])
        for i in key_item:new=insert(new,data[i],-1)
        write_v(fNew,new)
    return new_dir
    
def merge_data(data_dirs,data_name):
    #data_dir can be s single string or a list of strings.
    #data_name is the final name of data file and data directory.
    if(data_dirs==None or data_dirs=="" or len(data_dirs) is 0):return
    if not os.path.exists(data_name):
        os.makedirs(data_name)
    all_data=np.array([])
    if(type(data_dirs) is not list):
        data_dirs=[data_dirs]
    for data_dir in data_dirs:
        for data in os.listdir(data_dir):
            array=read_h(data_dir+"//"+data)
            print("merging data:",data,":",array.shape)
            old_shape=all_data.shape
            if(len(all_data)==0):all_data=array
            else:all_data=np.concatenate((all_data,array))
            new_shape=all_data.shape
            print(old_shape,"---->",new_shape)
    
    write_h(data_name+"//"+data_name+".txt",all_data)
    return data_name
def sort_data_by_row(data_dir,key_item=0):
    if(data_dir==None or data_dir=="" or key_item<0):return
    if not os.path.exists(data_dir+"-sorth"):
        os.makedirs(data_dir+"-sorth")
    for data in os.listdir(data_dir):
        if(os.path.isdir(data_dir+"\\"+data)):
            #print(data)
            continue
        array=read_h(data_dir+"\\"+data)
        print("sorting "+data+":"+str(array.shape))
        order=np.lexsort([array[key_item]],axis=0)
        array=(array.T[order]).T            
        if(len(array)==0):continue
        write_h(data_dir+"-sorth\\"+data,array)
    return data_dir+"-sorth"
def sort_data_by_col(data_dir,key_item=0):
    if(data_dir==None or data_dir=="" or key_item<0):return
    if not os.path.exists(data_dir+"-sortv"):
        os.makedirs(data_dir+"-sortv")
    for data in os.listdir(data_dir):
        if(os.path.isdir(data_dir+"\\"+data)):
            #print(data)
            continue
        array=read_v(data_dir+"\\"+data)
        print("sorting "+data+":"+str(array.shape))
        #for array in [arrays[i:i+41] for i in range(0,len(arrays),41)]:
        order=np.lexsort([array[key_item]],axis=0)
        array=(array.T[order]).T            
        #print(array)
        try:
            if(len(array)==0):continue
        except:
            continue
        write_v(data_dir+"-sortv\\"+data,array)
    return data_dir+"-sortv"
def peek_data(data_dir,N=10):
    if(data_dir==None or data_dir==""):return
    array=read_h(data_dir+"\\"+os.listdir(data_dir)[0])
    print("==========",data_dir,"==========")
    print(array[:N+1])
    print("========",data_dir,"end=========")
    return
def multiply_data(data_dir,key_item=0,factor=1,fname_filter=''):
    '''
    data_dir have to be a path string, not a list of strings
    '''
    if(data_dir is None or len(data_dir) is 0):return
    new_dir=data_dir+"-mul"
    if not os.path.exists(new_dir):os.makedirs(new_dir)
    

    for f in os.listdir(data_dir):
        fOld=os.path.join(data_dir,f)
        fNew=os.path.join(new_dir,f)
        if(os.path.isdir(fOld)):continue
        if(f.find(fname_filter)<0):
            shutil.copy(fOld,fNew)
            continue
        array=read_v(fOld)
        print("multiplying data:",f,":",array.shape)
        array[key_item]=array[key_item]*factor
        print("done multiplying:",f,":",array.shape)
        write_v(fNew,array)
    
    return new_dir
def magnet_len(H):
    #H in Tesla
    #return value in nm
    return 12.83744/np.sqrt(H)
def magnet_tau(H,D):
    #H in Tesla
    #return value in s
    return 1.65*10**(-16)/(H*D)


def calcMR(i):
    data=["E:\\Data\\20170208_BLG\\#"+str(i)]
    data+=[symm_data_col(data[-1],col=1,fname_filter='sweepB',operator='+')]
    data+=[symm_data_col(data[-1],col=3,fname_filter='sweepB',operator='-')]
    data+=[sigma(data[-1],Vxx_col=-2,Vxy_col=-1,aspect_ratio=4)]
    data+=[MR(data[-1],H_col=0,Vxx_col=5,Vxy_col=6,to_calc="%",fname_filter='sweepB')]
    for d in data[1:-1]:shutil.rmtree(d)
    
def collectMR(offset,maxround):
    vd=[-13,-20,-35,-38,-38,-40,-49,-50.0,-50.5,-53,-70,-72,-83,-102]
    new_dir="E:\\Data\\20170208_BLG\\nH_Vg="+str(offset)+'\\'
    for i in range(maxround+1):
        data="E:\\Data\\20170208_BLG\\#"+str(i)+'-symm-symm-sigma-MR'
        fil='Vg='+str(vd[i]+offset)
        prefix=str(i)+'H'
        #print(data,fil,new_dir)
        select_data(data,[0,-1],new_dir,fil,prefix)
    rm_str_in_name(new_dir,["_20nA","_Vsd=0.200000","9'.D.10'.11'.10.10'"])


#BLG in-situ hydrogenation
maxround=13
for nH in range(maxround+1):calcMR(nH)
for vg in [0,20,40,60]:collectMR(vg,maxround)

'''
## 0==T  1==H  2==Vg  3==Ig  4==Vxx  6==Vxy
#TLG data processing
data=[r"E:\Data_Cloud\TLG\all\sweepH"]
data+=[split_data(data[-1],key_item=2,prefix="Vg=",name_col=2,zero_col=1,threshold=1.0)]
data+=[MR(data[-1],H_col=1,Vxx_col=4,Vxy_col=6,to_calc='%',bench=0.0,fname_filter='50nA',curr_in_nA=50,carrier_factor=1)]
data+=[MR(data[-1],H_col=1,Vxx_col=4,Vxy_col=6,to_calc='%',bench=0.0,fname_filter='20nA',curr_in_nA=20,carrier_factor=1)]
data+=[symm_data_col(data[-1],col=-1,operator='+')]

#TLG data processing
data=[r"E:\Data_Cloud\TLG\all\sweepVg"]
data+=[split_data(data[-1],key_item=1,prefix="H=",name_col=1,zero_col=2,threshold=500)]
data+=[multiply_data(data[-1],key_item=4,factor=1/(5*10**(-8)),fname_filter='50nA')]
data+=[multiply_data(data[-1],key_item=4,factor=1/(2*10**(-8)),fname_filter='20nA')]
'''