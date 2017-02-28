# -*- coding: utf-8 -*-
"""
Data collection procedures, and instruments driver.
suits for gpib instruments.

TODO:support for visa and TCP/IP instruments(e.g., PPMS and MPMS)

v0.1 @Feb 22 2017, by shmcao.

"""
"""
Created on Tue Feb 21 13:19:02 2017

@author: admin
"""
import ctypes
import os
class GPIBInstr:
    def __init__(self,address,libpath=""):
        if(len(libpath) is 0):
            _dir=os.path.split(os.path.realpath(__file__))[0]
            self.lib_p=os.path.join(_dir,"PyInstr.dll")
            #print(self.lib_p,type(self.lib_p))
        else:self.lib_p=libpath
        self.addr=address
        self.lib=ctypes.cdll.LoadLibrary( self.lib_p )
        self._read=self.lib.GPIBRead
        self._write=self.lib.GPIBWrite
        self._pread=self.lib.PseudoRead
        self._pwrite=self.lib.PseudoWrite
        self._read.argtypes = [ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int]
        self._write.argtypes = [ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int]
        self._pread.argtypes = [ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int]
        self._pwrite.argtypes = [ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int]
        self._read.restypes = ctypes.c_void_p
        self._write.restypes = ctypes.c_void_p
        self._pread.restypes = ctypes.c_void_p
        self._pwrite.restypes = ctypes.c_void_p
        self.inner_message=""
        return None
    def read(self,maxLen=3200):
        a=ctypes.create_string_buffer(str(self.addr).encode())
        out=ctypes.create_string_buffer(maxLen)
        self._read(a,out,maxLen)
        out=str(out.value.decode()).rstrip(str(b'\x00'))
        return out
    def write(self,message):
        a=ctypes.create_string_buffer(str(self.addr).encode())
        msg=ctypes.create_string_buffer(str(message).encode())
        self._write(a,msg,len(msg))
        out=str(msg.value.decode()).rstrip(str(b'\x00'))
        return out
    def query(self,message):
        self.write(message)
        return self.read()
    def pread(self,maxLen=3200):
        a=ctypes.create_string_buffer(str(self.addr).encode())
        out=ctypes.create_string_buffer(maxLen)
        self._pread(a,out,maxLen)
        # This is in fact the address of the instrument
        addr=str(out.value.decode()).rstrip(str(b'\x00'))
        out="address="+addr+", message="+self.inner_message+"(end)"
        return out
    def pwrite(self,message):
        a=ctypes.create_string_buffer(str(self.addr).encode())
        msg=ctypes.create_string_buffer(str(message).encode())
        self._pwrite(a,msg,len(msg))
        out=str(msg.value.decode()).rstrip(str(b'\x00'))
        self.inner_message=out
        return out
    
test=GPIBInstr(1)
print(test.pread())
print(test.pwrite("hello"))
print(test.pread())
print(test.pwrite("hello again"))
print(test.pread())
print(test.pwrite("hi there"))
print(test.pread())

k=GPIBInstr(17)
print(k.write("*IDN?"))
print(k.read())
print(k.read())
print(k.read())