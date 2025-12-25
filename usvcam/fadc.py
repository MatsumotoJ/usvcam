'''
MIT License

Copyright (c) 2024 FRIENTECH INC.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
# ===========================================================================
# DLL定義
# ===========================================================================
# https://docs.python.org/ja/3.11/library/ctypes.html#fundamental-data-types
# https://docs.python.org/ja/3/library/ctypes.html#ctypes-fundamental-data-types-2
import sys
import ctypes
import os
from enum import IntEnum

script_dir = os.path.dirname(__file__)
os.add_dll_directory(script_dir)
#lib=ctypes.WinDLL('fadcd3xx.dll')
script_dir = os.path.dirname(__file__)
lib=ctypes.WinDLL(script_dir + '/fadcd3xx.dll')

class FT_STATUS(IntEnum):
    FT_OK = 0
    FT_INVALID_HANDLE = 1
    FT_DEVICE_NOT_FOUND = 2
    FT_DEVICE_NOT_OPENED = 3
    FT_IO_ERROR = 4
    FT_INSUFFICIENT_RESOURCES = 5
    FT_INVALID_PARAMETER = 6
    FT_INVALID_BAUD_RATE = 7
    FT_DEVICE_NOT_OPENED_FOR_ERASE = 8
    FT_DEVICE_NOT_OPENED_FOR_WRITE = 9
    FT_FAILED_TO_WRITE_DEVICE = 10
    FT_EEPROM_READ_FAILED = 11
    FT_EEPROM_WRITE_FAILED = 12
    FT_EEPROM_ERASE_FAILED = 13
    FT_EEPROM_NOT_PRESENT = 14
    FT_EEPROM_NOT_PROGRAMMED = 15
    FT_INVALID_ARGS = 16
    FT_NOT_SUPPORTED = 17
    FT_NO_MORE_ITEMS = 18
    FT_TIMEOUT = 19
    FT_OPERATION_ABORTED = 20
    FT_RESERVED_PIPE = 21
    FT_INVALID_CONTROL_REQUEST_DIRECTION = 22
    FT_INVALID_CONTROL_REQUEST_TYPE = 23
    FT_IO_PENDING = 24
    FT_IO_INCOMPLETE = 25
    FT_HANDLE_EOF = 26
    FT_BUSY = 27
    FT_NO_SYSTEM_RESOURCES = 28
    FT_DEVICE_LIST_NOT_READY = 29
    FT_DEVICE_NOT_CONNECTED = 30
    FT_INCORRECT_DEVICE_PATH = 31
    FT_OTHER_ERROR = 32

#------------------------------------------------
# fadcGetVersion
#------------------------------------------------
fadcGetVersion=lib.fadcGetVersion # 関数名
fadcGetVersion.restype=ctypes.c_uint32 #戻り値の型

#------------------------------------------------
# fadcOpenByIndex
#------------------------------------------------
fadcOpenByIndex=lib.fadcOpenByIndex # 関数名
fadcOpenByIndex.argtype=[ctypes.POINTER(ctypes.c_void_p),ctypes.c_uint32] #引数の型
fadcOpenByIndex.restype=ctypes.c_uint32 #戻り値の型

#------------------------------------------------
# fadcClose
#------------------------------------------------
fadcClose=lib.fadcClose # 関数名
fadcOpenByIndex.argtype=[ctypes.c_void_p] #引数の型
fadcClose.restype=ctypes.c_uint32 #戻り値の型

#------------------------------------------------
# fadcGetFpgaVersion
#------------------------------------------------
fadcGetFpgaVersion=lib.fadcGetFpgaVersion # 関数名
fadcGetFpgaVersion.argtype=[ctypes.c_void_p,ctypes.c_void_p] #引数の型
fadcGetFpgaVersion.restype=ctypes.c_uint32 #戻り値の型

#------------------------------------------------
# fadcStart
#------------------------------------------------
fadcStart=lib.fadcStart # 関数名
fadcStart.argtype=[ctypes.c_void_p,ctypes.c_uint32,ctypes.c_uint32] #引数の型
fadcStart.restype=ctypes.c_uint32 #戻り値の型

#------------------------------------------------
# fadcStop
#------------------------------------------------
fadcStop=lib.fadcStop # 関数名
fadcStop.argtype=[ctypes.c_void_p] #引数の型
fadcStop.restype=ctypes.c_uint32 #戻り値の型

#------------------------------------------------
# fadcRead
#------------------------------------------------
fadcRead=lib.fadcRead #関数
fadcRead.argtype=[ctypes.c_void_p,ctypes.POINTER(ctypes.c_ubyte),ctypes.c_uint32,ctypes.POINTER(ctypes.c_uint32)] #引数の型
fadcRead.restype=ctypes.c_int32 #戻り値の型
