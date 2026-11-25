#
# FadcDeviceクラス : fadcd3xx.dllをラッピングしたクラス
#
# Windows Python 3.12.2 AMD64
# Windows11 Pro 64bit
#
# フレンテック株式会社
#
# ■履歴
# 2024-06-03 1st release
#
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

import ctypes
from enum import IntEnum
import os
import time

script_dir = os.path.dirname(__file__)

os.add_dll_directory(script_dir)
#lib = ctypes.WinDLL('fadcd3xx.dll')
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

FPGA_9CH_SUPPORT_VERSION = 6 # FPAG version6から9CHサポート

# ctypes 関数定義
fadcGetVersion = lib.fadcGetVersion
fadcGetVersion.restype = ctypes.c_uint32

fadcOpenByIndex = lib.fadcOpenByIndex
fadcOpenByIndex.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_uint32]
fadcOpenByIndex.restype = ctypes.c_uint32

fadcClose = lib.fadcClose
fadcClose.argtypes = [ctypes.c_void_p]
fadcClose.restype = ctypes.c_uint32

fadcGetFpgaVersion = lib.fadcGetFpgaVersion
fadcGetFpgaVersion.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint16)]
fadcGetFpgaVersion.restype = ctypes.c_uint32

fadcSetChannelCount = lib.fadcSetChannelCount
fadcSetChannelCount.argtypes = [ctypes.c_void_p, ctypes.c_int32]
fadcSetChannelCount.restype = ctypes.c_uint32

fadcGetChannelCount = lib.fadcGetChannelCount
fadcGetChannelCount.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int32)]
fadcGetChannelCount.restype = ctypes.c_uint32

fadcStart = lib.fadcStart
fadcStart.argtypes = [ctypes.c_void_p, ctypes.c_uint32, ctypes.c_uint32]
fadcStart.restype = ctypes.c_uint32

fadcStop = lib.fadcStop
fadcStop.argtypes = [ctypes.c_void_p]
fadcStop.restype = ctypes.c_uint32

fadcRead = lib.fadcRead
fadcRead.argtypes = [
    ctypes.c_void_p,
    ctypes.POINTER(ctypes.c_ubyte),
    ctypes.c_uint32,
    ctypes.POINTER(ctypes.c_uint32)
]
fadcRead.restype = ctypes.c_int32

WAIT_AFTER_READ = 0.01  # 連続リード時読み出した値が不定になる対策

# ---------------------------------------------------------
#  ラッパークラス
# ---------------------------------------------------------
class FadcDevice:
    def __init__(self, index=0):
        self.handle = ctypes.c_void_p()
        ret = fadcOpenByIndex(ctypes.byref(self.handle), index)
        self._check(ret, "fadcOpenByIndex")
        self.fpga_version = self.get_fpga_version()

    def close(self):
        if self.handle:
            fadcStop(self.handle)
            fadcClose(self.handle)
            self.handle = None

    def __del__(self):
        self.close()

    # 共通エラーチェック
    def _check(self, ret, func_name):
        if ret != FT_STATUS.FT_OK:
            raise RuntimeError(f"{func_name} failed: {FT_STATUS(ret).name}")

    # -------------------------
    # API ラッパー
    # -------------------------
    def get_dll_version(self):
        ver = fadcGetVersion()
        return ver
    
    def get_fpga_version(self):
        ver = ctypes.c_uint16()
        ret = fadcGetFpgaVersion(self.handle, ctypes.byref(ver))
        self._check(ret, "fadcGetFpgaVersion")
        time.sleep(WAIT_AFTER_READ) # 連続リード時読み出した値が不定になる対策
        return ver.value

    def is_channel_count_supported(self):
        if self.fpga_version >= FPGA_9CH_SUPPORT_VERSION: # チャンネル数指定APIを使用可能か？
          return True
        else:
          return False

    def set_channel_count(self, count: int):
        if self.fpga_version >= FPGA_9CH_SUPPORT_VERSION: # 使用可能な場合のみ処理続行
          ret = fadcSetChannelCount(self.handle, int(count))
          self._check(ret, "fadcSetChannelCount")

    def get_channel_count(self):
        if self.fpga_version >= FPGA_9CH_SUPPORT_VERSION: # 使用可能な場合のみ処理続行
          count = ctypes.c_int32()
          ret = fadcGetChannelCount(self.handle, ctypes.byref(count))
          time.sleep(WAIT_AFTER_READ) # 連続リード時読み出した値が不定になる対策
          self._check(ret, "fadcGetChannelCount")
          return count.value
        else:
          return 0

    def start(self, sample_rate: int, range: int = 0):
        ret = fadcStart(self.handle, sample_rate, range)
        self._check(ret, "fadcStart")

    def stop(self):
        ret = fadcStop(self.handle)
        self._check(ret, "fadcStop")

    def read(self, size: int) -> bytes:
        buf = (ctypes.c_ubyte * size)()
        read_size = ctypes.c_uint32()

        ret = fadcRead(self.handle, buf, size, ctypes.byref(read_size))
        if ret < 0:
            raise RuntimeError(f"fadcRead failed: {ret}")

        return bytes(buf[:read_size.value])


# ---------------------------------------------------------
#  使い方
# ---------------------------------------------------------
if __name__ == "__main__":
    dev = FadcDevice(0)

    print("FPGA Version:", dev.get_fpga_version())

    dev.set_channel_count(8)
    print("Channel Count:", dev.get_channel_count())

    dev.start(0)
    data = dev.read(4096)
    print("Read:", len(data), "bytes", type(data))

    dev.stop()
    dev.close()
