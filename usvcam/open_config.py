import os
import subprocess
import platform
import sys

def open_config(config_path):

    if platform.system() == 'Darwin':       # macOS
        subprocess.call(('open', config_path))
    elif platform.system() == 'Windows':    # Windows
        os.startfile(config_path)
    else:                                   # linux variants
        subprocess.call(('xdg-open', config_path))

def main():
    
    script_dir = os.path.dirname(__file__)
    config_path = script_dir + '/config_fadc.yaml' 

    open_config(config_path)

def main_legacy():
    
    script_dir = os.path.dirname(__file__)
    config_path = script_dir + '/config_legacy.yaml' 

    open_config(config_path)


if __name__ == '__main__':
    main()
    