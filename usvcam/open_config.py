import os
import subprocess
import platform
import sys

def main():
    
    config_path = sys.prefix + '/etc/usvcam/config.yaml'

    if platform.system() == 'Darwin':       # macOS
        subprocess.call(('open', config_path))
    elif platform.system() == 'Windows':    # Windows
        os.startfile(config_path)
    else:                                   # linux variants
        subprocess.call(('xdg-open', config_path))

if __name__ == '__main__':
    main()
    