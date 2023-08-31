import sys
import shutil

file = sys.argv[-1]
dest_folder = "/mnt/c/Users/akb110/Desktop/HDFdump/"
shutil.copy2(file, dest_folder)