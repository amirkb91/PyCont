#!/bin/bash

rm -f /home/akb110/.local/lib/python3.8/site-packages/PyCont-lib.egg-link
rm -f /home/akb110/.local/lib/python3.8/site-packages/easy-install.pth
rm -rf PyCont_lib.egg-info
pip install --user -e .
