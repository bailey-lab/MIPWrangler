#!/usr/bin/env python3

import shutil, os, argparse, sys, stat
sys.path.append("scripts/pyUtils")
sys.path.append("scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper
def main():
    name = "MIPWrangler"
    libs = "seekdeep:v3.0.0"
    args = genHelper.parseNjhConfigureArgs()
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv, private = True)
    Utils.run(cmd)
    
main()

