#!/usr/bin/env python
import shutil, os, argparse, sys, stat

def parse_args_genServerProxyPassConfig():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--start',   type=int, default = 0)
    parser.add_argument('--stop',    type=int, default = 9)
    parser.add_argument('--outfile', type=str, required = True)
    parser.add_argument('--overWrite', action = "store_true")
    
    return parser.parse_args()

def genServerProxyPassConfig():
    args = parse_args_genServerProxyPassConfig()
    
    if args.stop < args.start:
        raise Exception("Stop should be greater than stop, start: " + str(args.start) + ", stop: " + str(args.stop))
    
    if not str(args.outfile).endswith(".conf"):
        args.outfile = args.outfile + ".conf"
    
    if os.path.exists(args.outfile) and not args.overWrite:
        raise Exception("Error, " + args.outfile + " already exits, use --overWrite to overWrite it")
    
    with open(args.outfile, "w") as outfile:
        for serverNum in xrange(args.start, args.stop  + 1):
            outfile.write("ProxyPass          /mip" + str(serverNum) + "  http://127.0.0.1:" + str(10000 + serverNum) + "/mip" + str(serverNum) + "\n")
            outfile.write("ProxyPassReverse   /mip" + str(serverNum) + "  http://127.0.0.1:" + str(10000 + serverNum) + "/mip" + str(serverNum) + "\n")
    
    
    
    
if __name__ == "__main__":
    genServerProxyPassConfig()




