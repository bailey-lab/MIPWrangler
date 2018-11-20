#!/usr/bin/env python
import shutil, os, argparse, sys, stat


class MipServerConfCreator:
    def __init__(self, outfileFnp):
        self.outfileFnp_ = outfileFnp
        if not self.outfileFnp_.endswith(".conf"):
            self.outfileFnp_ = self.outfileFnp_ + ".conf"        
        self.portsTaken_ = []
        self.namesTaken_ = []
        if os.path.exists(self.outfileFnp_):
            with open(self.outfileFnp_, 'r') as serverConfFile:
                for line in serverConfFile:
                    toks = line.strip().split()
                    if len(toks) == 3 and "ProxyPass" == toks[0]:
                        self.namesTaken_.append(toks[1][1:]);
                        colonPos = toks[2].find("127.0.0.1:")
                        slashPos = toks[2].find("/", colonPos)
                        if -1 != colonPos and -1 != slashPos:
                            self.portsTaken_.append(int(toks[2][(colonPos + len("127.0.0.1:")):slashPos]))
                        
    def addServer(self, portName, portNumber):
        maxPort = 65,535
        minPort = 1001
        
        if type(portName) is not str:
            raise Exception("Error in MipServerConfCreator, portName should be str, not " + type(portName))
        if type(portNumber) is not int:
            raise Exception("Error in MipServerConfCreator, portNumber should be int, not " + type(portNumber))
        if portNumber < minPort or portNumber > maxPort:
            raise Exception("Error in MipServerConfCreator, portNumber should be between 1000 and 65,536, not " + str(portNumber))
        if portName in self.namesTaken_:
            raise Exception("Error in MipServerConfCreator, portName, " + str(portName) + ", already taken")
        if portNumber in self.portsTaken_:
            raise Exception("Error in MipServerConfCreator, portNumber, " + str(portNumber) + ", already taken")
        if "/" in portName:
             raise Exception("Error in MipServerConfCreator, / not allowed in portName, " + str(portName))
        
        with open(self.outfileFnp_, "a") as outfile:
            outfile.write("ProxyPass          /" + str(portName) + "  http://127.0.0.1:" + str(portNumber) + "/" + str(portName) + "\n")
            outfile.write("ProxyPassReverse   /" + str(portName) + "  http://127.0.0.1:" + str(portNumber) + "/" + str(portName) + "\n")
            self.namesTaken_.append(portName)
            self.portsTaken_.append(portNumber)
    
    def addMipServers(self,serverBase, start, stop):
        if stop < start:
            raise Exception("Error in MipServerConfCreator::addMipServers, stop should be greater than stop, start: " + str(start) + ", stop: " + str(stop))
        for serverNumber in xrange(start, stop):
            self.addServer("mip" + str(serverNumber), serverBase + serverNumber);

    
    def removePort(self, portNumber):
        tempFnp = "/tmp/temp_" + os.path.basename(self.outfileFnp_)
        if os.path.exists(self.outfileFnp_):
            foundPort = False;
            with open(self.outfileFnp_, 'r') as serverConfFile:
                with open(tempFnp, "w") as tempFile:
                    for line in serverConfFile:
                        toks = line.strip().split()
                        if len(toks) == 3:
                            name = toks[1][1:];
                            colonPos = toks[2].find("127.0.0.1:")
                            slashPos = toks[2].find("/", colonPos)
                            if -1 != colonPos and -1 != slashPos:
                                port = int(toks[2][(colonPos + len("127.0.0.1:")):slashPos])
                                if port == portNumber:
                                    foundPort = True
                                    continue
                                else:
                                    tempFile.write(line)
            if foundPort:
                shutil.copy(tempFnp, self.outfileFnp_)
            else:
                sys.stderr.write("Didn't find port number: " + str(portNumber) + " in " + self.outfileFnp_  + " nothing done\n")
    def removeName(self, rmName):
        tempFnp = "/tmp/temp_" + os.path.basename(self.outfileFnp_)
        if os.path.exists(self.outfileFnp_):
            foundPort = False;
            with open(self.outfileFnp_, 'r') as serverConfFile:
                with open(tempFnp, "w") as tempFile:
                    for line in serverConfFile:
                        toks = line.strip().split()
                        if len(toks) == 3:
                            name = toks[1][1:];
                            colonPos = toks[2].find("127.0.0.1:")
                            slashPos = toks[2].find("/", colonPos)
                            if -1 != colonPos and -1 != slashPos:
                                port = int(toks[2][(colonPos + len("127.0.0.1:")):slashPos])
                                if name == rmName:
                                    foundPort = True
                                    continue
                                else:
                                    tempFile.write(line)
            if foundPort:
                shutil.copy(tempFnp, self.outfileFnp_)
            else:
                sys.stderr.write("Didn't find name: " + str(rmName) + " in " + self.outfileFnp_  + " nothing done\n")
                                      



def parse_args_genServerProxyPassConfig():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--removePort',type=int)
    parser.add_argument('--removeName',type=str)
    parser.add_argument('--port',      type=int)
    parser.add_argument('--name',      type=str)
    parser.add_argument('--start',     type=int, default = 0)
    parser.add_argument('--stop',      type=int, default = 9)
    parser.add_argument('--serverBase',type=int, default = 10000)
    
    parser.add_argument('--outfile',   type=str, required = True)
    parser.add_argument('--append',    action = "store_true")
    parser.add_argument('--overWrite', action = "store_true")
    
    return parser.parse_args()

def genServerProxyPassConfig():
    args = parse_args_genServerProxyPassConfig()
    if args.removePort:
        creator = MipServerConfCreator(args.outfile)
        creator.removePort(args.removePort)
        sys.exit(0)
    if args.removeName:
        creator = MipServerConfCreator(args.outfile)
        creator.removeName(args.removeName)
        sys.exit(0)
    mode = "write"
    if args.append:
        mode = "append"
    if args.overWrite:
        mode = "overWrite"
    append = False;
    overWrite = False;
    if mode == "append":
        append = True;
    if mode == "overWrite":
        overWrite = True;
    if os.path.exists(args.outfile) and not overWrite and not append:
        raise Exception("Error, " + args.outfile + " already exits, use --overWrite to overWrite it or --append to add to it")
    if os.path.exists(args.outfile) and overWrite:
        os.remove(args.outfile)
    
    creator = MipServerConfCreator(args.outfile)
    if args.port and args.name:
        creator.addServer(args.name, args.port)
    else:
        creator.addMipServers(args.serverBase, args.start, args.stop + 1)
    


if __name__ == "__main__":
    genServerProxyPassConfig()




