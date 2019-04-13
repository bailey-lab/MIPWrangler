MIPWrangler
========
Version 1.1.0

Bioinformatics tools to analyze molecular inversion probe (MIP) sequencing results Bailey Lab

Checkout the website bellow for more details  
[http://baileylab.brown.edu/MIPWrangler/](http://baileylab.brown.edu/MIPWrangler/)

Publications:  Aydemir O, et al. (2018) Drug-Resistance and Population Structure of Plasmodium falciparum Across the Democratic Republic of Congo Using High-Throughput Molecular Inversion Probes. J Infect Dis 218(6):946–955.

# Installing  
 
 See installing tab on [http://baileylab.brown.edu/MIPWrangler/](http://baileylab.brown.edu/MIPWrangler/) for full details for installing for each operating system. 
 
## Dependecnies
Need to have a c++ compiler, by default g++-7 (Ubuntu,RedHat) and clang++ (MAC OS). The install scripts use python3 which will also be needed for the install steps.  

Also though MIPWrangler does not use cmake, several of the libraries it uses do depend on cmake so it needs to be present.  

## To Install Latest Version    

The install.sh script will download other c++ libraries and compile them and then compile MIPWrangler. 

```bash
git clone https://github.com/bailey-lab/MIPWrangler.git   
cd MIPWrangler  
./install.sh
make   
```

To add MIPWrangler to path

```bash
export PATH=$HOME/MIPWrangler/bin/:$PATH
```




# Bash Completion  

MIPWrangler tends to have long flags so that they can be clear what they do but it's somewhat annoying to type them out so bash completion has been added.  Put the content of the file at bashCompletion/MIPWrangler into a file ~/.bash_completion and it will be source on your next login or use the bellow command while in the MIPWrangler directory  

```bash
./setup.py --addBashCompletion  
```

Which will actually do exactly described above, afterwards while typing flags use the tab key to complete them  


# Tutorials

Tutorials and detailed usages will be located at [http://baileylab.brown.edu/MIPWrangler](http://baileylab.brown.edu/MIPWrangler) or email nickjhathaway@gmail.com for more information  


