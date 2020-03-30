The Problem was consuming a lot of time and due to lots of server crashes and many kills on remote server, some of the values might be reduandant and while others are missing, detailed table is present in the main pdf file do consider that.
adaptive.csv contains all run for adaptive distribution and uniform.csv for uniform distribution
Problem6.csv contains a small amount of samples over which all the calculation are done which are present in the main pdf file 

Please install

numpy, tqdm, 
pip install numpy tqdm


TO RUN the files 

> python .\Adaptive.py

input file :Grids_15.uai
evidance file :Grids_15.uai.evid
w cutset: 4
N(number of itteration) :100
number of samples :3
write file name to store all the samples(.csv preferable) :temp.csv

> python .\UniformMarkov.py

input file :Grids_15.uai
evidance file :Grids_15.uai.evid
w cutset: 4
N(number of itteration) :100
number of samples :3
write file name to store all the samples(.csv preferable) :temp.csv



this is going to create a temp.csv(or any other file you want) file and store all the data in it for further analysis or computations.