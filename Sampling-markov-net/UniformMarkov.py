#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
from tqdm import tqdm
import numpy as np
import copy
import math
import random
import time


# In[2]:


def readEvidance(file, graph):
    if type(file) == str:
        with open(file) as f:
            ev = [int(x) for x in f.readlines()[0].strip().split()]
#             print(ev)
    else:
        ev = file
#         print(ev)
    if ev[0] == 0:
        return None
    for evidance_index in range(ev[0]):
        var = ev[(2*evidance_index)+1]
        val = ev[(2*evidance_index)+2]
        for func in graph.CPT.keys():
            if var in func:
                var_index = func.index(var)
                new_size = 1
                for v_index in range(len(func)):
                    if not v_index == var_index:
                        new_size *= graph.cardianlities[func[v_index]] 
                new_CPT = np.zeros(new_size)
                stride = 1
                for v_index in range(len(func)-2,var_index,-1):                    
                    stride *= graph.cardianlities[func[v_index+1]]
                
                new_CPT_point = 0
                old_CPT = graph.CPT[func]
                start_point = stride*val
                end_point = len(graph.CPT[func])
                step = stride*graph.cardianlities[func[var_index]]
                
                for i in range(start_point,end_point,step):
                    new_CPT[new_CPT_point:new_CPT_point+stride] = old_CPT[i:i+stride]
#                     new_CPT[new_CPT_point:new_CPT_point+stride+1] = old_CPT[i:i+stride+1]
                    new_CPT_point+=stride
                
                graph.CPT[func] = new_CPT
#                 print(var,val,func,old_CPT,new_CPT)
        graph.cardianlities[var]=1


# In[3]:


def readfile(file):
    with open(file) as f:
        lines = f.readlines()    

    index = nextIndex(lines, -1)
        
    if lines[index].strip() =='MARKOV':
        # storing variable_count 
        markov = Graph() 
        index = nextIndex(lines, index)
        
        markov.variable_count = int(lines[index].strip())

        # storing cardinality 
        index = nextIndex(lines, index)
        markov.cardianlities = [int(x) for x in lines[index].strip().split()]
        
        # storing function count
        index = nextIndex(lines, index)
        markov.function_count = int(lines[index].strip())
        
        for function_index in range(markov.function_count):
            index = nextIndex(lines, index)
            assert index,"invalid file format"
            clique = [int(x) for x in lines[index].strip().split()]
            markov.add(clique[1:])
            
        for function_index in range(markov.function_count):
            index = nextIndex(lines, index)
            assert index,"invalid file format"
            count = int(lines[index].strip())
            index = nextIndex(lines, index)
            assert index,"invalid file format"
            clique = tuple(markov.functions[function_index])
#             clique = '-'.join(str(x) for x in clique)
            
            markov.CPT[clique]= np.log(np.asarray([float(x) for x in lines[index].strip().split()]))
#             markov.CPT[clique]= np.asarray([float(x) for x in lines[index].strip().split()])
            
        return markov
    else:
        print(lines[index].strip()+" yet to write")


# In[4]:


def nextIndex(lines, index=-1):
    # to handell blank lines or comments
    index+=1
    while index< len(lines):
        content = lines[index].strip().split()
        if len(content)==0 or content[0]=='c': 
            index+=1
        else: 
            return index
    if index >= len(lines):
        print("Reached end of file")
        return None


# In[5]:


class Graph(object):
    def __init__(self):
        self.variable_count = 0
        self.cardianlities = []
        self.variables = []
        self.node = {}
        self.CPT = {}
        self.functions = []
        self.function_count = 0
        self.VE_order = []
        self.bucket = []
        self.bucket_var = []

    def __isEmpty__(self):
        if self.variables == []:
            return True
        return False
    
    def min_degree(self):
        VE_Order = []
        node = copy.deepcopy(self.node)
        while not self.__isEmpty__():
            min_var = None
            min_val = float('inf')
            for v in self.variables:
                v_count = len(node[v].intersection(self.variables))
                if v_count < min_val:
                    min_val = v_count
                    min_var = v
            VE_Order.append(min_var)
            self.remove_node_var(min_var, node)
        self.VE_order = VE_Order
        return VE_Order
    
    def add(self,var): #var is list of variables forming a clique 
        self.functions.append(var)
        for i in range(len(var)):
            if var[i] in self.variables:
                self.node[var[i]] = self.node[var[i]].union(var[0:i]+var[i+1:])
            else:
                self.variables.append(var[i])
#                 print(var,var[:i]+var[i+1:])
                self.node[var[i]] = set(var[0:i]+var[i+1:])
    
    def remove_node_var(self, var, node):
        if self.__isEmpty__():
            return False
        else:
            if var in self.variables:
                self.variables.remove(var)
                for adj_node in node[var]:
                    node[adj_node] = node[var].union(node[adj_node])-{adj_node}
            else:
                print(var+" does not exists")
                return False
            return True

    def add_to_bucket(self, func):
        for v in self.VE_order:
            if v in func:
                if func in self.bucket[v]:
                    break
                self.bucket[v].append(func)
                self.bucket_var[v]=self.bucket_var[v].union(func) - {v}
                break
    
    def remove_var_functions(self, var):
        var_functions = [] 
        if len(self.functions) > 0:
            for func in self.functions:
                if var in func:
                    var_functions.append(func)
        for func in var_functions:
            self.functions.remove(func)
        self.function_count = len(self.functions)
        
    def get_fuctions_with_var(self, var):
        var_functions = [] 
        if len(self.functions) > 0:
            for i in range(len(self.functions)):
                if var in self.functions[i]:
                    new_clique = tuple(self.functions[i])
#                     new_clique = "-".join(str(x) for x in self.functions[i])
                    var_functions.append(new_clique)
            if len(var_functions)>0:
                return var_functions
        print("No function for the variable", var)
        return False
    
    def create_bucket_functions(self):
        self.bucket = [[] for _ in range(self.variable_count)]
        self.bucket_var = [set() for _ in range(self.variable_count)]
        for func in self.CPT.keys():
            self.add_to_bucket(func)
        for min_var in self.VE_order:
            var_functions = self.bucket[min_var]
            new_clique, _ = self.get_new_clique(min_var,var_functions)
            self.add_to_bucket(new_clique)
    
    
    def get_new_clique(self,var, func):
        new_clique = set({})
        for f in func:
            new_clique = new_clique.union(set(f))
        new_clique -= {var} 
        new_clique = list(new_clique)
        list(new_clique).sort()
        new_clique_cardinality = [self.cardianlities[x] for x in new_clique]
        new_clique = tuple(new_clique)
        return new_clique,new_clique_cardinality
    
    
    def get_value(self,func , state):#        
#         state={}
#         func = (0,1,2)
#         state["order"] = (0,1,2,3)
#         state["value"] = (1,1,0,2)
#         cpt = [4.481689, 14.481689, 21.      , 40.481689]
        
#         state["cardinality"] = [2,2,1,2]
        cpt = self.CPT[func]
        index = 0 
        pre_card = 1
        func_len = len(func)
        for f_index in range(len(func)-1,-1,-1):
            f_order = state["order"].index(func[f_index])
            index += pre_card * state["value"][f_order] 
            pre_card *= state["cardinality"][f_order]  
        return cpt[index]

    
    def sum_product(self, var):
        if not self.bucket:
            self.create_bucket_functions()
        var_functions = self.bucket[var]
        assert var_functions,"Invalid VE Order"
        new_clique,new_clique_cardinality = self.get_new_clique(var,var_functions)
        new_clique_cardinality_size = 1 
        for i in new_clique_cardinality:
            new_clique_cardinality_size*=i 
#         clique_CPT = np.ones((self.cardianlities[var], new_clique_cardinality_size))
        temp_clique = [var]+list(new_clique)
        temp_clique_cardinalities = []
        for i in temp_clique:
            temp_clique_cardinalities.append(self.cardianlities[i])
        state = {}
        state["order"] = temp_clique
        state["cardinality"] = temp_clique_cardinalities
#         for i in range(len(temp_clique)):
#             state["order"][temp_clique[]]
        state["value"] = [0]*len(temp_clique) 
        state["value"] = [0]*len(temp_clique)
        index = len(temp_clique)
        var_sum = None #np.zeros(new_clique_cardinality_size,dtype=np.float64)
        for i in range(state["cardinality"][0]):
            val_product = np.zeros(new_clique_cardinality_size,dtype=np.float64)
            state["value"][0]=i
            for j in range(new_clique_cardinality_size):
                j_copy = copy.deepcopy(j)
                for k in range(len(temp_clique)-1,0,-1):
                    state["value"][k] = j_copy%temp_clique_cardinalities[k]
                    j_copy=j_copy//temp_clique_cardinalities[k]
                for func in var_functions:              
                    val_product[j] = val_product[j] + self.get_value(func , state)                    
#                     val_product[j] = np.multiply(val_product[j], self.get_value(func , state),dtype=np.float64)                    
            if not type(var_sum) == np.ndarray:
                var_sum = val_product
            else:
                var_sum = np.logaddexp(var_sum,val_product)
#             var_sum += val_product
        return new_clique,var_sum
    
    def do_VE(self):
        for min_var in self.VE_order:
            new_clique, new_clique_CPT = self.sum_product(min_var)
            if new_clique in self.CPT.keys():
                old_clique_CPT = self.CPT[new_clique]
                new_clique_CPT = new_clique_CPT + old_clique_CPT
#                 new_clique_CPT = np.multiply(new_clique_CPT,old_clique_CPT,dtype=np.float64)
            else:
                pass
            self.CPT[new_clique] = new_clique_CPT
#             print(new_clique,new_clique_CPT)
        return (self.CPT[()]/np.log(10))[0]
    
    def max_width(self):
        return max([len(x) for x in self.bucket_var])
    
    def max_occuring_var(self):
        VC = np.zeros(self.variable_count,dtype=np.float64)
        for buc_vr in self.bucket_var:
            for v in buc_vr:
                VC[v] += 1
        max_Count = max(VC)
        all_max_count_var = [i for i in range(self.variable_count) if VC[i]==max_Count]
        return random.choice(all_max_count_var)
    
    def remove_from_cluster(self,Xi):
        for i in range(self.variable_count):
            if Xi in self.bucket_var[i]:
                self.bucket_var[i].remove(Xi)
    
    def wCutset(self,w):
        X = []
        while(self.max_width() > w+1):
            Xi = self.max_occuring_var()
            self.remove_from_cluster(Xi)
            X.append(Xi)
        return X


# In[6]:


def multiplyList(myList) : 
      
    # Multiply elements one by one 
    result = 1
    for x in myList: 
         result = result * x  
    return result 


# In[7]:

'''
mar = readfile('1.uai')
readEvidance(file="1.uai.evid", graph=mar)
mar.min_degree()
mar.create_bucket_functions()
ans = mar.do_VE()
print(ans)
# print(math.log10(ans))
print(14.889866514255774)
'''

# In[8]:


def uniformQ(G,X):
    Q = []
    for x in X:
        Q.append(np.ones(G.cardianlities[x],dtype=np.float64)/G.cardianlities[x])
    return Q


# In[ ]:


def UniformAlgorithm1(G,w,N):
    start = time.time()
    Z = 0
    X = G.wCutset(w)
    Q = uniformQ(G,X)
#     for x in X:
#         Q.append(np.ones(G.cardianlities[x],dtype=np.float64)/G.cardianlities[x])
#     print(Q)
    vals = []
    total_time = []
    for i in tqdm(range(max(N))):
        G_i = copy.deepcopy(G)
        evidance,Q_val = generateSample(G,X,Q)
#         print(evidance)

        readEvidance(file = evidance, graph = G_i )        
        numerator = G_i.do_VE()
        denominator = Q_val
#         print(numerator,denominator)
        w = numerator - np.log10(denominator)
        Z = Z + w
        if i+1 in N:
            vals.append(Z/i+1)
            total_time.append(time.time()-start)
    return vals,total_time

# In[ ]:


def generateSample(G,X,Q):
    Q_val = 1
    evidance = [len(Q)]
    for i in range(len(Q)):
        evidance.append(X[i])
        rnd = random.random()
        insert = None
        for bound in range(len(Q[i])):
            if rnd < Q[i][bound]:
                insert = bound
                Q_val *= Q[i][bound]
                break
            else:
                rnd-=Q[i][bound]
        evidance.append(insert)
    return evidance,Q_val


# In[ ]:

'''
for w_cut_size in [1,2,3,4,5]:
    for N in [100,1000,10000,20000]:
        for itteration in range(10):
            mar = readfile('hw3-files/Grids_15.uai')
            readEvidance(file="hw3-files/Grids_15.uai.evid", graph=mar)

            mar.min_degree()
            mar.create_bucket_functions()

            start = time.time()
            G = UniformAlgorithm1(mar,w=w_cut_size,N=N)
            end = time.time()
            total_time = end-start
            with open("uniform.csv",'a+') as u:
                state = "Grids_15,"+str(w_cut_size)+","+str(N)+","+str(itteration)+","+str(G)+","+str(total_time)+"\n" 
                print(state)
                u.write(state)


# In[ ]:


for w_cut_size in [1,2,3,4,5]:
    for N in [100,1000,10000,20000]:
        for itteration in range(10):
            mar = readfile('hw3-files/Grids_14.uai')
            readEvidance(file="hw3-files/Grids_14.uai.evid", graph=mar)

            mar.min_degree()
            mar.create_bucket_functions()

            start = time.time()
            G = UniformAlgorithm1(mar,w=w_cut_size,N=N)
            end = time.time()
            total_time = end-start
            with open("uniform.csv",'a+') as u:
                state = "Grids_14,"+str(w_cut_size)+","+str(N)+","+str(itteration)+","+str(G)+","+str(total_time)+"\n" 
                print(state)
                u.write(state)


# In[ ]:


for w_cut_size in [1,2,3,4,5]:
    for N in [100,1000,10000,20000]:
        for itteration in range(10):
            mar = readfile('hw3-files/Grids_16.uai')
            readEvidance(file="hw3-files/Grids_16.uai.evid", graph=mar)

            mar.min_degree()
            mar.create_bucket_functions()

            start = time.time()
            G = UniformAlgorithm1(mar,w=w_cut_size,N=N)
            end = time.time()
            total_time = end-start
            with open("uniform.csv",'a+') as u:
                state = "Grids_16,"+str(w_cut_size)+","+str(N)+","+str(itteration)+","+str(G)+","+str(total_time)+"\n" 
                print(state)
                u.write(state)


# In[ ]:


for w_cut_size in [1,2,3,4,5]:
    for N in [100,1000,10000,20000]:
        for itteration in range(10):
            mar = readfile('hw3-files/Grids_17.uai')
            readEvidance(file="hw3-files/Grids_17.uai.evid", graph=mar)

            mar.min_degree()
            mar.create_bucket_functions()

            start = time.time()
            G = UniformAlgorithm1(mar,w=w_cut_size,N=N)
            end = time.time()
            total_time = end-start
            with open("uniform.csv",'a+') as u:
                state = "Grids_17,"+str(w_cut_size)+","+str(N)+","+str(itteration)+","+str(G)+","+str(total_time)+"\n" 
                print(state)
                u.write(state)


# In[ ]:



for w_cut_size in [1,2,3,4,5]:
    for N in [100,1000,10000,20000]:
        for itteration in range(10):
            mar = readfile('hw3-files/Grids_18.uai')
            readEvidance(file="hw3-files/Grids_18.uai.evid", graph=mar)

            mar.min_degree()
            mar.create_bucket_functions()

            start = time.time()
            G = UniformAlgorithm1(mar,w=w_cut_size,N=N)
            end = time.time()
            total_time = end-start
            with open("uniform.csv",'a+') as u:
                state = "Grids_18,"+str(w_cut_size)+","+str(N)+","+str(itteration)+","+str(G)+","+str(total_time)+"\n" 
                print(state)
                u.write(state)
'''
def prog(w_cut_size, N, Grid, writefile,itt = 10):
    for itteration in range(itt):
        mar = readfile('hw3-files/'+Grid+'.uai')
        readEvidance(file="hw3-files/"+Grid+".uai.evid", graph=mar)
        
        mar.min_degree()
        mar.create_bucket_functions()

        Gs,total_time = UniformAlgorithm1(mar,w=w_cut_size,N=N)
        with open(writefile,'a+') as u:                
            for n in range(len(N)):
                state = Grid+","+str(w_cut_size)+","+str(N[n])+","+str(itteration)+","+str(Gs[n])+","+str(total_time[n])+"\n" 
                print(state)
                u.write(state)
'''            
for Grid in ['Grids_18','Grids_17','Grids_16','Grids_14','Grids_15']:
    for w_cut_size in [5,4,3,2,1]:
        N = [100,1000,10000,20000]
        prog(w_cut_size,N,Grid,'tmp.csv',5)
'''
def prog_generic(w_cut_size, N, inp_file, evid_file, writefile,itt = 10):
    for itteration in range(itt):
        mar = readfile(inp_file)
        readEvidance(file= evid_file, graph=mar)
        
        mar.min_degree()
        mar.create_bucket_functions()

        Gs,total_time = UniformAlgorithm1(mar,w=w_cut_size,N=N)
        with open(writefile,'a+') as u:                
            for n in range(len(N)):
                state = inp_file+","+str(w_cut_size)+","+str(N[n])+","+str(itteration)+","+str(Gs[n])+","+str(total_time[n])+",Uniform Distribution\n" 
                print("Grid, cutset, N, sample_num, log10(part), time in sec, type of algo")
                print(state)
                u.write(state)
inp = input("input file :")
evid = input("evidance file :")
# Grid = 'Grids_16'
w_cut_size = int(input("w cutset: "))
N = [int(input("N(number of itteration) :"))]
#N = [100,1000,10000,20000]
samples = int(input("number of samples :"))
#prog(w_cut_size,N,Grid,'adaptive_16.csv',5)
writefile = input("write file name to store all the samples(.csv preferable) :")
prog_generic(w_cut_size, N, inp_file = inp, evid_file= evid, writefile = writefile,itt = samples)