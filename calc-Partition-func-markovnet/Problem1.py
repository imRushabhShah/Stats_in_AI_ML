# @author: Rushabh Shah
# @netid: rms170003

import re
import numpy as np
import copy
import math
import sys


def readEvidance(file, graph):
    with open(file) as f:
        ev = [int(x) for x in f.readlines()[0].strip().split()]
    if ev[0]==0:
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
                    new_CPT[new_CPT_point:new_CPT_point+stride+1] = old_CPT[i:i+stride+1]
                    new_CPT_point+=stride
                graph.CPT[func] = new_CPT
        graph.cardianlities[var]=1
            
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
            markov.CPT[clique]= np.asarray([float(x) for x in lines[index].strip().split()])
            
        return markov
    else:
        print(lines[index].strip()+" yet to write")
    

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
    

    def __isEmpty__(self):
        if self.variables == []:
            return True
        return False
    
    def min_degree(self):
        VE_Order = []
        node = copy.copy(self.node)
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
                self.bucket[v].append(func)
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
                    var_functions.append(new_clique)
            if len(var_functions)>0:
                return var_functions
        print("No function for the variable", var)
        return False
    
    def create_bucket_functions(self):
        self.bucket = [[] for _ in range(self.variable_count)]
        for func in self.CPT.keys():
            self.add_to_bucket(func)
    
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
        # state is a tuple which has the order of the new vars and value assigned and cardinalities
        # func = (0,1)
        # state.order = (0,1,2)
        # state.value = (1,0,0)
        # cpt = [4.481689, 14.481689, 21.      , 40.481689]
        cpt = self.CPT[func]
        # state.cardinality = [2,2,2]
        index = 0 
        pre_card = 0
        func_len = len(func)
        for f_index in range(len(func)):
            f_order = state["order"].index(func[len(func) - f_index - 1])
            index += pre_card**f_index * state["value"][f_order] 
            pre_card = state["cardinality"][f_order]    
        #  print(func,state)
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
        #  clique_CPT = np.ones((self.cardianlities[var], new_clique_cardinality_size))
        temp_clique = [var]+list(new_clique)
        temp_clique_cardinalities = []
        for i in temp_clique:
            temp_clique_cardinalities.append(self.cardianlities[i])
        state = {}
        state["order"] = temp_clique
        state["cardinality"] = temp_clique_cardinalities
        state["value"] = [0]*len(temp_clique) 
        state["value"] = [0]*len(temp_clique)
        index = len(temp_clique)
        var_sum = np.zeros(new_clique_cardinality_size)
        for i in range(state["cardinality"][0]):
            val_product = np.ones(new_clique_cardinality_size)
            state["value"][0]=i
            for j in range(new_clique_cardinality_size):
                j_copy = copy.copy(j)
                for k in range(len(temp_clique)-1,0,-1):
                    state["value"][k] = j_copy%temp_clique_cardinalities[k]
                    j_copy=j_copy//temp_clique_cardinalities[k]
                for func in var_functions:
                    val_product[j] = np.multiply(val_product[j], self.get_value(func , state))                    
            var_sum += val_product
        return new_clique,var_sum
                
    def do_VE(self):
        for min_var in self.VE_order:
            new_clique, new_clique_CPT = self.sum_product(min_var)
            if new_clique in self.CPT.keys():
                old_clique_CPT = self.CPT[new_clique]
                new_clique_CPT = np.multiply(new_clique_CPT,old_clique_CPT)
            else:
                self.add_to_bucket(new_clique)
            self.CPT[new_clique] = new_clique_CPT
            # print(new_clique,new_clique_CPT)
        return self.CPT[()]
            


#     mar = readfile('1.uai')
if (len(sys.argv) < 2):
    print("read files missing")

elif (len(sys.argv) < 3):
    print("Computing without evidence file")
    mar = readfile(sys.argv[1])
#     readEvidance(file=sys.argv[2], graph=mar)
    mar.min_degree()
    print("VE order ",mar.VE_order)
    mar.create_bucket_functions()
    ans = mar.do_VE()
    print("Partition function value",math.log10(ans))

else:
    mar = readfile(sys.argv[1])
    readEvidance(file=sys.argv[2], graph=mar)
    mar.min_degree()
    print("VE order ",mar.VE_order)
    mar.create_bucket_functions()
    ans = mar.do_VE()
    print("log Partition function value",math.log10(ans))




# mar = readfile('2.uai')
# readEvidance(file="2.uai.evid", graph=mar)
# # mar.VE_order.sort()
# mar.min_degree()
# mar.create_bucket_functions()
# print(math.log10(mar.do_VE()))


# # In[9]:


# mar = readfile('3.uai')
# readEvidance(file="3.uai.evid", graph=mar)
# # mar.VE_order.sort()
# mar.min_degree()
# mar.create_bucket_functions()
# print(math.log10(mar.do_VE()))


