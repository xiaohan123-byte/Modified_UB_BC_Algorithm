import numpy as np
from gurobipy import *
from input import *
from UBnBC_Node import Node
from Heuristic import Heuristic,Heuristic_gurobi
import heapq
import time

'''
stations = ['s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11']
S=stations

xS=[6 for _ in range(len(stations))] # 给每个station 6 块标准续航电池
xL=[6 for _ in range(len(stations))] # 给每个station 6 块长续航电池

k=0
XS={}
XL={}
for s in S:
    XS[s]=xS[k]
    XL[s]=xL[k]
    k+=1

print(XL)
print(XS)
'''

'''
with open('data/data_case.json', 'r') as f:
    path_data = json.load(f)

for p in path_data:
    source[p] = path_data[p]['source']
    AS[p] = tuplelist([tuple(arc) for arc in path_data[p]['strategyS']])
    AL[p] = tuplelist([tuple(arc) for arc in path_data[p]['strategyL']])
    Sq[p] = path_data[p]['stations']

path = path_data.keys()

print(path)
'''