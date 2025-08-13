import json
import pandas as pd
import numpy as np
import random
import math
from gurobipy import *
import time


random.seed(2)
# input parameters

num_scenario = 20
prob = (0.05,0.05)*10  # 每个scenario的概率都是0.1

I=21  #电池数量上限
cS = 300 #初始200，300gurobi差 标准续航电池的配置费用
cL = 350 #初始250，350gurobi差 长续航电池的配置费用
alpha = 300 #初始300
beta = 25 #charging power
gamma = 15 
E = (0.5,0.5,0.5,0.5,0.5,0.5,0.9,0.9,1.4,1.4,1.4,1.4,1.4,1.4,1.4,0.9,0.9,0.9,1.4,1.4,1.4,1.4,0.9,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.9,0.9,1.4,1.4,1.4,1.4) # 每个period的电价

KS =3 #标准电池的电池容量
KL =4 #长电池的电池容量
TimePeriods = 24#本来是36  
max_tau=12 # ？
SwapPriceS = tuple((x+0.5)*75 for x in E) #标准续航电池在每个时间段的换电费
SwapPriceL = tuple((x+0.5)*100 for x in E) #长续航电池在每个时间段的换电费

#stations = ['s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11']
stations = ['s1','s2','s3','s4','s5','s6']

XSbar = {}
XLbar = {}

for i in stations:
    XSbar[i] = 4 # 每个电站的初始标准续航电池数量
    XLbar[i] = 2 # 每个电站的初始长续航电池数量

T_matrix = -np.eye(2*len(stations))
#print(T_matrix)

# 从文件读取
with open('efficiency_data/data_6e.json', 'r') as f:
    path_data = json.load(f)
#print(path_data)
dist = pd.read_csv('efficiency_data/dist_6e.csv',header=0,index_col=0)
dist = dist.applymap(np.abs)

# input demand data
#f_data = read_data('f_data.csv')

with open('efficiency_data/demand_flow_S_6e.json', 'r') as file:
    json_data_fS = json.load(file)
fS = {eval(k): v for k, v in json_data_fS.items()}

with open('efficiency_data/demand_flow_L_6e.json', 'r') as file:
    json_data_fL = json.load(file)
fL = {eval(k): v for k, v in json_data_fL.items()}

tau={}
for p in path_data:
    for i in stations:
        tau[p,i] = math.ceil(dist[path_data[p]['source']][i]/100)

# input network data
AS={}
AL={}
Sq={}
source={}

for p in path_data:
    source[p] = path_data[p]['source']
    AS[p] = tuplelist([tuple(arc) for arc in path_data[p]['strategyS']])
    AL[p] = tuplelist([tuple(arc) for arc in path_data[p]['strategyL']])
    Sq[p] = path_data[p]['stations']

'''
输出: AS
{'p1': <gurobi.tuplelist (3 tuples, 2 values each):
 ( e1 , s1 )
 ( s1 , s2 )
 ( s2 , e3 )
>, 'p2': <gurobi.tuplelist (5 tuples, 2 values each):
 ( e1 , s1 )
 ( s3 , s4 )
 ( s1 , s2 )
 ( s4 , e5 )
 ( s2 , s3 )
'''
'''
输出: Sq

{'p1': ['s1', 's2'], 'p2': ['s1', 's2', 's3', 's4'], 'p3': ['s1', 's2', 's3', 's4', 's5', 's6']}

'''

T=TimePeriods 
#36
S=stations
# ['s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11']
path = path_data.keys()
# 输出 ['p1', 'p2', 'p3', 'p4', 'p5', 'p6'...] # 这里的 P 就是需求流 q

#--------------------------------------------------------------------------------------
#Build Gurobi Model



def BuildMasterP(num_scenario,prob): #整数规划主问题
    masterP = Model('masterP')

    # Add Variables
    XS = masterP.addVars(S,vtype=GRB.INTEGER,lb=0,ub=I,name='XS')
    XL = masterP.addVars(S,vtype=GRB.INTEGER,lb=0,ub=I,name='XL')
    PXS = masterP.addVars(S,vtype=GRB.INTEGER,lb=0,ub=I,name='PXS') # 新分配电池的数量 = 总的-原有的
    PXL = masterP.addVars(S,vtype=GRB.INTEGER,lb=0,ub=I,name='PXL')
    theta = masterP.addVars(num_scenario,vtype=GRB.CONTINUOUS,lb=-999999999999,name='theta')
    masterP._XS = XS
    masterP._XL = XL
    masterP._theta = theta

    # Add Constraints
    masterP.addConstrs(XS[i]+XL[i]<=I for i in S)
    masterP.addConstrs(PXS[i] >= XS[i] - XSbar[i] for i in S)
    masterP.addConstrs(PXL[i] >= XL[i] - XLbar[i] for i in S)

    # Add objective function
    masterP.setObjective(cS*quicksum(PXS[i] for i in S)+cL*quicksum(PXL[i] for i in S)+quicksum(prob[scenario_i]*theta[scenario_i] for scenario_i in range(num_scenario)))

    masterP.update()

    '''# 获取所有约束
    all_constraints = masterP.getConstrs()

    # 打印约束信息
    for constr in all_constraints:
        print(f"Constraint '{constr.ConstrName}': {constr.sense} {constr.RHS}")
        print("Linear part:", masterP.getRow(constr))
        print()'''

    return masterP

# 松弛子问题 relax SP
def BuildRelaxP(scenario_i,xS,xL):
    relaxP = Model('relaxP')

    k=0
    XS={}
    XL={}
    for s in S:
        XS[s]=xS[k]
        XL[s]=xL[k]
        k+=1

    '''
    {'s1': 6, 's2': 6, 's3': 6, 's4': 6, 's5': 6, 's6': 6, 's7': 6, 's8': 6, 's9': 6, 's10': 6, 's11': 6}
    {'s1': 6, 's2': 6, 's3': 6, 's4': 6, 's5': 6, 's6': 6, 's7': 6, 's8': 6, 's9': 6, 's10': 6, 's11': 6}
    '''

    # Add variables (continuous)
    yS=tupledict() # gurobi的元组字典类型 VarName: yS[p,t,ij]
    yL=tupledict()
    yS_h=tupledict()
    yL_h=tupledict()
    for p in path:
        for t in range(-max_tau,T):
            for ij in AS[p]:
                yS[p,t,ij] = relaxP.addVar(lb=0,vtype=GRB.CONTINUOUS,name='yS')
                yL_h[p,t,ij] = relaxP.addVar(lb=0,vtype=GRB.CONTINUOUS,name='yL')
            for ij in AL[p]:
                yL[p,t,ij] = relaxP.addVar(lb=0,vtype=GRB.CONTINUOUS,name='yL')
                yS_h[p,t,ij] = relaxP.addVar(lb=0,vtype=GRB.CONTINUOUS,name='yS')
        '''for t in range(-max_tau,0,1):
            for ij in AS[p]:
                yS[p,t,ij] = 0
                yL_h[p,t,ij] = 0
            for ij in AL[p]:
                yL[p,t,ij] = 0
                yS_h[p,t,ij] = 0'''
        

    wS = relaxP.addVars(S,T,KS+1,lb=0,vtype=GRB.CONTINUOUS,name='wS')
    wL = relaxP.addVars(S,T,KL+1,lb=0,vtype=GRB.CONTINUOUS,name='wL')
    zS = relaxP.addVars(S,T,KS,lb=0,vtype=GRB.CONTINUOUS,name='zS')
    zL = relaxP.addVars(S,T,KL,lb=0,vtype=GRB.CONTINUOUS,name='zL')

    # Add contraints

    relaxP.addConstrs(yS[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AS[p])
    relaxP.addConstrs(yL_h[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AS[p])
    relaxP.addConstrs(yL[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AL[p])
    relaxP.addConstrs(yS_h[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AL[p])

    relaxP.addConstrs(quicksum(yS[q,t,ij] for ij in AS[q].select(source[q],'*')) <= fS[scenario_i,q,t] for q in path for t in range(T))
    relaxP.addConstrs(quicksum(yL[q,t,ij] for ij in AL[q].select(source[q],'*')) <= fL[scenario_i,q,t] for q in path for t in range(T))

    relaxP.addConstrs(quicksum(yS_h[q,t,ij] for ij in AL[q].select(source[q],'*')) == 0 for q in path for t in range(T))
    relaxP.addConstrs(quicksum(yL_h[q,t,ij] for ij in AS[q].select(source[q],'*')) == 0 for q in path for t in range(T))

    relaxP.addConstrs(quicksum(yS[q,t,ij] for ij in AS[q].select(i,'*')) + quicksum(yS_h[q,t,ij] for ij in AL[q].select(i,'*')) 
             == quicksum(yS[q,t,ij] for ij in AS[q].select('*',i)) + quicksum(yS_h[q,t,ij] for ij in AL[q].select('*',i)) for q in path for t in range(T) for i in Sq[q])
    relaxP.addConstrs(quicksum(yL[q,t,ij] for ij in AL[q].select(i,'*')) + quicksum(yL_h[q,t,ij] for ij in AS[q].select(i,'*')) 
             == quicksum(yL[q,t,ij] for ij in AL[q].select('*',i)) + quicksum(yL_h[q,t,ij] for ij in AS[q].select('*',i)) for q in path for t in range(T) for i in Sq[q])

    XconstrS=relaxP.addConstrs((wS[i,0,KS] == XS[i] for i in S), name = 'Xconstr')
    XconstrL=relaxP.addConstrs((wL[i,0,KL] == XL[i] for i in S), name = 'Xconstr')

    relaxP._XconstrS = XconstrS
    relaxP._XconstrL = XconstrL

    relaxP.addConstrs(wS[i,T-1,k] == 0 for i in S for k in range(KS))
    relaxP.addConstrs(wL[i,T-1,k] == 0 for i in S for k in range(KL))

    relaxP.addConstrs(wS[i,0,k]==0 for i in S for k in range(KS))
    relaxP.addConstrs(wS[i,t,0]==wS[i,t-1,0]-zS[i,t-1,0]+quicksum(yS[q,t-1-tau[q,i], ij]+yL_h[q,t-1-tau[q,i], ij] for q in path for ij in AS[q].select('*',i)) for i in S for t in range(1,T))
    relaxP.addConstrs(wS[i,t,n]==wS[i,t-1,n]+zS[i,t-1,n-1]-zS[i,t-1,n] for i in S for t in range(1,T) for n in range(1,KS))
    relaxP.addConstrs(wS[i,t,KS]==wS[i,t-1,KS]+zS[i,t-1,KS-1]-quicksum(yS[q,t-1-tau[q,i],ij]+yL_h[q,t-1-tau[q,i],ij] for q in path for ij in AS[q].select(i,'*')) for i in S for t in range(1,T))

    relaxP.addConstrs(wL[i,0,k]==0 for i in S for k in range(KL))
    relaxP.addConstrs(wL[i,t,0]==wL[i,t-1,0]-zL[i,t-1,0]+quicksum(yL[q,t-1-tau[q,i],ij]+yS_h[q,t-1-tau[q,i],ij] for q in path for ij in AL[q].select('*',i)) for i in S for t in range(1,T))
    relaxP.addConstrs(wL[i,t,n]==wL[i,t-1,n]+zL[i,t-1,n-1]-zL[i,t-1,n] for i in S for t in range(1,T) for n in range(1,KL))
    relaxP.addConstrs(wL[i,t,KL]==wL[i,t-1,KL]+zL[i,t-1,KL-1]-quicksum(yL[q,t-1-tau[q,i],ij]+yS_h[q,t-1-tau[q,i],ij] for q in path for ij in AL[q].select(i,'*')) for i in S for t in range(1,T))

    relaxP.addConstrs(zS[i,t,n] <= wS[i,t,n] for i in S for t in range(T) for n in range(KS))
    relaxP.addConstrs(zL[i,t,n] <= wL[i,t,n] for i in S for t in range(T) for n in range(KL))

    # Add objective
    cost_electric = LinExpr(beta*quicksum(E[t]*zS[i,t,n] for i in S for t in range(T) for n in range(KS)) + beta*quicksum(E[t]*zL[i,t,n] for i in S for t in range(T) for n in range(KL)))
    penalty_unsatisfied = LinExpr(alpha*quicksum(fS[scenario_i,q,t]-quicksum(yS[q,t,ij] for ij in AS[q].select(source[q],'*'))+fL[scenario_i,q,t]-quicksum(yL[q,t,ij] for ij in AL[q].select(source[q],'*')) for t in range(T-max_tau) for q in path))
    income_swap = LinExpr(sum(SwapPriceS[t + tau[q, ij[0]]] * yS[q, t, ij] + (SwapPriceS[t + tau[q, ij[0]]]-gamma)*yL_h[q, t, ij] for t in range(T-max_tau) for q in path for ij in AS[q].select(path_data[q]['stations'],'*'))
                            +sum(SwapPriceL[t + tau[q, ij[0]]] * yL[q, t, ij] + SwapPriceS[t + tau[q, ij[0]]] * yS_h[q, t, ij] for t in range(T-max_tau) for q in path for ij in AL[q].select(path_data[q]['stations'],'*')))
    relaxP.setObjective(cost_electric-income_swap+penalty_unsatisfied)

    relaxP.update()
    return relaxP#, yS, yL, yS_h, yL_h, wS, wL, zS, zL

def BuildMW(scenario_i,xS,xL,xS0,xL0,ySP):
    relaxP = Model('relaxP')

    k=0
    XS={}
    XL={}
    XS0={}
    XL0={}
    for s in S:
        XS[s]=xS[k]
        XL[s]=xL[k]
        XS0[s]=xS0[k]
        XL0[s]=xL0[k]
        k+=1

    # Add variables (continuous)
    yS=tupledict()
    yL=tupledict()
    yS_h=tupledict()
    yL_h=tupledict()
    for p in path:
        for t in range(-max_tau,T):
            for ij in AS[p]:
                yS[p,t,ij] = relaxP.addVar(lb=0,vtype=GRB.CONTINUOUS,name='yS')
                yL_h[p,t,ij] = relaxP.addVar(lb=0,vtype=GRB.CONTINUOUS,name='yL')
            for ij in AL[p]:
                yL[p,t,ij] = relaxP.addVar(lb=0,vtype=GRB.CONTINUOUS,name='yL')
                yS_h[p,t,ij] = relaxP.addVar(lb=0,vtype=GRB.CONTINUOUS,name='yS')
        '''for t in range(-max_tau,0,1):
            for ij in AS[p]:
                yS[p,t,ij] = 0
                yL_h[p,t,ij] = 0
            for ij in AL[p]:
                yL[p,t,ij] = 0
                yS_h[p,t,ij] = 0'''
        

    wS = relaxP.addVars(S,T,KS+1,lb=0,vtype=GRB.CONTINUOUS,name='wS')
    wL = relaxP.addVars(S,T,KL+1,lb=0,vtype=GRB.CONTINUOUS,name='wL')
    zS = relaxP.addVars(S,T,KS,lb=0,vtype=GRB.CONTINUOUS,name='zS')
    zL = relaxP.addVars(S,T,KL,lb=0,vtype=GRB.CONTINUOUS,name='zL')

    xi = relaxP.addVar(vtype=GRB.CONTINUOUS,name='xi')
    relaxP._xi = xi

    # Add contraints

    relaxP.addConstrs(yS[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AS[p])
    relaxP.addConstrs(yL_h[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AS[p])
    relaxP.addConstrs(yL[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AL[p])
    relaxP.addConstrs(yS_h[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AL[p])

    relaxP.addConstrs(quicksum(yS[q,t,ij] for ij in AS[q].select(source[q],'*')) <= fS[scenario_i,q,t]-fS[scenario_i,q,t]*xi for q in path for t in range(T))
    relaxP.addConstrs(quicksum(yL[q,t,ij] for ij in AL[q].select(source[q],'*')) <= fL[scenario_i,q,t]-fS[scenario_i,q,t]*xi for q in path for t in range(T))

    relaxP.addConstrs(quicksum(yS_h[q,t,ij] for ij in AL[q].select(source[q],'*')) == 0 for q in path for t in range(T))
    relaxP.addConstrs(quicksum(yL_h[q,t,ij] for ij in AS[q].select(source[q],'*')) == 0 for q in path for t in range(T))

    relaxP.addConstrs(quicksum(yS[q,t,ij] for ij in AS[q].select(i,'*')) + quicksum(yS_h[q,t,ij] for ij in AL[q].select(i,'*')) 
             == quicksum(yS[q,t,ij] for ij in AS[q].select('*',i)) + quicksum(yS_h[q,t,ij] for ij in AL[q].select('*',i)) for q in path for t in range(T) for i in Sq[q])
    relaxP.addConstrs(quicksum(yL[q,t,ij] for ij in AL[q].select(i,'*')) + quicksum(yL_h[q,t,ij] for ij in AS[q].select(i,'*')) 
             == quicksum(yL[q,t,ij] for ij in AL[q].select('*',i)) + quicksum(yL_h[q,t,ij] for ij in AS[q].select('*',i)) for q in path for t in range(T) for i in Sq[q])

    XconstrS=relaxP.addConstrs((wS[i,0,KS] == XS[i]*xi-XS0[i] for i in S), name = 'Xconstr')
    XconstrL=relaxP.addConstrs((wL[i,0,KL] == XL[i]*xi-XL0[i] for i in S), name = 'Xconstr')

    relaxP._XconstrS = XconstrS
    relaxP._XconstrL = XconstrL

    relaxP.addConstrs(wS[i,T-1,k] == 0 for i in S for k in range(KS))
    relaxP.addConstrs(wL[i,T-1,k] == 0 for i in S for k in range(KL))

    relaxP.addConstrs(wS[i,0,k]==0 for i in S for k in range(KS))
    relaxP.addConstrs(wS[i,t,0]==wS[i,t-1,0]-zS[i,t-1,0]+quicksum(yS[q,t-1-tau[q,i], ij]+yL_h[q,t-1-tau[q,i], ij] for q in path for ij in AS[q].select('*',i)) for i in S for t in range(1,T))
    relaxP.addConstrs(wS[i,t,n]==wS[i,t-1,n]+zS[i,t-1,n-1]-zS[i,t-1,n] for i in S for t in range(1,T) for n in range(1,KS))
    relaxP.addConstrs(wS[i,t,KS]==wS[i,t-1,KS]+zS[i,t-1,KS-1]-quicksum(yS[q,t-1-tau[q,i],ij]+yL_h[q,t-1-tau[q,i],ij] for q in path for ij in AS[q].select(i,'*')) for i in S for t in range(1,T))

    relaxP.addConstrs(wL[i,0,k]==0 for i in S for k in range(KL))
    relaxP.addConstrs(wL[i,t,0]==wL[i,t-1,0]-zL[i,t-1,0]+quicksum(yL[q,t-1-tau[q,i],ij]+yS_h[q,t-1-tau[q,i],ij] for q in path for ij in AL[q].select('*',i)) for i in S for t in range(1,T))
    relaxP.addConstrs(wL[i,t,n]==wL[i,t-1,n]+zL[i,t-1,n-1]-zL[i,t-1,n] for i in S for t in range(1,T) for n in range(1,KL))
    relaxP.addConstrs(wL[i,t,KL]==wL[i,t-1,KL]+zL[i,t-1,KL-1]-quicksum(yL[q,t-1-tau[q,i],ij]+yS_h[q,t-1-tau[q,i],ij] for q in path for ij in AL[q].select(i,'*')) for i in S for t in range(1,T))

    relaxP.addConstrs(zS[i,t,n] <= wS[i,t,n] for i in S for t in range(T) for n in range(KS))
    relaxP.addConstrs(zL[i,t,n] <= wL[i,t,n] for i in S for t in range(T) for n in range(KL))

    # Add objective
    cost_electric = LinExpr(beta*quicksum(E[t]*zS[i,t,n] for i in S for t in range(T) for n in range(KS)) + beta*quicksum(E[t]*zL[i,t,n] for i in S for t in range(T) for n in range(KL)))
    penalty_unsatisfied = LinExpr(alpha*quicksum(fS[scenario_i,q,t]-quicksum(yS[q,t,ij] for ij in AS[q].select(source[q],'*'))+fL[scenario_i,q,t]-quicksum(yL[q,t,ij] for ij in AL[q].select(source[q],'*')) for t in range(T-max_tau) for q in path))
    income_swap = LinExpr(sum(SwapPriceS[t + tau[q, ij[0]]] * yS[q, t, ij] + (SwapPriceS[t + tau[q, ij[0]]]-gamma)*yL_h[q, t, ij] for t in range(T-max_tau) for q in path for ij in AS[q].select(path_data[q]['stations'],'*'))
                            +sum(SwapPriceL[t + tau[q, ij[0]]] * yL[q, t, ij] + SwapPriceS[t + tau[q, ij[0]]] * yS_h[q, t, ij] for t in range(T-max_tau) for q in path for ij in AL[q].select(path_data[q]['stations'],'*')))
    relaxP.setObjective(cost_electric-income_swap + ySP*xi+penalty_unsatisfied)

    relaxP.update()
    return relaxP

# 整数子问题
def BuildSubP(scenario_i,xS,xL):
    subP = Model('subP')

    k=0
    XS={}
    XL={}
    for s in S:
        XS[s]=xS[k]
        XL[s]=xL[k]
        k+=1

    # Add variables
    yS=tupledict()
    yL=tupledict()
    yS_h=tupledict()
    yL_h=tupledict()
    for p in path:
        for t in range(-max_tau,T):
            for ij in AS[p]:
                yS[p,t,ij] = subP.addVar(lb=0,vtype=GRB.INTEGER)
                yL_h[p,t,ij] = subP.addVar(lb=0,vtype=GRB.INTEGER)
            for ij in AL[p]:
                yL[p,t,ij] = subP.addVar(lb=0,vtype=GRB.INTEGER)
                yS_h[p,t,ij] = subP.addVar(lb=0,vtype=GRB.INTEGER)
        '''for t in range(-max_tau,0,1):
            for ij in AS[p]:
                yS[p,t,ij] = 0
                yL_h[p,t,ij] = 0
            for ij in AL[p]:
                yL[p,t,ij] = 0
                yS_h[p,t,ij] = 0'''

    wS = subP.addVars(S,T,KS+1,lb=0,vtype=GRB.INTEGER)
    wL = subP.addVars(S,T,KL+1,lb=0,vtype=GRB.INTEGER)
    zS = subP.addVars(S,T,KS,lb=0,vtype=GRB.INTEGER)
    zL = subP.addVars(S,T,KL,lb=0,vtype=GRB.INTEGER)

    subP._yS = yS
    subP._yL = yL
    subP._yS_h = yS_h
    subP._yL_h = yL_h
    subP._wS = wS
    subP._wL = wL
    subP._zS = zS
    subP._zL = zL

    # Add contraints
    subP.addConstrs(yS[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AS[p])
    subP.addConstrs(yL_h[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AS[p])
    subP.addConstrs(yL[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AL[p])
    subP.addConstrs(yS_h[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AL[p])

    subP.addConstrs(quicksum(yS[q,t,ij] for ij in AS[q].select(source[q],'*')) <= fS[scenario_i,q,t] for q in path for t in range(T))
    subP.addConstrs(quicksum(yL[q,t,ij] for ij in AL[q].select(source[q],'*')) <= fL[scenario_i,q,t] for q in path for t in range(T))

    subP.addConstrs(quicksum(yS_h[q,t,ij] for ij in AL[q].select(source[q],'*')) == 0 for q in path for t in range(T))
    subP.addConstrs(quicksum(yL_h[q,t,ij] for ij in AS[q].select(source[q],'*')) == 0 for q in path for t in range(T))

    subP.addConstrs(quicksum(yS[q,t,ij] for ij in AS[q].select(i,'*')) + quicksum(yS_h[q,t,ij] for ij in AL[q].select(i,'*')) 
             == quicksum(yS[q,t,ij] for ij in AS[q].select('*',i)) + quicksum(yS_h[q,t,ij] for ij in AL[q].select('*',i)) for q in path for t in range(T) for i in Sq[q])
    subP.addConstrs(quicksum(yL[q,t,ij] for ij in AL[q].select(i,'*')) + quicksum(yL_h[q,t,ij] for ij in AS[q].select(i,'*')) 
             == quicksum(yL[q,t,ij] for ij in AL[q].select('*',i)) + quicksum(yL_h[q,t,ij] for ij in AS[q].select('*',i)) for q in path for t in range(T) for i in Sq[q])

    subP.addConstrs((wS[i,0,KS] == XS[i] for i in S), name = 'Xconstr')
    subP.addConstrs((wL[i,0,KL] == XL[i] for i in S), name = 'Xconstr')

    subP.addConstrs(wS[i,T-1,k] == 0 for i in S for k in range(KS))
    subP.addConstrs(wL[i,T-1,k] == 0 for i in S for k in range(KL))

    subP.addConstrs(wS[i,0,k]==0 for i in S for k in range(KS))
    subP.addConstrs(wS[i,t,0]==wS[i,t-1,0]-zS[i,t-1,0]+quicksum(yS[q,t-1-tau[q,i], ij]+yL_h[q,t-1-tau[q,i], ij] for q in path for ij in AS[q].select('*',i)) for i in S for t in range(1,T))
    subP.addConstrs(wS[i,t,n]==wS[i,t-1,n]+zS[i,t-1,n-1]-zS[i,t-1,n] for i in S for t in range(1,T) for n in range(1,KS))
    subP.addConstrs(wS[i,t,KS]==wS[i,t-1,KS]+zS[i,t-1,KS-1]-quicksum(yS[q,t-1-tau[q,i],ij]+yL_h[q,t-1-tau[q,i],ij] for q in path for ij in AS[q].select(i,'*')) for i in S for t in range(1,T))

    subP.addConstrs(wL[i,0,k]==0 for i in S for k in range(KL))
    subP.addConstrs(wL[i,t,0]==wL[i,t-1,0]-zL[i,t-1,0]+quicksum(yL[q,t-1-tau[q,i],ij]+yS_h[q,t-1-tau[q,i],ij] for q in path for ij in AL[q].select('*',i)) for i in S for t in range(1,T))
    subP.addConstrs(wL[i,t,n]==wL[i,t-1,n]+zL[i,t-1,n-1]-zL[i,t-1,n] for i in S for t in range(1,T) for n in range(1,KL))
    subP.addConstrs(wL[i,t,KL]==wL[i,t-1,KL]+zL[i,t-1,KL-1]-quicksum(yL[q,t-1-tau[q,i],ij]+yS_h[q,t-1-tau[q,i],ij] for q in path for ij in AL[q].select(i,'*')) for i in S for t in range(1,T))

    subP.addConstrs(zS[i,t,n] <= wS[i,t,n] for i in S for t in range(T) for n in range(KS))
    subP.addConstrs(zL[i,t,n] <= wL[i,t,n] for i in S for t in range(T) for n in range(KL))

    # Add objective
    cost_electric = LinExpr(beta*quicksum(E[t]*zS[i,t,n] for i in S for t in range(T) for n in range(KS)) + beta*quicksum(E[t]*zL[i,t,n] for i in S for t in range(T) for n in range(KL)))
    penalty_unsatisfied = LinExpr(alpha*quicksum(fS[scenario_i,q,t]-quicksum(yS[q,t,ij] for ij in AS[q].select(source[q],'*'))+fL[scenario_i,q,t]-quicksum(yL[q,t,ij] for ij in AL[q].select(source[q],'*')) for t in range(T-max_tau) for q in path))
    income_swap = LinExpr(sum(SwapPriceS[t + tau[q, ij[0]]] * yS[q, t, ij] + (SwapPriceS[t + tau[q, ij[0]]]-gamma)*yL_h[q, t, ij] for t in range(T-max_tau) for q in path for ij in AS[q].select(path_data[q]['stations'],'*'))
                            +sum(SwapPriceL[t + tau[q, ij[0]]] * yL[q, t, ij] + SwapPriceS[t + tau[q, ij[0]]] * yS_h[q, t, ij] for t in range(T-max_tau) for q in path for ij in AL[q].select(path_data[q]['stations'],'*')))
    subP.setObjective(cost_electric-income_swap+penalty_unsatisfied)

    subP.update()
    return subP#, yS, yL, yS_h, yL_h, wS, wL, zS, zL


def BuildBaseP(scenario_i,xS,xL):
    subP = Model('subP')

    k=0
    XS={}
    XL={}
    for s in S:
        XS[s]=xS[k]
        XL[s]=xL[k]
        k+=1

    # Add variables
    yS=tupledict()
    yL=tupledict()
    yS_h=tupledict()
    yL_h=tupledict()
    for p in path:
        for t in range(-max_tau,T):
            for ij in AS[p]:
                yS[p,t,ij] = subP.addVar(lb=0,vtype=GRB.INTEGER)
                yL_h[p,t,ij] = subP.addVar(lb=0,vtype=GRB.INTEGER)
            for ij in AL[p]:
                yL[p,t,ij] = subP.addVar(lb=0,vtype=GRB.INTEGER)
                yS_h[p,t,ij] = subP.addVar(lb=0,vtype=GRB.INTEGER)
        '''for t in range(-max_tau,0,1):
            for ij in AS[p]:
                yS[p,t,ij] = 0
                yL_h[p,t,ij] = 0
            for ij in AL[p]:
                yL[p,t,ij] = 0
                yS_h[p,t,ij] = 0'''

    wS = subP.addVars(S,T,KS+1,lb=0,vtype=GRB.INTEGER)
    wL = subP.addVars(S,T,KL+1,lb=0,vtype=GRB.INTEGER)
    zS = subP.addVars(S,T,KS,lb=0,vtype=GRB.INTEGER)
    zL = subP.addVars(S,T,KL,lb=0,vtype=GRB.INTEGER)

    # Add contraints
    subP.addConstrs(yS[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AS[p])
    subP.addConstrs(yL_h[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AS[p])
    subP.addConstrs(yL[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AL[p])
    subP.addConstrs(yS_h[p,t,ij] == 0 for p in path for t in range(-max_tau,0) for ij in AL[p])

    subP.addConstrs(quicksum(yS[q,t,ij] for ij in AS[q].select(source[q],'*')) <= fS[scenario_i,q,t] for q in path for t in range(T))
    subP.addConstrs(quicksum(yL[q,t,ij] for ij in AL[q].select(source[q],'*')) <= fL[scenario_i,q,t] for q in path for t in range(T))

    subP.addConstrs(quicksum(yS_h[q,t,ij] for ij in AL[q].select(source[q],'*')) == 0 for q in path for t in range(T))
    subP.addConstrs(quicksum(yL_h[q,t,ij] for ij in AS[q].select(source[q],'*')) == 0 for q in path for t in range(T))

    subP.addConstrs(quicksum(yS[q,t,ij] for ij in AS[q].select(i,'*')) + quicksum(yS_h[q,t,ij] for ij in AL[q].select(i,'*')) 
             == quicksum(yS[q,t,ij] for ij in AS[q].select('*',i)) + quicksum(yS_h[q,t,ij] for ij in AL[q].select('*',i)) for q in path for t in range(T) for i in Sq[q])
    subP.addConstrs(quicksum(yL[q,t,ij] for ij in AL[q].select(i,'*')) + quicksum(yL_h[q,t,ij] for ij in AS[q].select(i,'*')) 
             == quicksum(yL[q,t,ij] for ij in AL[q].select('*',i)) + quicksum(yL_h[q,t,ij] for ij in AS[q].select('*',i)) for q in path for t in range(T) for i in Sq[q])

    subP.addConstrs((wS[i,0,KS] == XS[i] for i in S), name = 'Xconstr')
    subP.addConstrs((wL[i,0,KL] == XL[i] for i in S), name = 'Xconstr')

    subP.addConstrs(wS[i,T-1,k] == 0 for i in S for k in range(KS))
    subP.addConstrs(wL[i,T-1,k] == 0 for i in S for k in range(KL))

    subP.addConstrs(wS[i,0,k]==0 for i in S for k in range(KS))
    subP.addConstrs(wS[i,t,0]==wS[i,t-1,0]-zS[i,t-1,0]+quicksum(yS[q,t-1-tau[q,i], ij]+yL_h[q,t-1-tau[q,i], ij] for q in path for ij in AS[q].select('*',i)) for i in S for t in range(1,T))
    subP.addConstrs(wS[i,t,n]==wS[i,t-1,n]+zS[i,t-1,n-1]-zS[i,t-1,n] for i in S for t in range(1,T) for n in range(1,KS))
    subP.addConstrs(wS[i,t,KS]==wS[i,t-1,KS]+zS[i,t-1,KS-1]-quicksum(yS[q,t-1-tau[q,i],ij]+yL_h[q,t-1-tau[q,i],ij] for q in path for ij in AS[q].select(i,'*')) for i in S for t in range(1,T))

    subP.addConstrs(wL[i,0,k]==0 for i in S for k in range(KL))
    subP.addConstrs(wL[i,t,0]==wL[i,t-1,0]-zL[i,t-1,0]+quicksum(yL[q,t-1-tau[q,i],ij]+yS_h[q,t-1-tau[q,i],ij] for q in path for ij in AL[q].select('*',i)) for i in S for t in range(1,T))
    subP.addConstrs(wL[i,t,n]==wL[i,t-1,n]+zL[i,t-1,n-1]-zL[i,t-1,n] for i in S for t in range(1,T) for n in range(1,KL))
    subP.addConstrs(wL[i,t,KL]==wL[i,t-1,KL]+zL[i,t-1,KL-1]-quicksum(yL[q,t-1-tau[q,i],ij]+yS_h[q,t-1-tau[q,i],ij] for q in path for ij in AL[q].select(i,'*')) for i in S for t in range(1,T))

    subP.addConstrs(zS[i,t,n] <= wS[i,t,n] for i in S for t in range(T) for n in range(KS))
    subP.addConstrs(zL[i,t,n] <= wL[i,t,n] for i in S for t in range(T) for n in range(KL))

    # Add objective
    cost_electric = LinExpr(beta*quicksum(E[t]*zS[i,t,n] for i in S for t in range(T) for n in range(KS)) + beta*quicksum(E[t]*zL[i,t,n] for i in S for t in range(T) for n in range(KL)))
    penalty_unsatisfied = LinExpr(alpha*quicksum(fS[scenario_i,q,t]-quicksum(yS[q,t,ij] for ij in AS[q].select(source[q],'*'))+fL[scenario_i,q,t]-quicksum(yL[q,t,ij] for ij in AL[q].select(source[q],'*')) for t in range(T-max_tau) for q in path))
    income_swap = LinExpr(sum(SwapPriceS[t + tau[q, ij[0]]] * yS[q, t, ij] + (SwapPriceS[t + tau[q, ij[0]]]-gamma)*yL_h[q, t, ij] for t in range(T-max_tau) for q in path for ij in AS[q].select(path_data[q]['stations'],'*'))
                            +sum(SwapPriceL[t + tau[q, ij[0]]] * yL[q, t, ij] + SwapPriceS[t + tau[q, ij[0]]] * yS_h[q, t, ij] for t in range(T-max_tau) for q in path for ij in AL[q].select(path_data[q]['stations'],'*')))
    subP.setObjective(cost_electric-income_swap+penalty_unsatisfied)

    subP.update()
    return subP#, yS, yL, yS_h, yL_h, wS, wL, zS, zL


def BuildDEPgpt(num_scenario,prob):
    dep = Model('DEP')

    # Add Variables for master problem
    XS = dep.addVars(S, vtype=GRB.INTEGER, lb=0,ub=I, name='XS')
    XL = dep.addVars(S, vtype=GRB.INTEGER, lb=0,ub=I, name='XL')
    PXS = dep.addVars(S,vtype=GRB.INTEGER,lb=0,ub=I,name='PXS')
    PXL = dep.addVars(S,vtype=GRB.INTEGER,lb=0,ub=I,name='PXL')
    dep._XS = XS
    dep._XL = XL
    dep._PXS = PXS
    dep._PXL = PXL

    dep.addConstrs(XS[i]+XL[i]<=I for i in S)
    dep.addConstrs(PXS[i] >= XS[i] - XSbar[i] for i in S)
    dep.addConstrs(PXL[i] >= XL[i] - XLbar[i] for i in S)

    # Add Variables for each scenario in the subproblem
    yS = {}
    yL = {}
    yS_h = {}
    yL_h = {}
    wS = {}
    wL = {}
    zS = {}
    zL = {}

    for scenario_i in range(num_scenario):
        yS[scenario_i] = tupledict()
        yL[scenario_i] = tupledict()
        yS_h[scenario_i] = tupledict()
        yL_h[scenario_i] = tupledict()

        for p in path:
            for t in range(-max_tau, T):
                for ij in AS[p]:
                    yS[scenario_i][p, t, ij] = dep.addVar(lb=0, vtype=GRB.INTEGER)
                    yL_h[scenario_i][p, t, ij] = dep.addVar(lb=0, vtype=GRB.INTEGER)
                for ij in AL[p]:
                    yL[scenario_i][p, t, ij] = dep.addVar(lb=0, vtype=GRB.INTEGER)
                    yS_h[scenario_i][p, t, ij] = dep.addVar(lb=0, vtype=GRB.INTEGER)

        wS[scenario_i] = dep.addVars(S, T, KS+1, lb=0, vtype=GRB.INTEGER)
        wL[scenario_i] = dep.addVars(S, T, KL+1, lb=0, vtype=GRB.INTEGER)
        zS[scenario_i] = dep.addVars(S, T, KS, lb=0, vtype=GRB.INTEGER)
        zL[scenario_i] = dep.addVars(S, T, KL, lb=0, vtype=GRB.INTEGER)

        # Add constraints for each scenario
        dep.addConstrs(yS[scenario_i][p, t, ij] == 0 for p in path for t in range(-max_tau, 0) for ij in AS[p])
        dep.addConstrs(yL_h[scenario_i][p, t, ij] == 0 for p in path for t in range(-max_tau, 0) for ij in AS[p])
        dep.addConstrs(yL[scenario_i][p, t, ij] == 0 for p in path for t in range(-max_tau, 0) for ij in AL[p])
        dep.addConstrs(yS_h[scenario_i][p, t, ij] == 0 for p in path for t in range(-max_tau, 0) for ij in AL[p])

        dep.addConstrs(quicksum(yS[scenario_i][q, t, ij] for ij in AS[q].select(source[q], '*')) <= fS[scenario_i, q, t] for q in path for t in range(T))
        dep.addConstrs(quicksum(yL[scenario_i][q, t, ij] for ij in AL[q].select(source[q], '*')) <= fL[scenario_i, q, t] for q in path for t in range(T))

        dep.addConstrs(quicksum(yS_h[scenario_i][q, t, ij] for ij in AL[q].select(source[q], '*')) == 0 for q in path for t in range(T))
        dep.addConstrs(quicksum(yL_h[scenario_i][q, t, ij] for ij in AS[q].select(source[q], '*')) == 0 for q in path for t in range(T))

        dep.addConstrs(quicksum(yS[scenario_i][q, t, ij] for ij in AS[q].select(i, '*')) + quicksum(yS_h[scenario_i][q, t, ij] for ij in AL[q].select(i, '*')) 
                == quicksum(yS[scenario_i][q, t, ij] for ij in AS[q].select('*', i)) + quicksum(yS_h[scenario_i][q, t, ij] for ij in AL[q].select('*', i)) for q in path for t in range(T) for i in Sq[q])
        dep.addConstrs(quicksum(yL[scenario_i][q, t, ij] for ij in AL[q].select(i, '*')) + quicksum(yL_h[scenario_i][q, t, ij] for ij in AS[q].select(i, '*')) 
                == quicksum(yL[scenario_i][q, t, ij] for ij in AL[q].select('*', i)) + quicksum(yL_h[scenario_i][q, t, ij] for ij in AS[q].select('*', i)) for q in path for t in range(T) for i in Sq[q])

        dep.addConstrs((wS[scenario_i][i, 0, KS] == XS[i] for i in S), name=f'Xconstr_{scenario_i}')
        dep.addConstrs((wL[scenario_i][i, 0, KL] == XL[i] for i in S), name=f'Xconstr_{scenario_i}')

        dep.addConstrs(wS[scenario_i][i,T-1,k] == 0 for i in S for k in range(KS))
        dep.addConstrs(wL[scenario_i][i,T-1,k] == 0 for i in S for k in range(KL))

        dep.addConstrs(wS[scenario_i][i, 0, k] == 0 for i in S for k in range(KS))
        dep.addConstrs(wS[scenario_i][i, t, 0] == wS[scenario_i][i, t-1, 0] - zS[scenario_i][i, t-1, 0] + quicksum(yS[scenario_i][q, t-1-tau[q, i], ij] + yL_h[scenario_i][q, t-1-tau[q, i], ij] for q in path for ij in AS[q].select('*', i)) for i in S for t in range(1, T))
        dep.addConstrs(wS[scenario_i][i, t, n] == wS[scenario_i][i, t-1, n] + zS[scenario_i][i, t-1, n-1] - zS[scenario_i][i, t-1, n] for i in S for t in range(1, T) for n in range(1, KS))
        dep.addConstrs(wS[scenario_i][i, t, KS] == wS[scenario_i][i, t-1, KS] + zS[scenario_i][i, t-1, KS-1] - quicksum(yS[scenario_i][q, t-1-tau[q, i], ij] + yL_h[scenario_i][q, t-1-tau[q, i], ij] for q in path for ij in AS[q].select(i, '*')) for i in S for t in range(1, T))

        dep.addConstrs(wL[scenario_i][i, 0, k] == 0 for i in S for k in range(KL))
        dep.addConstrs(wL[scenario_i][i, t, 0] == wL[scenario_i][i, t-1, 0] - zL[scenario_i][i, t-1, 0] + quicksum(yL[scenario_i][q, t-1-tau[q, i], ij] + yS_h[scenario_i][q, t-1-tau[q, i], ij] for q in path for ij in AL[q].select('*', i)) for i in S for t in range(1, T))
        dep.addConstrs(wL[scenario_i][i, t, n] == wL[scenario_i][i, t-1, n] + zL[scenario_i][i, t-1, n-1] - zL[scenario_i][i, t-1, n] for i in S for t in range(1, T) for n in range(1, KL))
        dep.addConstrs(wL[scenario_i][i, t, KL] == wL[scenario_i][i, t-1, KL] + zL[scenario_i][i, t-1, KL-1] - quicksum(yL[scenario_i][q, t-1-tau[q, i], ij] + yS_h[scenario_i][q, t-1-tau[q, i], ij] for q in path for ij in AL[q].select(i, '*')) for i in S for t in range(1, T))

        dep.addConstrs(zS[scenario_i][i, t, n] <= wS[scenario_i][i, t, n] for i in S for t in range(T) for n in range(KS))
        dep.addConstrs(zL[scenario_i][i, t, n] <= wL[scenario_i][i, t, n] for i in S for t in range(T) for n in range(KL))

    dep._yS = yS
    dep._yL = yL
    dep._yS_h = yS_h
    dep._yL_h = yL_h
    dep._wS = wS
    dep._wL = wL
    dep._zS = zS
    dep._zL = zL

    # Add objective function
    cost_electric = quicksum(prob[scenario_i]*(beta * quicksum(E[t] * zS[scenario_i][i, t, n] for i in S for t in range(T) for n in range(KS)) + beta * quicksum(E[t] * zL[scenario_i][i, t, n] for i in S for t in range(T) for n in range(KL))) for scenario_i in range(num_scenario))
    penalty_unsatisfied = quicksum(prob[scenario_i]*(alpha * quicksum(fS[scenario_i, q, t] - quicksum(yS[scenario_i][q, t, ij] for ij in AS[q].select(source[q], '*')) + fL[scenario_i, q, t] - quicksum(yL[scenario_i][q, t, ij] for ij in AL[q].select(source[q], '*')) for t in range(T-max_tau) for q in path)) for scenario_i in range(num_scenario))
    income_swap = quicksum(prob[scenario_i]*(
        sum(SwapPriceS[t + tau[q, ij[0]]] * yS[scenario_i][q, t, ij] + (SwapPriceS[t + tau[q, ij[0]]]-gamma)* yL_h[scenario_i][q, t, ij] for t in range(T-max_tau) for q in path for ij in AS[q].select(path_data[q]['stations'], '*'))
        + sum(SwapPriceL[t + tau[q, ij[0]]] * yL[scenario_i][q, t, ij] + SwapPriceS[t + tau[q, ij[0]]] * yS_h[scenario_i][q, t, ij] for t in range(T-max_tau) for q in path for ij in AL[q].select(path_data[q]['stations'], '*')))
        for scenario_i in range(num_scenario))

    dep.setObjective(cS * quicksum(PXS[i] for i in S) + cL * quicksum(PXL[i] for i in S) + cost_electric - income_swap + penalty_unsatisfied)

    dep.update()

    dep.setParam(GRB.Param.TimeLimit, 7200)
    return dep



def main():
    print('Begin Solving')
    #ucnbc sol:
    XS=[ 6., -0., 21., 19., 17., 21., 17., 15.,  1.,  6., -0.]
    XL=[-0., 21.,  0.,  2.,  4.,  0.,  4.,  6., 20., -0., -0.]

    #subP,yS, yL, yS_h, yL_h, wS, wL, zS, zL=BuildBaseP(1,XS,XL)
    #start_time = time.time()
    #testP=BuildSubP(0,XS,XL)
    #testP.update()
    #build_time = time.time()-start_time

    # 设置启发式参数
    '''subP.setParam('Heuristics', 1.0)  # 启发式的权重设置为最大
    subP.setParam('Cuts', 0)          # 禁用切割平面
    subP.setParam('Presolve', 0)      # 禁用预求解
    subP.setParam('MIPFocus', 1)      # 优先寻找可行解
    subP.setParam('Symmetry', 0)      # 禁用对称性检测
    subP.setParam('VarBranch', 0)     # 禁用变量分支
    subP.setParam('NodeLimit', 1)     # 只探索一个节点
    subP.setParam('SolutionLimit', 1)     # 只寻找一个可行解
    subP.setParam('DualReductions', 0)    # 禁用对偶约简'''
    model = BuildRelaxP(0,XS,XL)
    print('finish build')
    model.setParam('Method', 2)  # 对偶单纯形法
    model.setParam('Presolve', 2)  # 激进预处理
    #model.setParam('Threads', 4)  # 4 个线程
    model.setParam('OptimalityTol', 1e-8)
    #model.setParam('IterationLimit', 100000)
    #model.setParam('TimeLimit', 3600)  # 1 小时
    model.setParam('FeasibilityTol', 1e-5)
    model.setParam('ScaleFlag', 1)  # 保守缩放
    model.setParam('DualReductions', 0)  # 关闭对偶缩减
    model.setParam('OutputFlag', 1)  # 开启日志
    #model.setParam('Heuristics', 0)
    #model.setParam('MIPGap', 0.0005)  # 设置允许的最小相对MIP Gap为5%
    model.optimize()
    for n,station in enumerate(stations):
        model._XconstrS[station].RHS = XS1[n]
        model._XconstrL[station].RHS = XL1[n]
    model.reset()
    model.optimize()
    #solve_time = time.time()-start_time


                            

if __name__ == "__main__":
    '''num_scenario = 10
    prob = (0.1,0.1)*5'''
    model = BuildDEPgpt(num_scenario,prob)
    model.optimize()
    if model.status == GRB.OPTIMAL:
        # 提取最优解中的各部分
        optimal_XS = [model._XS[i].x for i in S]
        optimal_XL = [model._XL[i].x for i in S]

        optimal_PXS = sum(model._PXS[i].x for i in S)
        optimal_PXL = sum(model._PXL[i].x for i in S)

        # 计算电能消耗成本 (cost_electric)
        optimal_cost_electric = sum(
            prob[scenario_i] * (
                beta * sum(E[t] * model._zS[scenario_i][i, t, n].x for i in S for t in range(T) for n in range(KS)) +
                beta * sum(E[t] * model._zL[scenario_i][i, t, n].x for i in S for t in range(T) for n in range(KL))
            )
            for scenario_i in range(num_scenario)
        )

        # 计算未满足需求的惩罚 (penalty_unsatisfied)
        optimal_penalty_unsatisfied = sum(
            prob[scenario_i] * (
                alpha * sum(
                    fS[scenario_i, q, t] - sum(model._yS[scenario_i][q, t, ij].x for ij in AS[q].select(source[q], '*')) +
                    fL[scenario_i, q, t] - sum(model._yL[scenario_i][q, t, ij].x for ij in AL[q].select(source[q], '*'))
                    for t in range(T-max_tau) for q in path
                )
            )
            for scenario_i in range(num_scenario)
        )

        # 计算电池交换的收入 (income_swap)
        optimal_income_swap = sum(
            prob[scenario_i] * (
                sum(SwapPriceS[t + tau[q, ij[0]]] * model._yS[scenario_i][q, t, ij].x + (SwapPriceS[t + tau[q, ij[0]]] - gamma) * model._yL_h[scenario_i][q, t, ij].x
                    for t in range(T-max_tau) for q in path for ij in AS[q].select(path_data[q]['stations'], '*')) +
                sum(SwapPriceL[t + tau[q, ij[0]]] * model._yL[scenario_i][q, t, ij].x + SwapPriceS[t + tau[q, ij[0]]] * model._yS_h[scenario_i][q, t, ij].x
                    for t in range(T-max_tau) for q in path for ij in AL[q].select(path_data[q]['stations'], '*'))
            )
            for scenario_i in range(num_scenario)
        )

        # 计算换电次数和充电次数
        num_swap_S = sum(
            prob[scenario_i] * (
                sum( model._yS[scenario_i][q, t, ij].x + model._yL_h[scenario_i][q, t, ij].x
                    for t in range(T-max_tau) for q in path for ij in AS[q].select(path_data[q]['stations'], '*'))
            )
            for scenario_i in range(num_scenario)
        )

        num_swap_L = sum(
            prob[scenario_i] * (
                sum(model._yL[scenario_i][q, t, ij].x + model._yS_h[scenario_i][q, t, ij].x
                    for t in range(T-max_tau) for q in path for ij in AL[q].select(path_data[q]['stations'], '*'))
            )
            for scenario_i in range(num_scenario)
        )

        num_charge_S = sum(
            prob[scenario_i] * (
                sum(model._zS[scenario_i][i, t, n].x for i in S for t in range(T) for n in range(KS))
            )
            for scenario_i in range(num_scenario)
        )

        num_charge_L = sum(
            prob[scenario_i] * (
                sum(model._zL[scenario_i][i, t, n].x for i in S for t in range(T) for n in range(KL))
            )
            for scenario_i in range(num_scenario)
        )

        # 计算目标函数的各部分
        total_cost_battery = cS * optimal_PXS + cL * optimal_PXL
        total_cost_electric = optimal_cost_electric
        total_penalty_unsatisfied = optimal_penalty_unsatisfied
        total_income_swap = optimal_income_swap

        # 输出目标函数的各部分结果
        print(f"Optimal solution:{optimal_XS,optimal_XL}")
        print(f"Total battery cost: {total_cost_battery}")
        print(f"Total electric cost: {total_cost_electric}")
        print(f"Total penalty for unsatisfied demand: {total_penalty_unsatisfied}")
        print(f"Total income from swap: {total_income_swap}")
        print(num_swap_S,num_swap_L,num_charge_S,num_charge_L)

        # 计算并输出总的目标函数值
        final_obj_value = total_cost_battery + total_cost_electric - total_income_swap + total_penalty_unsatisfied
        print(f"Final objective value: {final_obj_value}")

    
    '''XS = [6.0, 17.0, 18.0, 6.0, 8.0, 17.0, 11.0, 18.0, 13.0, 2.0, -0.0]
    XL = [0.0, 3.0, 0.0, 7.0, 5.0, 4.0, 6.0, 3.0, 2.0, -0.0, 5.0]
    y = cS*sum(XS)+cL*sum(XL)
    for i in range(num_scenario):
        model = BuildRelaxP(i,XS,XL)
        model.optimize()
        y_i = model.ObjVal
        y += prob[i]*y_i
    print(y)'''

    '''XS = [6.0, 17.0, 18.0, 6.0, 8.0, 17.0, 11.0, 18.0, 13.0, 2.0, -0.0]
    XL = [0.0, 3.0, 0.0, 7.0, 5.0, 4.0, 6.0, 3.0, 2.0, -0.0, 5.0]
    model = BuildRelaxP(2,XS,XL)
    model.optimize()'''