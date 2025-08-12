from input import *

#random.seed(2)

def Heuristic(scenario_i,XS,XL):
    
    # Input demand flow data, REMAINING DEMAND
    demandfS={(j, k): value for (i, j, k), value in fS.items() if i == scenario_i}
    demandfL={(j, k): value for (i, j, k), value in fL.items() if i == scenario_i}
    #print(sum(demandfL.values())+sum(demandfS.values()))
    
    # Construct Remain Battery Matrix
    RemainBatteryS = {}
    RemainBatteryL = {}

    for k,i in enumerate(stations):
        for t in range(TimePeriods+KL):
            RemainBatteryS[i,t] = XS[k]
            RemainBatteryL[i,t] = XL[k]

    # sort paths from longest to shortest
    length = {}
    for p in path_data:
        length[p] = dist[path_data[p]['source']][path_data[p]['root']]

    sorted_paths = sorted(length, key=length.get, reverse=True)
    #print(sorted_paths)

    # initial obj values
    SwapIncome = 0
    UnsatisfyPenalty = 0
    ChargingCost = 0

    #flowS={}
    #flowL={}
    for cur_t in range(TimePeriods-max_tau):
        for p in sorted_paths:
            #flowS[cur_t,p]=0
            #flowL[cur_t,p]=0
            for sn,strategyS in enumerate(path_data[p]['PATHS']):
                cur_max_flowS = demandfS[p,cur_t]
                for arc in range(len(strategyS) - 1):
                    station = strategyS[arc][1]
                    for k in range(KS):
                        cur_max_flowS=min(cur_max_flowS,RemainBatteryS[station,cur_t+tau[p,station]+k])
                #print('path',p,' :',max_flowS)
                #flowS[cur_t,p]+=cur_max_flowS

                demandfS[p,cur_t] -= cur_max_flowS
                for arc in range(len(strategyS) - 1):
                    station = strategyS[arc][1]
                    SwapIncome += SwapPriceS[cur_t+tau[p,station]]*cur_max_flowS
                    for k in range(KS):
                        RemainBatteryS[station,cur_t+tau[p,station]+k]-=cur_max_flowS

            for sn, strategyL in enumerate(path_data[p]['PATHL']):
                cur_max_flowL = demandfL[p,cur_t]
                for arc in range(len(strategyL) - 1):
                    station = strategyL[arc][1]
                    for k in range(KL):
                        cur_max_flowL = min(cur_max_flowL, RemainBatteryL[station,cur_t + tau[p,station] + k])
                #print('path', p, ' :', max_flowL)
                #flowL[cur_t,p] += cur_max_flowL

                demandfL[p,cur_t] -= cur_max_flowL
                for arc in range(len(strategyL) - 1):
                    station = strategyL[arc][1]
                    SwapIncome += SwapPriceL[cur_t+tau[p,station]]*cur_max_flowL
                    for k in range(KL):
                        RemainBatteryL[station,cur_t + tau[p,station] + k] -= cur_max_flowL

    # calculate obj
    for t in range(T-max_tau):
        for p in path_data:
            UnsatisfyPenalty += demandfS[p,t] + demandfL[p,t] 
    

    for t in range(TimePeriods-max_tau):
        for k,i in enumerate(stations):
            #if XS[k]-RemainBatteryS[i,t]<0 or XL[k]-RemainBatteryL[i,t]<0:
            #    raise KeyError('Wrong!,num charging>num battery')
            ChargingCost += E[t]*beta*(XS[k]-RemainBatteryS[i,t]+XL[k]-RemainBatteryL[i,t])

    print(ChargingCost,alpha*UnsatisfyPenalty,SwapIncome)

    return ChargingCost-SwapIncome+alpha*UnsatisfyPenalty


def Heuristic_gurobi(scenatio_i,XS,XL):
    model = BuildSubP(scenatio_i,XS,XL)

    # 设置参数以使用启发式算法
    model.setParam('OutputFlag', 0)
    model.setParam('MIPFocus', 1)  # 更加注重寻找可行解
    model.setParam('Heuristics', 1)  # 增强启发式搜索的强度
    #model.setParam('TimeLimit', 10)  # 设置求解时间限制为60秒
    model.setParam('MIPGap', 0.1)  # 设置允许的最小相对MIP Gap为5%
    #model.setParam('NodeLimit', 2)  # 限制分支节点数量

    # 优化模型
    model.optimize()

    return model.objVal


def main():
    #XS=[ 0., 13., 13.,  2., 13., 13.,  0.,  0., 13.,  0., 13.]
    #XL=[ 0.,  0.,  0., 11.,  0.,  0., 13., 13.,  0.,  7.,  0.]
    XS=[ 8.,  6., 10., 11., 14., 17., 11., 14., 17., 13., -0.]
    XL=[ 7., 15., 11., 10.,  7.,  4., 10.,  7.,  4.,  3.,  5.]
    print('yH:',Heuristic(1,XS,XL))

if __name__ == "__main__":
    main()


                



