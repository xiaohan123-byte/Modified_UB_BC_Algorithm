import numpy as np
from gurobipy import *
from input import *
from UBnBC_Node import Node
from Heuristic import Heuristic_gurobi
import heapq
import time

class UBnBC:
    
    def __init__(self, branching_rule="frac", node_rule="dfs", ub=np.inf):
        """Unified Branch and Benders Cut

        Parameters
        ----------
        branching_rule : str, optional
            How to define in which variable to branch on, by default "min"
            Options are:
                - "min": branching on the smallest violation to floor
                - "max": branching on the largest violation to ceil
                - "frac": branching on the most fractional value
        
        node_rule : str, optional
            How to explore the search tree, by default "dfs"
            Options are:
            - "dfs": Depth first search
            - "bfs": Breadth first search
        """
        self.branching_rule = branching_rule  
        self.node_rule = node_rule # 应该没有被用到，待check
        self.ub=ub # 全局ub
        self.notbranched = [] # 还未被分支的node的变量名称列表
        self.xs=[6 for _ in range(len(stations))] # 初始解 给每个站都分配6块小电池
        self.xl=[6 for _ in range(len(stations))] # 初始解 给每个站都分配6块大电池

    def __call__(   # 在你将实例当作函数调用时才会被调用，对象示例.xx()
        self,
        int_tol: float = 1e-6,   # 容差
        max_iter: int = 10000,   # 最大迭代次数
        verbose: bool = False     # 控制是否输出中间过程
    ):
        
        # Start Timer
        start_time = time.time()
        
        # Initialize values
        explored = 0 #迭代次数
        fathomed = 0 #已经被查明（剪枝）的节点数
        Z = [] # candidate solution set

        min_improved_lb = np.inf # 集合Z中的最优下界值
        min_solution = None # 集合Z中的最优下界值对应的解
        modified3 = True

        # Select rootnode
        baseMP=BuildMasterP(num_scenario,prob)  # 这俩是全局变量
        base_notbranched=[]  # 初始化未分支的结点的变量名称列表，rootnode的变量是全没有分支的--这就是论文中的待分支变量集合V
        for i in stations:
            base_notbranched.append(baseMP._XS[i].VarName)
            base_notbranched.append(baseMP._XL[i].VarName)
        #print(base_notbranched)

        '''
        输出：
            [
        'XS[s1]', 'XL[s1]',
        'XS[s2]', 'XL[s2]',
        'XS[s3]', 'XL[s3]',
        'XS[s4]', 'XL[s4]',
        'XS[s5]', 'XL[s5]',
        'XS[s6]', 'XL[s6]',
        'XS[s7]', 'XL[s7]',
        'XS[s8]', 'XL[s8]',
        'XS[s9]', 'XL[s9]',
        'XS[s10]', 'XL[s10]',
        'XS[s11]', 'XL[s11]'
        ]
        '''

        '''for var in baseMP.getVars():
            if var.VarName.startswith('XS'):
                base_notbranched.append(var.VarName)
            elif var.VarName.startswith('XL'):
                base_notbranched.append(var.VarName)'''
        
        self.relaxP = [] # 松弛子问题的所有场景下的求解对象
        self.MW = [] # MW问题的所有场景下的求解对象
        for i in range(num_scenario):
            self.relaxP.append(BuildRelaxP(i,self.xs,self.xl)) 
            self.MW.append(BuildMW(i,self.xs,self.xl,self.xs,self.xl,10))

        # 创建根节点
        node = Node(BuildMasterP(num_scenario,prob),self.branching_rule,int_tol,base_notbranched)

        # Warm Start
        '''XS_ws,XL_ws = warm_start_H()
        #Compute SP_w(x) for each w
        OptCut={}
        for i in range(num_scenario):
            yrelaxSP,pi = self.solve_relaxP(i,XS_ws,XL_ws)
            OptCut[i]=(yrelaxSP,pi)

        node.n_masterP=self.add_opt_cut(node.n_masterP,XS_ws,XL_ws,OptCut)'''


        # queue 表示分支定界树上的pendant node 论文里面的N
        queue = []
        explored = explored + 1

        # 为什么这里没有把根节点加入queue？

        # 用来迭代M-W点
        XS0 = np.array([0 for _ in range(len(stations))]) #初始化为0列表
        XL0 = np.array([0 for _ in range(len(stations))])

        # Iterate until all good nodes are explored
        while explored < max_iter:
            # Solve the master problem
            node.solve()  # 求到整数最优
            print(f"current UB:{self.ub}")

            # Check feasibility
            if not node.feasible:
                fathomed = fathomed +1
                if verbose:
                    print(f"Fathom node {explored}: Infeasible")
                if len(queue) > 0:
                    _,_,node = heapq.heappop(queue) # 用于从堆队列（优先队列）中弹出并返回最小的元素
                    #node = self._pop_next(queue)
                    explored = explored + 1
                    continue
                else:
                    break

            # Check if cTx+theta>ub
            elif node.objVal>self.ub:
                fathomed = fathomed +1
                if verbose:
                    print(f"Fathom node {explored}: cTx+theta>UB")
                if len(queue) > 0:
                    _,_,node = heapq.heappop(queue)
                    #node = self._pop_next(queue)
                    explored = explored + 1
                    continue
                else:
                    break

            else:
                # Compute SP_w(x) for each w
                flag_optimality = 1  #用于判断是否添加pareto optimal cut  等于0就添加
                ySP_opt=[]   #记录每个子问题的线性松弛问题的目标函数值
                sum_yrelaxSP=0 #记录每个场景下子问题的加权目标和
                SPstart_time = time.time() #记录所有场景线性松弛子问题的求解时间
                if verbose:
                    print('Start Solve SP relaxation')
                OptCut={}  # 最优割集
                for i in range(num_scenario):
                    yrelaxSP,pi = self.solve_relaxP(i,node.XS,node.XL) #求解每个场景的线性松弛子问题
                    sum_yrelaxSP+=prob[i]*yrelaxSP
                    ySP_opt.append(yrelaxSP)
                    if node.theta[i]+int_tol<yrelaxSP:  #如果theta i < 该场景的子问题的最优值 
                        flag_optimality = 0
                        #OptCut[i]=[yrelaxSP,pi]
                SPend_time = time.time() 
                if verbose:
                    print('Solve SP relaxation time:',SPend_time-SPstart_time)
                
                if flag_optimality==0:
                    if verbose:
                            print(f'Node {explored}:Exist theta<QLP_w(x), Add continuous cut')
                    #node.n_masterP=self.add_opt_cut(node.n_masterP,node.XS,node.XL,OptCut) 

                    # Solve Magnanti-Wong
                    # 更新 core point
                    XS0 = node.XS*0.5 + XS0*0.5
                    XL0 = 0.5*node.XL + 0.5*XL0
                    
                    for i in range(num_scenario):
                        yrelaxSP_p,pi_p = self.solve_MW(i,node.XS,node.XL,XS0,XL0,ySP_opt[i])
                        OptCut[i]=[yrelaxSP_p,pi_p]

                    node.n_masterP=self.add_opt_cut(node.n_masterP,XS0,XL0,OptCut)
                    continue
                
                # flag_optimality==1 即该结点的松弛子问题已求解完，不能再添加最优割
                else:
                    improved_lb = node.cTx+sum_yrelaxSP  # 提升下界
                    if improved_lb+int_tol>self.ub:   # 如果 最小化问题的节点下界已大于全局UP，剪枝
                        fathomed = fathomed +1
                        if verbose:
                            print(f"Fathom node {explored}: cTx+yrelaxSP>UB")
                        if len(queue) > 0:
                            #node = self._pop_next(queue)
                            _,_,node = heapq.heappop(queue)
                            explored = explored + 1
                            continue
                        else:
                            break
                    else: # 存下这个整数解
                        Z.append((node.XS,node.XL,improved_lb,node.cTx))

                        if modified3:

                            # 如果最优下界更新了，记录最小的improved_lb解，求子问题到整数最优
                            if improved_lb < min_improved_lb:
                                min_improved_lb = improved_lb
                                print(f"improved_lb:{improved_lb}")

                                min_solution = Z.pop(0)
                                ysubP=solve_subP(min_solution[0],min_solution[1])
                                cur_opt=min_solution[3]+ysubP
                                if cur_opt<=self.ub:
                                    self.ub=cur_opt
                                    self.XS_opt,self.XL_opt=min_solution[0],min_solution[1]
                                    print(f"cur_opt:{cur_opt}")

                            else:
                                #Solve subproblem for all scenario
                                sum_yH=0
                                for i in range(num_scenario):
                                    yH = Heuristic_gurobi(i,node.XS,node.XL)
                                    sum_yH += prob[i]*yH
                                improved_ub = node.cTx+sum_yH
                                if improved_ub<self.ub:
                                    self.ub = improved_ub


                        else:
                            #Solve subproblem for all scenario
                            sum_yH=0
                            for i in range(num_scenario):
                                yH = Heuristic_gurobi(i,node.XS,node.XL)
                                sum_yH += prob[i]*yH
                            improved_ub = node.cTx+sum_yH
                            if improved_ub<self.ub:
                                self.ub = improved_ub

                        # branch on a new variable
                        if node.notbranched:
                            node.branch_integer()
                            #queue.extend(node.children)
                            heapq.heappush(queue,(-node.objVal,id(node.children[0]),node.children[0]))
                            heapq.heappush(queue,(-node.objVal,id(node.children[1]),node.children[1]))
                            

                            if verbose:
                                print(f"Node {explored}: Heuristic solved, branch on node {explored}, remain {len(queue)} pendant nodes")
                        if len(queue) > 0:
                            _,_,node = heapq.heappop(queue)
                            #node = self._pop_next(queue)
                            explored = explored + 1
                            continue
                        else:
                            break
                        #node = self._pop_next(queue)

        print(f"candidate number:{len(Z)}")
        # Rank Z by ascending order
        if Z:
            Z_sorted = sorted(Z, key=lambda x: x[2])
            print(Z_sorted)
            while len(Z_sorted)>0:
                cur_Z = Z_sorted.pop(0)
                if cur_Z[2]<=self.ub:
                    ysubP=solve_subP(cur_Z[0],cur_Z[1])
                    cur_opt=cur_Z[3]+ysubP
                    if cur_opt<=self.ub:
                        self.ub=cur_opt
                        self.XS_opt,self.XL_opt=cur_Z[0],cur_Z[1]
        
        end_time = time.time()
        print('UCnBC run time: ',end_time-start_time)
        
        return self.ub,self.XS_opt,self.XL_opt,end_time-start_time

    # 用不到
    def _pop_next(self,queue):
        if self.node_rule == "dfs":
            node = queue.pop(-1)
        elif self.node_rule == "bfs":
            node = queue.pop(0)
        else:
            raise ValueError("node_rule must be either 'dfs' or 'bfs'")
        return node

    def add_opt_cut(self,n_masterP,xS,xL,OptCut):
        x_opt = np.hstack((xS,xL))
        theta = []
        x = []

        for var in n_masterP.getVars():
            if var.VarName.startswith('X'):
                x.append(var)
            elif var.VarName.startswith('theta'):
                theta.append(var)

        for i in OptCut:
            y_opt = OptCut[i][0]
            pi = OptCut[i][1]
            piT = np.dot(pi,T_matrix)
            piTx = np.dot(piT,x_opt.transpose())

            expr = LinExpr(piT,x)
            n_masterP.addConstr(theta[i]+expr>=piTx+y_opt)

        n_masterP.update()

        return n_masterP

        
    
    def solve_relaxP(self,i,xS,xL):
        pi_opt = np.empty((0,2*len(S)))
        #relaxedP = BuildRelaxP(i,xS,xL)
        relaxedP = self.relaxP[i]
        for n,station in enumerate(stations):
            relaxedP._XconstrS[station].RHS = xS[n]
            relaxedP._XconstrL[station].RHS = xL[n]
        relaxedP.reset()
        # 选择内点法（Barrier Method）求解
        relaxedP.setParam('Method', 2)
        relaxedP.setParam('Presolve', 2)  # 激进预处理
        #relaxedP.setParam('Heuristics', 0.1)
        relaxedP.setParam('OutputFlag', 0)
        relaxedP.optimize()
        if relaxedP.status!=GRB.Status.OPTIMAL:
            raise ValueError("Relaxed Problem Infeasible!")
        y_opt = relaxedP.ObjVal
        pi_opt = []
        for station in stations:
            pi_opt.append(relaxedP._XconstrS[station].Pi)
        for station in stations:
            pi_opt.append(relaxedP._XconstrL[station].Pi)

        '''for constr in relaxedP.getConstrs():
            if constr.ConstrName.startswith('Xconstr'):
                pi_opt.append(constr.Pi)'''

        return y_opt,pi_opt

    def solve_MW(self,i,xS,xL,xS0,xL0,ySP):
        pi_opt = np.empty((0,2*len(S)))
        #relaxedP = BuildRelaxP(i,xS,xL)
        relaxedP = self.MW[i]
        for n,station in enumerate(stations):
            relaxedP._XconstrS[station].RHS = xS0[n]
            relaxedP._XconstrL[station].RHS = xL0[n]

            relaxedP.chgCoeff(relaxedP._XconstrS[station], relaxedP._xi, xS[n])
            relaxedP.chgCoeff(relaxedP._XconstrL[station], relaxedP._xi, xL[n])
        relaxedP._xi.setAttr('Obj', ySP)
        relaxedP.reset()
        # 选择内点法（Barrier Method）求解
        relaxedP.setParam('Method', 2)
        relaxedP.setParam('Presolve', 2)  # 激进预处理
        #relaxedP.setParam('Heuristics', 0.1)
        relaxedP.setParam('OutputFlag', 0)
        relaxedP.optimize()
        if relaxedP.status!=GRB.Status.OPTIMAL:
            raise ValueError("Relaxed Problem Infeasible!")
        y_opt = relaxedP.ObjVal
        pi_opt = []
        for station in stations:
            pi_opt.append(relaxedP._XconstrS[station].Pi)
        for station in stations:
            pi_opt.append(relaxedP._XconstrL[station].Pi)

        return y_opt,pi_opt

# 直接求解子问题到最优
def solve_subP(xS,xL):
    y_opt = 0
    for i in range(num_scenario):
        subP = BuildSubP(i,xS,xL)
        subP.setParam('OutputFlag', 0)
        subP.optimize()
        if subP.status!=GRB.Status.OPTIMAL:
            raise ValueError("SubProblem Infeasible!")
        y = subP.ObjVal
        y_opt += prob[i]*y
    return y_opt


def numerical_efficiency(num):
    global num_scenario
    num_scenario = num
    prob_lst = []
    for i in range(num):
        prob_lst.append(1/num)
    global prob
    prob = tuple(prob_lst)

    # UCnBC
    test = UBnBC(branching_rule='frac',node_rule='dfs')
    alg_result = test(max_iter=10000)

    # gurobi DEP
    starttime = time.time()
    model = BuildDEPgpt(num_scenario,prob)
    model.optimize()
    endtime = time.time()
    if model.status == GRB.OPTIMAL:
        y_opt = model.objVal
        XS_sol = [ model.getVarByName(f'XS[{i}]').x for i in S]
        XL_sol = [ model.getVarByName(f'XL[{i}]').x for i in S]
        print("gurobi optim")
        GAP = model.MIPGap

    elif model.status == GRB.TIME_LIMIT:
        y_opt = model.objVal
        XS_sol = [ model.getVarByName(f'XS[{i}]').x for i in S]
        XL_sol = [ model.getVarByName(f'XL[{i}]').x for i in S]
        print("gurobi time limit")
        GAP = model.MIPGap
    gurobi_result = y_opt,XS_sol,XL_sol,GAP,endtime-starttime

    alogorithm_gap = (gurobi_result[0]-alg_result[0])/(gurobi_result[0])

    return gurobi_result,alg_result,alogorithm_gap


def numerical_compare(num):
    global num_scenario
    num_scenario = num
    prob_lst = []
    for i in range(num):
        prob_lst.append(1/num)
    global prob
    prob = tuple(prob_lst)

    # UCnBC
    test = UBnBC(branching_rule='frac',node_rule='dfs')
    y_opt,xs_opt,xl_opt,runtime = test(max_iter=10000)

    total_cost_battery = cS * sum(max(x - 4, 0) for x in xs_opt) + cL * sum(max(x - 2, 0) for x in xl_opt)
    total_cost_electric = 0
    total_penalty_unsatisfied = 0
    total_income_swap = 0

    for scenario_i in range(num_scenario):
        subP = BuildSubP(scenario_i,xs_opt,xl_opt)
        subP.setParam('OutputFlag', 0)
        subP.optimize()

        cost_electric_value = beta * (sum(E[t] * subP._zS[i, t, n].x for i in S for t in range(T) for n in range(KS)) + 
                                    sum(E[t] * subP._zL[i, t, n].x for i in S for t in range(T) for n in range(KL)))

        penalty_unsatisfied_value = alpha * sum(fS[scenario_i, q, t] - sum(subP._yS[q, t, ij].x for ij in AS[q].select(source[q], '*')) + 
                                                fL[scenario_i, q, t] - sum(subP._yL[q, t, ij].x for ij in AL[q].select(source[q], '*')) 
                                                for t in range(T-max_tau) for q in path)

        income_swap_value = sum(SwapPriceS[t + tau[q, ij[0]]] * subP._yS[q, t, ij].x + 
                                (SwapPriceS[t + tau[q, ij[0]]] - gamma) * subP._yL_h[q, t, ij].x 
                                for t in range(T-max_tau) for q in path for ij in AS[q].select(path_data[q]['stations'], '*')) +  sum(SwapPriceL[t + tau[q, ij[0]]] * subP._yL[q, t, ij].x + 
                            SwapPriceS[t + tau[q, ij[0]]] * subP._yS_h[q, t, ij].x 
                            for t in range(T-max_tau) for q in path for ij in AL[q].select(path_data[q]['stations'], '*'))
        
        total_cost_electric += cost_electric_value
        total_penalty_unsatisfied += penalty_unsatisfied_value
        total_income_swap += income_swap_value

    print(f"Optimal solution:{xs_opt,xl_opt,y_opt}")
    print(f"Total battery cost: {total_cost_battery}")
    print(f"Electric Cost: {total_cost_electric/num_scenario}")
    print(f"Penalty for Unsatisfied Demand: {total_penalty_unsatisfied/num_scenario}")
    print(f"Swap Income: {total_income_swap/num_scenario}")



if __name__ == "__main__":
    '''global num_scenario
    num_scenario = 20
    prob_lst = []
    for i in range(20):
        prob_lst.append(1/20)
    global prob
    prob = tuple(prob_lst)
    test1 = UBnBC(branching_rule='frac',node_rule='dfs')
    print(test1(max_iter=10000))'''

    test1 = UBnBC(branching_rule='frac',node_rule='dfs')
    print(test1(max_iter=10000))
    
   #  print(numerical_efficiency(20))

    
