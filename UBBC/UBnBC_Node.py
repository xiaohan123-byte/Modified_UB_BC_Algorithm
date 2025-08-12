from typing import List, Any, Tuple
import numpy as np
from gurobipy import *
from input import *
import math
import copy

class Node:
    def __init__(self,n_masterP,branching_rule,int_tol,notbranched):
        self.n_masterP=n_masterP
        self.branching_rule=branching_rule
        self.int_tol=int_tol
        self.notbranched=notbranched
        self.notfixed = notbranched
        '''
        Node of the UBnBC method

        Parameters
        ----------
        self.n_masterP: Gurobi model
        
        self.branching_rule:

        self.int_tol:

        self.feasibility:

        self.objVal:

        self.sol:

        self.X:

        self.cTx:

        self.theta:


        '''

    
    @property

    def integer(self):
        # 如果问题可行，返回所有整数变量名称
        # 如果不可行，返回False
        if self.feasible:
            residuals = self._get_residuals()
            return all(residual <= self.int_tol for residual in residuals.values())
        else:
            return False

    def solve(self):

        cur_masterP = self.n_masterP.copy() #当前主问题是传入主问题的copy（不改变传入主问题本身）
        cur_masterP.setParam('OutputFlag', 0) # 禁用gurobi日志输出
        cur_masterP.optimize() #gurobi求解主问题到整数最优

        if cur_masterP.status==GRB.Status.OPTIMAL: #如果主问题已求解到最优解
            self.feasible = True 
            self.sol = {} # 解字典
            XS_vec = [] # XS的值列表
            XL_vec = [] 
            theta_vec = [] 
            PXS_vec = []
            PXL_vec = []
            for var in cur_masterP.getVars():
                self.sol[var.VarName] = var.x
                if var.VarName.startswith('XS'):
                    XS_vec.append(var.x)
                elif var.VarName.startswith('XL'):
                    XL_vec.append(var.x)
                elif var.VarName.startswith('theta'):
                    theta_vec.append(var.x)
                elif var.VarName.startswith('PXS'):
                    PXS_vec.append(var.x)
                elif var.VarName.startswith('PXL'):
                    PXL_vec.append(var.x)
            self.XS = np.array(XS_vec)
            self.XL = np.array(XL_vec)
            self.theta = np.array(theta_vec)
            self.objVal = cur_masterP.objVal # 目标值还得算上第二阶段的theta
            self.cTx = cS*sum(PXS_vec)+cL*sum(PXL_vec)  # 不算第二阶段的第一阶段目标值
            
        else: # 主问题不可行
            self.feasible = False

    '''def solve_integer(self):
        cur_masterP = self.n_masterP.copy()
        for var in cur_masterP.getVars():
            if var.varName.startswith('X'):
                var.vtype = GRB.INTEGER
        cur_masterP.setParam('OutputFlag', 0)
        cur_masterP.optimize()
        return cur_masterP.objVal'''


    def _get_residuals(self):  
        # 主问题的变量的的小数部分的字典
        residual_dict = {}
        for var_name, var_value in self.sol.items():
            if not var_name.startswith("theta"):
                residual = var_value - math.floor(var_value)
                residual_dict[var_name] = residual
        return residual_dict

    # 定义三种分支规则  实际上没有用到
    def _minimum_violation(self, residuals):
        frac_var = ((var_name,var_value) for var_name,var_value in residuals.items() if var_value>self.int_tol)
        min_var = min(frac_var, key=lambda x: x[1])
        return min_var
    
    def _maximum_violation(self, residuals: np.ndarray):
        frac_var = ((var_name,var_value) for var_name,var_value in residuals.items() if var_value>self.int_tol)
        max_var = max(frac_var, key=lambda x: x[1])
        return max_var
    
    # 绝对值？
    def _most_fractional(self, residuals: np.ndarray):
        frac_var = ((var_name,var_value) for var_name,var_value in residuals.items() if var_value>self.int_tol)
        most_var = min(frac_var, key=lambda x: x[1]-0.5)
        return most_var

    def find_branch_var(self):
        if self.branching_rule == "min":
            branch_var = self._minimum_violation(self._get_residuals())
        if self.branching_rule == "max":
            branch_var = self._maximum_violation(self._get_residuals())
        if self.branching_rule == "frac":
            branch_var = self._most_fractional(self._get_residuals())
        
        return branch_var[0],math.floor(self.sol[branch_var[0]]),math.floor(self.sol[branch_var[0]])+1

    # 老的分支方法 (没有固定整数求解)
    def branch(self):
        var_name,var_floor,var_ceil = self.find_branch_var()
        try:
            self.notbranched.remove(var_name)
        except ValueError:
            pass  # 忽略异常
        n_masterP_c1 = self.n_masterP.copy()
        var_c1 = n_masterP_c1.getVarByName(var_name)
        n_masterP_c1.addConstr(var_c1<=var_floor)
        n_masterP_c1.update()
        n_masterP_c2 = self.n_masterP.copy()
        var_c2 = n_masterP_c2.getVarByName(var_name)
        n_masterP_c2.addConstr(var_c2>=var_ceil)
        n_masterP_c2.update()
        self.children = [
            Node(n_masterP_c1,self.branching_rule,self.int_tol,self.notbranched),
            Node(n_masterP_c2,self.branching_rule,self.int_tol,self.notbranched)]
        return var_name
    
    # 现在的分支方法 
    def branch_integer(self):
        var_name=self.notbranched.pop(0) 
        var_floor = self.sol[var_name]-1
        var_ceil = self.sol[var_name]+1
        n_masterP_c1 = self.n_masterP.copy()
        var_c1 = n_masterP_c1.getVarByName(var_name)
        n_masterP_c1.addConstr(var_c1<=var_floor)
        n_masterP_c1.update()
        n_masterP_c2 = self.n_masterP.copy()
        var_c2 = n_masterP_c2.getVarByName(var_name)
        n_masterP_c2.addConstr(var_c2>=var_ceil)
        n_masterP_c2.update()
        self.children = [
            Node(n_masterP_c1,self.branching_rule,self.int_tol,self.notbranched),
            Node(n_masterP_c2,self.branching_rule,self.int_tol,self.notbranched)]

    # 不用 not_branched，用 not_fixed  深度优先搜索
    def branch_optim1(self):
        var_name = copy.deepcopy(self.notfixed[0])
        var_floor = self.sol[var_name]-1
        var_ceil = self.sol[var_name]+1

        n_masterP_c1 = self.n_masterP.copy()
        var_c1 = n_masterP_c1.getVarByName(var_name)
        n_masterP_c1.addConstr(var_c1<=var_floor)
        n_masterP_c1.update()

        n_masterP_c2 = self.n_masterP.copy()
        var_c2 = n_masterP_c2.getVarByName(var_name)
        n_masterP_c2.addConstr(var_c2>=var_ceil)
        n_masterP_c2.update()

        n_masterP_c3 = self.n_masterP.copy()
        var_c3 = n_masterP_c3.getVarByName(var_name)
        n_masterP_c3.addConstr(var_c3==self.sol[var_name])
        n_masterP_c3.update()
        notfixed_new = self.notfixed.copy()
        notfixed_new.remove(var_name)

        self.children = [
            Node(n_masterP_c1,self.branching_rule,self.int_tol,self.notfixed),
            Node(n_masterP_c2,self.branching_rule,self.int_tol,self.notfixed),
            Node(n_masterP_c3,self.branching_rule,self.int_tol,notfixed_new)]

    # 
    def branch_optim(self):
        var_name=self.notbranched.pop(0) 
        var_floor = self.sol[var_name]-1
        var_ceil = self.sol[var_name]+1

        n_masterP_c1 = self.n_masterP.copy()
        var_c1 = n_masterP_c1.getVarByName(var_name)
        n_masterP_c1.addConstr(var_c1<=var_floor)
        n_masterP_c1.update()

        n_masterP_c2 = self.n_masterP.copy()
        var_c2 = n_masterP_c2.getVarByName(var_name)
        n_masterP_c2.addConstr(var_c2>=var_ceil)
        n_masterP_c2.update()

        n_masterP_c3 = self.n_masterP.copy()
        var_c3 = n_masterP_c3.getVarByName(var_name)
        n_masterP_c3.addConstr(var_c3==self.sol[var_name])
        n_masterP_c3.update()

        self.children = [
            Node(n_masterP_c1,self.branching_rule,self.int_tol,self.notbranched),
            Node(n_masterP_c2,self.branching_rule,self.int_tol,self.notbranched),
            Node(n_masterP_c3,self.branching_rule,self.int_tol,self.notbranched)]


    def remove_element(self,A, n):
        # 确保n是有效的索引
        if n < 0 or n >= len(A):
            raise IndexError("Index out of range")
        # 使用列表切片来去掉第n个元素
        B = A[:n] + A[n+1:]
        return B

    # 根节点
    def branch_root(self):
        self.children = []
        for n,var_name in enumerate(self.notbranched):
            notbranched_c = self.remove_element(self.notbranched,n)
            var_floor = self.sol[var_name]-1
            var_ceil = self.sol[var_name]+1

            n_masterP_c1 = self.n_masterP.copy()
            var_c1 = n_masterP_c1.getVarByName(var_name)
            n_masterP_c1.addConstr(var_c1<=var_floor)
            n_masterP_c1.update()
            self.children.append(Node(n_masterP_c1,self.branching_rule,self.int_tol,notbranched_c))

            n_masterP_c2 = self.n_masterP.copy()
            var_c2 = n_masterP_c2.getVarByName(var_name)
            n_masterP_c2.addConstr(var_c2>=var_ceil)
            n_masterP_c2.update()
            self.children.append(Node(n_masterP_c2,self.branching_rule,self.int_tol,notbranched_c))


