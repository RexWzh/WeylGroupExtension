def ExpOfNilpotentMat(mat, max_depth=36, base_ring=ZZ):
    """
    1. 返回幂零矩阵的指数矩阵
    2. max_depth 指定 x^n=0 中，n 的最大值
    3. 默认基环为整数环
    4. 输入必须为 Sage 矩阵类型
    """
    dim = mat.dimensions()[0]
    exp_mat = matrix.identity(dim) # 单位矩阵
    new,i = mat,1
    while not new.is_zero(): # 运算直到 new 为零矩阵
        exp_mat += 1/factorial(i) * new
        i += 1
        new *= mat
    return matrix(ZZ,exp_mat)

# 必须用伴随表示，不能是自然表示
def Thetas(s,n):
    """返回给定型的生成元集合 {\theta_i}_{i=1}^n"""
    exp = ExpOfNilpotentMat # 简写函数名
    L = LieAlgebra(ZZ,cartan_type=[s,n])
    # 一组基元
    basis = L.basis()
    # 上三角 e_i
    es = [L.e(i+1).adjoint_matrix() for i in range(n)]
    # 下三角 f_i
    fs = [L.f(i+1).adjoint_matrix() for i in range(n)]
    thetas = [exp(es[i])*exp(-fs[i])*exp(es[i]) for i in range(n)]
    return thetas

def OrderMatrixOfGens(gens):
    """
    1. 计算 gens 的群阶矩阵
    2. gens 的基环必须为整数环（或支持 .multiplicative_order 方法）
    """
    orders = [[(a*b).multiplicative_order() for a in gens] for b in gens]
    return matrix(ZZ,orders)

def GroupByOrderMatrix(order_mat):
    """群阶矩阵 -> 自由群商群"""
    n = order_mat.dimensions()[0]
    FG = FreeGroup(n) # 自由群
    gens = FG.gens() # 生成元
    # 生成关系
    relations = [(gens[i]*gens[j])^order_mat[i,j] for i in range(n) for j in range(n)]
    return FG/relations

def GroupTreeOfMaxDepth(G,depth=-1):
    """
    1. BFS 获取群元素，按长度分层
    2. 获取前 depth 层，首层记0，默认获取整树（有限群）
    3. 允许生成元重复
    """
    # 群树，初始为零层
    tree = [[G.one()]] 
    if depth == 0: # 只有零层
        return tree
    gens = G.gens() # 生成元
    S = [G.one()] # 已生成的群元素
    tree.append([]) # 添加一层
    # 将 gens 化简（去重复元素，以及单位元）
    for a in gens:
        if a not in S:
            S.append(a)
            tree[-1].append(a)
    if depth == 1: # 只有一层
        return tree
    gens = tree[-1] # 简化的生成元
    
    # 循环开始，直到深度足够，或者群生成完毕
    while depth-1 and len(tree[-1]):
        depth -= 1
        tree.append([])
        for a in tree[-2]:
            for b in gens:
                c = a * b
                if c not in S:
                    S.append(c)
                    tree[-1].append(c)
    return tree