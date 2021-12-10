r"""
## 主要函数
1. ExpOfNilpotentMat(mat, max_depth=36, base_ring=ZZ) # 返回幂零阵的指数矩阵
2. Thetas(s,n) 返回生成元集合 \{\theta_i\}_{i=1}^n
3. OrderMatrixOfGens(gens) 返回生成元的群阶矩阵
4. GroupByOrderMathix(order_mat) 群阶矩阵作为泛性质，生成群
5. GroupTreeOfMaxDepth(G,depth) 生成群树（BFS）
6. UniversalPropertyOfGroup(G,max_depth=-1) 求群生成元的泛性质
7. relations2elements(relations,gens) 生成元通过 gens 转为群的子集
"""


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
        assert i < max_depth,"矩阵的%d次幂非零，计算终止"%max_depth
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


def UniversalPropertyOfGroup(G,max_depth=-1):
    """
    1. 求群生成元的泛性质
    2. G 的生成元不能存在相等，或者等于1
    3. max_depth 获取层数，默认遍历整树
    """
    if max_depth == 0:
        return [[G.one()]],[[""]],{}
    # 检验输入
    gens = list(G.gens()) # 生成元
    S = [G.one()] + gens # 已生成群元素
    assert len(S)==len(set(S)),"G 的生成元不能存在相等，或者等于1"
    
    # 群树：元素形式和字符串形式
    tree = [[G.one()],gens] 
    tree_str = [[""],["%d"%i for i in range(len(gens))]]
    # 生成关系
    relations = {}
    while max_depth-1:
        max_depth -= 1
        new_layer = []
        new_layer_str = []
        for s,a in zip(tree_str[-1],tree[-1]):
            for si,b in enumerate(gens):
                c = a * b
                c_str = s + "-%d"%si
                if c in new_layer: # 化简到同一层
                    if not any(s in c_str for s in relations):
                        relations[c_str] = new_layer_str[new_layer.index(c)]
                    continue
                for i,layer in enumerate(tree):
                    if c in layer: # 化简到前边层
                        if not any(s in c_str for s in relations):
                            relations[c_str] = tree_str[i][layer.index(c)]
                        break
                else:
                    new_layer.append(c)
                    new_layer_str.append(c_str)
        if len(new_layer):
            tree.append(new_layer)
            tree_str.append(new_layer_str)
        else:
            break
    return tree,tree_str,relations


def relations2elements(relations,gens):
    """将生成元转为群元素"""
    one = gens[0]/gens[0]
    elements = []
    for key in relations:
        new = one
        for i in key.split("-"):
            new *= gens[int(i)]
        if not relations[key]: # 指向单位元
            elements.append(new)
            continue
        for i in relations[key][::-1].split("-"):
            new /= gens[int(i)]
        elements.append(new)
    return elements