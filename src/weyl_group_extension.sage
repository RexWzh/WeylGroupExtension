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


def NilExp(mat,base_ring=ZZ, max_order=None):
    """
    返回幂零矩阵的指数矩阵，默认基环为整数环
    """
    n = mat.dimensions()[0] if max_order is None else max_order
    assert (mat^n).is_zero(),"输入的矩阵非幂零！"
    return matrix(base_ring,sum(mat^i/factorial(i) for i in range(n)))

# 必须用伴随表示，不能是自然表示
def Thetas(s, n, reduced=False):
    """
    生成元集合 {\theta_i}_{i=1}^n
    默认略去 Cartan 子代数部分的作用
    """
    max_order = dict(zip("ABCDEFG",[3]*6+[4]))[s]
    Exp = lambda mat:NilExp(mat,max_order=max_order)
    L = LieAlgebra(ZZ,cartan_type=[s,n])
    es = [L.e(i+1).adjoint_matrix() for i in range(n)]
    fs = [L.f(i+1).adjoint_matrix() for i in range(n)]
    thetas = [Exp(es[i])*Exp(-fs[i])*Exp(es[i]) for i in range(n)]
    if reduced: # omit the action on the Cartan subalgebra
        pos_num = (len(L.basis())-n)/2 # number of positive roots
        ind = list(range(pos_num)) + list(range(-pos_num,0))
        thetas = [mat[ind,ind] for mat in thetas]
    return thetas

# replace the default hash function which collapsed
myhash = lambda x: hash(tuple(tuple(line) for line in matrix(x)))

def GroupTree(G, depth=-1, s="s", check=True):
    """
    1. BFS 求群树，按长度分层
    2. depth 为求解层数，首层记0，默认取整树（有限群）
    3. s 设置群元素字符
    返回字符形式的群树，以及哈希字典
    """
    # initialize
    tree,tree_str = [[G.one()]],[[""]]
    hash2str = {myhash(G.one()):""} # hash value to string form
    if depth == 0:
        return tree,[[""]],hash2str
    # layer 1
    gens = list(G.gens()) # generators of the set
    genstr = [s+str(i) for i in range(len(gens))]
    tree_str.append(genstr)
    tree.append(gens)
    hash2str.update({myhash(x):str(i) for x,i in zip(gens,tree_str[-1])})
    assert not check or len(hash2str)==len(gens)+1, "G 的生成元存在相等或者包含1"
    
    # 循环开始，直到深度达到，或者群生成完毕
    while depth-1:
        depth -= 1
        new_str,new_layer = [],[]
        for expr,a in zip(tree_str[-1],tree[-1]):
            for i,b in zip(genstr,gens):
                c = a*b
                hash_c = myhash(c)
                if hash_c not in hash2str: # update new layer
                    hash2str[hash_c] = expr+i
                    new_layer.append(c)
                    new_str.append(expr+i)
        if not len(new_str):
            break
        tree_str.append(new_str)
        tree.append(new_layer)
    return tree,tree_str,hash2str

def UniversalPropertyOfGroup(G,depth=-1, s="s"):
    """
    求群的生成元泛性质
    depth 设置层数，返回 depth-1 层关系
    """
    tree,tree_str,hash2str = GroupTree(G,depth=depth,s="s")
    gens = list(G.gens())
    genstr = [s+str(i) for i in range(len(gens))]
    relations,depth = [],len(tree)
    for i in range(1,depth-1):
        for expr,a in zip(tree_str[i],tree[i]):
            for si,b in zip(genstr,gens):
                hc = myhash(a*b)
                if expr+si not in tree_str[i+1]:
                    relations.append((expr+si,hash2str[hc]))
    # remove extra relations
    i,m = 0,len(relations)
    while i<m-1:
        rel = relations[i][0]
        for j in range(m-1,i):
            if rel in relations[j][0]:
                del relations[j]
        i += 1
        m = len(relations)
    return relations

# def relation2element(relation,gens):
#     """将生成元转为群元素"""
#     one = gens[0]/gens[0]
#     if not relation:return one
#     new = one
#     for i in relation.split("-"):
#         new *= gens[int(i)]
#     return new

# relations2elements = lambda relations,gens: [relation2element(a,gens)/relation2element(b,gens) for a,b in relations]

## 其他函数
def OrderMatrixOfGens(gens):
    """
    1. 计算生成元的群阶矩阵
    2. gens 的基环必须为整数环（或支持 .multiplicative_order 方法）
    """
    # 自带函数的 bug：单位元无法判断阶数
    orders =[[1 if (a*b).is_one() else (a*b).multiplicative_order() for a in gens] for b in gens]
    return matrix(ZZ,orders)

def GroupByOrderMatrix(order_mat):
    """群阶矩阵 -> 自由群商群"""
    n = order_mat.dimensions()[0]
    FG = FreeGroup(n)
    gens = FG.gens()
    relations = [(gens[i]*gens[j])^order_mat[i,j] for i in range(n) for j in range(n)]
    return FG/relations