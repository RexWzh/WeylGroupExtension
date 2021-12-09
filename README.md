# Weyl 扩张群及泛性质相关工具
## 主要函数

1. Thetas(s,n) # 返回生成元集合 $\{\theta_i\}_{i=1}^n$
2. OrderMatrixOfGens(gens) # 返回生成元的群阶矩阵
3. GroupByOrderMathix(order_mat) # 群阶矩阵作为泛性质，生成群
4. GroupTreeOfMaxDepth(G,depth) # 生成群树（BFS）
5. UniversalPropertyOfGroup(G,max_depth) # 返回生成元的泛性质关系 
