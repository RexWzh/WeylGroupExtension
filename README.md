# Weyl 扩张群及泛性质相关工具
## 主要函数

1. Thetas(s,n) # 返回生成元集合 $\{\theta_i\}_{i=1}^n$
2. OrderMatrixOfGens(gens) # 返回生成元的群阶矩阵
3. GroupByOrderMathix(order_mat) # 群阶矩阵作为泛性质，生成群
4. GroupTreeOfMaxDepth(G,depth) # 生成群树（BFS）
5. UniversalPropertyOfGroup(G,max_depth) # 返回群的泛性质关系

## 关于求泛性质函数
1. 这里使用的算法存在缺陷，导出的商集并不是极小的！
2. 利用 Humphreys 反射群与 Coxeter 群 1.9 的技巧，可以给出 Coxeter 的极小关系
3. 但对于一般的情况，极小性很难得到（文献？）
4. 求一般群作为自由群商群的极小商集，这应该是个复杂的组合问题（大概率不存在一般的可行算法，否则 GAP 在处理某些情况时，也不至于卡死）
5. 调用这里的函数时，一些情形可以通过设置 max_depth，导出极小商集
