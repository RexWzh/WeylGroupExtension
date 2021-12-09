### 考虑计算效率，使用gap进行交互 ###
LieExp := function(mat) 
    ## 求导子的指数函数,第二参数默认取0 ##
    local i,exp_mat,new_mat;
    i := 1;
    exp_mat := IdentityMat(Size(mat),Rationals);
    new_mat := mat;
    while not IsZero(new_mat) do #循环至添加矩阵为0
        exp_mat := exp_mat + 1/Factorial(i)*new_mat;
        i := i + 1;
        new_mat := mat^i;
    od;
    return exp_mat;
end;;

LieThetas := function(s,n)
    ## 求thetai=exp(adei)*exp(-adfi)*exp(adei) 的矩阵表达 ##
    local L,b,m,ade,adf,theta;
    L := SimpleLieAlgebra(s,n,Rationals); #新建李代数
    b := Basis(L); #获取李代数一组基
    ade := List([1..n],i->AdjointMatrix(b,b[i])); #初始化ei
    m := (Size(b)-n)/2; #正根总数
    adf := List([1..n],i->AdjointMatrix(b,b[m+i])); #初始化fi
    theta := List([1..n],i->LieExp(ade[i])*LieExp(-adf[i])*LieExp(ade[i])); #初始化thetai
    return theta;
end;


QuotientGroup := function(M) 
    # M 为生成关系矩阵，只用到下三角部分
    local G,K,gen,i,j,n;
    n := Size(M);
    G := FreeGroup(n); # n 元自由群
    gen := GeneratorsOfGroup(G); # 生成元
    K := []; # 生成关系
    for i in [1..n] do  # 第i行
        for j in [1..i] do # 第j列
            Add(K,(gen[i]*gen[j])^M[i][j]);
        od;
    od;
    return G/K;
end;;

# 群生成树
GroupTree := function(G,k)
# 生成 k+1行，k>0
local gens,L,x,y,i;
gens := GeneratorsOfGroup(G); # 生成元
L := List([1..k+1],x->[]); # 初始化树
L[1] := [gens[1]^2]; # 第1行
L[2] := List(gens); # 第2行
for i in [3..k+1] do # 开始生成第i行
    for x in L[i-1] do # 对第i-1行操作
        if IsList(x) then # 被标记则跳过
            continue;
        fi;
        for y in gens do #开始生成
            y := x*y; #右乘生分支
            if y in L[i-2] then # 一类断点
                Add(L[i],[y,1]);
            elif y in L[i] then # 二类断点
                Add(L[i],[y,2]);
            else
                Add(L[i],y); # 新分支
            fi;
        od;
    od;
od;
return L;;
end;;


#coxeter群生成函数
CoxeterGroup:=function(k,L) 
#生成元个数为k，L为对应mij，下三角从上到下读
local G,K,gen,i,j,m;
G := FreeGroup(k);#初始化自由群
gen := GeneratorsOfGroup(G);#生成元
K := List(gen,x->x^2);#kernal
m := 0;
for i in [1..k-1] do  #第i行
for j in [1..i] do #第j列
Add(K,(gen[i+1]*gen[j])^L[m+j]);
od;
m := m + i;
od;
#PPrint(K); #调试，检验生成元
return G/K;
end;;