################ test tools #################
from pprint import pprint
class Time():
    import time
    def tic(self):
        self.time = time.time()
    def toc(self,text='用时：'):
        t = self.time
        self.tic()
        print(text+'%.3f'%(self.time-t))
tt = Time()
tic = tt.tic
toc = tt.toc

# 保存数据
def save_vari(v,name):
    '''文件存储'''
    import pickle
    fout = open(name,'wb')
    pickle.dump(v,fout)
    fout.close()
    return

# 读取数据
def read_vari(name):
    '''读取文件'''
    import pickle
    fin = open(name,'rb')
    a = pickle.load(fin)
    fin.close()
    return a

###########################################

### 检验泛性质 ###
s,k = "A",4
for n in range(1,6):
    print(s,n)
    gens = Thetas(s,n)
    G = MatrixGroup(gens)
    # 计算结果
    tic()
    _,_,relations = UniversalPropertyOfGroup(G,k)
    toc()
    # 定义商群
    FG = FreeGroup(n)
    elements = relations2elements(relations,FG.gens())
    newG = FG/elements
    tic()
    len_newG = len(newG)
    print(len(G)==len_newG)
    toc()
    

### 泛性质实验 ###
# 初始化
s,n,k = "A",5,4
print(s,n)
gens = Thetas(s,n)
G = MatrixGroup(gens)
# 计时
tic()
tree = GroupTreeOfMaxDepth(G,k)
toc("获取群树用时")
record_len = [len(i) for i in tree] # 记录各层长度数据
print("各层数目",record_len)

# 测试
FG = FreeGroup(n)
gens = FG.gens()
# type_1 = [(gens[i]^4,FG.one()) for i in range(n)]
type_2 = [(gens[i]*gens[i+1]*gens[i],gens[i+1]*gens[i]*gens[i+1]) for i in range(n-1)]
type_3 = [(gens[i]*gens[j],gens[j]*gens[i]) for i in range(n) for j in range(n) if i-j>1]
type_4 = [(gens[i]*gens[i+1]^2*gens[i],gens[i+1]^2) for i in range(n-1)]
odd = FG.one()
for i in range((n+1)//2):
    odd *= gens[2*i]^2
type_5 = [(odd,FG.one())] if n%2 else []
elements = type_1 + type_2 + type_3 + type_4 + type_5
QG = FG / [a/b for a,b in elements] # 商群
tic()
tree = GroupTreeOfMaxDepth(QG,k)
toc()
tree_len = [len(i) for i in tree]
print("修改后的各层数目",tree_len,"\n是否相等",tree_len==record_len)