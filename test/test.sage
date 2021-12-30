################ initialize #################
load("../src/weyl_group_extension.sage")

## lie algebras of simple type
LieDatas = [["A",range(1,9)],["B",range(2,9)],["C",range(3,9)],["D",range(4,9)],
            ["E",range(6,9)],["F",[4]],["G",[2]]]

################ test tools #################
from pprint import pprint
class Time():
    import time
    def tic(self):
        self.time = time.time()
    def toc(self,text='time cost:'):
        t = self.time
        self.tic()
        print(text+'%.3f'%(self.time-t))
tt = Time()
tic = tt.tic
toc = tt.toc

# 保存数据
def save_vari(v,name):
    '''save variables'''
    import pickle
    fout = open(name,'wb')
    pickle.dump(v,fout)
    fout.close()
    return

# 读取数据
def read_vari(name):
    '''read variables'''
    import pickle
    fin = open(name,'rb')
    a = pickle.load(fin)
    fin.close()
    return a

###########################################

### Group structure ###
# e.g. Type A
res = [['level','|G|','|W|','|K|','Structure of G','Structure of W','Structure of W']] 
s = "A" # lie algebra of Type A
for n in range(1,8):
    print(n,end="\t")
    thetas = [gap(mat) for mat in Thetas(s,n,reduced=True)]
    G = gap.Group(thetas)
    K = gap.Group([x^2 for x in thetas])
    W = WeylGroup([s,n])
    res.append([n, gap.Size(G), W.order(), gap.Size(K),
                gap.StructureDescription(G), W.structure_description(), gap.StructureDescription(K)])
table(res) # pprint(res)


### Chevalley basis ###
def check(s,n):
    "are all coefficients lie in 1,-1,0"
    pos_num = (len(L.basis())-n)/2
    thetas = Thetas(s,n,reduced=True)
    for mat in thetas:
        if max(max(mat)) > 1 or min(min(mat)) < -1:
            return False
    return True

for s,l in LieDatas:
    print(s,end="\t")
    for n in l:
        print(n,end="\t")
        assert check(s,n),"Error guessing!"
    print()

    
### universal property ###
s,n = "A",3
G = WeylGroup([s,n])
relations = UniversalPropertyOfGroup(G)
table([rel for rel in relations if len(rel[0])<2*4])


### equalities by Coxeter relations ###
# e.g. Type A
# initialize
s,n = "A",6
print(DynkinDiagram([s,n]))
M = matrix(CoxeterMatrix([s,n]))
thetas = Thetas(s,n,reduced=True)
thetas_square = [t^2 for t in thetas] # generators of K
K = MatrixGroup(thetas_square)
_,_,hash2str = GroupTree(K,s="t",check=False)
# experiment
res = []
for i,t1 in enumerate(thetas):
    res.append([])
    for j,t2 in enumerate(thetas):
        t = (t1*t2)^M[i,j]
        res[-1].append(hash2str[myhash(t)])