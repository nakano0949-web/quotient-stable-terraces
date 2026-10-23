import sys, time
from itertools import product
from collections import Counter

def is_directed(t,n):
    return len(set(t))==n and len({(t[i+1]-t[i])%n for i in range(n-1)})==n-1 and 0 not in {(t[i+1]-t[i])%n for i in range(n-1)}

def is_terrace(t,m):
    """BOP definition (cyclic group Z_m): permutation with t_0 = 0, derived 2-sequencing condition."""
    if len(t)!=m or len(set(t))!=m: return False
    cnt=Counter()
    for i in range(m-1):
        d=(t[i+1]-t[i])%m
        cnt[min(d,(-d)%m)]+=1      # class of {d,-d}
    for x in range(1,m):
        c=min(x,(-x)%m)
        if x!=(-x)%m:
            if x==c and cnt[c]!=2: return False
        else:
            if cnt[c]!=1: return False
    return True

def count_basic_terraces(m):
    """number of terraces of Z_m (BOP sense) with t_0 = 0, by DFS."""
    cnt=0
    t=[0]; used=[False]*m; used[0]=True; cls=Counter()
    def cap(c): return 1 if c==(-c)%m else 2
    def rec():
        nonlocal cnt
        if len(t)==m:
            cnt+=1; return
        for x in range(1,m):
            if used[x]: continue
            d=(x-t[-1])%m; c=min(d,(-d)%m)
            if cls[c]>=cap(c): continue
            used[x]=True; cls[c]+=1; t.append(x)
            rec()
            t.pop(); cls[c]-=1; used[x]=False
    rec(); return cnt

def symmetric_directed_terraces(n):
    """all directed terraces of Z_n with t_0=0 and symmetric derived sequencing (d_{n-2-k} = -d_k)."""
    m=n//2
    out=[]
    d=[]; t=[0]; used=[False]*n; res=[False]*m; res[0]=True
    usedd=[False]*n
    def rec(k):
        if k==m-1:
            T=t[:]          # t_0..t_{m-1}
            T_full=T+[0]*m
            for j in range(m):   # T_{n-1-j} = T_j + m
                T_full[n-1-j]=(T[j]+m)%n
            if is_directed(T_full,n): out.append(tuple(T_full))
            return
        for x in range(1,n):
            if x==m or usedd[x] or usedd[(-x)%n]: continue
            y=(t[-1]+x)%n
            if res[y%m]: continue
            usedd[x]=True; res[y%m]=True; t.append(y)
            rec(k+1)
            t.pop(); res[y%m]=False; usedd[x]=False
    rec(0)
    return out

def projection_info(T,n):
    m=n//2
    P=[x%m for x in T[:m]]
    return P

t0=time.time()
print("Z_m basic terraces (BOP sense) :", {m:count_basic_terraces(m) for m in range(2,9)})
print()
print(f"{'n':>3} {'#symm. directed':>16} {'BOP l*#terr(Z_m)':>17} {'half-proj is terrace':>21} {'half-proj directed':>19}")
for n in [4,6,8,10,12,14,16]:
    m=n//2
    S=symmetric_directed_terraces(n)
    f=1 if m%2==0 else 0                 # elements of order 2 in Z_m
    l=2**((m+f-1)//2)
    pred=l*count_basic_terraces(m)
    allterr=all(is_terrace(projection_info(T,n),m) for T in S)
    ndir=sum(1 for T in S if is_directed(projection_info(T,n),m)) if m>1 else len(S)
    print(f"{n:>3} {len(S):>16} {pred:>17} {str(allterr):>21} {ndir:>19}   ({time.time()-t0:.0f}s)",flush=True)
