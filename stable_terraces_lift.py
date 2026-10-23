Lifting directed terraces from Z_m to Z_{2m} with the structure theorem (m even).

lift(P', R', m): given basic directed terraces P', R' of Z_m, solve the GF(2) system
Sigma(P', R') (Gaussian elimination, polynomial time) and return one directed terrace T of Z_{2m}
whose halves reduce to P' and R' + m/2, or None if the pair is not liftable.
sym_R(P', m) gives the partner R' for the symmetric lift (Q' = rev(P')), which is always liftable
with a fibre of dimension m/2 (Bailey-Ollis-Preece, Theorem 3).

Running this file builds the tower Z_4 -> Z_8 -> ... -> Z_256 by symmetric lifts and checks every step,
and tests mixed lifts between the Walecki-type and triangular terrace orbits at m = 8.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from quotient_stable_terraces import all_terraces, is_directed_terrace, is_quotient_stable
from math import gcd

def system(Pp,Rp,m):
    Qp=[(x+m//2)%m for x in Rp]
    A=lambda k:1<<k; B=lambda j:1<<(m+j)
    rows=[(A(0),0)]
    posP={r:k for k,r in enumerate(Pp)}; posQ={r:j for j,r in enumerate(Qp)}
    for r in range(m): rows.append((A(posP[r])|B(posQ[r]),1))
    dP=[(Pp[k+1]-Pp[k])%m for k in range(m-1)]; dQ=[(Qp[j+1]-Qp[j])%m for j in range(m-1)]
    wP=[1 if Pp[k+1]<Pp[k] else 0 for k in range(m-1)]; wQ=[1 if Qp[j+1]<Qp[j] else 0 for j in range(m-1)]
    kP={d:k for k,d in enumerate(dP)}; kQ={d:j for j,d in enumerate(dQ)}
    for rho in range(1,m):
        k=kP[rho]; j=kQ[rho]
        rows.append((A(k)^A(k+1)^B(j)^B(j+1),1^wP[k]^wQ[j]))
    rows.append((B(0)|A(m-1),1))
    return rows,Qp

def solve(rows,nv,rng=None):
    """Gaussian elimination over GF(2); return (consistent, one random solution, dim)."""
    piv={}
    for mask,rhs in rows:
        while mask:
            hb=mask.bit_length()-1
            if hb in piv: pm,pr=piv[hb]; mask^=pm; rhs^=pr
            else: piv[hb]=(mask,rhs); break
        else:
            if rhs: return False,None,None
    free=[i for i in range(nv) if i not in piv]
    x=[0]*nv
    for i in free: x[i]=rng.randint(0,1) if rng else 0
    for hb in sorted(piv):               # ascending pivot bit: lower bits already known
        mask,rhs=piv[hb]; v=rhs
        for i in range(hb):
            if (mask>>i)&1: v^=x[i]
        x[hb]=v
    return True,x,len(free)

def lift(Pp,Rp,m,rng=None):
    rows,Qp=system(Pp,Rp,m)
    ok,x,dim=solve(rows,2*m,rng)
    if not ok: return None,None
    a=x[:m]; b=x[m:]
    T=[Pp[k]+m*a[k] for k in range(m)]+[Qp[j]+m*b[j] for j in range(m)]
    return tuple(T),dim

def norm(t,n): return tuple((x-t[0])%n for x in t)
def sym_R(P,m):   # R' giving the symmetric pair Q' = rev(P')
    return tuple((x-m//2)%m for x in P[::-1])

if __name__=="__main__":
    rng=random.Random(1)
    # 1) tower: symmetric lifts Z_4 -> Z_8 -> ... -> Z_256
    c=(0,3,1,2); m=4
    print("Tower of lifts (random solution at each step; checks: directed terrace of Z_n, quotient-stable):")
    while m<=128:
        T,dim=lift(c,sym_R(c,m),m,rng); n=2*m
        print(f"  Z_{m} -> Z_{n}: fibre dim {dim} (=m/2: {dim==m//2}), directed terrace: {is_directed_terrace(list(T),n)}, stable: {is_quotient_stable(T)}")
        c=T; m=n
    # 2) mixing: which pairs (P',R') of Z_8 mixing Walecki-type and triangular terraces are liftable?
    m=8
    W=tuple((i//2) if i%2==0 else m-1-i//2 for i in range(m)); tri=tuple((k*(k+1)//2)%m for k in range(m))
    units=[1,3,5,7]
    def orbit(t):
        out=set()
        for u in units:
            for rev in (False,True):
                base=t[::-1] if rev else t
                out.add(norm(tuple(u*x%m for x in base),m))
        return out
    OW,OT=orbit(W),orbit(tri)
    print("\nm=8: |Walecki orbit|=%d, |triangular orbit|=%d, disjoint: %s"%(len(OW),len(OT),OW.isdisjoint(OT)))
    cnt={}
    for name1,O1 in (("W",OW),("T",OT)):
        for name2,O2 in (("W",OW),("T",OT)):
            c_ok=0; tot=0; mixsym=0
            for P in O1:
                for R in O2:
                    T,dim=lift(P,R,m)
                    tot+=1
                    if T is not None: c_ok+=1
            cnt[(name1,name2)]=(c_ok,tot)
    for k,v in cnt.items(): print(f"  P' in {k[0]}-orbit, R' in {k[1]}-orbit: liftable pairs {v[0]} of {v[1]}")
    # 3) one explicit mixed lift: P' Walecki-orbit, R' triangular-orbit
    done=False
    for P in sorted(OW):
        for R in sorted(OT):
            T,dim=lift(P,R,m)
            if T is not None and not done:
                n=2*m
                print("\nExample mixed lift Z_8 -> Z_16: P' =",P," R' =",R," fibre dim",dim)
                print("  T =",T," directed:",is_directed_terrace(list(T),n)," stable:",is_quotient_stable(T))
                print("  pi(P) =",tuple(x%m for x in T[:m])," pi(Q) =",tuple(x%m for x in T[m:]))
                d=[(T[i+1]-T[i])%n for i in range(n-1)]
                print("  symmetric:",all((d[k]+d[n-2-k])%n==0 for k in range(n-1)))
                done=True
