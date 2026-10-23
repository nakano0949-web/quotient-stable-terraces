Fiber structure of quotient-stable terraces of Z_{2m} (m even).

A quotient-stable terrace T=(P,Q) (t_0 = 0) satisfies:
  * pi(P) = P' is a basic directed terrace of Z_m and ends at m/2;
  * pi(Q) = R' + m/2 with R' a basic directed terrace of Z_m  (the junction difference is m);
  * given (P', R'), the lifts of T are the solutions of an affine system over GF(2)
    (coset condition, split of the two lifts of each difference, junction = m),
    so the fiber over (P', R') is empty or has 2^k elements.
Total number of stable terraces with t_0 = 0  =  sum over pairs (P', R') of 2^k.

Usage: python3 stable_terraces_fibers.py [m ...]     e.g.  python3 stable_terraces_fibers.py 2 4 6 8 10
(m = 12, i.e. n = 24, has 3856^2 ~ 1.5e7 pairs: expect several minutes in pure Python.)
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from quotient_stable_terraces import all_terraces, is_directed_terrace, is_quotient_stable
from collections import Counter

def solve_gf2(rows, nvars):
    """rows: list of (mask, rhs). Return (consistent, rank)."""
    piv={}  # pivot bit -> (mask,rhs)
    for mask,rhs in rows:
        while mask:
            hb=mask.bit_length()-1
            if hb in piv:
                pm,pr=piv[hb]; mask^=pm; rhs^=pr
            else:
                piv[hb]=(mask,rhs); break
        else:
            if rhs: return False,0
    return True,len(piv)

def fiber(Pp,Rp,m):
    """Pp, Rp basic directed terraces of Z_m. Returns (consistent, dim) of the set of stable terraces over (Pp, Qp = Rp+m/2)."""
    Qp=[(x+m//2)%m for x in Rp]
    # variables: a_0..a_{m-1} -> bits 0..m-1 ; b_0..b_{m-1} -> bits m..2m-1
    A=lambda k:1<<k
    B=lambda j:1<<(m+j)
    rows=[(A(0),0)]                       # a_0 = 0
    posP={r:k for k,r in enumerate(Pp)}; posQ={r:j for j,r in enumerate(Qp)}
    for r in range(m):                    # (1) cosets: a(r) xor b(r) = 1
        rows.append((A(posP[r])|B(posQ[r]),1))
    dP=[(Pp[k+1]-Pp[k])%m for k in range(m-1)]; dQ=[(Qp[j+1]-Qp[j])%m for j in range(m-1)]
    wP=[1 if Pp[k+1]<Pp[k] else 0 for k in range(m-1)]; wQ=[1 if Qp[j+1]<Qp[j] else 0 for j in range(m-1)]
    kP={d:k for k,d in enumerate(dP)}; kQ={d:j for j,d in enumerate(dQ)}
    for rho in range(1,m):                # (2) lifts of difference rho split between P and Q
        k=kP[rho]; j=kQ[rho]
        mask=A(k)^A(k+1)^B(j)^B(j+1)
        rows.append((mask,1^wP[k]^wQ[j]))
    rows.append((B(0)|A(m-1),1))          # (3) junction = m
    ok,rank=solve_gf2(rows,2*m)
    return ok,(2*m-rank) if ok else None

def build_all(Pp,Rp,m):
    """brute-force the fiber (for validation, small m): enumerate bit assignments."""
    n=2*m; Qp=[(x+m//2)%m for x in Rp]; out=[]
    for abits in range(1<<(m-1)):
        a=[0]+[(abits>>i)&1 for i in range(m-1)]
        P=[Pp[k]+m*a[k] for k in range(m)]
        for bb in range(1<<m):
            b=[(bb>>j)&1 for j in range(m)]
            Q=[Qp[j]+m*b[j] for j in range(m)]
            T=tuple(P+Q)
            if is_directed_terrace(list(T),n): out.append(T)
    return out

if __name__=="__main__":
    totals={4:2,8:8,12:56,16:896,20:106160}
    ms=[int(a) for a in sys.argv[1:] if a.isdigit()] or [2,4,6,8,10]
    for m in ms:
        n=2*m; t0=time.time()
        Bm=all_terraces(m)
        assert all(t[-1]==m//2 for t in Bm)
        tot=0; dist=Counter(); consistent=0
        for P in Bm:
            for R in Bm:
                ok,dim=fiber(P,R,m)
                if ok:
                    consistent+=1; tot+=2**dim; dist[dim]+=1
        print(f"n={n} (m={m}): basic terraces of Z_m D={len(Bm)}, pairs={len(Bm)**2}, consistent pairs={consistent}, "
              f"sum of fibers={tot}, expected stable={totals.get(n,'?')}, match={tot==totals.get(n)}  ({time.time()-t0:.0f}s)")
        print("   fiber dimension distribution (dim: #pairs):",dict(sorted(dist.items())),flush=True)
    # validate fiber sizes against brute force for m=4,6
    for m in (4,6):
        Bm=all_terraces(m); bad=0
        for P in Bm:
            for R in Bm:
                ok,dim=fiber(P,R,m); bf=len(build_all(P,R,m))
                if bf!=(2**dim if ok else 0): bad+=1
        print(f"validation m={m}: fibers differing from brute force: {bad}")
