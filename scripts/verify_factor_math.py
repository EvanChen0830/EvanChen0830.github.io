# Numerical checks for every identity in the Factor note and the factor sections of the
# Universe, Predictive Geometry, Selection and Label-Alpha notes. Run: python3 scripts/verify_factor_math.py
import numpy as np
rng=np.random.default_rng(42); ok=[]
def check(name,cond): ok.append(cond); print(("PASS " if cond else "FAIL ")+name)
cen=lambda v: v-v.mean(0)
T,K,m,G=12,2,3,4
Ns=[int(rng.integers(40,80)) for _ in range(T)]
# ---------- site pipeline: labels -> Z -> Phi (D,V,T), W = 1/(T N_d)
Fs,Ys=[],[]
for N in Ns:
    grp=np.r_[np.arange(G),rng.integers(0,G,N-G)]; H=np.eye(G)[grp]
    beta=rng.normal(1,.3,N)
    F=cen(np.column_stack([beta,H[:,1:],rng.normal(size=N)])); Fs.append(F)
    Ys.append(cen(np.outer(beta,rng.normal(size=K))+H@rng.normal(size=(G,K))*.5+rng.normal(size=(N,K))))
w=np.concatenate([np.full(N,1/(T*N)) for N in Ns]); dates=np.concatenate([[d]*N for d,N in enumerate(Ns)])
Y=np.vstack(Ys); SY=Y.T@(w[:,None]*Y); ev,E=np.linalg.eigh(SY); Z=Y@E@np.diag(ev**-.5)
check("exact whitening Z'WZ=I",np.allclose(Z.T@(w[:,None]*Z),np.eye(K)))
def centd(M):
    out=M.copy()
    for d in range(T): idx=dates==d; out[idx]-=out[idx].mean(0)
    return out
iu=np.triu_indices(K); sv=lambda z: np.array([ (z[i]*z[j])*(1 if i==j else np.sqrt(2)) for i,j in zip(*iu)])
GV=centd(np.array([sv(z) for z in Z])); GT=centd(np.linalg.norm(Z,axis=1)[:,None]*Z)
def orthW(R):
    Xi=R.T@(w[:,None]*R); l,U=np.linalg.eigh(Xi); keep=l>1e-10*l.max(); return R@U[:,keep]@np.diag(l[keep]**-.5)
PD=Z; RV=GV-PD@(PD.T@(w[:,None]*GV)); PV_=orthW(RV)
RT=GT-PD@(PD.T@(w[:,None]*GT))-PV_@(PV_.T@(w[:,None]*GT)); PT=orthW(RT)
Phi=np.column_stack([PD,PV_,PT]); blocks=[range(0,K),range(K,K+PV_.shape[1]),range(K+PV_.shape[1],Phi.shape[1])]
check("basis orthonormal Phi'WPhi=I",np.allclose(Phi.T@(w[:,None]*Phi),np.eye(Phi.shape[1])))
# ---------- features with missing values (zero-filled) and global normalization
def make_feature(load):
    x=np.zeros(len(w)); valid=np.zeros(len(w),bool)
    for d,N in enumerate(Ns):
        idx=np.where(dates==d)[0]; v=rng.random(N)>0.12
        raw=Fs[d]@load+0.7*Ys[d][:,0]+rng.normal(size=N)
        r=raw[v]; x[idx[v]]=(r-r.mean())/r.std(); valid[idx[v]]=True
    return x/np.sqrt(x@(w*x)),valid
def split(x,valid):
    xP=np.zeros_like(x)
    for d in range(T):
        idx=np.where((dates==d)&valid)[0]; Fv=cen(Fs[d][valid[dates==d]])
        xP[idx]=Fv@np.linalg.pinv(Fv)@x[idx]
    return xP,x-xP
c=lambda s: s@(w[:,None]*Phi); nW=lambda s: s@(w*s)
x,valid=make_feature(np.array([1.,.5,.5,.5,.3])); xP,xM=split(x,valid)
check("35.1 missing names stay zero in x^P and x^M",np.allclose(xP[~valid],0) and np.allclose(xM[~valid],0))
check("35.1 both pieces within-date centered",all(abs(xP[dates==d].sum())<1e-10 and abs(xM[dates==d].sum())<1e-10 for d in range(T)))
check("35.1 vector identity c=c_P+c_M",np.allclose(c(x),c(xP)+c(xM)))
rho2=nW(xP); check("35.1 norm identity ||xP||^2+||xM||^2=1",np.isclose(rho2+nW(xM),1))
for b,blk in zip("DVT",blocks):
    cb=lambda s: c(s)[list(blk)]
    check(f"35.1 energy identity block {b}",np.isclose(cb(x)@cb(x),cb(xP)@cb(xP)+cb(xM)@cb(xM)+2*cb(xP)@cb(xM)))
# library C, SVD, P_V projection
X=np.column_stack([make_feature(rng.normal(size=5)*s)[0] for s in (1,0,.5,2,1,0)])
Cm=X.T@(w[:,None]*Phi); U,S_,Vt=np.linalg.svd(Cm,full_matrices=False); Vr=Vt[:2].T; PVr=Vr@Vr.T
check("35.1 identity after projection onto V_r",np.isclose((c(x)@PVr)@(c(x)@PVr),(c(xP)@PVr)@(c(xP)@PVr)+(c(xM)@PVr)@(c(xM)@PVr)+2*(c(xP)@PVr)@(c(xM)@PVr)))
xMt=xM/np.sqrt(nW(xM)); Q=lambda s: c(s)@c(s)
check("35.3 Q(x~M)=Q(xM)/(1-rho2)",np.isclose(Q(xMt),Q(xM)/(1-rho2)))
check("Cor16.3-type bound Q(xM)<=1-rho2",Q(xM)<=1-rho2+1e-12)
# per-date restandardization would break the relation
xMd=xM.copy()
for d in range(T):
    idx=(dates==d)&valid; xMd[idx]/=xM[idx].std()
xMd/=np.sqrt(nW(xMd)); check("35.3 per-date restandardization differs (so global is required)",not np.isclose(Q(xMd),Q(xM)/(1-rho2)))
# library level C=C^P+C^M
XP=np.column_stack([split(X[:,j],X[:,j]!=0)[0] for j in range(X.shape[1])]); XM=X-XP
check("35.6 C = C^P + C^M",np.allclose(Cm,XP.T@(w[:,None]*Phi)+XM.T@(w[:,None]*Phi)))
# circular shift null: (x_d,F_d) move together == shifting Phi; equivalence of shifting X vs Phi when N_d equal
# matched thinning: F must be recentered on the subset
F0=Fs[0]; sub=rng.choice(Ns[0],Ns[0]-8,replace=False); xs=cen(rng.normal(size=len(sub)))
Pn=lambda F:F@np.linalg.pinv(F)
check("35.4 thinning needs recentring of F",(not np.isclose((xs-Pn(F0[sub])@xs).sum(),0)) and np.isclose((xs-Pn(cen(F0[sub]))@xs).sum(),0))
# ---------- alpha four cells
a,_=make_feature(np.array([.2,1,0,0,.4])); lin=0.5*x+0.3*X[:,3]; eps=a-lin
parts=[split(v,np.ones(len(w),bool)) for v in (lin,eps)]
check("28.2 four-cell sum",np.allclose(sum(c(p) for pr in parts for p in pr),c(a)))
# ---------- Factor page: Lemma F.0, F.1, projector, F.3, Sec 6, F.4, attribution
N=60; Cn=np.eye(N)-np.ones((N,N))/N; one=np.ones(N); g=1.5
Fr=np.column_stack([one,rng.normal(1,.3,N),rng.normal(size=(N,2))]); F=Cn@Fr
Sp=np.diag(rng.uniform(.5,2,4)); Su=np.diag(rng.uniform(.2,4,N)); A=np.linalg.inv(Su)
om=Cn@rng.normal(size=N)
check("F.0 exposures equal for DN portfolio",np.allclose(om@Fr,om@F))
check("F.0 risk equal for DN portfolio",np.isclose(om@(Fr@Sp@Fr.T+Su)@om,om@(F@Sp@F.T+Su)@om))
check("F.0 intercept column vanishes",np.allclose(F[:,0],0))
F=F[:,1:]; Sp=Sp[1:,1:]; Sr=F@Sp@F.T+Su; Fb=np.column_stack([F,one]); S=np.linalg.inv(Sp)
def Pw(G,Pi): return G@np.linalg.pinv(G.T@Pi@G)@G.T@Pi
y=Cn@rng.normal(size=N)
for Pi,nm in [(np.eye(N),"I"),(A,"GLS")]:
    P=Pw(F,Pi); M=np.eye(N)-P
    check(f"F.1 centered residual (Pi={nm})",np.isclose(one@M@y,0))
    check(f"projector idempotent & Pi-selfadjoint (Pi={nm})",np.allclose(P@P,P) and np.allclose(Pi@P,(Pi@P).T))
    u_,z_=rng.normal(size=N),rng.normal(size=N)
    check(f"inner-product identity (Pi={nm})",np.isclose(u_@Pi@z_,(P@u_)@Pi@(P@z_)+(M@u_)@Pi@(M@z_)))
    check(f"rotation invariance (Pi={nm})",np.allclose(Pw(F@rng.normal(size=(3,3)),Pi),P))
Fspan=np.column_stack([F,one*0]); 
Pi1=np.column_stack([F,one]); Mi=np.eye(N)-Pw(Pi1,A)
check("F.1 remark: 1 in span gives only 1'Pi M y=0",np.isclose(one@A@Mi@y,0))
def kkt(G,mu):
    k=G.shape[1]; Kk=np.block([[g*Sr,G],[G.T,np.zeros((k,k))]]); return np.linalg.solve(Kk,np.r_[mu,np.zeros(k)])[:N]
mu=rng.normal(size=N); MA=np.eye(N)-Pw(Fb,A)
w_star=kkt(Fb,mu)
check("F.3a closed form",np.allclose(w_star,A@MA@mu/g))
check("F.3a dollar-neutral & factor-neutral",np.isclose(one@w_star,0) and np.allclose(F.T@w_star,0))
check("F.3b invariance to mu+c1, M^I_F mu, M^I_Fb mu, mu+Fb b",all(np.allclose(kkt(Fb,v),w_star) for v in (mu+2*one,(np.eye(N)-Pw(F,np.eye(N)))@mu,(np.eye(N)-Pw(Fb,np.eye(N)))@mu,mu+Fb@rng.normal(size=4))))
w_dn=kkt(one[:,None],mu); Si=np.linalg.inv(Sr); nu=(one@Si@mu)/(one@Si@one)
check("F.3c exposure formula",np.allclose(F.T@w_dn,S@np.linalg.inv(S+F.T@A@F)@F.T@A@(mu-nu*one)/g))
check("F.3c zero for A-orthogonal signal",np.allclose(F.T@kkt(one[:,None],MA@mu),0))
check("F.3c nonzero for OLS-neutral signal (Su not prop. I)",np.linalg.norm(F.T@kkt(one[:,None],(np.eye(N)-Pw(Fb,np.eye(N)))@mu))>1e-6)
Sr_save=Sr; Sr=F@Sp@F.T+0.7*np.eye(N)
check("F.3c zero for OLS-neutral signal when Su prop. I",np.allclose(F.T@kkt(one[:,None],(np.eye(N)-Pw(Fb,np.eye(N)))@mu),0)); Sr=Sr_save
eta=3.0; Sr2=F@(Sp+eta/g*np.eye(3))@F.T+Su
wpen=np.linalg.solve(np.block([[g*Sr2,one[:,None]],[one[None,:],np.zeros((1,1))]]),np.r_[mu,0])[:N]
Kp=np.block([[g*Sr+eta*F@F.T,one[:,None]],[one[None,:],np.zeros((1,1))]]); wpen2=np.linalg.solve(Kp,np.r_[mu,0])[:N]
check("quadratic exposure penalty == larger Sigma_psi",np.allclose(wpen,wpen2))
yr=rng.normal(size=N)+3
for Pi,nm,expect in [(np.eye(N),"OLS",True),(A,"GLS",False)]:
    bi=np.linalg.solve(Fb.T@Pi@Fb,Fb.T@Pi@yr)[:3]; bc=np.linalg.solve(F.T@Pi@F,F.T@Pi@(Cn@yr)); ba=np.linalg.solve(Fb.T@Pi@Fb,Fb.T@Pi@(Cn@yr))[:3]
    check(f"Sec6 centered==intercept is {expect} for {nm}",np.allclose(bi,bc)==expect)
    check(f"Sec6 augmented invariant to demeaning y ({nm})",np.allclose(bi,ba))
R=np.linalg.solve(Fb.T@A@Fb,Fb.T@A)
check("Sec6 pure factor portfolios: R Fb=I and dollar-neutral",np.allclose(R@Fb,np.eye(4)) and np.allclose(R[:3]@one,0))
grp=np.r_[np.arange(G),rng.integers(0,G,N-G)]; Hc=Cn@np.eye(G)[grp]
check("F.4 centered dummies rank G-1",np.linalg.matrix_rank(Hc)==G-1)
check("F.4 pinv = equal-weight sum-to-zero",np.isclose((np.linalg.pinv(Hc)@yr).sum(),0))
Fh=np.column_stack([F,Hc[:,1:]]); Fh2=np.column_stack([F,Hc[:,:-1]])
check("F.4 restriction choice leaves fitted values & projector unchanged",np.allclose(Pw(Fh,np.eye(N)),Pw(Fh2,np.eye(N))))
psi=rng.normal(size=3); u=rng.normal(size=N); yy=Cn@(F@psi+u); uh=yy-F@np.linalg.lstsq(F,yy,rcond=None)[0]
check("Sec7 return identity holds by construction",np.isclose(om@yy,(F.T@om)@np.linalg.lstsq(F,yy,rcond=None)[0]+om@uh))
wn=rng.normal(size=N); check("Sec7 variance split needs DN (fails otherwise)",np.isclose(om@Cn@Su@Cn@om,om@Su@om) and not np.isclose(wn@Cn@Su@Cn@wn,wn@Su@wn))
# alpha 28.3 coincidence condition
at=Cn@rng.normal(size=N); MI_F=np.eye(N)-Pw(F,np.eye(N)); MI_Fb=np.eye(N)-Pw(Fb,np.eye(N))
check("28.3 equal Su & mu=c*alpha: alpha_dec = c * alpha^M",np.allclose(MI_Fb@(2.5*at),2.5*MI_F@at))
vol=rng.uniform(.5,2,N); Su_h=np.diag(rng.uniform(.2,4,N)); Ah=np.linalg.inv(Su_h); MG=np.eye(N)-Pw(Fb,Ah)
dec=MG@(vol*at); geo=MI_F@at; cpar=(dec@geo)/(geo@geo)
check("28.3 unequal Su or vol-scaled mu: alpha_dec not a multiple of alpha^M",not np.allclose(dec,cpar*geo))
# ---------- Universe membership share
for d in range(2):
    F=Fs[d]; N=Ns[d]; grp_mask=np.isclose(F[:,1],F[:,1].max())  # members of sector 1 (in F)
    mt=cen(grp_mask.astype(float)); P=Pn(F)
    check(f"21.3 rho2=1 for sector in F (date {d})",np.isclose((P@mt)@(P@mt)/(mt@mt),1))
    # dropped sector 0: rows where all sector dummies are at their minimum
    m0=np.all(np.isclose(F[:,1:4],F[:,1:4].min(0)),axis=1); mt0=cen(m0.astype(float))
    check(f"21.3 rho2=1 for dropped reference sector (date {d})",np.isclose((P@mt0)@(P@mt0)/(mt0@mt0),1))
    r=cen((rng.random(N)<.3).astype(float)); v=(P@r)@(P@r)/(r@r)
    check(f"21.3 rho2 in [0,1], equals R^2 of regression on [1,F] (date {d})",0<=v<=1 and np.isclose(v,1-np.sum((r-np.column_stack([np.ones(N),F])@np.linalg.lstsq(np.column_stack([np.ones(N),F]),r,rcond=None)[0])**2)/(r@r)))
print(f"\n{sum(ok)}/{len(ok)} checks passed")
