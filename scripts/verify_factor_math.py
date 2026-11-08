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
print(f"\nidentity block: {sum(ok)}/{len(ok)} checks passed")
_ok_first=list(ok)

# ===== Per-factor attribution (Predictive Geometry 35.7, Label-Alpha 28.4-28.5) and PCA alignment (Factor 5.3-5.4) =====
def _attribution_block():
    rng=np.random.default_rng(11); ok=[]
    def check(n,c): ok.append(bool(c)); print(("PASS " if c else "FAIL ")+n)
    cen=lambda v:v-v.mean(0)
    T,K,m=60,3,4; Ns=[int(rng.integers(50,90)) for _ in range(T)]
    Fs,Ys,Xs,valids=[],[],[],[]
    for d,N in enumerate(Ns):
        F=cen(rng.normal(size=(N,m))); psi=rng.normal(size=(m,K))*.4+np.array([.3,0,0,0])[:,None]
        u=rng.normal(size=(N,K)); Y=cen(F@psi+u)
        x=np.sign(psi[1,0])*F[:,1]+0.6*F[:,0]+0.5*u[:,0]+rng.normal(size=N)   # timing f1, static f0, specific
        v=rng.random(N)>0.1; xs=np.zeros(N); xs[v]=(x[v]-x[v].mean())/x[v].std()
        Fs.append(F); Ys.append(Y); Xs.append(xs); valids.append(v)
    W=np.concatenate([np.full(N,1/(T*N)) for N in Ns]); Y=np.vstack(Ys)
    SY=Y.T@(W[:,None]*Y); ev,E=np.linalg.eigh(SY); B=E@np.diag(ev**-.5)
    Z=[Yd@B for Yd in Ys]
    nrm=np.sqrt(sum(xd@xd/(T*N) for xd,N in zip(Xs,Ns))); Xs=[xd/nrm for xd in Xs]
    cD=sum(Xs[d]@Z[d]/(T*Ns[d]) for d in range(T))
    # per-date label factor regression on the feature's analysis cross-section
    ebar=np.zeros((T,m)); Psi=np.zeros((T,m,K)); spec=np.zeros(K); cP=np.zeros(K); cM=np.zeros(K)
    for d in range(T):
        v=valids[d]; Fv=cen(Fs[d][v]); Yv=Ys[d][v]; xv=Xs[d][v]
        Ps=np.linalg.pinv(Fv)@Yv; U=Yv-Fv@Ps
        ebar[d]=Fv.T@xv/Ns[d]; Psi[d]=Ps
        spec+=xv@U@B/(T*Ns[d])
        xP=Fv@np.linalg.pinv(Fv)@xv; xM=xv-xP
        cP+=xP@Yv@B/(T*Ns[d]); cM+=xM@Yv@B/(T*Ns[d])
    a=np.array([(1/T)*np.einsum('d,dk->k',ebar[:,k],Psi[:,k,:])@B for k in range(m)])
    check("D block = sum_k factor terms + specific",np.allclose(cD,a.sum(0)+spec))
    check("factor terms == D-representation of x^P (same cross-section)",np.allclose(a.sum(0),cP))
    check("specific term == D-representation of x^M",np.allclose(spec,cM))
    stat=np.array([ebar[:,k].mean()*Psi[:,k,:].mean(0)@B for k in range(m)])
    tim=np.array([((ebar[:,k]-ebar[:,k].mean())[:,None]*(Psi[:,k,:]-Psi[:,k,:].mean(0))).mean(0)@B for k in range(m)])
    check("static + timing = factor term (each k)",np.allclose(stat+tim,a))
    ystar=sum([0])  # target direction: first horizon
    rD=np.zeros(K); rD[0]=1; rD=sum(Xs[d]*0 for d in range(1)) if False else rD
    r=np.linalg.solve(B,np.eye(K))[:,0]*0+ (B.T@np.array([1,0,0]))  # direction of target y* = Y[:,0] in whitened coords
    r=r/np.linalg.norm(r)
    A=cD@r; parts=a@r
    check("target alignment additive: A = sum_k A_k + A_spec",np.isclose(A,parts.sum()+spec@r))
    print("   A_k static/timing (horizon target):",np.round(stat@r,4),np.round(tim@r,4)," spec:",round(spec@r,4))
    # circular shift: static term invariant, timing term changes
    def shifted(s):
        e2=np.roll(ebar,s,axis=0); return np.array([e2[:,k].mean()*Psi[:,k,:].mean(0)@B for k in range(m)]),np.array([((e2[:,k]-e2[:,k].mean())[:,None]*(Psi[:,k,:]-Psi[:,k,:].mean(0))).mean(0)@B for k in range(m)])
    s2,t2=shifted(7)
    check("static term invariant under circular date shift",np.allclose(s2,stat))
    check("timing term changes under shift",not np.allclose(t2,tim))
    # basis dependence: rotate F -> F A
    Arot=rng.normal(size=(m,m)); Ablk=np.zeros((m,m)); Ablk[:2,:2]=rng.normal(size=(2,2)); Ablk[2:,2:]=rng.normal(size=(2,2))
    def terms(Amat):
        out=np.zeros((m,K)); tot=np.zeros(K)
        for d in range(T):
            v=valids[d]; Fv=cen(Fs[d][v])@Amat; Ps=np.linalg.pinv(Fv)@Ys[d][v]; e=Fv.T@Xs[d][v]/Ns[d]
            out+=(e[:,None]*Ps)@B/T
        return out
    tr=terms(Arot); tb=terms(Ablk)
    check("total factor term invariant to any invertible basis change",np.allclose(tr.sum(0),a.sum(0)))
    check("per-factor terms change under general rotation",not np.allclose(tr,a))
    check("group sums invariant under block-diagonal change",np.allclose(tb[:2].sum(0),a[:2].sum(0)) and np.allclose(tb[2:].sum(0),a[2:].sum(0)))
    # PCA factors: sign flip / rotation inside near-degenerate eigenspace; Procrustes alignment recovers reference basis
    N=200; L=rng.normal(size=(N,3)); L[:,2]*=0.999
    R=cen(np.column_stack([rng.normal(size=(300,3))@np.diag([3,1.5,1.49])@L.T]).reshape(300,N)) if False else None
    Ret=rng.normal(size=(300,3))@np.diag([3,1.5,1.49])@cen(L).T+0.5*rng.normal(size=(300,N))
    def pca_load(Rt,k=3):
        Rt=Rt-Rt.mean(0); U,S,Vt=np.linalg.svd(Rt,full_matrices=False); return Vt[:k].T
    V1=pca_load(Ret[:150]); V2=pca_load(Ret[150:])
    ang=np.degrees(np.arccos(np.clip(np.linalg.svd(V1.T@V2,compute_uv=False),-1,1)))
    print("   principal angles between windows (deg):",np.round(ang,2))
    col_cos=np.abs(np.sum(V1*V2,axis=0)); print("   column-wise |cos| (unaligned):",np.round(col_cos,3))
    Uq,_,Vq=np.linalg.svd(V2.T@V1); Q=Uq@Vq; V2a=V2@Q
    print("   column-wise cos after Procrustes:",np.round(np.sum(V1*V2a,axis=0),3))
    check("Procrustes alignment does not change the projector",np.allclose(V2a@V2a.T,V2@V2.T))
    check("subspace stable (largest angle small) while columns 2-3 may be unstable",ang.max()<25)
    # Label-Alpha 28.5: attribution is linear in the signal -> A_k(sum beta_p x_p) = sum beta_p A_k(x_p) (same cross-section)
    X2=[np.where(valids[d],rng.normal(size=Ns[d]),0) for d in range(T)]
    X2=[np.where(valids[d],x-x[valids[d]].mean(),0) for d,x in enumerate(X2)]
    def Ak(sig):
        out=np.zeros(m); st=np.zeros(m); ti=np.zeros(m); es=np.zeros((T,m)); gs=np.zeros((T,m))
        for d in range(T):
            v=valids[d]; Fv=cen(Fs[d][v]); Ps=np.linalg.pinv(Fv)@Ys[d][v]
            es[d]=Fv.T@sig[d][v]/Ns[d]; gs[d]=Ps@B@r
        return (es*gs).mean(0), es.mean(0)*gs.mean(0), ((es-es.mean(0))*(gs-gs.mean(0))).mean(0)
    b1,b2=0.7,-1.3
    comb=[b1*Xs[d]+b2*X2[d] for d in range(T)]
    ok_lin=all(np.allclose(Ak(comb)[i],b1*Ak(Xs)[i]+b2*Ak(X2)[i]) for i in range(3))
    check("28.5 A_k, static and timing are linear in the signal (feature tracing)",ok_lin)
    print(f"\nattribution block: {sum(ok)}/{len(ok)} passed")
    
    return ok
_ok_second=_attribution_block()
print(f"\nALL: {sum(_ok_first)+sum(_ok_second)}/{len(_ok_first)+len(_ok_second)} checks passed")
