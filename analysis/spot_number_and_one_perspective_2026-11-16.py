"""Rule 4's kept spot in the 3x3 octonion relationship matrix: Peirce grading of the Albert algebra (Jordan 1934; Freudenthal 1954)."""
import numpy as np, itertools
def cd_mul(a,b):
    n=len(a)
    if n==1: return a*b
    h=n//2; a1,a2=a[:h],a[h:]; b1,b2=b[:h],b[h:]
    conj=lambda x: np.concatenate([x[:1],-x[1:]])
    return np.concatenate([cd_mul(a1,b1)-cd_mul(conj(b2),a2), cd_mul(b2,a1)+cd_mul(a2,conj(b1))])
conj=lambda x: np.concatenate([x[:1],-x[1:]])
# basis of 3x3 Hermitian octonion matrices: 3 diagonal reals + 3 off-diagonal octonions (8 each) = 27
def basis():
    B=[]
    for i in range(3):
        M=np.zeros((3,3,8)); M[i,i,0]=1; B.append(M)
    for (i,j) in ((0,1),(0,2),(1,2)):
        for k in range(8):
            M=np.zeros((3,3,8)); e=np.zeros(8); e[k]=1; M[i,j]=e; M[j,i]=conj(e); B.append(M)
    return B
def mm(A,B):
    C=np.zeros((3,3,8))
    for i in range(3):
        for j in range(3):
            for k in range(3): C[i,j]+=cd_mul(A[i,k],B[k,j])
    return C
jordan=lambda A,B: 0.5*(mm(A,B)+mm(B,A))
Bs=basis(); flat=np.array([b.reshape(-1) for b in Bs]); P=np.linalg.pinv(flat.T)
e=Bs[0]   # the observer's spot: diag(1,0,0)
L=np.array([P@jordan(e,b).reshape(-1) for b in Bs]).T
ev=np.round(np.linalg.eigvals(L).real,6)
vals,counts=np.unique(ev,return_counts=True)
print("eigenvalues of 'multiply by the spot' and how many states each:",dict(zip(vals,counts)))
print("spot number = 6*(eigenvalue - 1/3):",{round(6*(v-1/3)):c for v,c in zip(vals,counts)})
# which states are in the eigenvalue-0 block (the far side, 10 states)
idx0=[n for n,b in enumerate(Bs) if np.allclose(jordan(e,b),0)]
lab=["d1","d2","d3"]+[f"x12_{k}" for k in range(8)]+[f"x13_{k}" for k in range(8)]+[f"x23_{k}" for k in range(8)]
print("far side (spot number -2):",[lab[n] for n in idx0])
print("  -> with facing u = e7: colourless = d2, d3, x23_0, x23_7 (4 = the weak pair);  coloured = x23_1..6 (6 = D, Dbar). Same spot number.")
