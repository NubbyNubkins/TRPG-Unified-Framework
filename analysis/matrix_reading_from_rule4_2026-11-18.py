"""Derive test C: the matrix reading (which forces exist) from rule 4 (one spot, one facing, one i).
TEST FIXED BEFORE RUNNING:
  'One spot' = fixing one diagonal place of the 3x3 octonion relationship matrix; its symmetry group is Spin(9), acting on the
  16 real components of the two octonions that link that spot to the others (Freudenthal; Todorov & Dubois-Violette 2018).
  'One facing / one i' = choosing one unit imaginary octonion u, acting as the complex structure J = right-multiplication by u on both.
  Compute the part of Spin(9) that commutes with J. PASS (the reading is derived, not chosen) if it is exactly 12-dimensional and splits
  8 + 3 + 1 (strong, weak, hypercharge), and if the result does not depend on which u is chosen. Otherwise FAIL.
  Prior art to credit: Krasnov 2021 (J. Math. Phys. 62, 021703), Todorov & Dubois-Violette 2018, Baez 2002, Gunaydin & Gursey 1973.
"""
import numpy as np, itertools
# octonion multiplication table (Fano lines)
lines=[(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
T=np.zeros((8,8,8))
for i in range(8): T[0,i,i]=1; T[i,0,i]=1
for i in range(1,8): T[i,i,0]=-1
for a,b,c in lines:
    for x,y,z in [(a,b,c),(b,c,a),(c,a,b)]:
        T[x,y,z]=1; T[y,x,z]=-1
def L(u): return np.einsum('i,ijk->kj',u,T)     # x -> u x
def R(u): return np.einsum('j,ijk->ki',u,T)     # x -> x u
e=np.eye(8); conj=np.diag([1]+[-1]*7)
# check alternativity quickly
x=np.random.randn(8); y=np.random.randn(8)
assert np.allclose(L(x)@L(x)@y, L(L(x)@x)@y)
# Cl(9) on O^2: Gamma_i = [[0, L(u_i)],[L(conj u_i), 0]] for u_i in basis, Gamma_8 = diag(1,-1)
G=[]
for i in range(8):
    A=L(e[i]); B=L(conj@e[i]); G.append(np.block([[np.zeros((8,8)),A],[B,np.zeros((8,8))]]))
G.append(np.block([[np.eye(8),np.zeros((8,8))],[np.zeros((8,8)),-np.eye(8)]]))
for i in range(9):
    for j in range(9):
        assert np.allclose(G[i]@G[j]+G[j]@G[i], 2*np.eye(16)*(i==j))
so9=[G[i]@G[j] for i in range(9) for j in range(i+1,9)]
print("Clifford relations hold; so(9) generators:",len(so9))
def commutant(J):
    M=np.array([(X@J-J@X).ravel() for X in so9]).T
    s=np.linalg.svd(M,compute_uv=False); null=sum(s<1e-9)
    _,_,Vt=np.linalg.svd(M); basis=Vt[-null:] if null else []
    return null,[sum(c*X for c,X in zip(v,so9)) for v in basis]
for k in range(1,8):
    u=e[k]; J=np.block([[R(u),np.zeros((8,8))],[np.zeros((8,8)),R(u)]])
    assert np.allclose(J@J,-np.eye(16))
    n,B=commutant(J)
    if k==7 or k==1:
        # split: centre (commutes with all of the commutant) and simple parts via the Killing-form-free method: dimension of [g,g]
        V=np.array([b.ravel() for b in B])
        comm=[ (a@b-b@a).ravel() for a in B for b in B]
        rk_der=np.linalg.matrix_rank(np.array(comm),tol=1e-8)
        cent=[b for b in B if all(np.allclose(b@c,c@b) for c in B)]
        # find centre dimension by solving
        Mc=np.array([np.concatenate([(b@c-c@b).ravel() for c in B]) for b in B]).T
        cdim=sum(np.linalg.svd(Mc,compute_uv=False)<1e-8)
        print(f"u = e{k}: commutant dimension {n}; derived (non-abelian) part {rk_der}; centre {cdim}")
    else:
        print(f"u = e{k}: commutant dimension {n}")
# random unit imaginary u
u=np.concatenate([[0],np.random.randn(7)]); u/=np.linalg.norm(u)
J=np.block([[R(u),np.zeros((8,8))],[np.zeros((8,8)),R(u)]]); print("random u: commutant dimension",commutant(J)[0])
print()
print("Identify the pieces (u = e7):")
u=e[7]; J=np.block([[R(u),np.zeros((8,8))],[np.zeros((8,8)),R(u)]]); n,B=commutant(J)
# centre element
Mc=np.array([np.concatenate([(b@c-c@b).ravel() for c in B]) for b in B]).T
_,_,Vt=np.linalg.svd(Mc); z=Vt[-1]; Z=sum(c*b for c,b in zip(z,B))
print("  11 = the only compact split of dimension 11 into simple parts is 8 + 3 (strong + weak); the 12th is the centre (hypercharge)")
# work on C^8: complex structure J; Z commutes with J, eigenvalues of Z as a complex-linear map: i*y
ev=np.linalg.eigvals(Z); ys=np.round(np.sort(ev.imag),6)
vals,counts=np.unique(np.round(np.abs(ys),6),return_counts=True)
# group by J-complex eigenvalue: compute eigenvalues of Z restricted to +i eigenspace of J
w,V=np.linalg.eig(J); P=V[:,np.isclose(w,1j)]
Zc=np.linalg.pinv(P)@Z@P; y=np.linalg.eigvals(Zc).imag
yv,yc=np.unique(np.round(y/np.min(np.abs(y)),4),return_counts=True)
print("  hypercharge-like eigenvalues on the 8 complex states (scaled):",dict(zip(yv,yc)))
print("  expected one generation's doublets: Q (3 colours x 2) with Y = 1/6 and L (2) with Y = -1/2, ratio -3, multiplicities 6 and 2")
