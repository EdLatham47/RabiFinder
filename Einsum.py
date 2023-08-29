import numpy as np
from Consts import hbar, Planck, c, MHzmT, SpinOp4, SpinOp2

SpinOp = SpinOp2
STranspose = SpinOp.transpose(2,1,0)

D = np.diag([0, 0, 1]) 
H = np.einsum('iab,ij,jbc -> ac', SpinOp, D, SpinOp)
#sum through i and j to get D mat components multiplied by the entire spin operator, then mat mult SpinOps with ab bc

print("------------------------------------------------------------")

exit("Check if we did it right.")

a = np.array([[2, -1], [-1, 2]])
Val,Vec = np.linalg.eig(a)
print("Vectors",Vec, "Values",Val)
print("------------------------------------------------------------")
Inverse = np.linalg.inv(Vec)
print("Inverse",Inverse)
b = Inverse @ a @ Vec
print("b",b)
print("------------------------------------------------------------")

exit()
x = np.array([[0,1],
            [1,0]])
y = np.array([[0,-1j],
            [1j,0]])
z = np.array([[1,0],
            [0,-1]])
S = np.array([x,y,z])
print("Starting....")
B = np.array([0, 0, 1])
#g = np.diag([1., 2., 2.25])
g = np.array([[1., 0.1, 0.01],
            [0.1, 2.0, 0.1],
            [0.01, 0.1, 2.25]])
STranspose = S.transpose(0,2,1)
# print(STranspose.shape, "STranspose", STranspose)
# B and g first, then spin op.
Zeeman = np.einsum('i...,ij,j...',STranspose, g, B)
EigenvaluesZeeman = np.linalg.eigvals(Zeeman)
print(Zeeman.shape, "Zeeman", Zeeman)
print("EigenvaluesZeeman:",EigenvaluesZeeman)
D = np.array([[1, 1, 1],
            [0.05, 0.02, 0.01],
            [100, 100, 100] ])
ZFS = np.einsum('i...,ij,j...',STranspose, D, S)
print(ZFS.shape, "ZFS", ZFS)
EigenvaluesZFS = np.linalg.eigvals(ZFS)

MicrowaveB = np.array([1, 1, 0])
Microwave = np.einsum('i...,ij,j...',STranspose, g, MicrowaveB)

print(Zeeman, "Zeeman", Zeeman.shape)
exit("In1")
Zeeman = np.einsum('jkl,ji->ikl',STranspose, g)
print(Zeeman, "Zeeman", Zeeman.shape)
Zeeman = np.einsum('i...,i',Zeeman, B)
print(Zeeman, "Zeeman", Zeeman.shape)
#alternative = np.einsum('jkl,ij,a->ia', STranspose, g, B)
print("------------------------------------------------------------")
# Tensor contraction of spin operator and its transpose
ZFS = np.einsum('ij,jkl->ikl',D,S)
print(ZFS.shape, "expectation", ZFS)
print("------------------------LAST------------------------------")
ZFS = np.einsum('i...,i...', STranspose, ZFS)
print(ZFS.shape, "expectation", ZFS)
exit()
ZFS = np.einsum('akl,ikl->ia',STranspose, ZFS)
print(ZFS.shape, "expectation", ZFS)
exit()
print("---------------------------WELL DONE---------------------------------")
alternative = np.einsum('akl,ij,jkl->ia', STranspose, D, S)
print(alternative.shape, "alternative", alternative)
EigenvaluesZFS = np.linalg.eigvals(ZFS)
print(EigenvaluesZFS)