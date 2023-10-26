import numpy as np
from Consts import hbar, Planck, c, MHzmT, SpinOp4, SpinOp2, SpinOp3
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
print("Starting....     ----------------------------------------------------------")

# Add an orientation mixer to average over it all? 
# B field with sinesoidal patterns for different oriantations, return the collection of diagonalised matrices' energy differences, 
# Add Spin Orbit Coupling first from the two G centres????

class H0:
    def __init__(self, 
                 #D = np.diag([0, 0, 0]),
                 D = np.array([[0, 0, 0],       #In Mhz
                               [0,0, 0],
                               [0, 0, 0]]),
                 g = np.diag([0, 0, 0]), 
                 coil = np.array([0,0,320]),        # In mT
                 SpinOp = SpinOp3,
                 MicrowaveField = np.array([0.001,0.001,0.001])): # in MHz of the Frequency of the Microwaves, 
        self.D = D
        #print ("D = ", self.D)
        self.g = g
        #print ("g = ", self.g)
        self.coil = coil
        #print ("coil = ", self.coil)
        self.SpinOp = SpinOp
        self.SOT = np.transpose(self.SpinOp, (0,2,1))
        #print ("SpinOp = ", self.SpinOp)
        #print ("SOT = ", self.SOT)
        self.MicrowaveField = MicrowaveField
        #print ("MicrowaveField = ", self.MicrowaveField)
        self.ZFS = np.einsum('iab,ij,jbc -> ac', self.SOT, self.D, self.SpinOp) # D matrix is expressed in MHz already. 
        #print ("ZFS Frequency = ", np.around(self.ZFS, decimals=2))
        self.Zeeman = np.einsum('iab,ij,j -> ab', self.SOT, self.g, self.coil) * MHzmT  # Convert from mT to MHz 
        #print ("Zeeman Frequency = ", np.around(self.Zeeman, decimals=2))
        self.Microwave = np.einsum('iab,ij,j -> ab', self.SOT, self.g, self.MicrowaveField) * MHzmT # Convert from mT to MHz
        #https://encyclopedia.pub/entry/9965 Equation sources. 
        #print ("Microwave Frequency = ", np.around(self.Microwave, decimals=2)  )
        self.H = self.ZFS + self.Zeeman
        self.EigenValues, self.EigenVectors = np.linalg.eig(self.H)
        self.Diagonaliser = np.linalg.inv(self.EigenVectors)
        self.H_Diagonalized = self.Diagonaliser @ self.H @ self.EigenVectors
        self.M_Diagonalized = self.Diagonaliser @ self.Microwave @ self.EigenVectors
        self.H_Combined = self.H_Diagonalized + self.M_Diagonalized


    def plotMatrix(self, Matrix):
        fig = plt.figure(figsize=(3,2))
        ax_real = fig.add_subplot(311, projection='3d')
        ax_imag = fig.add_subplot(313, projection='3d')
        ax_absolute = fig.add_subplot(131, projection='3d')
        #Plot the element position with the real value
        for i in range(len(Matrix)):
            for j in range(len(Matrix[i])):
                ax_real.bar3d(i, j, 0, color = '#0000FF', dx=0.5, dy=0.5, dz=np.real(Matrix[i][j]))
                ax_imag.bar3d(i, j, 0, color = '#FF0000', dx=0.5, dy=0.5, dz=np.imag(Matrix[i][j]))
                ax_absolute.bar3d(i, j, 0, color = '#00FF00', dx=0.5, dy=0.5, dz=np.absolute(Matrix[i][j]))
        ax_real.title.set_text('Real'), ax_imag.title.set_text('Imaginary'), ax_absolute.title.set_text('Absolute')
        #ax.set_zscale('log')
        plt.show()
    def MicorwaveFrequency(self):
        for i in range(len(self.H)):
            for j in range(len(self.H[i])):
                if i != j:
                    print("Transition between states ", i, " and ", j, " is ", (self.H[i][i] - self.H[j][j])/10**6, "MHz")
        
    def updateMicrowave(self, x, y):
        self.MicrowaveField = np.array([x,y,0])
        self.Microwave = np.einsum('iab,ij,j -> ab', self.SOT, self.g, self.MicrowaveField)
        self.M_Diagonalized = self.Diagonaliser @ self.Microwave @ self.EigenVectors
        self.H_Combined = self.H_Diagonalized + self.M_Diagonalized
        return self.H_Combined, self.M_Diagonalized


a = H0(D=np.diag([6,8,13])*10**6, 
       g=np.diag([2.00, 2.00, 2.25]),
       coil=np.array([0, 0, 320]), 
       SpinOp=SpinOp4, 
       MicrowaveField=np.array([0.0000,0.0000,0.0000]))

a.MicorwaveFrequency()



#a.plotMatrix(a.H_Combined)

"""
a = H0( D=np.diag([-500, -500, 1000]), 
        g=np.diag([2.00, 2.00, 2.25]), 
        coil=np.array([0,0,320]), 
        SpinOp = SpinOp2, SOT = SOT2, 
        MicrowaveField = np.array([.3,.3,0]))

a makes a rabi frequency, 
compare that to energy gap between two states. 
adjust the rabi frequency to be the same as the energy gap.
then see what happens to the other states.
loop. 

print("Ending....       ----------------------------------------------------------")

    def PerturbedEnergyGap(self):
        PerturbedEnergyGap = np.diag([0.0 + 0.0j, 0.0+ 0.0j, 0.0+0.0j, 0.0+0.0j])
        for i in range(len(self.H_Evolve)-1):
            PerturbedEnergyGap[i][i] = (self.H_Evolve[i][i] - self.H_Evolve[i+1][i+1])/10**(9)
        #https://courses.lumenlearning.com/suny-physics/chapter/24-4-energy-in-electromagnetic-waves/
        #https://iopscience.iop.org/article/10.1088/1367-2630/15/9/093016/pdf
        return PerturbedEnergyGap """