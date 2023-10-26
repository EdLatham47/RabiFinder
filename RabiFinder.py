import scipy.optimize as opt
from GenerateFrequencies import H0
import numpy as np
from Consts import hbar, Planck, c, MHzmT, SpinOp4, SpinOp2, SpinOp3

#All Frequencies in MHz for Spinnach Convenience. Will also get the mT value of Microwaves.

x,y = 100, 100 #mT of the Microwave field

#Hamiltonian of the system
H = H0(D = np.array([[-100, 0.1, 0],
                    [0.1, -100, 0.1],
                    [0, 0.1, 200]]), #MHz
        g = np.diag([2.00, 2.00, 2.25]), #Unitless
        coil = np.array([0,0,320]), #mT
        MicrowaveField = np.array([x,0,0]), #mT
        SpinOp = SpinOp3)

initial_guess = [x, y]
RabiFreq = np.array([])

for i in range(len(H.H_Combined)-1):
    H.updateMicrowave(x,y)
    FreqDiff = 100
    while np.real(FreqDiff) > 1: #MHz
        EnFreq = H.H_Combined[i][i] - H.H_Combined[i+1][i+1]
        MicFreq = H.M_Diagonalized[i][i+1]
        FreqDiff = EnFreq - MicFreq
        if FreqDiff > 0:
            H.H_Combined, H.M_Diagonalized = H.updateMicrowave(H.MicrowaveField[0] + FreqDiff/2, H.MicrowaveField[1] + FreqDiff/2)
        else:
            H.H_Combined, H.M_Diagonalized = H.updateMicrowave(H.MicrowaveField[0] - FreqDiff/2, H.MicrowaveField[1] - FreqDiff/2)
        print("FreqDiff: ", FreqDiff)
        #H.plotMatrix(H.H_Combined)
    np.append(RabiFreq, H.MicrowaveField)
    print("Mic Field x,y,z: ", H.MicrowaveField)

print("RabiFreq: ", RabiFreq)
#print("Microwave Field Strength:", H.MicrowaveField/MHzmT, "mT")
print("H_Combined: ", np.abs(H.H_Combined))


print("EigenVec, Val:", np.abs(H.EigenVectors), H.EigenValues)


"""
for i in range(len(H.H_Combined)-1):
    result = opt.least_squares( CalcDiff, initial_guess, 
                                ftol=10, 
                                bounds = ([0,0],[0.001,0.001]),
                                cost = 0)
    if result.success:
        fitted_params = result.x
        print(fitted_params)
    else:
        raise ValueError(result.message)
"""