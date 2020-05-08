#!/usr/bin/python3
import numpy as np
import time
import math
import random
def Writer_Split(Matrix,Number_1,Number_2):
    if Number_2 != 0 :
        for i in range(Number_1):
            for j in range(Number_2):
                print("%.3f" % Matrix[i][j],end=" \t")
            print ("")
    else :
        for i in range (Number_1):
            print("%.3f" % Matrix[i],end="  ")
    return 0

Nx = int(input("Введите кол-во неизвестных: "))
N = 6  
R1 = 10*10**(-3) #Внутренний радиус кольца м
R2 = 50*10**(-3) #Внешний радиус кольца м
Lx = R2-R1 #длина расчетной области м
delta = 0.1*(N+4)*10**(-3) #Толщина м
Bi_value = N/100 #Число Био
m_value = ((2*Bi_value)**0.5)/delta #Число m
Theta_0 = 5*N #Разность температур *C
Lambda = 101 #коэфф теплоп-ти Вт/м*К
T_liquid = 0 #*C
Qv = 0	

#Nx = 300

dx = (R2-R1)/Nx
XFace = np.linspace(R1,R2,Nx+1)
AFace = np.zeros(Nx+1)
for j in range(Nx+1):
        AFace[j] = 2*math.pi*(R1 + dx*j)*delta
XP = np.zeros(Nx)
for j in range(Nx):
        XP[j] = (XFace[j] + XFace[j+1])/2.0
VolumP = np.zeros(Nx) 
for j in range(Nx):
        VolumP[j] = math.pi*(XFace[j+1]**2-XFace[j]**2)*delta
Gamma = np.zeros(Nx)
Gamma[:] = Lambda
Alpha = np.zeros(Nx,dtype=float)
Alpha = 2.0*Gamma/dx
AW = np.zeros(Nx,dtype=float)
AW[1:Nx] = 1.0/(1.0/Alpha[0:Nx-1]+1.0/Alpha[1:Nx])*AFace[1:Nx]
AE = np.zeros(Nx,dtype=float)
AE[0:Nx-1] = 1.0/(1.0/Alpha[0:Nx-1]+1.0/Alpha[1:Nx])*AFace[1:Nx]
CSource = np.zeros(Nx,dtype=float)
CSource[:] = 2*Bi_value*Lambda/(delta**2) #10**7
VSource = np.zeros(Nx,dtype=float)
VSource[:] = T_liquid 
AP = np.zeros(Nx,dtype=float)
AP = AW + AE + CSource*VolumP
b = np.zeros(Nx,dtype=float)
b = CSource*VSource*VolumP
AP[0] += Alpha[0]*AFace[0] 
AP[Nx-1] += (Alpha[Nx-1])*AFace[Nx]
b[0] = Theta_0+T_liquid
b[Nx-1] = (Bi_value*Lambda/delta)*AFace[Nx]*T_liquid
A = np.zeros((Nx,Nx),dtype=float)
#for i in range(Nx):
#    for j in range(Nx):
#        A[i][j] = random.uniform(1.0,10.0)
for j in range(1,Nx-1):
        A[j,j-1] = -AW[j] #коэфф aw
        A[j,j+1] = -AE[j] #коэфф ae
        A[j,j] = AP[j] #коэфф аp
A[(Nx-1),(Nx-1)] = AP[Nx-1]
A[Nx-1,Nx-2] = -AW[Nx-1]
A[0,0] = 1
A[0,1] = -0

Number_Vector = Nx
#Matrix_A = [[2,1,3],[1,-2,1],[3,2,2]]
#Vector_B = [9,-2,7]
#Number_Vector = 3
Matrix_A = A
Vector_B = b

#print('')
#print('System:')
#Writer_Split(Matrix_A,Number_Vector,Number_Vector)
#print(Vector_B)

Time_NP = time.clock()

Solve_NP = np.linalg.solve(Matrix_A,Vector_B)
print('')
#print('Numpy solve system:',Solve_NP)
#print('NP 1,2:',Solve_NP[0],' , ',Solve_NP[1])
print('Numpy calculating time:',time.clock()-Time_NP)
print('')
#print('Приведенные к диагональному виду:')

#Приведение к треугольному виду:

Time_Solve = time.clock()

for j in range(Number_Vector-1):
    for i in range(j+1,Number_Vector):
        Mu = Matrix_A[i][j]/Matrix_A[j][j]
        Vector_B[i] -= Mu*Vector_B[j]
        for k in range(Number_Vector):
            Matrix_A[i][k] -= Mu*Matrix_A[j][k]
#print('')
#Writer_Split(Matrix_A,Number_Vector,Number_Vector)
#print(Vector_B)
#print('')

for i in range(Number_Vector):
    Mu = 1.0/Matrix_A[i][i]
    Vector_B[i] /= Matrix_A[i][i]
    for j in range(Number_Vector):
        Matrix_A[i][j] *= Mu 
#print('Приведенные к стандартному виду:')
#Writer_Split(Matrix_A,Number_Vector,Number_Vector)
#print(Vector_B)
#print('')

#Решение СЛАУ:
Solve_Vector = [0 for i in range(Number_Vector)]
i = Number_Vector - 1
while(i !=-1):
    Sum = 0
    j = Number_Vector - 1
    while(j != i):
        Sum += Solve_Vector[j]*Matrix_A[i][j]
        j -=1
    #print(Sum)
    Solve_Vector[i] = Vector_B[i] - Sum
    i -=1
#print('')
#print('Gauss solve:',Solve_Vector)
#print('Gauss solve 1,2:',Solve_Vector[0],' , ',Solve_Vector[1])
print('Gauss calculating time:',time.clock()-Time_Solve)
