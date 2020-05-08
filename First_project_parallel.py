#!/usr/bin/python3
import numpy as np
import time
import math

#time_mpi = time.clock()
from mpi4py import MPI

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

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if rank == 0:
        Number_proc = size    
        Nx = int(input("Введите кол-во неизвестных: "))
#######################################################################
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
        for j in range(1,Nx-1):
                A[j,j-1] = -AW[j] #коэфф aw
                A[j,j+1] = -AE[j] #коэфф ae
                A[j,j] = AP[j] #коэфф аp
        A[(Nx-1),(Nx-1)] = AP[Nx-1]
        A[Nx-1,Nx-2] = -AW[Nx-1]
        A[0,0] = 1
        A[0,1] = -0

        Number_X =Nx 
        Vector = b  #Правая часть
        Matrix = A  # Матрица
        data_Vector = [0 for i in range(Number_X+1)]
        data_Matrix = [[0 for j in range(Number_X+1)]for i in range(2)]
########################################################################       
        time_mpi = time.clock()
        print('Calculating...')
        
        k = 0
        for p in range(1,Number_proc):
            comm.send(Number_X+1,dest = p)

        while (k != Number_X-1):

            for p in range(1,Number_proc):
                if p > k:
                    comm.send(True,dest = p)
                else:
                    comm.send(False,dest = k)
                

            for i in range(k+1,Number_X):
                for j in range(Number_X):
                    data_Matrix[0][j] = Matrix[k][j] 
                    data_Matrix[1][j] = Matrix[i][j]
                data_Matrix[0][Number_X] = Vector[k]
                data_Matrix[1][Number_X] = Vector[i]
                comm.send(data_Matrix,dest = i)
                
            for p in range(k+1,Number_proc):
                data_Vector = comm.recv(source = p)
                for j in range(Number_X):
                    Matrix[p][j] = data_Vector[j]
                Vector[p] = data_Vector[Number_X]
                
            k +=1
        comm.send(False,dest = Number_proc-1)
        
        #Обратный ход

#        print('Диагональный вид:')
#        Writer_Split(Matrix,Number_X,Number_X)
#        print(Vector)
#        print('')
        
        Solve_Vector = [0 for i in range(Number_X)]
        i = Number_X - 1
        while(i !=-1):
            Sum = 0
            j = Number_X - 1
            while(j != i):
                Sum += Solve_Vector[j]*Matrix[i][j]
                j -=1
            Solve_Vector[i] = Vector[i] - Sum
            Solve_Vector[i] /= Matrix[i][i]
            i -=1
#        print('Solve system:',Solve_Vector)
        print('Calculating time:',time.clock()-time_mpi)

else:
    Number_X = comm.recv(source = 0)
    key = True
    while(key):
        key = comm.recv(source = 0)
        #print('key = ',key,',rank = ',rank)
        if key==False:
            break
        #print('Calc,rank =',rank)
        data_old = comm.recv(source = 0) #Получение данных , source - номер процесса
        for k in range(Number_X):
            if data_old[0][k] != 0:
                parameter = k
                break
        Mu = data_old[1][parameter]/data_old[0][parameter]
        for i in range(Number_X):
            data_old[1][i] -=Mu*data_old[0][i]
        data = [0 for i in range(Number_X)]
        for i in range(Number_X):
            data[i] = data_old[1][i]
        
        comm.send(data,dest = 0)  #Отправка данных , dest - номер процесса
    
        del data_old
        del data
        del Mu
    del key 
