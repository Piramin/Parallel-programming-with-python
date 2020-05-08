#!/usr/bin/python3
#import random
import numpy as np
import time
import math
#import sys


#time_mpi = time.clock()
from mpi4py import MPI

def Writer_Split(Matrix):
    if len(Matrix[0]) != 0 :
        for i in range(len(Matrix)):
            for j in range(len(Matrix[0])):
                print("%.3f" % Matrix[i][j],end=" \t")
            print ("")
    else :
        for i in range (len(Matrix)):
            print("%.3f" % Matrix[i],end="  ")
    return 0

def Cannon(Matrix):
    core = 0
    result = False
    for i in range(len(Matrix[0])-1):
        core += abs(Matrix[i+1][i])
    if core < 10**(-12):
        result = True
    return result

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
        Nx = int(input("Введите кол-во контрольных объемов: "))
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
        #A = [[1 for j in range(Nx)]for i in range(Nx)]
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

        Vector = b  #Правая часть
        Matrix = A  # Матрица
        #print('Исходная матрица:')
        #Writer_Split(Matrix)
        #print('')
        #print('Исходный вектор:')
        #print(Vector)
        Number_proc = size
        data_Matrix = [[0 for j in range(len(Matrix[0])+1)]for i in range(len(Matrix[0])//Number_proc)]
        Sort_Vector = [0 for i in range(len(Matrix))]
########################################################################       
        time_mpi = time.clock()
        print('Calculating...')
        
        #print('Matrix:')
        #Writer_Split(Matrix)
        #print('')
        #print('Vector: ',Vector)
        
        #Первоначальная сборка
        for p in range(1,Number_proc):
            k = 0
            for i in range(len(data_Matrix)*p,len(data_Matrix)*(p+1)):
                if k ==len(data_Matrix):
                    k = 0
                for j in range(len(Matrix[0])):
                    data_Matrix[k][j] = Matrix[i][j]
                data_Matrix[k][len(data_Matrix[0])-1] = Vector[i]
                k +=1
            '''
            print('')
            print('data_Matrix TO p=',p)
            Writer_Split(data_Matrix)
            print('')
            '''
            comm.send(data_Matrix, dest = p)

        N = 0 #Начальное значение флага при зацикливании
        Line = 0
        while(Cannon(Matrix)!=True):

            #Вычисление верхней мини-матрицы(равномерная нагрузка)
            #print('Вычисление верхней части')
            for j in range(len(Matrix[0])):
                for i in range(j+1,len(data_Matrix)):
                    Mu = Matrix[i][j]/Matrix[j][j]
                    Vector[i] -= Mu*Vector[j]
                    for k in range(len(Matrix[0])):
                        Matrix[i][k] -= Mu*Matrix[j][k]
                        if abs(Matrix[i][k]) < 10**(-6):
                            Matrix[i][k] = 0
             
            #Перезаписывание матрицы
            for p in range(1,Number_proc):
                data_Matrix = comm.recv( source = p )
                k = 0
                for i in range(len(data_Matrix)*p,len(data_Matrix)*(p+1)):
                    if k == len(data_Matrix):
                        k = 0
                    for j in range(len(Matrix[0])):
                        Matrix[i][j] = data_Matrix[k][j]
                    Vector[i] = data_Matrix[k][len(data_Matrix[0])-1]
                    k += 1
            
            #Пересборка матрицы
            for i in range(len(Matrix)):
                Sort_Vector[i] = 0
                for j in range(i):
                    if Matrix[i][j] == 0:
                        Sort_Vector[i] += 1
            for j in range(len(Matrix)):
                for i in range(len(Matrix)-j-1):
                    if Sort_Vector[i] > Sort_Vector[i+1]:
                        Sort_Vector[i] , Sort_Vector[i+1] = Sort_Vector[i+1] , Sort_Vector[i]
                        Vector[i] , Vector[i+1] = Vector[i+1] , Vector[i]
                        for k in range(len(Matrix[0])):
                            Matrix[i][k] , Matrix[i+1][k] = Matrix[i+1][k] , Matrix[i][k]
            
            #Отправка мини-матриц по процессам
            for p in range(1,Number_proc):
                k = 0
                for i in range(len(data_Matrix)*p,len(data_Matrix)*(p+1)):
                    if k ==len(data_Matrix):
                        k = 0
                    for j in range(len(Matrix[0])):
                        data_Matrix[k][j] = Matrix[i][j]
                    data_Matrix[k][len(data_Matrix[0])-1] = Vector[i]
                    k +=1
                comm.send(data_Matrix , dest = p)
            
            N +=1
            if N>(len(Matrix[0]))*2:
                print('Error! Circle algoritm!')
                break
            
            for i in range(2,len(Matrix)):
                    if Matrix[i][i-2] == 0:
                        Line += 1
            if Line == len(Matrix[0])-2:
                Line = 1
            else:
                Line = 0

            if Line == 1:
                for k in range(len(data_Matrix)-1,len(Matrix)-1):
                    Mu = Matrix[k+1][k]/Matrix[k][k]
                    Vector[k+1] -= Mu*Vector[k]
                    for j in range(len(data_Matrix)-1,len(Matrix[0])):
                        Matrix[k+1][j] -= Mu*Matrix[k][j]
                break
        
        #print(Number_proc)
        for p in range(1,Number_proc):
            comm.send(False,dest = p)
        #if Number_proc == 2:
        #    comm.send(False, dest = 1)
        if N>(len(Matrix[0]))*2:
            quit()
    
        Solve_Vector = [0 for i in range(len(Matrix[0]))]
        i = len(Matrix[0]) - 1
        while(i !=-1):
            Sum = 0
            j = len(A[0]) - 1
            while(j != i):
                Sum += Solve_Vector[j]*Matrix[i][j]
                j -=1
            Solve_Vector[i] = Vector[i] - Sum
            Solve_Vector[i] /= Matrix[i][i]
            i -=1
        #print('Solve system:',Solve_Vector)
        #print('1: ',Solve_Vector[0],' 2: ',Solve_Vector[1])
        print('Calculating time:',time.clock()-time_mpi)
       # sys.exit()
        
else:
    while(True): 
        data_old = comm.recv( source = 0 )
        '''
        if rank == 1:
            if data_old != False:
                print('')
                print('get data PORT',rank)
                Writer_Split(data_old)
                print('')
            else:
                print('False')
        '''
        #Проверка на *смысловую нагрузку*
        if data_old == False :
            quit()
        
        Zero = 0
        for j in range(len(data_old[0])-1):
            if data_old[0][j]!=0:
                Zero = j
                break
        for k in range(len(data_old)-1):
            for i in range(k+1,len(data_old)):
                Mu = data_old[i][Zero]/data_old[k][Zero]
                for j in range(Zero,len(data_old[0])):
                    data_old[i][j] -= Mu*data_old[k][j]
            '''
            if rank == 1:
                print('')
                print('it old')
                Writer_Split(data_old)
                print('')
            '''
            Zero += 1
        '''
        if rank == 1:
            print('')
            print('send data PORT',rank)
            Writer_Split(data_old)
            print('')
        '''
        comm.send(data_old,dest = 0)
        del data_old
        del Zero
    #quit()

