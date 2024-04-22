from math import sin
import time
from matplotlib import pyplot as plt
import methods as meth

# A task 
index_number = 201267
d = index_number % 10 #7
c = int(index_number/10)%10 #6
e = int(index_number/100)%10 #2
f = int(index_number/1000)%10 #1
a1 = 5 + e
a2 = a3 = -1
N = 9 * c *  d

def createMatrix(a1,a2,a3,N):
    A = []
    for i in range(N):
        row = []
        for j in range(N):
            if i == j:
                row.append(a1)
            elif abs(i - j) == 1:
                row.append(a2)
            elif abs(i - j) == 2:
                row.append(a3)
            else:
                row.append(0)
        A.append(row)
    return A

def createVector(N,f):
    b = []
    for i in range(N):
            b.append(sin(i * (f + 1)))
    return b
            
# B task
def jacobiMethod(A, b, tol=1e-9, max_iter=3000):
    startTime = time.time()
    matA = meth.matcopy(A)
    vecB = meth.veccopy(b)
    k = 0 
    n = len(matA)
    x = [0] * n  # Initial guess for the solution vector
    while k < max_iter:
        x_new = [0] * n  # Temporary storage for updated solution vector
        for i in range(n):
            summation = 0
            for j in range(n):
                if j != i:
                    summation += matA[i][j] * x[j]
            x_new[i] = (vecB[i] - summation) / matA[i][i]
        
        # Calculate the residual vector
        residual = [vecB[i] - sum(matA[i][j] * x_new[j] for j in range(n)) for i in range(n)]
        
        # Check convergence
        if all(abs(element) < tol for element in residual):
            break
        
        # Update the solution vector
        x = x_new
        k += 1
    
    print("Jacobi's method")
    print('time:', time.time() - startTime)
    print('iterations:', k)
    print()
    return time.time() - startTime

        
def gaussSeidelMethod(A, b):
    startTime = time.time()
    matA = meth.matcopy(A)
    vecB = meth.veccopy(b)
    k = 0
    n = len(matA)
    vecX = [0] * n
    while True:
        for i in range(n):
            tmp = vecB[i]
            for j in range(n):
                if i != j:
                    tmp -= matA[i][j] * vecX[j] 
            tmp /= matA[i][i]
            vecX[i] = tmp
        
        res = [vecB[i] - sum(matA[i][j] * vecX[j] for j in range(n)) for i in range(n)]
        
        if meth.normalization(res) < 1e-9:
            break
        k += 1
    
    print("Gauss-Seidel's method")
    print("time:", time.time() - startTime)
    print("iterations:", k)
    print()
    return time.time() - startTime
# C task 
a1_D = 3
a2_D = a3_D = -1 
# D task 

def luFactoryzationMethod(A, b):
    # Initialize variables and matrices
    startTime = time.time()
    n = len(A)

    matA = meth.matcopy(A)
    mat_l = meth.diagToSquare(meth.vecones(n))
    mat_u = meth.matzeros(n, n)

    vecB = meth.veccopy(b)
    vecX = meth.vecones(n)
    vecY = meth.veczeros(n)

    # 1. Perform LU decomposition: LUx = b
    for j in range(n):
        # Calculate U elements
        for i in range(j + 1):
            mat_u[i][j] += matA[i][j]
            for k in range(i):
                mat_u[i][j] -= mat_l[i][k] * mat_u[k][j]
        
        # Calculate L elements
        for i in range(j + 1, n):
            for k in range(j):
                mat_l[i][j] -= mat_l[i][k] * mat_u[k][j]
            mat_l[i][j] += matA[i][j]
            mat_l[i][j] /= mat_u[j][j]

    # 3. Solve Ly = b
    for i in range(n):
        value = vecB[i]
        for j in range(i):
            value -= mat_l[i][j] * vecY[j]
        vecY[i] = value / mat_l[i][i]

    # 4. Solve Ux = y
    for i in range(n - 1, -1, -1):
        value = vecY[i]
        for j in range(i + 1, n):
            value -= mat_u[i][j] * vecX[j]
        vecX[i] = value / mat_u[i][i]

    # Calculate residual
    res = [vecB[i] - sum(matA[i][j] * vecX[j] for j in range(n)) for i in range(n)]

    # Print results
    print("LU method")
    print('Time taken:', time.time() - startTime)
    print("Residual norm:", meth.normalization(res))
    print()

    return time.time() - startTime

# Using:
# A task
A = createMatrix(a1,a2,a3,N)
b = createVector(N,f)
# B task
jacobiMethod(A,b)
gaussSeidelMethod(A,b)
# C task
C = createMatrix(a1_D,a2_D,a3_D,N)
#jacobiMethod(C,b) #error
#gaussSeidelMethod(C,b) #error

# D task
luFactoryzationMethod(C,b)

# E task
N = [100, 500, 1000, 1200]
timeJacobi = []
timeGauss = []
timeLU = []

for n in N:
    print("Size:", n)
    matrix_A = createMatrix(a1,a2,a3,n)
    vector_b = createVector(n,f)

    timeJacobi.append(jacobiMethod(matrix_A, vector_b))
    timeGauss.append(gaussSeidelMethod(matrix_A, vector_b))
    timeLU.append(luFactoryzationMethod(matrix_A, vector_b))

plt.plot(N, timeJacobi, label="Jacobi", color="red")
plt.plot(N, timeGauss, label="Gauss-Seidel", color="blue")
plt.plot(N, timeLU, label="LU", color="green")
plt.legend()
plt.grid(True)
plt.ylabel('Czas (s)')
plt.xlabel('Liczba niewiadomych')
plt.title('Zależność czasu od liczby niewiadomych')
plt.show()
