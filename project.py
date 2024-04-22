from math import sin
import time
from matplotlib import pyplot
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
def jacobiMethod(A, b, tol=1e-9, max_iter=1000):
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

# A task
A = createMatrix(a1,a2,a3,N)
b = createVector(N,f)
# B task
jacobiMethod(A,b)
gaussSeidelMethod(A,b)









