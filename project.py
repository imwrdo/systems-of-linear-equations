
from matplotlib import pyplot as plt
import methods as meth
import initMethods as inmeth 

# Task A 
index_number = 201267
d = index_number % 10 #7
c = int(index_number/10)%10 #6
e = int(index_number/100)%10 #2
f = int(index_number/1000)%10 #1
a1 = 5 + e
a2 = a3 = -1
N = 967 #9cd
            
# Task C 
a1_D = 3
a2_D = a3_D = -1 

# Using:
# Task A
A = inmeth.createMatrix(a1,a2,a3,N)
b = inmeth.createVector(N,f)

# Task B
jacobi_time_B, jacobi_residuals_B = inmeth.jacobiMethod(A, b)
gauss_time_B, gauss_residuals_B = inmeth.gaussSeidelMethod(A, b)

plt.figure(figsize=(10, 6))
plt.semilogy(range(len(jacobi_residuals_B)), jacobi_residuals_B, label="Jacobi")
plt.semilogy(range(len(gauss_residuals_B)), gauss_residuals_B, label="Gauss-Seidel")
plt.axhline(y=1e-9,color='g',linestyle='--',label="Residual border")
plt.xlabel("Iterations")
plt.ylabel("Residual norm")
plt.title("Jacobi and Gauss-Seidel")
plt.legend()
plt.grid(True)
plt.show()

# Task C

C = inmeth.createMatrix(a1_D,a2_D,a3_D,N)
jacobi_time_C, jacobi_residuals_C = inmeth.jacobiMethod(C, b,500)
gauss_time_C, gauss_residuals_C = inmeth.gaussSeidelMethod(C, b,500)

plt.figure(figsize=(10, 6))
plt.semilogy(range(len(jacobi_residuals_C)), jacobi_residuals_C, label="Jacobi method")
plt.semilogy(range(len(gauss_residuals_C)), gauss_residuals_C, label="Gauss-Seidel's method")
plt.xlabel("Iterations")
plt.ylabel("Residual norm")
plt.title("Jacobi and Gauss-Seidel's methods")
plt.legend()
plt.grid(True)
plt.show()


# Task D
inmeth.luFactoryzationMethod(C,b)

# Task E         
N = [100, 500, 1000, 1200,1500]
timeJacobi = []
timeGauss = []
timeLU = []

for n in N:
    print("Size:", n)
    matrix_A = inmeth.createMatrix(a1,a2,a3,n)
    vector_b = inmeth.createVector(n,f)
    
    jacobi_time, _ = inmeth.jacobiMethod(matrix_A, vector_b)
    gauss_time, _ = inmeth.gaussSeidelMethod(matrix_A, vector_b)   
    LU_time = inmeth.luFactoryzationMethod(matrix_A, vector_b)
    
    timeJacobi.append(jacobi_time)
    timeGauss.append(gauss_time)
    timeLU.append(LU_time)

plt.plot(N, timeJacobi, label="Jacobi", color="red")
plt.plot(N, timeGauss, label="Gauss-Seidel", color="blue")
plt.plot(N, timeLU, label="LU", color="green")
plt.legend()
plt.grid(True)
plt.ylabel('Time taken(s)')
plt.xlabel('Matrix size')
plt.title('Zależność czasu od liczby niewiadomych')
plt.show()
