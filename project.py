import methods

# 1 task 
index_number = 201267
e = d = index_number % 10
c = int(index_number/10)%10
d = int(index_number/100)%10
a1 = 5 + e
a2 = a3 = -1
N = 9 * c *  d

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
print(A)