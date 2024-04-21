def matcopy(matrix):
    num_rows = len(matrix)
    num_cols = len(matrix[0]) if matrix else 0
    new_matrix = [[0] * num_cols for _ in range(num_rows)]
    for i in range(num_rows):
        for j in range(num_cols):
            new_matrix[i][j] = matrix[i][j]
    
    return new_matrix

def diag(matrix):
    if len(matrix) != len(matrix[0]):
        raise ValueError("Matrix must be square to extract diagonal elements.")
    diagonal = [matrix[i][i] for i in range(len(matrix))]
    return diagonal

def matadd(matrix1, matrix2):
    if len(matrix1) != len(matrix2) or len(matrix1[0]) != len(matrix2[0]):
        raise ValueError("Matrices must have the same dimensions for addition.")
    
    num_rows = len(matrix1)
    num_cols = len(matrix1[0])
    
    result_matrix = [[0] * num_cols for _ in range(num_rows)]
    
    for i in range(num_rows):
        for j in range(num_cols):
            result_matrix[i][j] = matrix1[i][j] + matrix2[i][j]
    
    return result_matrix

def matsub(matrix1, matrix2):
   
    if len(matrix1) != len(matrix2) or len(matrix1[0]) != len(matrix2[0]):
        raise ValueError("Matrices must have the same dimensions for subtraction.")
    
    num_rows = len(matrix1)
    num_cols = len(matrix1[0])
    
    result_matrix = [[0] * num_cols for _ in range(num_rows)]
    
    for i in range(num_rows):
        for j in range(num_cols):
            result_matrix[i][j] = matrix1[i][j] - matrix2[i][j]
    
    return result_matrix

def matzeros(rows, cols):
    zero_matrix = [[0] * cols for _ in range(rows)]
    return zero_matrix
