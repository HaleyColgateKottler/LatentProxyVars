import numpy as np
import scipy as sp

def sigmaSim(p, k, communality):
    # p number of variables
    # k number of factors
    
    A = np.zeros((p,k))
    
    # need row entries to sum to k-1
    A[:,0] = np.random.randint(0, k, p)
    
    for i in range(p):
        for j in range(1, k-1):
            current_sum = sum(A[i,:])
            if current_sum == k-1:
                A[i,j] = 0
            else:
                a = np.random.randint(0, k-current_sum)
                A[i,j] = a
    A[:,k-1] = (k-1)*np.ones((1,p)) - np.sum(A,1)
    
    # Add normal deviation
    c = .1*np.random.randint(7, 10, k)
    x1 = np.random.normal(size = (p, k))
    d = 1/np.linalg.norm(x1, ord=2, axis=1)
    
    Y = A*c + d.reshape(p,1)*x1*np.sqrt(1-c**2)
    
    # Apply skewing function
    Y2 = Y + abs(Y) + 0.2
    Y3 = abs(Y) + 0.2
    Z = (1.2/2.2) * (Y*Y2) / Y3
    
    g = 1/np.linalg.norm(Z, ord=2, axis=1)
    
    # Scale to set communality
    if communality == 1:
        B1 = np.diag(.1*np.random.randint(2,5,p))
    elif communality == 2:
        B1 = np.diag(.1*np.random.randint(2,9,p))
    elif communality == 3:
        B1 = np.diag(.1*np.random.randint(6,9,p))
    else:
        B1 = np.zeros(p,p)
        
    B2 = np.identity(p) - B1
    
    # Final factor loading matrix
    lambda_common = np.dot(np.sqrt(B1), g).reshape(p,1)*Z
    lambda_unique = np.sqrt(B2)
    
    sigma = np.dot(lambda_common, np.transpose(lambda_common)) + np.dot(lambda_unique, np.transpose(lambda_unique))
    
    return(sigma, lambda_common, lambda_unique)

p = 3
k = 5
communality = 1
print(sigmaSim(p, k, communality))
