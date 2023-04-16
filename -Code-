
# Importing libraries

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Loading condition check

print('Enter combined loading case: ')

# initiation of all the values and intake

H = float(input('Enter H = '))              # Lateral load
p = float(input('Enter P = '))              # Axial load
M = float(input('Enter M = '))              # Moment
L = int(input('Enter L =') )                # Length of pile
d = float(input('Enter d = '))              # Diameter of pile
n = int(input('Enter n =') )                # no of Sections the pile is divided into
f = float(input('Enter f =') )              # Grade of concrete
ei = 5000*(f**(1/2))                        # As per IS 456 code
v = float(input('Enter v = '))              # Poisson's Ratio of soil 
n = n+1

# Initial calculations for required values

de = (L/(n-1))                              # Length of pile Section
ip = (math.pi*d**4)/64                      # Moment of inertia of pile section

# Inputs for the conditions

K = float(input('Enter 0 if k is constant or 1 for varying k = ')) 
G = float(input('Enter 0 for fixed headed pile and 1 for free headed pile '))

# Declaration of embedded list to get a nxn matrix to store our values

A = []                                      # Main matrix storing the coefficients of deflections from the equations
B = []                                      # constant matrix storing the constant values of equations

# initiation of A and B

A = np.zeros((n,n))
B = np.zeros(n)

# Calculations to simplify the assignment

M1 = (M*(de**2))/(ei*ip)
H1 = (H*(de**3))/(ei*ip)

# calculating the 'a' value from Es

w = (ei*ip)/(de**2)
o = p-4*w
if K!=0:
    a=[]
    k=[]
    for i in range (n):
        C = float(input('C = '))        # Cohesion of soil
        es = 100*C
        k.append((0.65/d)*((((es*d**4)/(ei*ip))**(1/12))*(es/(1-v**2))))
        a.append((6*w) - (2*p) + (k[i]*(de**2)))
else:
    C = float(input('C = '))            # Cohesion of clay
    es = 100*C 
    k = (0.65/d)*((((es*d**4)/(ei*ip))**(1/12))*(es/(1-v**2)))
    a = (6*w) - (2*p) + (k*(de**2))

# Filling the coefficients of all intermediate node displacement equations 

for i in range(2,n-2):
    A[i][i-2] = A[i][i+2] = w
    A[i][i-1] = A[i][i+1] = o
    if K!=0:
        A[i][i]=a[i]
    else:
        A[i][i]=a


# Fill values based on head condition

if(G == 1):
    if(K != 0):
        A[0][0] = 2+(k[0]*(de**2))/w
        A[1][1] = (a[1]/w)-1
    else:
        A[0][0] = 2+(k*(de**2))/w
        A[1][1] = (a/w)-1
    A[0][1] = -4
    A[1][0] = A[1][2] = (p/w)-2
else:
    if(K != 0):
        A[0][0] = a[0]/w
        A[1][1] = a[1]/w+1
    else:
        A[0][0] = a/w
        A[1][1] = a/w+1
    A[0][1] = 2*o
    A[1][2] = A[1][0] = o
A[0][2] = 2
A[1][3] = 1
if(K != 0):
    A[n-2][n-2] = (a[n-2]/w) -1
    A[n-1][n-1] = (a[n-1]/w) + 2*o +4
else:
    A[n-2][n-2] = (a/w)-1
    A[n-1][n-1] = (a/w) + 2*o +4
A[n-1][n-3] = 2
A[n-2][n-1] = o+2
A[n-2][n-3] = o
A[n-2][n-4] = 1
A[n-1][n-2] =-4

# Calculate X inverse

z=np.linalg.inv(A)

# Filling Constant matrix



if G == 0:
    B[0] = H1
    B[1] = 0
else:
    B[0] = H1 + 2*M1 -p*M1
    B[1] = -M1

# Multiplying the constant matrix and inverse of coefficient matrix

res = z @ B

# To plot the results on the graph

if res[0]<0 and G == 0:
    res = res*(-1)
res_df=pd.DataFrame(res)
res_df.to_csv("res.csv")
# The Y values corresponding to the nodes are being stored in y

y = []
g = L                                           # To get the depth values of nodes
while(g>0):
    y.append(g)
    g = g - (L/n)
y=y[:-1]
print(y) 



# Printing the displacement values

print('Resultant lateral displacement of pile due to applied load = ',res_df)

# Plotting the displacement profile

plt.plot(res,y)
plt.xlabel("Displacemnt")
plt.ylabel("Depth")

if G == 1:
    plt.title("Displacement due to combined load on free pile")
else:
    plt.title("Displacement due to combined load on fixed head pile")
plt.show()
