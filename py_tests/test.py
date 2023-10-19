import ASCsoft.bla as bla
import numpy as np

m = bla.Matrix(3,3)
for p in [(i,j,i+j) for i in range(3) for j in range(3)]:
  (r,c,v) = p
  m[r,c] = v

print(m)
print(3*m)
print(m*m)

print(np.array(m))