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

import pickle
f = open("file.txt", 'wb')
pickle.dump(m, f)
del f
f2 = open("file.txt", 'rb')
val = pickle.load(f2)
print (val)
