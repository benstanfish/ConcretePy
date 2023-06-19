import numpy as np
import matplotlib.pyplot as plt

a = np.array([2,1,5,7,4,6,8,14,10,9,18,20,22])
# plt.plot(a)
# plt.show()

b = np.array([[-2,-2],[-4,3],[1,6],[3,-1]])
c = np.copy(b)
# print(c.shape[0])
# print(c.size)
x = b[:,0]
y = b[:,1]

# plt.plot(x,y)
# plt.show()
print(c)
d = np.array([b[0]])
print(d)
