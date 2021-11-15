import numpy as np
import matplotlib.pyplot as plt

def f(x1, x2):
  return x1**2 + x2**2

x1 = np.arange(-5.0, 5.0, 0.1)
x2 = np.arange(-5.0, 5.0, 0.1)

X1, X2 = np.meshgrid(x1, x2)
F = f(X1, X2)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('F')
ax.set_title('F(x1, x2)=x1**2+x2**2')

ax.plot_wireframe(X1, X2, F)

plt.savefig('3Dgraph.jpg')