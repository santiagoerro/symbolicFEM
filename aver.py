import numpy as np
import matplotlib.pyplot as plt

a = np.linspace(0.01, 300, 500)
massGeneral = (2*a**3*np.exp(4*a) + 4*a**3*np.exp(3*a) + 24*a**3*np.exp(2*a) + 4*a**3*np.exp(a) + 2*a**3 - 15*a**2*np.exp(4*a) - 24*a**2*np.exp(3*a) + 24*a**2*np.exp(a) + 15*a**2 + 36*a*np.exp(4*a) - 36*a*np.exp(3*a) - 36*a*np.exp(a) + 36*a - 18*np.exp(4*a) + 36*np.exp(3*a) - 36*np.exp(a) + 18)/(6*a**3*(a**2*np.exp(4*a) - 2*a**2*np.exp(2*a) + a**2 - 4*a*np.exp(4*a) + 8*a*np.exp(3*a) - 8*a*np.exp(a) + 4*a + 4*np.exp(4*a) - 16*np.exp(3*a) + 24*np.exp(2*a) - 16*np.exp(a) + 4))
massBigA = (2*a**3 - 15*a**2 + 36*a - 18)/(6*a**3*(a**2- 4*a + 4))
massSmallA = 1/105 - a**2/3150 + 149*a**4/14553000 - 361*a**6/1135134000


plt.figure()
plt.plot(a,massGeneral)
plt.plot(a,massSmallA)
plt.plot(a,massBigA)
plt.ylim([-0.002,0.015])

plt.show()