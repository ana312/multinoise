#Returns Graph of First Moment against 1st-moment parameter - x_0 in Shiino  
#

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


#a=1=theta
def R(m,nu,beta,a=1,theta=1):
	z_eval = lambda x : (x**2/nu+1/beta)**((a-theta-nu**-1-(nu/beta))/(2*nu))*np.exp(-nu*x**2+(m*theta)*np.sqrt(beta*nu)*np.arctan(np.sqrt(beta/nu)*x))
	z = quad(z_eval,-5.,5.)
	r_eval = lambda x : x*(x**2/nu+1/beta)**((a-theta-nu**-1-(nu/beta))/(2*nu))*np.exp(-nu*x**2+(m*theta)*np.sqrt(beta*nu)*np.arctan(np.sqrt(beta/nu)*x))*(z[0])**-1
	r = quad(r_eval,-5.,5.)
	return r[0]

mrange=100
jrange=6
madmis=np.zeros(mrange)
mstore=np.zeros((jrange,mrange))

for k in range(jrange):
	nu=.75+k
	nu=round(nu,1)
	for i in range(mrange):
		m=-5+0.1*i
		madmis[i]=m
		mstore[k,i]=R(m,nu,10.)
	plt.plot(madmis,mstore[k,:], label=nu)



plt.plot(madmis,madmis, '--')


plt.xlabel('First Moment Parameter')
plt.ylabel('Actual First Moment')
plt.xlim([-0.75,0.75])
plt.ylim([-0.75,0.75])
plt.title('Self-Consistency')
plt.legend(loc='lower right')
plt.show()
