#
# test iminuit with gradients
# also test iminuit pickle
#

import numpy as np
import pickle
from iminuit import Minuit

class testing(object):

    def __init__(self):
        self.x = np.arange(0.,100.,1.)
        self.yerr = np.random.normal(size=100)

    def yfunc(self,par):
        y = par[0] + par[1]*self.x + par[2]*self.x**2
        return y

    def mkdata(self):
        truepar = np.array([4.0,-0.4,2.0])
        self.ydata = self.yfunc(truepar) + 1.*self.yerr

    def chisq(self,par):

        # calculate chi2
        dy = self.ydata-self.yfunc(par)
        chisquared = np.sum(dy**2)
        #print("chisq = ",chisquared)
        return chisquared

    def calcgrad(self,par):

        gradArr = np.zeros(par.shape[0])
        gradArr[0] = -2. * np.sum( (self.ydata-self.yfunc(par)) )
        gradArr[1] = -2. * np.sum( (self.ydata-self.yfunc(par)) * self.x )
        gradArr[2] = -2. * np.sum( (self.ydata-self.yfunc(par)) * self.x**2 )
        return gradArr


#chisq.errordef = Minuit.LEAST_SQUARES

startpar = [3.,0.,1.]
names = ["a","b","c"]

atest = testing()
atest.mkdata()

gMinuit = Minuit(atest.chisq,startpar,name=names,grad=atest.calcgrad)
#gMinuit = Minuit(chisq,startpar,name=names)
gMinuit.print_level = 1 # 3 has lots more details
gMinuit.strategy = 1

# set limits...
for aname in names:
    gMinuit.limits[aname] = (-10,10.)
    gMinuit.errors[aname] = 0.01

# fix a parameter
gMinuit.fixed['a'] = True

# print setup
gMinuit.init_params

# do fit
gMinuit.migrad()

# print out results
print(gMinuit.values)
print(gMinuit.fmin)
print(gMinuit.errors)
print(gMinuit.covariance)

# release a
gMinuit.fixed['a'] = False

# do fit
gMinuit.migrad()

# print out results
print(gMinuit)
print(gMinuit.values)
print(gMinuit.fmin)
print(gMinuit.errors)
print(gMinuit.covariance)

# pickle
pickle.dump(gMinuit,open('output/test-iminuit.pkl','wb'))

pMinuit = pickle.load(open('output/test-iminuit.pkl','rb'))

print(pMinuit.values)
print(pMinuit)
