from  QuantLib import *
import numpy as np
from math import *


#indexes definitions

def A(t,T, crvToday,sigma):
    evaldate=Settings.instance().evaluationDate
    forward = crvToday.forwardRate(t, t,Continuous, NoFrequency).rate()
    value = B(t,T)*forward - 0.25*sigma*B(t,T)*sigma*B(t,T)*B(0.0,2.0*t);
    return exp(value)*crvToday.discount(T)/crvToday.discount(t);

def B(t,T):
    a=0.376739
    return (1.0-exp(-a*(T-t)))/a;

def gamma(t, crvToday, sigma):
    a=0.376739
    forwardRate =crvToday.forwardRate(t, t, Continuous, NoFrequency).rate()
    temp = sigma*(1.0 - exp(-a*t))/a
    return (forwardRate + 0.5*temp*temp)

def gamma_v(t, crvToday, sigma): #vectorized version of gamma(t), this is not vectorized
    res=np.zeros(len(t))
    for i in xrange(len(t)):
        res[i]=gamma(t[i], crvToday, sigma)
    return res


def make_index(term="6m"):
    forecastTermStructure = RelinkableYieldTermStructureHandle()
    index = Euribor(Period(term),forecastTermStructure)
    return (index,forecastTermStructure)

#swap 1 definition
def make_swap(maturity, index, startDate):
    fixedSchedule = Schedule(startDate, maturity,Period("6m"), TARGET(),ModifiedFollowing,ModifiedFollowing,DateGeneration.Forward, False)
    floatingSchedule = Schedule(startDate, maturity,Period("6m"),TARGET() ,ModifiedFollowing,ModifiedFollowing,DateGeneration.Forward, False)
    swap1 = VanillaSwap(VanillaSwap.Receiver, 10000000,fixedSchedule,0.02, Actual360(),floatingSchedule, index, 0,Actual360())  #0.01215
    return (swap1, floatingSchedule)

def calc_cva(swap, T, floatingSchedule, index, Nsim, crvToday, crvMat, rmean, npvMat, Dates, forecastTermStructure):
    for iT in xrange(len(T)):
        Settings.instance().evaluationDate=Dates[iT]
        allDates= list(floatingSchedule)
        fixingdates=[index.fixingDate(floatingSchedule[iDate]) for iDate in xrange(len(allDates)) if index.fixingDate(floatingSchedule[iDate])<=Dates[iT]]
        if fixingdates:
            for date in fixingdates[:-1]:
                try:index.addFixing(date,0.0)
                except:pass
            try:index.addFixing(fixingdates[-1],rmean[iT])
            except:pass
        discountTermStructure = RelinkableYieldTermStructureHandle()
        swapEngine = DiscountingSwapEngine(discountTermStructure)
        swap.setPricingEngine(swapEngine)

        for iSim in xrange(Nsim):
            crv=crvMat[iSim][iT]
            discountTermStructure.linkTo(crv)
            forecastTermStructure.linkTo(crv)
            npvMat[iSim][iT]=swap.NPV()

    npvMat=np.array(npvMat)
    npv=npvMat[0,0]
    #replace negative values with 0
    npvMat[npvMat<0]=0
    EE=np.mean(npvMat,axis=0)
    S=0.05
    R=0.4

    ## Calc CVA
    cva_sum=0
    for i in xrange(len(T)-1):
        cva_sum +=0.5*crvToday.discount(T[i+1])*(EE[i]+EE[i+1])*(exp(-S*T[i]/(1.0-R))-exp(-S*T[i+1]/(1.0-R)))
    CVA=(1.0-R)*cva_sum
    return CVA
