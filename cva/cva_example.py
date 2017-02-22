from  QuantLib import *
import numpy as np
from math import *
from CVA import *
from contexttimer import Timer

Nsim=1000
a=0.376739
sigma=0.0209835
todaysDate=Date(26,12,2013);
Settings.instance().evaluationDate=todaysDate;
crvTodaydates=[Date(26,12,2013),
               Date(30,6,2014),
               Date(30,7,2014),
               Date(29,8,2014),
               Date(30,9,2014),
               Date(30,10,2014),
               Date(28,11,2014),
               Date(30,12,2014),
               Date(30,1,2015),
               Date(27,2,2015),
               Date(30,3,2015),
               Date(30,4,2015),
               Date(29,5,2015),
               Date(30,6,2015),
               Date(30,12,2015),
               Date(30,12,2016),
               Date(29,12,2017),
               Date(31,12,2018),
               Date(30,12,2019),
               Date(30,12,2020),
               Date(30,12,2021),
               Date(30,12,2022),
               Date(29,12,2023),
               Date(30,12,2024),
               Date(30,12,2025),
               Date(29,12,2028),
               Date(30,12,2033),
               Date(30,12,2038),
               Date(30,12,2043),
               Date(30,12,2048),
               Date(30,12,2053),
               Date(30,12,2058),
               Date(31,12,2063)]
crvTodaydf=[1.0,
            0.998022,
            0.99771,
            0.99739,
            0.997017,
            0.996671,
            0.996337,
            0.995921,
            0.995522,
            0.995157,
            0.994706,
            0.994248,
            0.993805,
            0.993285,
            0.989614,
            0.978541,
            0.961973,
            0.940868,
            0.916831,
            0.890805,
            0.863413,
            0.834987,
            0.807111,
            0.778332,
            0.750525,
            0.674707,
            0.575192,
            0.501258,
            0.44131,
            0.384733,
            0.340425,
            0.294694,
            0.260792
            ]

crvToday=DiscountCurve(crvTodaydates,crvTodaydf,Actual360(),TARGET())
#crvToday=FlatForward(todaysDate,0.0121,Actual360())

r0=forwardRate =crvToday.forwardRate(0,0, Continuous, NoFrequency).rate()
months=range(3,12*5+1,3)
sPeriods=[str(month)+"m" for month in months]
print sPeriods
Dates=[todaysDate]+[todaysDate+Period(s) for s in sPeriods]
T=[0]+[Actual360().yearFraction(todaysDate,Dates[i]) for i in xrange(1,len(Dates))]
T=np.array(T)
rmean=r0*np.exp(-a*T)+ gamma_v(T, crvToday,sigma) -gamma(0,crvToday, sigma)*np.exp(-a*T)
# rvar=sigma*sigma/2.0/a*(1.0-np.exp(-2.0*a*T))
# rstd=np.sqrt(rvar)
np.random.seed(1)
stdnorm = np.random.standard_normal(size=(Nsim,len(T)-1))

rmat=np.zeros(shape=(Nsim,len(T)))
rmat[:,0]=r0
for iSim in xrange(Nsim):
    for iT in xrange(1,len(T)):
        mean=rmat[iSim,iT-1]*exp(-a*(T[iT]-T[iT-1]))+gamma(T[iT], crvToday,sigma)-gamma(T[iT-1],crvToday,sigma)*exp(-a*(T[iT]-T[iT-1]))
        var=0.5*sigma*sigma/a*(1-exp(-2*a*(T[iT]-T[iT-1])))
        rmat[iSim,iT]=mean+stdnorm[iSim,iT-1]*sqrt(var)

# bonds=np.zeros(shape=rmat.shape)

# #E(E(exp(rt)|ti) test
# for iSim in xrange(Nsim):
#     for iT in xrange(1,len(T)):
#         bonds[iSim,iT]=bonds[iSim,iT-1]+rmat[iSim,iT]*(T[iT]-T[iT-1])

# bonds=-bonds;
# bonds=np.exp(bonds)

#bondsmean=np.mean(bonds,axis=0)
#plot(T,bondsmean)
#plot(T,[crvToday.discount(T[iT]) for iT in xrange(len(T))])
#show()


with Timer() as t:
    crvMat= [ [ 0 for i in xrange(len(T)) ] for iSim in range(Nsim) ]
    npvMat= [ [ 0 for i in xrange(len(T)) ] for iSim in range(Nsim) ]

    for row in crvMat:
        row[0]=crvToday

    for iT in xrange(1,len(T)):
        for iSim in xrange(Nsim):
            crvDate=Dates[iT];
            crvDates=[crvDate]+[crvDate+Period(k,Years) for k in xrange(1,21)]
            rt=rmat[iSim,iT]
            #if (rt<0): rt=0
            crvDiscounts=[1.0]+[A(T[iT],T[iT]+k, crvToday, sigma)*exp(-B(T[iT],T[iT]+k)*rt) for k in xrange(1,21)]
            crvMat[iSim][iT]=DiscountCurve(crvDates,crvDiscounts,Actual360(),TARGET())
    print "time for curve creation: ", t.elapsed

with Timer() as t:
    maturity = Date(26,12,2018);
    startDate=Date(26,12,2013);
    (index,forecastTermStructure)=make_index()
    (swap1, floatingSchedule) = make_swap(maturity=maturity,
                                      startDate=startDate,
                                      index=index)
print "time for swap creation: ", t.elapsed

with Timer() as t:
    CVA = calc_cva(swap = swap1,
                   T = T,
                   floatingSchedule=floatingSchedule,
                   index=index,
                   Nsim = Nsim,
                   crvToday=crvToday,
                   crvMat=crvMat,
                   rmean= rmean,
                   npvMat = npvMat,
                   Dates = Dates,
                   forecastTermStructure=forecastTermStructure)
    print "CVA is: ", CVA
    print "time for CVA: ", t.elapsed

#print '\nEE:\n',EE
#print "\nnpv=",npv
#print "\nCVA=",CVA
#print "\nnpv",npvMat
#print '\nrmat:\n',rmat
#print '\nrmean:\n',rmean
#print '\nrstd:\n',rstd
#print '\n95% are in \n',zip(rmean-2*rstd,rmean+2*rstd)


def plot_cva(T, EE, rmat, rmean):
    #plot(T,EE)
    nn#title('Expected Exposure')
    #pxlabel('Time in years')
    #plot(T,np.mean(rmat,axis=0))
    #plot(T,rmean)
    #plot(T,[npvMat[0,0]]*len(T))
    #show()
