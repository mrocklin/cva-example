from QuantLib import *
import numpy as np
from math import *
from contexttimer import Timer


# indexes definitions

def A(t, T, crvToday, sigma):
    evaldate = Settings.instance().evaluationDate
    forward = crvToday.forwardRate(t, t, Continuous, NoFrequency).rate()
    value = B(t, T) * forward - 0.25 * sigma * \
        B(t, T) * sigma * B(t, T) * B(0.0, 2.0 * t)
    return exp(value) * crvToday.discount(T) / crvToday.discount(t)


def B(t, T):
    a = 0.376739
    return (1.0 - exp(-a * (T - t))) / a


def gamma(t, crvToday, sigma):
    a = 0.376739
    forwardRate = crvToday.forwardRate(t, t, Continuous, NoFrequency).rate()
    temp = sigma * (1.0 - exp(-a * t)) / a
    return (forwardRate + 0.5 * temp * temp)


def gamma_v(t, crvToday, sigma):  # vectorized version of gamma(t), this is not vectorized
    res = np.zeros(len(t))
    for i in range(len(t)):
        res[i] = gamma(t[i], crvToday, sigma)
    return res


def make_index(term="6m"):
    forecastTermStructure = RelinkableYieldTermStructureHandle()
    index = Euribor(Period(term), forecastTermStructure)
    return (index, forecastTermStructure)

# swap 1 definition


def make_swap(maturity, index, startDate):
    fixedSchedule = Schedule(startDate, maturity, Period("6m"), TARGET(),
                             ModifiedFollowing, ModifiedFollowing,
                             DateGeneration.Forward, False)
    floatingSchedule = Schedule(startDate, maturity, Period("6m"), TARGET(),
                                ModifiedFollowing, ModifiedFollowing,
                                DateGeneration.Forward, False)
    swap = VanillaSwap(VanillaSwap.Receiver, 10000000, fixedSchedule, 0.02,
                       Actual360(), floatingSchedule, index, 0, Actual360())
    return (swap, floatingSchedule)


def calc_cva(swap, floatingSchedule, index, Nsim, forecastTermStructure,
             crvTodaydates, crvTodaydf, todaysDate, sigma):

    # Insert dask for make_curves
    (crvToday, npvMat, crvMat, rmean, Dates, T) = make_curves(
            crvTodaydates=crvTodaydates, crvTodaydf=crvTodaydf,
            todaysDate=todaysDate, sigma=sigma, Nsim=Nsim)

    for iT in range(len(T)):
        Settings.instance().evaluationDate = Dates[iT]
        allDates = list(floatingSchedule)
        fixingdates = [index.fixingDate(floatingSchedule[iDate])
                       for iDate in range(len(allDates))
                       if str(index.fixingDate(floatingSchedule[iDate]))
                       <= str(Dates[iT])]
        if fixingdates:
            for date in fixingdates[:-1]:
                try:
                    index.addFixing(date, 0.0)
                except:
                    pass
            try:
                index.addFixing(fixingdates[-1], rmean[iT])
            except:
                pass
        discountTermStructure = RelinkableYieldTermStructureHandle()
        swapEngine = DiscountingSwapEngine(discountTermStructure)
        swap.setPricingEngine(swapEngine)

        for iSim in range(Nsim):
            crv = crvMat[iSim][iT]
            discountTermStructure.linkTo(crv)
            forecastTermStructure.linkTo(crv)
            npvMat[iSim][iT] = swap.NPV()

    npvMat = np.array(npvMat)
    npv = npvMat[0, 0]
    # replace negative values with 0
    npvMat[npvMat < 0] = 0
    EE = np.mean(npvMat, axis=0)
    S = 0.05
    R = 0.4

    # Calc CVA
    cva_sum = 0
    for i in range(len(T) - 1):
        cva_sum += 0.5 * crvToday.discount(T[i + 1]) * (EE[i] + EE[i + 1]) * (
            exp(-S * T[i] / (1.0 - R)) - exp(-S * T[i + 1] / (1.0 - R)))
    CVA = (1.0 - R) * cva_sum
    return CVA


def make_curves(crvTodaydates, crvTodaydf, todaysDate, sigma, Nsim):
    crvToday = DiscountCurve(crvTodaydates, crvTodaydf, Actual360(), TARGET())
    # crvToday=FlatForward(todaysDate,0.0121,Actual360())
    a = 0.376739

    r0 = forwardRate = crvToday.forwardRate(0, 0, Continuous, NoFrequency).rate()
    months = range(3, 12 * 5 + 1, 3)
    sPeriods = [str(month) + "m" for month in months]
    print(sPeriods)
    Dates = [todaysDate] + [todaysDate + Period(s) for s in sPeriods]
    T = [0] + [Actual360().yearFraction(todaysDate, Dates[i])
               for i in range(1, len(Dates))]
    T = np.array(T)
    rmean = (r0 * np.exp(-a * T) + gamma_v(T, crvToday, sigma) -
             gamma(0, crvToday, sigma) * np.exp(-a * T))

    np.random.seed(1)
    stdnorm = np.random.standard_normal(size=(Nsim, len(T) - 1))

    rmat = np.zeros(shape=(Nsim, len(T)))
    rmat[:, 0] = r0
    for iSim in range(Nsim):
        for iT in range(1, len(T)):
            mean = (rmat[iSim, iT - 1]
                  * exp(-a * (T[iT] - T[iT - 1]))
                  + gamma(T[iT], crvToday, sigma)
                  - gamma(T[iT - 1], crvToday, sigma)
                  * exp(-a * (T[iT] - T[iT - 1]))
            var = (0.5 * sigma * sigma / a *
                  (1 - exp(-2 * a * (T[iT] - T[iT - 1]))))
            rmat[iSim, iT] = mean + stdnorm[iSim, iT - 1] * sqrt(var)

    with Timer() as t:
        crvMat = [[0 for i in range(len(T))] for iSim in range(Nsim)]
        npvMat = [[0 for i in range(len(T))] for iSim in range(Nsim)]

        for row in crvMat:
            row[0] = crvToday

        for iT in range(1, len(T)):
            for iSim in range(Nsim):
                crvDate = Dates[iT]
                crvDates = [crvDate] + [crvDate + Period(k, Years)
                                        for k in range(1, 21)]
                rt = rmat[iSim, iT]
                #if (rt<0): rt=0
                crvDiscounts = [1.0] + [A(T[iT], T[iT] + k, crvToday, sigma)
                                        * exp(-B(T[iT], T[iT] + k) * rt) for k in range(1, 21)]
                crvMat[iSim][iT] = DiscountCurve(crvDates, crvDiscounts, Actual360(), TARGET())
        # print("time for curve creation: ", t.elapsed)
    return (crvToday, npvMat, crvMat, rmean, Dates, T)


def plot_cva(T, EE, rmat, rmean):
    # plot(T,EE)
    nn  # title('Expected Exposure')
    #pxlabel('Time in years')
    # plot(T,np.mean(rmat,axis=0))
    # plot(T,rmean)
    # plot(T,[npvMat[0,0]]*len(T))
    # show()
