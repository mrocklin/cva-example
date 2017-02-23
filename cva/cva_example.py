from  QuantLib import *
import numpy as np
from math import *
from CVA import *
from contexttimer import Timer
import locale
locale.setlocale( locale.LC_ALL, '' )


a=0.376739
sigma=0.0209835
todaysDate=Date(26,12,2016);
Settings.instance().evaluationDate=todaysDate;
crvTodaydates=[Date(26,12,2016),
               Date(30,6,2017),
               Date(30,7,2017),
               Date(29,8,2017),
               Date(30,9,2017),
               Date(30,10,2017),
               Date(28,11,2017),
               Date(30,12,2017),
               Date(30,1,2018),
               Date(27,2,2018),
               Date(30,3,2018),
               Date(30,4,2018),
               Date(29,5,2018),
               Date(30,6,2018),
               Date(30,12,2018),
               Date(30,12,2019),
               Date(29,12,2020),
               Date(31,12,2021),
               Date(30,12,2022),
               Date(30,12,2023),
               Date(30,12,2024),
               Date(30,12,2025),
               Date(29,12,2026),
               Date(30,12,2027),
               Date(30,12,2028),
               Date(29,12,2031),
               Date(30,12,2036),
               Date(30,12,2041),
               Date(30,12,2046),
               Date(30,12,2051),
               Date(30,12,2056),
               Date(30,12,2061),
               Date(31,12,2066)]
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




def make_swap_portfolio(years = None, end_year=2023, months = None, end_month=None, days=None):
    if years is None:
        years = range(2021,end_year+1)

    if months is None:
        if end_month is None:
            months = range(1,4)
        if type(end_month) is int:
            months = range(1,end_month+1)

    if days is None:
        days = [1,5,10,20]
    swaps = []
    for year in years:
        for month in months:
            for day in days:
                (index,forecastTermStructure)=make_index()
                maturity = Date(26,12,year)
                startDate = Date(26,12,2016)
                (swap, floating_schedule)  = make_swap(maturity  = maturity,
                                                       startDate = startDate,
                                                       index     = index)
                swaps.append((swap, floating_schedule, forecastTermStructure,index))

    print "Created portfolio of ", len(swaps), " swaps"
    return swaps



one_swap_portfolio   = make_swap_portfolio(end_year = 2021, months=[12], days = [26])

## Starting Example (1)
swap_portfolio  = large_swap_portfolio

## Small Portfolio (2)
# small_swap_portfolio = make_swap_portfolio(end_year = 2023, end_month = 3, days = [1,5])
# swap_portfolio  = small_swap_portfolio

## Large Portfolio  (3)
large_swap_portfolio = make_swap_portfolio(end_year = 2031, end_month = 12, days = [1,5,10,15,20])
swap_portfolio = large_swap_portfolio

Nsim=1000
with Timer() as t:
    ## Insert Dask here for calc_cva
    CVA_results = [(calc_cva(swap = swap,
                             floatingSchedule=floatingSchedule,
                             index=index,
                             Nsim = Nsim,
                             forecastTermStructure=forecastTermStructure,
                             crvTodaydates=crvTodaydates,
                             crvTodaydf = crvTodaydf,
                             todaysDate = todaysDate,
                             sigma = sigma),swap) for (swap,floatingSchedule, forecastTermStructure,index) in swap_portfolio]

    for (CVA,swap) in CVA_results:
             print "CVA is: {CVA} for swap maturiting {year}/{month}/{day}".format(CVA=locale.currency(CVA, grouping=True),
                                                                                   year=swap.maturityDate().year(),
                                                                                   month=swap.maturityDate().month(),
                                                                                   day=swap.maturityDate().dayOfMonth())
