from CVA import calc_cva, make_swap_portfolio
from contexttimer import Timer
from datetime import date
import locale
locale.setlocale(locale.LC_ALL, '')

a = 0.376739
sigma = 0.0209835
todaysDate_dt=date(day=26,month=12,year=2016);

from datetime import date
a=0.376739
sigma=0.0209835
todaysDate_dt=date(day=26,month=12,year=2016);
crvTodaydates_dt=[date(day=26,month=12,year=2016),
               date(day=30,month=6,year=2017),
               date(day=30,month=7,year=2017),
               date(day=29,month=8,year=2017),
               date(day=30,month=9,year=2017),
               date(day=30,month=10,year=2017),
               date(day=28,month=11,year=2017),
               date(day=30,month=12,year=2017),
               date(day=30,month=1,year=2018),
               date(day=27,month=2,year=2018),
               date(day=30,month=3,year=2018),
               date(day=30,month=4,year=2018),
               date(day=29,month=5,year=2018),
               date(day=30,month=6,year=2018),
               date(day=30,month=12,year=2018),
               date(day=30,month=12,year=2019),
               date(day=29,month=12,year=2020),
               date(day=31,month=12,year=2021),
               date(day=30,month=12,year=2022),
               date(day=30,month=12,year=2023),
               date(day=30,month=12,year=2024),
               date(day=30,month=12,year=2025),
               date(day=29,month=12,year=2026),
               date(day=30,month=12,year=2027),
               date(day=30,month=12,year=2028),
               date(day=29,month=12,year=2031),
               date(day=30,month=12,year=2036),
               date(day=30,month=12,year=2041),
               date(day=30,month=12,year=2046),
               date(day=30,month=12,year=2051),
               date(day=30,month=12,year=2056),
               date(day=30,month=12,year=2061),
               date(day=31,month=12,year=2066)]
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

one_swap_portfolio = make_swap_portfolio(end_year=2021, months=[12], days=[26])

# Starting Example (1)
# swap_portfolio  = large_swap_portfolio

# Small Portfolio (2)
small_swap_portfolio = make_swap_portfolio(
    end_year = 2023, end_month = 3, days = [20])
swap_portfolio  = small_swap_portfolio

# Large Portfolio  (3)
large_swap_portfolio = make_swap_portfolio(
    end_year=2031, end_month=12, days=[1, 5, 10, 15, 20])

# swap_portfolio = large_swap_portfolio

Nsim = 1000
with Timer() as t:
    # Insert Dask here for calc_cva
    CVA_results = [(calc_cva(dt_startdate=start_date,
                             dt_maturitydate=maturity_date,
                             Nsim = Nsim,
                             crvTodaydates_dt=crvTodaydates_dt,
                             crvTodaydf = crvTodaydf,
                             todaysDate_dt = todaysDate_dt,
                             sigma = sigma),maturity_date) for (start_date, maturity_date) in swap_portfolio]
    for (CVA,maturity_date) in CVA_results:
        print(
            "CVA: {CVA} for swap maturing {year}/{month}/{day}".format(
                CVA=locale.currency(
                    CVA,
                    grouping=True),
                year  = maturity_date.year,
                month = maturity_date.month,
                day   = maturity_date.day))
