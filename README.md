# Indoor air infiltration rate estimator

This module estimates the air mixing/air change rate from AQI measurements. Data must be given
in a CSV file, with three data columns: time, indoor, and outdoor air quality observations.  The
timestep must be sufficiently small for the implied air change rate.  In addition, the estimator
assumes no filtering on the measurements.

The estimate computes three values, the initial value `c`, the time-constant rate `r`, and the
implied air-changes per hour `ach`.  These values correspond to the continous model

![equation](http://www.sciweavers.org/tex2img.php?im=jpg&fs=12&eq=y=c~\exp(rt))

and the discrete model

![equation](http://www.sciweavers.org/tex2img.php?im=jpg&fs=12&eq=y_t=y_{t-1}~(1-t_s~ach)+{1\over2}(x_{t-1}-x_{t})~t_s~ach)

where `t` is the time in hours, `x` is the outdoor AQI, and `y` is the indoor AQI.

# Example data

File: `test1.csv`
~~~
time,outdoor,indoor
0.000,40,133
0.167,40,128
0.333,40,86
0.500,40,76
0.667,40,67
0.833,40,61
1.000,40,56
1.167,40,53
1.333,40,51
1.500,40,48
1.667,40,47
1.833,40,45
2.000,40,42
2.167,40,43
2.333,40,41
~~~

# Example estimation

~~~
host% python3 -q
>>> import AQI
>>> A = AQI.estimate("test1.csv")
>>> print(repr(A),'\n-->',str(A))
AQI.estimate(csvfile='test1.csv',tlabel='time',xlabel='indoor',ylabel='outdoor',constrain=None,prec=0.001): 
--> {'r': -1.7933034171096698, 'c': 133.00692175930675, 'ach': 1.546985672682819}
~~~
