import warnings
import csv
import numpy as np
import scipy.stats as st
import scipy.optimize as opt

class estimate:

	"""Estimate the air mixing/air change rate from AQI measurements

	Properties:

		csvfile (str)     the pathname of the CSV file containing the AQI measurements
		constrain (str)   type of constaint on the fit
		prec (float)      precision of timestep required
		tlabel (str)      CSV time column label
		xlabel (str)      CSV indoor air column label
		ylabel (str)      CSV outdoor air column label
		t (list)          time data vector
		x (list)          indoor AQI data vector
		y (list)          outdoor AQI data vector
		r (float)         mixing rate linear regression result
		c (float)         initial value linear regression result
		rvalue (float)    R-value of linear regression of mixing rate
		pvalue (float)    P-value of linear regression of mixing rate 
		stderr (float)    standard error of linear regression of mixing rate
		ach (float)       inferred air-changes per hour

	Description:

		This class contains the method for estimate the implied indoor air mixing rate from outdoor air.
		Two methods are used depending on whether the first observation is considered the initial value.
		If the initial value is known, then the air-change rate is estimated using a simple optimization
		to minimize the error.  If the initial value is not known, then a simple linear regression is
		performed.

	"""

	def __init__(self,csvfile,tlabel='time',xlabel='indoor',ylabel='outdoor',constrain=None,prec=1e-3):

		"""Create an estimate from AQI data

		Parameters:

			csvfile (str)     pathname of the CSV file containing the AQI measurements
			tlabel (str)      label of the time column (default is 'time')
			xlabel (str)      label of the indoor AQI column (default is 'indoor')
			ylabel (str)      label of the outdoor AQI column (default is 'outdoor')
			constrain (str)   constrain fit to initial condition (default None)

		"""		

		data = {}
		with open(csvfile,"r") as fh:
			reader = csv.reader(fh)
			header = []
			for row in reader:
				if not header:
					header = row
					for item in header:
						data[item] = []
				else:
					for n in range(len(header)):
						data[header[n]].append(float(row[n]))
		self.csvfile = csvfile
		self.constrain = constrain
		self.prec = prec
		self.tlabel = tlabel
		self.xlabel = xlabel
		self.ylabel = ylabel
		self.t = data[tlabel]
		self.x = data[xlabel]
		self.y = data[ylabel]
		t = np.array(self.t)
		x = np.array(self.x)
		y = np.array(self.y)
		if np.std(t[1:]-t[0:-1]) / np.mean(t) > prec:
			raise "timestep is not uniform to specified precision"
		ts = t[1]-t[0]
		def predict(t,x0,ach,ts,z=0):
			y = [x0]
			for n in range(len(t)-1):
				y.append(y[-1]*(1-ts*ach)+(z[n]+z[n+1])/2*ach*ts)
			return y
		if constrain == 'init': # initial value is contrained to observed value
			def objective(ach):
				z = predict(t,x[0],ach,ts,y) - x
				return sum(z*z)
		elif constrain:
			raise "invalid constraint"
		else: # initial value is not constrained to observed value
			dx = np.log(np.array(self.x)-np.array(self.y))
			self.r,self.c,self.rvalue,self.pvalue,self.stderr = st.linregress(t,dx)
			def objective(ach):
				z = predict(t,np.exp(self.c)+y[0],ach,ts,y) - x
				return sum(z*z)
		res = opt.minimize_scalar(objective,bounds=[0,10])
		if res.success:
			self.ach = res.x
			dx = np.log(np.array(predict(t,x[0],self.ach,ts,y))-np.array(self.y))
			self.r,self.c,self.rvalue,self.pvalue,self.stderr = st.linregress(t,dx)
			self.c = np.exp(self.c)+np.mean(self.y)
			if self.ach > 1/ts:
				warnings.warn("estimated ACH is too high for timestep")
			elif self.ach < 0:
				warnings.warn("estimated ACH is negative")
		else:
			self.r = self.c = self.ach = float('nan')

	def to_dict(self):
		return {
			'r': self.r,
			'c': self.c,
			'ach' : self.ach
		}

	def __repr__(self):
		"""Return the constructor for this object as a string"""
		if self.constrain:
			constrain = f"'self.constrain'"
		else:
			constrain = 'None'
		return f"{__name__}.estimate(csvfile='{self.csvfile}',tlabel='{self.tlabel}',xlabel='{self.xlabel}',ylabel='{self.ylabel}',constrain={constrain},prec={self.prec}):"

	def __str__(self):
		"""Return properties of this object as string of the dictionary"""
		return str(self.to_dict())

import unittest
class TestEstimate(unittest.TestCase):

	def test1_init(self):
		A = estimate("test1.csv",constrain='init')
		self.assertEqual(A.r,-1.7177908615058413)
		self.assertEqual(A.c,133.00663028695413)
		self.assertEqual(A.ach,1.4907421330551835)

	def test1_noinit(self):
		A = estimate("test1.csv").to_dict()
		self.assertEqual(A['r'],-1.7933034171096698)
		self.assertEqual(A['c'],133.00692175930675)
		self.assertEqual(A['ach'],1.546985672682819)

if __name__ == '__main__':
	unittest.main()

