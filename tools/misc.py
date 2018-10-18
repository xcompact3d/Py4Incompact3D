"""
.. module:: misc
    :synopsis: For tools which don't fit anywhere else

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

def moving_avg(x, nsample=16):
	""" Compute the moving average over nsamples. 

  :param x: The input array/list-like object to compute moving average of.
  :param nsample: How many samples to compute average over.

  :type x: list
  :type nsample: int

  :returns: xavg - the moving average of x.
  :rtype: list
  """

	n = len(x)
	xavg = []

	for i in range(n):
		xavg.append(0)

		j0 = max(i - nsample // 2, 0)
		j1 = min(i + nsample // 2, n)
		for j in range(j0, j1):
			xavg[-1] += x[j]
		xavg[-1] /= (j1 - j0)

	return xavg
