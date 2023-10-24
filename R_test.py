import rpy2.rinterface.tests
import unittest

# the verbosity level can be increased if needed
tr = unittest.TextTestRunner(verbosity = 1)
suite = rpy2.rinterface.tests.suite()
tr.run(suite)
