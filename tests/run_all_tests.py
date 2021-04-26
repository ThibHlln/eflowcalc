import unittest
import doctest

import eflowcalc


if __name__ == '__main__':
    test_loader = unittest.TestLoader()
    test_suite = unittest.TestSuite()

    all_tests = test_loader.discover('.', pattern='test_*.py')
    test_suite.addTests(all_tests)

    test_suite.addTests(doctest.DocTestSuite(eflowcalc.eflowcalc))

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(test_suite)
