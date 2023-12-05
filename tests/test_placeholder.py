import unittest
import random
import jsf

random.seed(1234)

class TestPlaceholderFunction(unittest.TestCase):
    def test_output(self):
        self.assertEqual(1, 1)


class TestSISExample(unittest.TestCase):
    # Set up the simulation

    def setUp(self):
        x0 = [1000-3, 3]
        def rates(x, time):
            s = x[0]
            i = x[1]
            return [2e-3*s*i, 1.0*i]
        nu_reactants = [[1, 1],
                        [0, 1]]
        nu_products = [[0, 2],
                       [1, 0]]
        stoich = {'nu': [[a - b for a, b in zip(r1, r2)] for r1, r2 in zip(nu_products, nu_reactants)],
                  'DoDisc': [1, 1],
                  'nuReactant': nu_reactants,
                  'nuProduct': nu_products}
        self.threshold = 40
        my_opts = {'EnforceDo': [0, 0],
                   'dt': 0.1,
                   'SwitchingThreshold': [self.threshold,
                                          self.threshold]}
        self.sim = jsf.JumpSwitchFlowSimulator(x0, rates, stoich, 10.0, my_opts)

        self.susceptible_timeseries = self.sim[0][0]
        self.infected_timeseries = self.sim[0][1]

    def test_infected_timeseries(self):

        # extract the elements from the infected timeseries that are
        # less than or equal to the threshold and check that they all
        # take values that are within 1e-6 of an integer. This checks
        # that the process if

        small_i_values = [i for i in self.infected_timeseries if i <= self.threshold]
        for i in small_i_values:
            self.assertTrue(abs(i - round(i)) < 1e-6)

        # extract all the values from the infected timeseries that are
        # strictly greater than the threshold plus one and check that
        # they are not within 1e-6 of an integer. This checks that the
        # process is evolving continuously.

        large_i_values = [i for i in self.infected_timeseries if i > self.threshold + 1]
        for i in large_i_values:
            self.assertTrue(abs(i - round(i)) > 1e-6)



if __name__ == '__main__':
    unittest.main()
