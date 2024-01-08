import unittest
import random
import math
import jsf

random.seed(1234)


class TestPlaceholderFunction(unittest.TestCase):
    def test_output(self):
        self.assertEqual(1, 1)


class AnEdgeCase(unittest.TestCase):
    def setUp(self):

        self.x0 = [5, 12.1]
        self.threshold = 12
        self.t_max = 5.0
        self.dt = 0.1
        self.rates = lambda x, t: [1.0 if min(x) >= 1.0 else 0.0]
        self.nu_reactants = [[1, 1]]
        self.nu_products = [[0, 0]]

        self.stoich = {
            "nu": [
                [a - b for a, b in zip(r1, r2)]
                for r1, r2 in zip(self.nu_products, self.nu_reactants)
            ],
            "DoDisc": [1, 1],
            "nuReactant": self.nu_reactants,
            "nuProduct": self.nu_products,
        }
        self.opts = {
            "EnforceDo": [0, 0],
            "dt": self.dt,
            "SwitchingThreshold": [self.threshold, self.threshold],
        }

    def test_values(self):
        try:
            self.sim = jsf.jsf(self.x0, self.rates, self.stoich, self.t_max, config=self.opts, method="operator-splitting")
        except ZeroDivisionError:
            self.assertTrue(False)


class TestBirthDeathExample(unittest.TestCase):
    def setUp(self):
        self.num_reps = 500

        self.x0 = 5
        self.x0s = [[self.x0]] * self.num_reps
        self.birth_rate = 1.0
        self.death_rate = 0.5
        self.t_max = 6.0
        self.dt = 0.01
        self.mean_field_soln = lambda t: self.x0 * math.exp(
            (self.birth_rate - self.death_rate) * t
        )
        self.threshold = 10

        def rates(x, t):
            return [self.birth_rate * x[0], self.death_rate * x[0]]

        nu_reactants = [[1], [1]]
        nu_products = [[2], [0]]

        # In the nu-matrix each row is a reaction and each column
        # describes the number items of that species used in the
        # reaction.
        stoich = {
            "nu": [
                [a - b for a, b in zip(r1, r2)]
                for r1, r2 in zip(nu_products, nu_reactants)
            ],
            "DoDisc": [1],
            "nuReactant": nu_reactants,
            "nuProduct": nu_products,
        }

        my_opts = {
            "EnforceDo": [0],
            "dt": self.dt,
            "SwitchingThreshold": [self.threshold],
        }

        self.sims = [
            jsf.jsf(x0, rates, stoich, self.t_max, config=my_opts, method="operator-splitting")
            for x0 in self.x0s
        ]
        self.x_timeseriess = [sim[0][0] for sim in self.sims]

        self.sims_exact = [
            jsf.jsf(x0, rates, stoich, self.t_max, config=my_opts, method="exact")
            for x0 in self.x0s
        ]
        self.x_timeseriess_exact = [sim[0][0] for sim in self.sims_exact]

    def test_output_shape(self):
        self.assertTrue(len(self.sims) == self.num_reps)
        self.assertTrue(len(self.sims_exact) == self.num_reps)

    def test_op_split_mean_field_soln(self):
        # import pdb; pdb.set_trace()
        self.assertTrue(all([foo[0] == self.x0 for foo in self.x_timeseriess]))
        final_x = [xs[-1] for xs in self.x_timeseriess]
        mean_final_x = sum(final_x) / len(final_x)
        # remember to use the unbiased estimator of the variance
        std_final_x = math.sqrt(
            sum([(x - mean_final_x) ** 2 for x in final_x]) / (len(final_x) - 1)
        )
        std_err_final_x = std_final_x / math.sqrt(self.num_reps)
        thry_mean_final_x = self.mean_field_soln(self.t_max)
        self.assertTrue(abs(mean_final_x - thry_mean_final_x) < 2 * std_err_final_x)

    def test_exact_mean_field_soln(self):
        self.assertTrue(all([foo[0] == self.x0 for foo in self.x_timeseriess_exact]))
        final_x = [xs[-1] for xs in self.x_timeseriess_exact]
        mean_final_x = sum(final_x) / len(final_x)
        # remember to use the unbiased estimator of the variance
        std_final_x = math.sqrt(
            sum([(x - mean_final_x) ** 2 for x in final_x]) / (len(final_x) - 1)
        )
        std_err_final_x = std_final_x / math.sqrt(self.num_reps)
        thry_mean_final_x = self.mean_field_soln(self.t_max)
        self.assertTrue(abs(mean_final_x - thry_mean_final_x) < 2 * std_err_final_x)

class TestSISExample(unittest.TestCase):

    def setUp(self):
        # Because this is a stochastic process we need to set the
        # random seed so that we get reproducible results. It is not
        # sufficient to set the seed at the start of the test module
        # because we want to be able to run the tests in any order.
        random.seed(1234)
        x0 = [1000 - 3, 3]

        def rates(x, time):
            s = x[0]
            i = x[1]
            return [2e-3 * s * i, 1.0 * i]

        # In the nu-matrix each row is a reaction and each column
        # describes the number items of that species used in the
        # reaction.
        nu_reactants = [[1, 1], [0, 1]]
        nu_products = [[0, 2], [1, 0]]
        stoich = {
            "nu": [
                [a - b for a, b in zip(r1, r2)]
                for r1, r2 in zip(nu_products, nu_reactants)
            ],
            "DoDisc": [1, 1],
            "nuReactant": nu_reactants,
            "nuProduct": nu_products,
        }
        self.threshold = 40
        my_opts = {
            "EnforceDo": [0, 0],
            "dt": 0.1,
            "SwitchingThreshold": [self.threshold, self.threshold],
        }
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


if __name__ == "__main__":
    unittest.main()
