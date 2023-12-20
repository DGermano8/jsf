import unittest
import random
import jsf


class TestCython(unittest.TestCase):

    def setUp(self):
        self.fib_seq = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
        self.fib_calc = [jsf.fib(i) for i in range(10)]

    def test_output(self):
        for i in range(10):
            self.assertEqual(self.fib_seq[i], self.fib_calc[i])
