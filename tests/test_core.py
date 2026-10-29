import io
import unittest
from contextlib import redirect_stdout

from effect_size_calculator.cli import main
from effect_size_calculator.core import (cohens_d, d_from_t, d_to_r, from_raw, hedges_g, interpret_d,
                                         number_needed_to_treat, odds_ratio, r_to_d, risk_difference, risk_ratio)


class MeanTests(unittest.TestCase):
    def test_cohens_d(self):
        r = cohens_d(10, 2, 20, 8, 2, 20)
        self.assertAlmostEqual(r["d"], 1.0)
        self.assertAlmostEqual(r["se"], (40 / 400 + 1 / 80) ** 0.5)
        self.assertAlmostEqual(r["ci"][0], 1 - 1.959964 * r["se"], places=5)

    def test_hedges_correction(self):
        g = hedges_g(10, 2, 20, 8, 2, 20)
        self.assertAlmostEqual(g["J"], 1 - 3 / 151)
        self.assertLess(g["g"], 1.0)

    def test_from_raw_and_t(self):
        r = from_raw([1, 2, 3, 4, 5], [3, 4, 5, 6, 7], kind="cohen")
        self.assertAlmostEqual(r["d"], -2 / 2.5 ** 0.5)
        self.assertAlmostEqual(d_from_t(2.0, 25), 0.4)
        self.assertAlmostEqual(d_from_t(3.0, 10, 10), 3 * (0.2 ** 0.5))

    def test_conversions_round_trip(self):
        self.assertAlmostEqual(r_to_d(d_to_r(0.8)), 0.8)
        self.assertEqual([interpret_d(x) for x in (0.1, -0.3, 0.6, 1.2)], ["negligible", "small", "medium", "large"])


class TableTests(unittest.TestCase):
    def test_odds_ratio_woolf(self):
        o = odds_ratio(20, 80, 10, 90)
        self.assertAlmostEqual(o["or"], 2.25)
        self.assertAlmostEqual(o["ci"][0], 0.9943, places=3)
        self.assertAlmostEqual(o["ci"][1], 5.092, places=2)

    def test_risk_ratio_difference_nnt(self):
        r = risk_ratio(20, 80, 10, 90)
        self.assertAlmostEqual(r["rr"], 2.0)
        self.assertAlmostEqual(risk_difference(20, 80, 10, 90)["rd"], 0.1)
        self.assertAlmostEqual(number_needed_to_treat(20, 80, 10, 90), 10.0)

    def test_zero_cell_correction(self):
        self.assertAlmostEqual(odds_ratio(0, 10, 5, 5)["or"], (0.5 * 5.5) / (10.5 * 5.5))

    def test_cli(self):
        with redirect_stdout(io.StringIO()) as buf:
            main(["table", "20", "80", "10", "90"])
        self.assertIn("OR = 2.250", buf.getvalue())
        with redirect_stdout(io.StringIO()) as buf:
            main(["means", "10", "2", "20", "8", "2", "20"])
        self.assertIn("Cohen's d = 1.000", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
