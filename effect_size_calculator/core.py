"""Standardized mean differences and 2x2 table effect sizes, with large-sample CIs."""
import math
import statistics

Z = statistics.NormalDist().inv_cdf


def _z(level):
    return Z(1 - (1 - level) / 2)


def _summ(x):
    return statistics.fmean(x), statistics.stdev(x), len(x)


def cohens_d(m1, s1, n1, m2, s2, n2, level=0.95):
    """Cohen's d with pooled SD; CI from the large-sample SE (Hedges & Olkin)."""
    sp = math.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    d = (m1 - m2) / sp
    se = math.sqrt((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2)))
    z = _z(level)
    return {"d": d, "se": se, "ci": (d - z * se, d + z * se), "pooled_sd": sp}


def hedges_g(m1, s1, n1, m2, s2, n2, level=0.95):
    """Hedges' g: d multiplied by the small-sample correction J = 1 - 3 / (4 df - 1)."""
    r = cohens_d(m1, s1, n1, m2, s2, n2, level)
    j = 1 - 3 / (4 * (n1 + n2 - 2) - 1)
    g, se = r["d"] * j, r["se"] * j
    z = _z(level)
    return {"g": g, "se": se, "ci": (g - z * se, g + z * se), "J": j}


def glass_delta(m1, m2, s_control):
    """Glass's delta: difference scaled by the control group's SD."""
    return (m1 - m2) / s_control


def from_raw(x, y, kind="hedges", level=0.95):
    a, b = _summ(x), _summ(y)
    return (hedges_g if kind == "hedges" else cohens_d)(*a, *b, level=level)


def d_from_t(t, n1, n2=None):
    """d from an independent (n1, n2) or paired/one-sample (n1 only) t statistic."""
    return t * math.sqrt(1 / n1 + 1 / n2) if n2 else t / math.sqrt(n1)


def d_to_r(d, n1=None, n2=None):
    a = 4 if not (n1 and n2) else (n1 + n2) ** 2 / (n1 * n2)
    return d / math.sqrt(d * d + a)


def r_to_d(r):
    return 2 * r / math.sqrt(1 - r * r)


def interpret_d(d):
    """Cohen's (1988) conventional benchmarks; context matters more than labels."""
    a = abs(d)
    return "negligible" if a < 0.2 else "small" if a < 0.5 else "medium" if a < 0.8 else "large"


def _table(a, b, c, d, correction):
    if correction and 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5    # Haldane-Anscombe correction
    return a, b, c, d


def odds_ratio(a, b, c, d, level=0.95, correction=True):
    """2x2 table: rows exposed/unexposed, columns event/no event:  [[a, b], [c, d]]. Woolf log CI."""
    a, b, c, d = _table(a, b, c, d, correction)
    orr = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    z = _z(level)
    return {"or": orr, "log_se": se, "ci": (math.exp(math.log(orr) - z * se), math.exp(math.log(orr) + z * se))}


def risk_ratio(a, b, c, d, level=0.95, correction=True):
    a, b, c, d = _table(a, b, c, d, correction)
    r1, r2 = a / (a + b), c / (c + d)
    rr = r1 / r2
    se = math.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
    z = _z(level)
    return {"rr": rr, "log_se": se, "ci": (math.exp(math.log(rr) - z * se), math.exp(math.log(rr) + z * se)),
            "risk_exposed": r1, "risk_unexposed": r2}


def risk_difference(a, b, c, d, level=0.95):
    r1, r2 = a / (a + b), c / (c + d)
    se = math.sqrt(r1 * (1 - r1) / (a + b) + r2 * (1 - r2) / (c + d))
    z = _z(level)
    rd = r1 - r2
    return {"rd": rd, "se": se, "ci": (rd - z * se, rd + z * se)}


def number_needed_to_treat(a, b, c, d):
    rd = risk_difference(a, b, c, d)["rd"]
    return math.inf if rd == 0 else 1 / abs(rd)
