"""Command-line interface for effect-size-calculator."""
import argparse
import json

from .core import (cohens_d, d_from_t, d_to_r, glass_delta, hedges_g, interpret_d, number_needed_to_treat, odds_ratio,
                   risk_difference, risk_ratio)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="effect-size-calculator", description="Compute effect sizes with CIs.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("means", help="from group means, SDs and sizes")
    for a in ("m1", "s1", "n1", "m2", "s2", "n2"):
        m.add_argument(a, type=float)
    m.add_argument("--level", type=float, default=0.95)
    t = sub.add_parser("t", help="from a t statistic")
    t.add_argument("t", type=float)
    t.add_argument("n1", type=int)
    t.add_argument("n2", type=int, nargs="?")
    tb = sub.add_parser("table", help="from a 2x2 table a b c d (exposed: a events, b none; unexposed: c, d)")
    for a in "abcd":
        tb.add_argument(a, type=float)
    tb.add_argument("--level", type=float, default=0.95)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    if args.cmd == "means":
        n1, n2 = int(args.n1), int(args.n2)
        d = cohens_d(args.m1, args.s1, n1, args.m2, args.s2, n2, args.level)
        g = hedges_g(args.m1, args.s1, n1, args.m2, args.s2, n2, args.level)
        res = {"cohens_d": d, "hedges_g": g, "glass_delta": glass_delta(args.m1, args.m2, args.s2),
               "r": d_to_r(d["d"], n1, n2), "interpretation": interpret_d(g["g"])}
        text = ("Cohen's d = %.3f [%.3f, %.3f]\nHedges' g = %.3f [%.3f, %.3f] (%s)\nGlass's delta = %.3f\nr = %.3f" %
                (d["d"], *d["ci"], g["g"], *g["ci"], res["interpretation"], res["glass_delta"], res["r"]))
    elif args.cmd == "t":
        d = d_from_t(args.t, args.n1, args.n2)
        res = {"d": d, "interpretation": interpret_d(d)}
        text = "d = %.3f (%s)" % (d, res["interpretation"])
    else:
        o = odds_ratio(args.a, args.b, args.c, args.d, args.level)
        r = risk_ratio(args.a, args.b, args.c, args.d, args.level)
        rd = risk_difference(args.a, args.b, args.c, args.d, args.level)
        res = {"odds_ratio": o, "risk_ratio": r, "risk_difference": rd,
               "nnt": number_needed_to_treat(args.a, args.b, args.c, args.d)}
        text = ("OR = %.3f [%.3f, %.3f]\nRR = %.3f [%.3f, %.3f]\nRD = %.4f [%.4f, %.4f]\nNNT = %.1f" %
                (o["or"], *o["ci"], r["rr"], *r["ci"], rd["rd"], *rd["ci"], res["nnt"]))
    print(json.dumps(res, indent=2, default=str) if args.json else text)
    return 0
