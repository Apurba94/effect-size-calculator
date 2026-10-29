# effect-size-calculator

**Report effect sizes, not just p-values: Cohen's d, Hedges' g, Glass's delta, odds ratios,
risk ratios, risk differences and NNT, each with a confidence interval.**

Part of the *Research Automation Toolkit*: 100 open tools for daily research work.

## Why

APA, CONSORT and most journals require effect sizes with confidence intervals; meta-analysts need
them computed consistently from what papers report (means and SDs, t statistics, 2x2 tables).
Online calculators are not reproducible and rarely state their formulas. This tool does, in code.

## Features

- From means/SDs/ns: Cohen's d (pooled SD), Hedges' g (small-sample correction J), Glass's delta, r
- From a t statistic (independent or paired); d <-> r conversions
- From 2x2 tables: odds ratio (Woolf log CI), risk ratio (log CI), risk difference (Wald CI), NNT
- Haldane-Anscombe 0.5 correction for zero cells; any confidence level; JSON output
- Cohen's benchmarks as a rough label (interpret in context)

## Usage

```bash
pip install .
effect-size-calculator means 24.1 5.2 30 20.3 4.8 32
effect-size-calculator t 2.45 18 20
effect-size-calculator table 20 80 10 90 --level 0.99
```

## Testing

```bash
python -m unittest discover -s tests -t .
```

## License

MIT (c) JANIN A APURBA
