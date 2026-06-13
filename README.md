# Constellation Time

### Using Predictable Orbits to Synchronize Clocks in Orbital Data Centers

**Ming Leong Tsui** · Boston University · 2026

Low-Earth-orbit constellations are turning from communication relays into computers, and a constellation that *computes* inherits the oldest problem in distributed systems: agreement on time. This repository holds the paper and the single script that reproduces every number in it.

## TL;DR

GPS manufactures time on the ground and broadcasts it one way to passive receivers. An orbital data center inverts every assumption behind that design — commodity oscillators instead of atomic standards, a ground station that sees each satellite roughly 0.5% of the time, time consumed *inside* the constellation by neighboring satellites, and a relativistic rate offset that changes sign. The paper:

- makes the GPS timekeeping assumptions explicit and shows each one fails in orbit;
- derives that **ordering correctness between satellites needs agreement only at the inter-satellite light time** — about 70 µs at the closest approach of a dispersed shell, and a third of a microsecond in close formation — a bar GNSS clears by orders of magnitude;
- shows that for an operator unwilling to depend on GNSS, **two-way time transfer over inter-satellite links** works in the commodity regime: the dominant error (satellite motion during the exchange, hundreds of nanoseconds) collapses to **sub-picosecond** once corrected with predicted orbits — three orders of magnitude below the nanosecond hardware-timestamping floor;
- proposes **constellation time**, a layered architecture that serves time as an interval with an explicit uncertainty bound (in the spirit of Spanner's TrueTime), models the orbit's relativistic and thermal clock dynamics, and maintains the bound through scheduled exchanges computed ahead from geometry.

## Paper

[`main.pdf`](main.pdf) — full text. LaTeX source in [`main.tex`](main.tex).

## Reproducing the numbers

Every computed value in the paper comes from `derivations.py` — orbital geometry, relativistic rates, and error budgets, all from physical constants and the stated assumptions.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python derivations.py
```

It prints a sectioned report:

1. Proper-time rates, circular orbit vs. geoid clock
2. Differential rates inside / between shells
3. Walker shell geometry (72 planes × 22 sats, 550 km, 53°)
   - 3b. Closest approach between any same-shell pair (the Lamport light-time floor)
   - 3c. Like-for-like visibility and rotation-aided pass bound
4. Two-way time-transfer error budget over an ISL
5. Holdover: time to exceed an uncertainty target

and regenerates the figures in `figs/` (`fig_rates.pdf`, `fig_isl.pdf`, `fig_holdover.pdf`).

## Layout

```
derivations.py     reproduces all numbers and figures
main.pdf           the paper
main.tex           LaTeX source
figs/              generated figures + the architecture diagram (arch.tex)
requirements.txt   numpy, scipy, matplotlib
```

## Citation

```bibtex
@misc{tsui2026constellationtime,
  author       = {Tsui, Ming Leong},
  title        = {Using Predictable Orbits to Synchronize Clocks in Orbital Data Centers},
  year         = {2026},
  howpublished = {\url{https://github.com/lucastsui/constellation-time}}
}
```

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). The paper text and figures (`main.pdf`, `main.tex`, `figs/`) are © 2026 Ming Leong Tsui.
