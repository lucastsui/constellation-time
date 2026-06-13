"""
Derivations for "Predictable Orbits as a Synchronization Substrate for
Orbital Data Centers". Every number quoted in the paper is computed here
from physical constants and stated assumptions, so the error budget is
reproducible. Run with the playground venv python.
"""

import numpy as np

# ----------------------------------------------------------------------
# Constants (IERS / WGS-84)
# ----------------------------------------------------------------------
GM = 3.986004418e14        # m^3/s^2, Earth gravitational parameter
C = 299792458.0            # m/s
C2 = C * C
R_E = 6.371e6              # m, mean Earth radius
R_EQ = 6.378137e6          # m, equatorial radius
W0 = 62636856.0            # m^2/s^2, geoid potential (gravity + rotation)
J2 = 1.0826e-3
OMEGA_E = 7.2921150e-5     # rad/s, Earth rotation
SEC_PER_DAY = 86400.0

out = []
def rec(label, value, unit=""):
    out.append((label, value, unit))
    print(f"{label:62s} {value:>14} {unit}")

print("=" * 95)
print("1. PROPER-TIME RATES, CIRCULAR ORBIT vs GEOID CLOCK")
print("=" * 95)
# For a circular orbit of radius r: d(tau_sat)/dt = 1 - (GM/r + v^2/2)/c^2
#                                                 = 1 - (3GM / 2r)/c^2
# For a geoid clock:                d(tau_g)/dt   = 1 - W0/c^2
# Fractional rate of satellite clock vs geoid clock:
#     y(r) = (W0 - 3GM/(2r)) / c^2
def y_rate(r):
    return (W0 - 1.5 * GM / r) / C2

def us_per_day(y):
    return y * SEC_PER_DAY * 1e6

r_gps = 26561.75e3
r_leo = R_E + 550e3
r_iss = R_E + 420e3

rec("GPS (r=26561.75 km) rate vs geoid", f"{us_per_day(y_rate(r_gps)):+.2f}", "us/day")
rec("LEO 550 km (r=6921.0 km) rate vs geoid", f"{us_per_day(y_rate(r_leo)):+.2f}", "us/day")
rec("ISS 420 km cross-check (literature ~ -25)", f"{us_per_day(y_rate(r_iss)):+.2f}", "us/day")
r_cross = 1.5 * GM / W0
rec("Null-shift radius", f"{r_cross/1e3:.0f}", "km")
rec("Null-shift altitude", f"{(r_cross - R_E)/1e3:.0f}", "km")

v_leo = np.sqrt(GM / r_leo)
T_orb = 2 * np.pi * np.sqrt(r_leo**3 / GM)
rec("Orbital speed at 550 km", f"{v_leo:.1f}", "m/s")
rec("Orbital period at 550 km", f"{T_orb:.0f} s = {T_orb/60:.1f}", "min")

print()
print("=" * 95)
print("2. DIFFERENTIAL RATES INSIDE / BETWEEN SHELLS")
print("=" * 95)
# dy/dr = (3GM / 2r^2) / c^2  : rate difference per metre of altitude
dydr = 1.5 * GM / (r_leo**2 * C2)
rec("Rate sensitivity at 550 km", f"{dydr:.3e}", "frac per m")
rec("  = per 1 km altitude separation", f"{dydr*1e3*SEC_PER_DAY*1e9:.1f}", "ns/day")
rec("  = shells 30 km apart (e.g. 540 vs 570)", f"{dydr*30e3*SEC_PER_DAY*1e9:.0f}", "ns/day")
# Same-shell circular orbits: identical r and v magnitude -> identical rate.
rec("Same-shell circular orbits, rate difference", "0 (first order)", "")

# Eccentricity periodic term: dt = (2/c^2) sqrt(GM a) e sin E
ecc_coef = 2.0 * np.sqrt(GM * r_leo) / C2
rec("Eccentricity periodic amplitude at a=6921 km", f"{ecc_coef*1e9:.2f}", "ns per unit e")
rec("  at e = 1e-3", f"{ecc_coef*1e-3*1e9:.2f}", "ns")
rec("  at e = 1e-4 (typical Starlink)", f"{ecc_coef*1e-4*1e9:.3f}", "ns")
# GPS cross-check vs Ashby (~46 ns at e=0.02 half-amplitude... he quotes ~46 ns for e=0.02)
ecc_gps = 2.0 * np.sqrt(GM * r_gps) / C2
rec("GPS cross-check, e=0.02 amplitude (Ashby ~46 ns)", f"{ecc_gps*0.02*1e9:.1f}", "ns")

# J2 periodic effect on the clock at LEO: fractional amplitude ~ GM J2 R^2 / (r^3 c^2) at 2u
yJ2 = GM * J2 * R_EQ**2 / (r_leo**3 * C2)
tJ2 = yJ2 / (2 * (2*np.pi/T_orb))   # integrate sinusoid at twice orbital frequency
rec("J2 periodic fractional amplitude at 550 km", f"{yJ2:.2e}", "")
rec("J2 periodic time amplitude", f"{tJ2*1e9:.2f}", "ns")

print()
print("=" * 95)
print("3. WALKER SHELL GEOMETRY (72 planes x 22 sats, 550 km, 53 deg)")
print("=" * 95)
INC = np.radians(53.0)
P_PLANES, S_SATS = 72, 22
def sat_pos_vel(plane, slot, t, F=0):
    """ECI position/velocity, circular orbit. Walker delta phasing F."""
    raan = 2*np.pi * plane / P_PLANES
    u0 = 2*np.pi * (slot / S_SATS + F * plane / (P_PLANES*S_SATS))
    n = 2*np.pi / T_orb
    u = u0 + n*t
    # orbital plane -> ECI
    cu, su = np.cos(u), np.sin(u)
    co, so = np.cos(raan), np.sin(raan)
    ci, si = np.cos(INC), np.sin(INC)
    r_vec = r_leo * np.array([cu*co - su*ci*so, cu*so + su*ci*co, su*si])
    v_vec = v_leo * np.array([-su*co - cu*ci*so, -su*so + cu*ci*co, cu*si])
    return r_vec, v_vec

t = np.linspace(0, T_orb, 4000)
# intra-plane neighbour: same plane, adjacent slot
rng_intra = []
for tt in t[::40]:
    r1, _ = sat_pos_vel(0, 0, tt)
    r2, _ = sat_pos_vel(0, 1, tt)
    rng_intra.append(np.linalg.norm(r2 - r1))
rng_intra = np.array(rng_intra)
rec("Intra-plane neighbour range (constant)", f"{rng_intra.mean()/1e3:.0f}", "km")
rec("Intra-plane range variation over orbit", f"{(rng_intra.max()-rng_intra.min()):.2e}", "m")

def crossplane_stats(F):
    rng, rr = [], []
    for tt in t:
        r1, v1 = sat_pos_vel(0, 0, tt, F)
        r2, v2 = sat_pos_vel(1, 0, tt, F)
        d = r2 - r1
        rho = np.linalg.norm(d)
        rng.append(rho)
        rr.append(np.dot(d, v2 - v1) / rho)
    return np.array(rng), np.array(rr)

for F in (0, 11, 39):
    rng_x, rr_x = crossplane_stats(F)
    rec(f"Cross-plane neighbour, F={F:2d}  range min-max",
        f"{rng_x.min()/1e3:.0f} - {rng_x.max()/1e3:.0f}", "km")
    rec(f"Cross-plane neighbour, F={F:2d}  |range rate| max",
        f"{np.abs(rr_x).max():.0f}", "m/s")

# crossing link (ascending vs descending, worst case closing speed)
v_cross_max = 2 * v_leo * np.sin(INC)
rec("Crossing-link closing speed bound 2 v sin(i)", f"{v_cross_max:.0f}", "m/s")

print()
print("=" * 95)
print("3b. CLOSEST APPROACH BETWEEN ANY SAME-SHELL PAIR (Lamport floor)")
print("=" * 95)
# The anomalous-behavior floor is the minimum light time over ALL causally
# connectable pairs, not just grid neighbours. Crossing-track satellites pass
# close together once per orbit; sample all pairs at fine time steps.
def min_approach(F, nt=2000):
    ts = np.linspace(0, T_orb, nt)
    best = np.inf
    # all satellites, positions per time step (vectorised)
    planes = np.arange(P_PLANES)
    slots = np.arange(S_SATS)
    raan = 2*np.pi*planes/P_PLANES
    n = 2*np.pi/T_orb
    ci, si = np.cos(INC), np.sin(INC)
    for tt in ts:
        u = 2*np.pi*(slots[None, :]/S_SATS + F*planes[:, None]/(P_PLANES*S_SATS)) + n*tt
        cu, su = np.cos(u), np.sin(u)
        co, so = np.cos(raan)[:, None], np.sin(raan)[:, None]
        x = r_leo*(cu*co - su*ci*so)
        y_ = r_leo*(cu*so + su*ci*co)
        z = r_leo*su*si
        pts = np.stack([x.ravel(), y_.ravel(), z.ravel()], axis=1)
        # nearest non-identical pair via sorted sweep is overkill; use cdist on
        # a random subsample? exact: chunked pairwise min
        from scipy.spatial import cKDTree
        tree = cKDTree(pts)
        d, _ = tree.query(pts, k=2)
        m = d[:, 1].min()
        if m < best:
            best = m
    return best

try:
    from scipy.spatial import cKDTree  # noqa: F401
    for F in (11, 17, 39):
        m = min_approach(F)
        rec(f"Min same-shell pair distance, F={F:2d}", f"{m/1e3:.1f} km -> light time {m/C*1e6:.0f}", "us")
except ImportError:
    # fallback without scipy: brute force on coarser sampling
    def min_approach_brute(F, nt=600):
        ts = np.linspace(0, T_orb, nt)
        best = np.inf
        planes = np.arange(P_PLANES)
        slots = np.arange(S_SATS)
        raan = 2*np.pi*planes/P_PLANES
        n = 2*np.pi/T_orb
        ci, si = np.cos(INC), np.sin(INC)
        for tt in ts:
            u = 2*np.pi*(slots[None, :]/S_SATS + F*planes[:, None]/(P_PLANES*S_SATS)) + n*tt
            cu, su = np.cos(u), np.sin(u)
            co, so = np.cos(raan)[:, None], np.sin(raan)[:, None]
            pts = np.stack([(r_leo*(cu*co - su*ci*so)).ravel(),
                            (r_leo*(cu*so + su*ci*co)).ravel(),
                            (r_leo*su*si).ravel()], axis=1)
            dif = pts[:, None, :2]  # placeholder, replaced below
            # chunked exact pairwise minimum
            mbest = np.inf
            for i in range(0, len(pts), 264):
                d2 = ((pts[i:i+264, None, :] - pts[None, :, :])**2).sum(-1)
                np.fill_diagonal(d2[:, i:i+264], np.inf)
                mbest = min(mbest, np.sqrt(d2.min()))
            best = min(best, mbest)
        return best
    for F in (11, 17, 39):
        m = min_approach_brute(F)
        rec(f"Min same-shell pair distance, F={F:2d}", f"{m/1e3:.1f} km -> light time {m/C*1e6:.0f}", "us")

print()
print("=" * 95)
print("3c. LIKE-FOR-LIKE VISIBILITY AND ROTATION-AIDED PASS BOUND")
print("=" * 95)
def vis_fraction(r, el_deg):
    el = np.radians(el_deg)
    lam = np.arccos((R_E / r) * np.cos(el)) - el
    return (1 - np.cos(lam)) / 2, lam
for el in (0.0, 25.0):
    f_leo, _ = vis_fraction(r_leo, el)
    f_gps, _ = vis_fraction(r_gps, el)
    rec(f"Visibility fraction at {el:.0f} deg mask, LEO vs GPS",
        f"{f_leo*100:.2f}% vs {f_gps*100:.1f}%", "")
# co-rotating ground track: relative angular rate reduced by Earth rotation
_, lam25 = vis_fraction(r_leo, 25.0)
w_orb = 2*np.pi/T_orb
pass_corot = 2*lam25 / (w_orb - OMEGA_E)
rec("Pass bound with co-rotating Earth (25 deg mask)", f"{pass_corot:.0f} s = {pass_corot/60:.1f}", "min")

print()
print("=" * 95)
print("4. TWO-WAY TIME TRANSFER ERROR BUDGET OVER AN ISL")
print("=" * 95)
# Estimator bias when geometry changes during the exchange:
#   bias = -rho_dot (d + tau_turn) / (2c),  d = rho/c one-way delay
def twtt_bias(rho, rho_dot, tau_turn):
    d = rho / C
    return rho_dot * (d + tau_turn) / (2 * C)

cases = [
    ("Intra-plane ISL  rho=1972 km, rdot~1 m/s,    tau=5 ms",  1972e3, 1.0,   5e-3),
    ("Cross-plane ISL  rho=1300 km, rdot=350 m/s,  tau=5 ms",  1300e3, 350.0, 5e-3),
    ("Cross-plane ISL  rho=1300 km, rdot=350 m/s,  tau=1 ms",  1300e3, 350.0, 1e-3),
    ("Crossing link    rho=2000 km, rdot=12 km/s,  tau=5 ms",  2000e3, 12e3,  5e-3),
]
for label, rho, rd, tau in cases:
    rec(f"Uncorrected bias  {label}", f"{abs(twtt_bias(rho, rd, tau))*1e9:.1f}", "ns")

# After ephemeris correction the residual is the same expression with
# rho_dot replaced by the relative-velocity knowledge error.
for dv in (1e-3, 1e-2, 1.0):
    rec(f"Residual bias, velocity known to {dv*1e3:6.1f} mm/s (rho=3000 km, tau=5 ms)",
        f"{abs(twtt_bias(3000e3, dv, 5e-3))*1e12:.3f}", "ps")

# Shapiro delay on a 3000 km ISL between r1 = r2 = r_leo
rho_s = 3000e3
shapiro = (2*GM/C**3) * np.log((2*r_leo + rho_s) / (2*r_leo - rho_s))
rec("Shapiro delay, 3000 km ISL (cancels in two-way)", f"{shapiro*1e12:.1f}", "ps")

print()
print("=" * 95)
print("5. HOLDOVER: TIME TO EXCEED AN UNCERTAINTY TARGET, eps = y * t")
print("=" * 95)
classes = [
    ("TCXO, thermally disturbed      y=1e-7",  1e-7),
    ("OCXO, calibrated               y=1e-9",  1e-9),
    ("CSAC / steered OCXO            y=1e-11", 1e-11),
    ("USO, GRACE class, rate-modeled y=1e-13", 1e-13),
]
targets = [(100e-6, "100 us"), (1e-6, "1 us"), (100e-9, "100 ns"), (10e-9, "10 ns")]
hdr = " " * 36 + "".join(f"{lbl:>12}" for _, lbl in targets)
print(hdr)
def fmt_t(s):
    if s < 60: return f"{s:.0f} s"
    if s < 3600: return f"{s/60:.0f} min"
    if s < SEC_PER_DAY: return f"{s/3600:.1f} h"
    return f"{s/SEC_PER_DAY:.0f} d"
for label, yy in classes:
    row = "".join(f"{fmt_t(eps/yy):>12}" for eps, _ in targets)
    print(f"{label:36s}{row}")

print()
print("=" * 95)
print("5b. CSAC AGING AND FLEET STALENESS TERMS")
print("=" * 95)
# CSAC datasheet aging <9e-10/month integrates quadratically.
a_csac = 9e-10 / (30 * SEC_PER_DAY)        # fractional frequency per second
def csac_err(t, y0=1e-11):
    return y0 * t + 0.5 * a_csac * t * t
for label, tt in (("1 day", SEC_PER_DAY), ("1 week", 7*SEC_PER_DAY),
                  ("1 month", 30*SEC_PER_DAY), ("3 months", 90*SEC_PER_DAY)):
    rec(f"CSAC holdover incl. aging at {label}", f"{csac_err(tt)*1e6:.1f}", "us")
t_100us = (-1e-11 + np.sqrt(1e-22 + 2*a_csac*100e-6)) / a_csac
rec("CSAC time to exceed 100 us (aging model)", f"{t_100us/SEC_PER_DAY:.1f}", "days")
t_1ms = (-1e-11 + np.sqrt(1e-22 + 2*a_csac*1e-3)) / a_csac
rec("CSAC time to exceed 1 ms (aging model)", f"{t_1ms/SEC_PER_DAY:.1f}", "days")
# fleet bound staleness: eps_fleet ~ h * (r_link + y * T_slot)
h_depth = (P_PLANES + S_SATS) // 2
for y_osc, T_slot in ((1e-9, 1.0), (1e-9, 0.1), (1e-11, 1.0)):
    rec(f"Fleet bound h=47, r=0.5 ns, y={y_osc:.0e}, T={T_slot:g} s",
        f"{h_depth*(0.5e-9 + y_osc*T_slot)*1e9:.0f}", "ns")

print()
print("=" * 95)
print("6. GROUND CONTACT AND ECLIPSE GEOMETRY AT 550 km")
print("=" * 95)
el = np.radians(25.0)
lam = np.arccos((R_E / r_leo) * np.cos(el)) - el          # Earth-central half angle
pass_max = 2 * lam / (2*np.pi / T_orb)
frac_area = (1 - np.cos(lam)) / 2
rec("Earth-central half-angle of visibility cone (25 deg mask)", f"{np.degrees(lam):.1f}", "deg")
rec("Longest possible single-station pass", f"{pass_max:.0f} s = {pass_max/60:.1f}", "min")
rec("Fraction of sphere where one station is visible", f"{frac_area*100:.2f}", "%")
# GPS contrast: visibility cone at 0 deg mask from MEO
lam_gps = np.arccos(R_E / r_gps)
rec("GPS MEO visibility half-angle (0 deg mask)", f"{np.degrees(lam_gps):.1f}", "deg")
rec("Fraction of sphere where one station sees a GPS sat", f"{(1-np.cos(lam_gps))/2*100:.0f}", "%")
# eclipse, worst case beta = 0
half_shadow = np.arcsin(R_E / r_leo)
ecl_frac = half_shadow / np.pi
rec("Worst-case eclipse fraction per orbit", f"{ecl_frac*100:.1f}", "%")
rec("Worst-case eclipse duration", f"{ecl_frac*T_orb/60:.1f}", "min")

print()
print("=" * 95)
print("7. ORDERING COST: COMMIT-WAIT vs ISL PROPAGATION")
print("=" * 95)
rec("One ISL hop, 1972 km, propagation", f"{1972e3/C*1e3:.2f}", "ms")
rec("10-hop constellation path ~20000 km", f"{20000e3/C*1e3:.0f}", "ms")
rec("TrueTime production average eps (Spanner paper)", "~4", "ms")
rec("Fabric eps target (this proposal)", "0.1-1", "us")
# time for eps to reach the 1-hop propagation delay (6.6 ms) per class
for label, yy in classes:
    rec(f"  eps growth to 6.6 ms 1-hop RTT-scale: {label.split()[0]}", fmt_t(6.6e-3/yy), "")

# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "figure.dpi": 150, "lines.linewidth": 1.3,
})
FIGDIR = "/Users/tsuimingleong/Desktop/leo-sync-paper/figs"

# Figure A: rate offset vs altitude
alts = np.linspace(200e3, 36000e3, 800)
ys = np.array([us_per_day(y_rate(R_E + a)) for a in alts])
fig, ax = plt.subplots(figsize=(3.4, 2.1))
ax.plot(alts/1e3, ys, color="#1a4d8f")
ax.axhline(0, color="0.6", lw=0.7)
for a_km, lbl, dy in ((550, "Starlink shell", -6), (3175, "null shift", 4),
                      (20190, "GPS", 4), (35786, "GEO", -6)):
    yv = us_per_day(y_rate(R_E + a_km*1e3))
    ax.plot([a_km], [yv], "o", ms=3.5, color="#b3322e")
    ax.annotate(lbl, (a_km, yv), textcoords="offset points", xytext=(4, dy), fontsize=7)
ax.set_xscale("log")
ax.set_xlabel("altitude / km")
ax.set_ylabel("clock rate vs geoid\nmicroseconds per day")
fig.tight_layout()
fig.savefig(f"{FIGDIR}/fig_rates.pdf", bbox_inches="tight")

# Figure B: cross-plane ISL range and range rate over one orbit (F=11)
rng_x, rr_x = crossplane_stats(11)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.4, 2.9), sharex=True)
ax1.plot(t/60, rng_x/1e3, color="#1a4d8f")
ax1.set_ylabel("range / km")
ax2.plot(t/60, rr_x, color="#b3322e")
ax2.set_ylabel("range rate / m s$^{-1}$")
ax2.set_xlabel("time / min, one orbital period")
rho_mid = rng_x.mean()   # use this link's own mean range for the conversion
sec = ax2.secondary_yaxis("right",
        functions=(lambda v: np.abs(v)*(rho_mid/C + 5e-3)/(2*C)*1e9*np.sign(v),
                   lambda b: b*2*C/((rho_mid/C + 5e-3))/1e9))
sec.set_ylabel("uncorrected bias / ns")
for ax in (ax1, ax2):
    ax.grid(alpha=0.25, lw=0.4)
fig.tight_layout()
fig.savefig(f"{FIGDIR}/fig_isl.pdf", bbox_inches="tight")

# Figure C: holdover growth. CSAC carries its datasheet aging term so the
# curve bends up beyond a day instead of understating long holdover.
fig, ax = plt.subplots(figsize=(3.4, 2.3))
tt = np.logspace(0, 6, 400)
curves = [
    (r"TCXO  $y=10^{-7}$",            1e-7 * tt,                       "#b3322e"),
    (r"OCXO  $y=10^{-9}$",            1e-9 * tt,                       "#d98032"),
    (r"CSAC  $y=10^{-11}$ + aging",   1e-11*tt + 0.5*a_csac*tt*tt,     "#1a4d8f"),
    (r"USO  $y=10^{-13}$, rate-modeled", 1e-13 * tt,                   "#3d8f5f"),
]
for lbl, err, col in curves:
    ax.loglog(tt, err*1e9, label=lbl, color=col)
ax.axhline(100, color="0.5", lw=0.7, ls="--")
ax.text(1.4, 130, "100 ns target", fontsize=6.5, color="0.35")
ax.axhline(6.6e6, color="0.5", lw=0.7, ls=":")
ax.text(1.4, 9e6, "one ISL hop, 6.6 ms", fontsize=6.5, color="0.35")
ax.axhline(7e4, color="0.5", lw=0.7, ls="-.")
ax.text(1.4, 1.0e5, "shell light-time floor, 70 us", fontsize=6.5, color="0.35")
ax.axvline(T_orb, color="0.7", lw=0.7)
ax.text(T_orb*1.2, 0.03, "one orbit", rotation=90, fontsize=6.5, color="0.35")
ax.axvline(SEC_PER_DAY, color="0.7", lw=0.7)
ax.text(SEC_PER_DAY*1.2, 0.03, "one day", rotation=90, fontsize=6.5, color="0.35")
ax.set_xlabel("holdover duration / s")
ax.set_ylabel("uncertainty growth / ns")
ax.set_ylim(1e-2, 1e10)
ax.legend(loc="upper left", framealpha=0.9)
ax.grid(alpha=0.25, lw=0.4, which="both")
fig.tight_layout()
fig.savefig(f"{FIGDIR}/fig_holdover.pdf", bbox_inches="tight")

print("\nFigures written to", FIGDIR)
