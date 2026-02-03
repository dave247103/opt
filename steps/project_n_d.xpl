<plask loglevel="result">

<defines>
  <define name="lam_min" value="550"/>
  <define name="n_ref" value="4.0"/>
  <define name="n_layer" value="n_ref**0.5"/>
</defines>

<materials>
  <material name="mySi" base="Si">
    <nr>n_ref</nr>
    <absp>0</absp>
  </material>
  <material name="mySiO2" base="SiO2">
    <Nr>n_layer</Nr>
  </material>
  <module name="si_nk"/>
</materials>

<geometry>
  <cartesian2d name="simple" axes="x,y " left="periodic" right="periodic" bottom="Si_Shinke">
    <stack>
      <rectangle name="layer_SiO2" material="SiO2" dx="1" dy="1"/>
      <rectangle name="layer_Si3N4" material="Si3N4" dx="1" dy="1"/>
    </stack>
  </cartesian2d>
</geometry>

<solvers>
  <optical name="OPTICAL" solver="Fourier2D" lib="modal">
    <geometry ref="simple"/>
    <expansion size="0"/>
  </optical>
</solvers>

<script><![CDATA[
import numpy as np
import matplotlib.pyplot as plt
import plask.material as pm

L1 = GEO.layer_SiO2
L2 = GEO.layer_Si3N4

SIDE = "top"
EPS_NM = 5e-4

D_BOUNDS = (10.0, 200.0)     # nm
TOTAL_MAX = 500.0           # nm
N_BOUNDS = (1.10, 3.40)     # (-)

LAMS_OPT   = np.arange(400.0, 700.1, 10.0)
THETAS_OPT = np.arange(0.0, 90.1, 5.0)

LAMS_FINE   = np.arange(400.0, 700.1, 5.0)
THETAS_FINE = np.arange(0.0, 90.1, 2.0)

LAMS_MAP   = np.arange(400.0, 700.1, 5.0)
THETAS_MAP = np.arange(0.0, 90.1, 2.0)

POP, GENS, SEED = 16, 8, 1
F, CR = 0.8, 0.9


def set_layers(d1_nm, d2_nm, n1, n2):
    L1.height = max(float(d1_nm), EPS_NM) * 1e-3
    L2.height = max(float(d2_nm), EPS_NM) * 1e-3
    L1.nr = n1
    L2.nr = n2

def set_bare():
    set_layers(EPS_NM, EPS_NM, 1.0, 1.0)

def ktran_um(lam_nm, theta_deg):
    return float(2e3 * np.pi / lam_nm * np.sin(np.deg2rad(theta_deg)))

def R_unpol(lam_nm, theta_deg):
    OPTICAL.ktran = ktran_um(lam_nm, theta_deg)
    r_te = float(OPTICAL.compute_reflectivity(float(lam_nm), SIDE, "TE"))
    r_tm = float(OPTICAL.compute_reflectivity(float(lam_nm), SIDE, "TM"))
    return 0.5 * (r_te + r_tm)

def mean_R(d1_nm, d2_nm, n1, n2, lams=LAMS_OPT, thetas=THETAS_OPT):
    set_layers(d1_nm, d2_nm, n1, n2)
    acc = 0.0
    for th in thetas:
        for lam in lams:
            acc += R_unpol(lam, th)
    return acc / (len(lams) * len(thetas))

def objective(x):
    d1, d2, n1, n2 = map(float, x)

    if not (D_BOUNDS[0] <= d1 <= D_BOUNDS[1] and D_BOUNDS[0] <= d2 <= D_BOUNDS[1]):
        return np.inf
    if TOTAL_MAX is not None and (d1 + d2) > TOTAL_MAX:
        return np.inf

    if not (N_BOUNDS[0] <= n1 <= N_BOUNDS[1] and N_BOUNDS[0] <= n2 <= N_BOUNDS[1]):
        return np.inf
    if n1 > n2:
        return np.inf

    return mean_R(d1, d2, n1, n2)

def differential_evolution(bounds, pop=POP, gens=GENS, F=F, CR=CR, seed=SEED):
    rng = np.random.default_rng(seed)
    lo = np.array([b[0] for b in bounds], float)
    hi = np.array([b[1] for b in bounds], float)
    dim = len(bounds)

    X = rng.uniform(lo, hi, size=(pop, dim))
    fX = np.array([objective(x) for x in X], float)

    best_i = int(np.argmin(fX))
    best, fbest = X[best_i].copy(), float(fX[best_i])

    for g in range(gens):
        for i in range(pop):
            idx = [j for j in range(pop) if j != i]
            a, b, c = X[rng.choice(idx, 3, replace=False)]
            v = np.clip(a + F * (b - c), lo, hi)

            jrand = rng.integers(dim)
            mask = (rng.random(dim) < CR)
            mask[jrand] = True
            u = np.where(mask, v, X[i])

            fu = float(objective(u))
            if fu < fX[i]:
                X[i], fX[i] = u, fu
                if fu < fbest:
                    best, fbest = u.copy(), fu

        print_log("info",
                  f"gen {g+1}/{gens}: d1={best[0]:.2f} nm, d2={best[1]:.2f} nm, "
                  f"n1={best[2]:.3f}, n2={best[3]:.3f}, meanR={fbest:.6f}%")
    return best, fbest

def spectrum(lams, theta_deg, d1, d2, n1, n2):
    set_layers(d1, d2, n1, n2)
    return np.array([R_unpol(lam, theta_deg) for lam in lams], float)

def map_R(lams, thetas, d1, d2, n1, n2):
    set_layers(d1, d2, n1, n2)
    Z = np.empty((thetas.size, lams.size), float)
    for i, th in enumerate(thetas):
        for j, lam in enumerate(lams):
            Z[i, j] = R_unpol(lam, th)
    return Z

def heatmap(Z, title, vmax):
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(THETAS_MAP, LAMS_MAP, Z.T, shading="auto", vmin=0, vmax=vmax)
    plt.xlabel("angle [deg]")
    plt.ylabel("wavelength [nm]")
    plt.title(title)
    plt.colorbar(label="reflectance [%]")
    plt.tight_layout()

# baseline
J_bare = mean_R(EPS_NM, EPS_NM, 1.0, 1.0)
print_log("result", f"bare meanR (opt grid) = {J_bare:.6f}%")

# optimize (d1, d2, n1, n2)
bounds = [D_BOUNDS, D_BOUNDS, N_BOUNDS, N_BOUNDS]
best, Jopt = differential_evolution(bounds)
d1_star, d2_star, n1_star, n2_star = map(float, best)

J_bare_chk = mean_R(EPS_NM, EPS_NM, 1.0, 1.0, lams=LAMS_FINE[::10], thetas=THETAS_OPT)
J_opt_chk  = mean_R(d1_star, d2_star, n1_star, n2_star, lams=LAMS_FINE[::10], thetas=THETAS_OPT)

print_log("result", f"opt: d1={d1_star:.3f} nm, d2={d2_star:.3f} nm, n1={n1_star:.4f}, n2={n2_star:.4f}")
print_log("result", f"meanR bare={J_bare_chk:.6f}%  opt={J_opt_chk:.6f}%  (check grid)")


# heatmaps
Rb_map = map_R(LAMS_MAP, THETAS_MAP, EPS_NM, EPS_NM, 1.0, 1.0)
Ro_map = map_R(LAMS_MAP, THETAS_MAP, d1_star, d2_star, n1_star, n2_star)
heatmap(Rb_map, "Bare Si: unpolarized R(λ,θ)", vmax=50)
heatmap(Ro_map, f"2-layer AR: d1={d1_star:.1f} nm, d2={d2_star:.1f} nm, n1={n1_star:.3f}, n2={n2_star:.3f}", vmax=15)

plt.show()

]]></script>

</plask>
