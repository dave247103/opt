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
L1.material = pm.SiO2()
L2.material = pm.Si3N4()

SIDE = "top"
EPS_NM = 5e-4

D_BOUNDS = (10.0, 200.0)   # nm
TOTAL_MAX = 300.0          # nm

LAMS_OPT   = np.arange(400.0, 700.1, 15.0)   # coarse (objective)
THETAS_OPT = np.arange(0.0,   60.1, 15.0)

LAMS_FINE   = np.arange(400.0, 700.1, 2.0)   # spectra
THETAS_FINE = np.array([0.0, 30.0, 60.0])

LAMS_MAP   = np.arange(400.0, 700.1, 4.0)    # heatmap
THETAS_MAP = np.arange(0.0,   70.1, 4.0)

POP, GENS, SEED = 10, 10, 1
F, CR = 0.8, 0.9

def set_stack(d1_nm, d2_nm):
    L1.height = max(float(d1_nm), EPS_NM) * 1e-3  # um
    L2.height = max(float(d2_nm), EPS_NM) * 1e-3

def ktran_um(lam_nm, theta_deg, n0=1.0):
    return float(2e3 * np.pi / lam_nm * n0 * np.sin(np.deg2rad(theta_deg)))  # 1/um

def R_unpol(lam_nm, theta_deg):
    OPTICAL.ktran = ktran_um(lam_nm, theta_deg)
    Rte = float(OPTICAL.compute_reflectivity(float(lam_nm), SIDE, "TE"))
    Rtm = float(OPTICAL.compute_reflectivity(float(lam_nm), SIDE, "TM"))
    return 0.5 * (Rte + Rtm)

def mean_R_unpol(d1_nm, d2_nm, lams=LAMS_OPT, thetas=THETAS_OPT):
    set_stack(d1_nm, d2_nm)
    acc = 0.0
    for th in thetas:
        for lam in lams:
            acc += R_unpol(lam, th)
    return acc / (len(lams) * len(thetas))

def objective(x):
    d1, d2 = float(x[0]), float(x[1])
    lo, hi = D_BOUNDS
    if not (lo <= d1 <= hi and lo <= d2 <= hi):
        return np.inf
    if TOTAL_MAX is not None and (d1 + d2) > TOTAL_MAX:
        return np.inf
    return mean_R_unpol(d1, d2)

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

        print_log("info", f"gen {g+1}/{gens}: best d1={best[0]:.2f} nm, d2={best[1]:.2f} nm, meanR={fbest:.6f}%")

    return best, fbest

def spectrum_unpol(lams, theta_deg, d1_nm, d2_nm):
    set_stack(d1_nm, d2_nm)
    return np.array([R_unpol(lam, theta_deg) for lam in lams], float)

def map_unpol(lams, thetas, d1_nm, d2_nm):
    set_stack(d1_nm, d2_nm)
    Z = np.empty((thetas.size, lams.size), float)
    for i, th in enumerate(thetas):
        for j, lam in enumerate(lams):
            Z[i, j] = R_unpol(lam, th)
    return Z

# baseline
J_bare = mean_R_unpol(EPS_NM, EPS_NM)
print_log("result", f"bare meanR (opt grid) = {J_bare:.6f}%")

# optimize
best, Jopt = differential_evolution([D_BOUNDS, D_BOUNDS])
d1_star, d2_star = float(best[0]), float(best[1])

J_bare_fine = mean_R_unpol(EPS_NM, EPS_NM, lams=LAMS_FINE[::5], thetas=THETAS_OPT)  # quick fine-ish
J_star_fine = mean_R_unpol(d1_star, d2_star, lams=LAMS_FINE[::5], thetas=THETAS_OPT)
print_log("result", f"opt: d1={d1_star:.3f} nm, d2={d2_star:.3f} nm")
print_log("result", f"meanR bare={J_bare_fine:.6f}%  opt={J_star_fine:.6f}%  (check grid)")

# spectra
plt.figure()
for th in THETAS_FINE:
    Rb = spectrum_unpol(LAMS_FINE, th, EPS_NM, EPS_NM)
    Ro = spectrum_unpol(LAMS_FINE, th, d1_star, d2_star)
    plt.plot(LAMS_FINE, Rb, label=f"bare θ={th:g}°")
    plt.plot(LAMS_FINE, Ro, label=f"opt  θ={th:g}°  (d1={d1_star:.0f}, d2={d2_star:.0f} nm)")
plt.xlabel("wavelength [nm]")
plt.ylabel("reflectance [%]")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()

# heatmaps
Rb_map = map_unpol(LAMS_MAP, THETAS_MAP, EPS_NM, EPS_NM)
Ro_map = map_unpol(LAMS_MAP, THETAS_MAP, d1_star, d2_star)

def heatmap(Z, title, vmax):
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(THETAS_MAP, LAMS_MAP, Z.T, shading="auto", vmin=0, vmax=vmax)
    plt.xlabel("angle [deg]")
    plt.ylabel("wavelength [nm]")
    plt.title(title)
    plt.colorbar(label="reflectance [%]")
    plt.tight_layout()

heatmap(Rb_map, "Bare Si: unpolarized R(λ,θ)", vmax=50)
heatmap(Ro_map, f"2-layer AR: unpolarized R(λ,θ)  d1={d1_star:.1f} nm, d2={d2_star:.1f} nm", vmax=15)

plt.show()
]]></script>

</plask>
