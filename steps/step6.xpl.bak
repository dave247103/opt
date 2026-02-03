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

OPT_STEP = 6.0          # nm; 2=accurate, 6-10=faster
POP = 16
GENS = 20
SEED = 1

LAMS_FINE = np.arange(400.0, 700.1, 2.0)
LAMS_OPT  = np.arange(400.0, 700.1, OPT_STEP)

D_BOUNDS  = (10.0, 200.0)   # nm
TOTAL_MAX = 300.0           # nm or None
POL = "TE"                 

L1 = GEO.layer_SiO2
L2 = GEO.layer_Si3N4

def set_thicknesses_nm(d1_nm, d2_nm):
    L1.height = float(d1_nm) * 1e-3  # um
    L2.height = float(d2_nm) * 1e-3

def R_mean_on(lams, d1_nm, d2_nm):
    set_thicknesses_nm(d1_nm, d2_nm)
    R = np.asarray(OPTICAL.compute_reflectivity(lams, "top", POL), dtype=float)
    return float(np.mean(R))

def objective(x):
    d1, d2 = float(x[0]), float(x[1])
    lo, hi = D_BOUNDS
    if d1 < lo or d1 > hi or d2 < lo or d2 > hi:
        return np.inf
    if TOTAL_MAX is not None and (d1 + d2) > TOTAL_MAX:
        return np.inf
    return R_mean_on(LAMS_OPT, d1, d2)

def differential_evolution(bounds, pop=POP, gens=GENS, F=0.8, CR=0.9, seed=SEED):
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

        if (g == 0) or ((g + 1) % 5 == 0):
            print_log("info", f"gen {g+1}/{gens}: best d1={best[0]:.2f} nm, d2={best[1]:.2f} nm, meanR_opt={fbest:.6f}%")

    return best, fbest

# baseline + optimize
set_thicknesses_nm(100.0, 100.0)
R_init = np.asarray(OPTICAL.compute_reflectivity(LAMS_FINE, "top", POL), dtype=float)
print_log("result", f"init meanR={float(np.mean(R_init)):.6f}%")

best, Jopt = differential_evolution([D_BOUNDS, D_BOUNDS])
d1_star, d2_star = float(best[0]), float(best[1])

# evaluate on fine grid
set_thicknesses_nm(d1_star, d2_star)
R_opt = np.asarray(OPTICAL.compute_reflectivity(LAMS_FINE, "top", POL), dtype=float)
Jfine = float(np.mean(R_opt))
print_log("result", f"best: d_SiO2={d1_star:.3f} nm, d_Si3N4={d2_star:.3f} nm, meanR_fine={Jfine:.6f}% (opt grid step {OPT_STEP:g} nm)")

plt.figure()
plt.plot(LAMS_FINE, R_init, label="init (100 nm, 100 nm)")
plt.plot(LAMS_FINE, R_opt,  label=f"opt (d1={d1_star:.1f} nm, d2={d2_star:.1f} nm)")
plt.xlabel("wavelength [nm]")
plt.ylabel("reflectance [%]")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.show()
]]></script>

</plask>
