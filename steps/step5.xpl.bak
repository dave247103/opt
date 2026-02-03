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
  <module name="si_nk"/>
  <material name="mySiO2" base="SiO2">
    <Nr>n_layer</Nr>
  </material>
</materials>

<geometry>
  <cartesian2d name="simple" axes="x,y " left="periodic" right="periodic" bottom="Si_Shinke">
    <rectangle name="coating" material="SiO2" dx="1" dy="{lam_min * 1e-3 / (4 * n_layer)}"/>
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

lams = np.arange(400.0, 700.1, 2.0)

coating = GEO.coating

def set_thickness_nm(d_nm: float):
    coating.height = d_nm * 1e-3

def R_spectrum(d_nm: float):
    set_thickness_nm(d_nm)
    R = OPTICAL.compute_reflectivity(lams, 'top', 'TE')  # normal incidence => TE==TM
    return R

def objective_mean_R(d_nm: float) -> float:
    return float(np.mean(R_spectrum(d_nm)))    # mean R


ds = np.linspace(10.0, 200.0, 191)
R_avg = np.array([objective_mean_R(d) for d in ds])
d0 = float(ds[np.argmin(R_avg)])
R_avg0 = float(R_avg.min())
print_log('result', f"coarse best: d={d0:.2f} nm, R_avg={R_avg0:.4f}%")



def section_min(f, a, b, tol=1e-3, max_iter=200):
    """Golden-section search on [a,b]."""
    gr = (np.sqrt(5) - 1) / 2  # 0.618...
    c = b - gr*(b-a)
    d = a + gr*(b-a)
    fc, fd = f(c), f(d)
    for _ in range(max_iter):
        if abs(b-a) < tol:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr*(b-a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr*(b-a)
            fd = f(d)
    x = (a+b)/2
    return float(x), float(f(x))


# bracket around the coarse best
a = max(10.0, d0 - 20.0)
b = min(200.0, d0 + 20.0)
d_star, R_star = section_min(objective_mean_R, a, b, tol=1e-2)
print_log('result', f"refined best: d={d_star:.3f} nm, R={R_star:.6f}%")


R_opt = R_spectrum(d_star)
d_qw = lam_min/(4*n_layer)
R_qw =  R_spectrum(d_qw)

plt.figure()
plt.plot(lams, R_qw,  label=f"quarter-wave @ {lam_min:.0f} nm (d={d_qw:.2f} nm)")
plt.plot(lams, R_opt, label=f"optimized mean (d={d_star:.2f} nm)")
plt.xlabel("wavelength [nm]")
plt.ylabel("reflectance [%]")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.show()

# overkill for 1d
# def differential_evolution_1d(f, dmin, dmax, pop=20, F=0.7, CR=0.9, iters=60, seed=0):
#     rng = np.random.default_rng(seed)
#     X = rng.uniform(dmin, dmax, size=pop)
#     FX = np.array([f(x) for x in X], dtype=float)
# 
#     for _ in range(iters):
#         for i in range(pop):
#             idx = [j for j in range(pop) if j != i]
#             r1, r2, r3 = rng.choice(idx, size=3, replace=False)
#             v = X[r1] + F*(X[r2] - X[r3])
#             v = np.clip(v, dmin, dmax)
# 
#             u = v if rng.random() < CR else X[i]
#             fu = f(u)
#             if fu < FX[i]:
#                 X[i], FX[i] = u, fu
# 
#     best = int(np.argmin(FX))
#     return float(X[best]), float(FX[best])

# d_star, R_star = differential_evolution_1d(objective_mean_R, 10.0, 200.0)
# print_log('result', f"DE best: d={d_star:.3f} nm, R={R_star:.6f}%")

]]></script>

</plask>
