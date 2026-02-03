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
import pygad

L1, L2 = GEO.layer_SiO2, GEO.layer_Si3N4
SIDE = "top"
EPS_NM = 5e-4

D_BOUNDS = (10.0, 200.0)     # nm
N_BOUNDS = (1.10, 3.40)      # refractive index
TOTAL_MAX = 500.0            # nm

LAMS_OPT   = np.arange(400.0, 700.1, 10.0)
THETAS_OPT = np.arange(0.0, 90.1, 5.0)

LAMS_MAP   = np.arange(400.0, 700.1, 5.0)
THETAS_MAP = np.arange(0.0, 90.1, 2.0)

POP, GENS, SEED = 16, 8, 1

def set_layers(d1_nm, d2_nm, n1, n2):
    L1.height = max(float(d1_nm), EPS_NM) * 1e-3  # um
    L2.height = max(float(d2_nm), EPS_NM) * 1e-3
    L1.nr = float(n1)
    L2.nr = float(n2)

def ktran_um(lam_nm, theta_deg):
    return float(2e3 * np.pi / lam_nm * np.sin(np.deg2rad(theta_deg)))  # 1/um

def R_unpol(lam_nm, theta_deg):
    OPTICAL.ktran = ktran_um(lam_nm, theta_deg)
    rte = float(OPTICAL.compute_reflectivity(float(lam_nm), SIDE, "TE"))
    rtm = float(OPTICAL.compute_reflectivity(float(lam_nm), SIDE, "TM"))
    return 0.5 * (rte + rtm)

def mean_R(d1_nm, d2_nm, n1, n2):
    set_layers(d1_nm, d2_nm, n1, n2)
    acc = 0.0
    for th in THETAS_OPT:
        for lam in LAMS_OPT:
            acc += R_unpol(lam, th)
    return acc / (len(LAMS_OPT) * len(THETAS_OPT))

def fitness_func(ga, sol, sol_idx):
    d1, d2, n1, n2 = map(float, sol)

    if (d1 + d2) > TOTAL_MAX:
        return -1e9

    n1, n2 = sorted((n1, n2))  # enforce n1 <= n2

    J = mean_R(d1, d2, n1, n2)  # minimization target (in %)
    return -J                   # PyGAD maximizes fitness

gene_space = [
    {"low": D_BOUNDS[0], "high": D_BOUNDS[1]},
    {"low": D_BOUNDS[0], "high": D_BOUNDS[1]},
    {"low": N_BOUNDS[0], "high": N_BOUNDS[1]},
    {"low": N_BOUNDS[0], "high": N_BOUNDS[1]},
]

ga = pygad.GA(
    num_generations=GENS,
    sol_per_pop=POP,
    num_parents_mating=POP // 2,
    num_genes=4,
    gene_space=gene_space,
    gene_type=float,
    fitness_func=fitness_func,          # (ga_instance, solution, solution_idx)
    parent_selection_type="sss",
    crossover_type="single_point",
    mutation_type="random",
    mutation_probability=0.02,
    keep_elitism=1,
    random_seed=SEED,
)

# baseline
J_bare = mean_R(EPS_NM, EPS_NM, 1.0, 1.0)
print_log("result", f"bare meanR = {J_bare:.6f}%")

ga.run()
sol, fit, _ = ga.best_solution()
d1, d2, n1, n2 = map(float, sol)
n1, n2 = sorted((n1, n2))
print_log("result", f"best: d1={d1:.3f} nm, d2={d2:.3f} nm, n1={n1:.4f}, n2={n2:.4f}, meanR={-fit:.6f}%")

# plots
def map_R(d1_nm, d2_nm, n1, n2):
    set_layers(d1_nm, d2_nm, n1, n2)
    Z = np.empty((THETAS_MAP.size, LAMS_MAP.size), float)
    for i, th in enumerate(THETAS_MAP):
        for j, lam in enumerate(LAMS_MAP):
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

Z_bare = map_R(EPS_NM, EPS_NM, 1.0, 1.0)
Z_opt  = map_R(d1, d2, n1, n2)
heatmap(Z_bare, "Bare Si: unpolarized R(λ,θ)", vmax=60)
heatmap(Z_opt, f"2-layer AR: d1={d1:.1f} nm, d2={d2:.1f} nm, n1={n1:.3f}, n2={n2:.3f}", vmax=25)
plt.show()

]]></script>

</plask>
