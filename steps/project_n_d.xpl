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

# plask
L1 = GEO.layer_SiO2
L2 = GEO.layer_Si3N4
eps = 5e-4
side = "top"
lams   = np.arange(400.0, 700.1, 20.0)  # nm
angles = np.arange(0.0, 90.1, 10.0)      # deg

# GA
d_bounds = (10.0, 200.0)     # nm
# TOTAL_MAX = 500.0           # nm
n_bounds = (1.10, 3.40)     # (-)
P_size= 4
generations = 8
patience = 5
tol = 0.1
seed = 8




def set_layers(d1, d2, n1, n2):
    L1.height = max(d1, eps) * 1e-3
    L2.height = max(d2, eps) * 1e-3
    L1.nr = n1
    L2.nr = n2

def set_bare():
    set_layers(eps, eps, 1.0, 1.0)

def refl(lam, angle):
    OPTICAL.ktran = 2e3 * np.pi / lam * np.sin(np.pi / 180. * angle)
    if side == 'bottom':
        OPTICAL.ktran *= material.get('Si_Shinke').nr(lam)
    r_TE = OPTICAL.compute_reflectivity(lam, side, "TE")
    r_TM = OPTICAL.compute_reflectivity(lam, side, "TM")
    return 0.5 * (r_TE + r_TM)

def mean_R(d1, d2, n1, n2, lams=lams, thetas=angles):
    set_layers(d1, d2, n1, n2)
    sum_R = 0.0
    for th in thetas:
        for lam in lams:
            sum_R += refl(lam, th)
    return sum_R / (len(lams) * len(thetas))

def fitness(x):
    d1, d2, n1, n2 = map(float, x)

    if not (d_bounds[0] <= d1 <= d_bounds[1] and d_bounds[0] <= d2 <= d_bounds[1]):
        return np.inf
    # if TOTAL_MAX is not None and (d1 + d2) > TOTAL_MAX:
    #     return np.inf

    if not (n_bounds[0] <= n1 <= n_bounds[1] and n_bounds[0] <= n2 <= n_bounds[1]):
        return np.inf
    # if n1 > n2:
    #     return np.inf

    return mean_R(d1, d2, n1, n2)



# GA

class GA:
    def __init__(self, P_size):
        rng = np.random.default_rng(seed)
        d_0 = rng.uniform(d_bounds[0], d_bounds[1], [P_size, 2])
        n_0 = rng.uniform(n_bounds[0], n_bounds[1], [P_size, 2])
        self.t = 0
        self.curr_patience = 1
        self.y_star = 50
        self.P_t = np.zeros((P_size, 4))
        for i in range(P_size):
            self.P_t[i] = d_0[i][0], d_0[i][1], n_0[i][0], n_0[i][1]
        self.fitness_t = [fitness(c) for c in self.P_t]
        # if self.y_star - sorted(y_0)[0] > tol:
        #    self.curr_patience = 1
        # else:
        #     self.curr_patience += 1
        #self.y_star = y_0[0]
    
    def MP(self, alg="roulette"):
        # roulette wheel selection
        if alg == "roulette":
            fitness_sum = np.sum(self.fitness_t)
            p = [f/fitness_sum for f in self.fitness_t]
            id_P = np.random.choice(len(self.fitness_t), size=len(p), p=p)
            return [self.P_t[i] for i in id_P]
        else if alg == "tournament":
        # TODO: tournament selection
            pass
    
    def recombine(self, mating_pool, alg="convec"):
        if algo = "convex":
            pass
        else if alg == "average":
        # TODO: Average Crossover
            pass
        else if alg = "discrete":
        # TODO: Uniform Crossover
            pass              
        
    
    def run_GA(self, n_gen, patience, tol):
        while (self.t<n_gen and self.curr_patience < patience):
            P1 = GA.MP(self)
            P2 = GA.recombine(P1)
            # P3 = mutate(P2)
            # select(P3)
            self.t+=1
    

g = GA(P_size)
#print(g.t, g.curr_patience, g.y_star, g.P_t)
g.run_GA(generations, patience, tol)
    # t = 0
    # curr_patiance = 1
    # y_star = 50
    # P_t = np.zeros((P_size, 4))
    # for i in range(P_size):
    #     P_t[i] = d_0[i][0], d_0[i][1], n_0[i][0], n_0[i][1]
    # y_0 = sorted([fitness(c) for c in P_0])
    # if y_star - y_0 > tol:
    #     curr_patiance = 1
    # else:
    #     patience += 1
    # y_star = y_0[0]
    # while (t<generations and curr_patiance < patiance):
    #     P1 = MP(P_t)
    #     P2 = recombine(P1)
    #     P3 = mutate(P2)
    #     P_t = select(P3, P_t)
    #     t+=1
    






# def spectrum(lams, theta, d1, d2, n1, n2):
#     set_layers(d1, d2, n1, n2)
#     return np.array([refl(lam, theta) for lam in lams], float)
# 
# def map_R(lams, thetas, d1, d2, n1, n2):
#     set_layers(d1, d2, n1, n2)
#     Z = np.empty((thetas.size, lams.size), float)
#     for i, th in enumerate(thetas):
#         for j, lam in enumerate(lams):
#             Z[i, j] = refl(lam, th)
#     return Z
# 
# def heatmap(Z, title, vmax):
#     plt.figure(figsize=(8, 6))
#     plt.pcolormesh(angles, lams, Z.T, shading="auto", vmin=0, vmax=vmax)
#     plt.xlabel("angle [deg]")
#     plt.ylabel("wavelength [nm]")
#     plt.title(title)
#     plt.colorbar(label="reflectance [%]")
#     plt.tight_layout()
# 
# # baseline
# J_bare = mean_R(eps, eps, 1.0, 1.0)
# print_log("result", f"bare meanR (opt grid) = {J_bare:.6f}%")
# 
# J_bare_chk = mean_R(eps, eps, 1.0, 1.0, lams=lams, thetas=angles)
# print_log("result", f"meanR bare={J_bare_chk:.6f}% (check grid)")
# 
# 
# # heatmaps
# Rb_map = map_R(lams, angles, eps, eps, 1.0, 1.0)
# heatmap(Rb_map, "Bare Si: unpolarized R(λ,θ)", vmax=100)
# 
# plt.show()






# def differential_evolution(bounds, pop=POP, gens=GENS, F=F, CR=CR, seed=SEED):
#     rng = np.random.default_rng(seed)
#     lo = np.array([b[0] for b in bounds], float)
#     hi = np.array([b[1] for b in bounds], float)
#     dim = len(bounds)
# 
#     X = rng.uniform(lo, hi, size=(pop, dim))
#     fX = np.array([objective(x) for x in X], float)
# 
#     best_i = int(np.argmin(fX))
#     best, fbest = X[best_i].copy(), float(fX[best_i])
# 
#     for g in range(gens):
#         for i in range(pop):
#             idx = [j for j in range(pop) if j != i]
#             a, b, c = X[rng.choice(idx, 3, replace=False)]
#             v = np.clip(a + F * (b - c), lo, hi)
# 
#             jrand = rng.integers(dim)
#             mask = (rng.random(dim) < CR)
#             mask[jrand] = True
#             u = np.where(mask, v, X[i])
# 
#             fu = float(objective(u))
#             if fu < fX[i]:
#                 X[i], fX[i] = u, fu
#                 if fu < fbest:
#                     best, fbest = u.copy(), fu
# 
#         print_log("info",
#                   f"gen {g+1}/{gens}: d1={best[0]:.2f} nm, d2={best[1]:.2f} nm, "
#                   f"n1={best[2]:.3f}, n2={best[3]:.3f}, meanR={fbest:.6f}%")
#     return best, fbest

]]></script>

</plask>
