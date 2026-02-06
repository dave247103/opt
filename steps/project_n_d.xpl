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
Si = pm.get("Si_Shinke")

# plask
L1 = GEO.layer_SiO2
L2 = GEO.layer_Si3N4
eps = 5e-4
side = "top"
lams   = np.arange(400.0, 700.1, 10.0)  # nm
angles = np.arange(0.0, 90.1, 15.0)      # deg

lams_plot   = np.arange(400.0, 700.1, 2.0)  # nm
angles_plot = np.arange(0.0, 90.1, 2.0)      # deg


# GA
d_bounds = (10.0, 200.0)     # nm
# TOTAL_MAX = 500.0           # nm
n_bounds = (1.10, 3.40)     # (-)
P_size= 64
generations = 16
patience = 8
tol = 0.1
seed = 8
rng = np.random.default_rng(seed)





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
        OPTICAL.ktran *= float(np.real(Si.nr(lam)))
    r_TE = OPTICAL.compute_reflectivity(lam, side, "TE")
    r_TM = OPTICAL.compute_reflectivity(lam, side, "TM")
    return 0.5 * (r_TE + r_TM)

# standart mean_R
def s_mean_R(d1, d2, n1, n2, lams=lams, thetas=angles, w_lam=None, w_th=None):
    set_layers(d1, d2, n1, n2)
    sum_R = 0.0
    for th in thetas:
        for lam in lams:
            sum_R += refl(lam, th)
    return sum_R / (len(lams) * len(thetas))

# weigthed mean_R
def mean_R(d1, d2, n1, n2, lams=lams, thetas=angles, w_lam=None, w_th=None):
    set_layers(d1, d2, n1, n2)

    w_lam = np.ones_like(lams, float) if w_lam is None else np.asarray(w_lam, float)
    w_th  = np.ones_like(thetas, float) if w_th  is None else np.asarray(w_th,  float)

    w_lam = w_lam / w_lam.sum()
    w_th  = w_th  / w_th.sum()

    acc = 0.0
    for i, th in enumerate(thetas):
        for j, lam in enumerate(lams):
            acc += w_th[i] * w_lam[j] * refl(lam, th)
    return acc


w_lam = np.ones_like(lams, float)
w_th  = np.cos(np.deg2rad(angles)).clip(min=0.0)


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

    return mean_R(d1, d2, n1, n2, w_lam=w_lam, w_th=w_th)

lo = np.array([d_bounds[0], d_bounds[0], n_bounds[0], n_bounds[0]], float)
hi = np.array([d_bounds[1], d_bounds[1], n_bounds[1], n_bounds[1]], float)

def clamp(P):
    return np.clip(P, lo, hi)


# GA

class GA:
    def __init__(self, P_size, tol=0.1):
        d_0 = rng.uniform(d_bounds[0], d_bounds[1], [P_size, 2])
        n_0 = rng.uniform(n_bounds[0], n_bounds[1], [P_size, 2])
        self.t = 0
        self.curr_patience = 1
        self.y_star = 50
        self.tol = float(tol)
        self.P_t = np.column_stack([d_0, n_0]).astype(float)
        self.fitness_t = np.array([fitness(c) for c in self.P_t], float)
        self.y_star = float(self.fitness_t.min())
    
    def MP(self, alg="roulette"):
        # roulette wheel selection
        if alg == "roulette":
            f = np.asarray(self.fitness_t, float)
            w = 1.0 / (f + 1e-12)
            w[~np.isfinite(w)] = 0.0
            w /= w.sum()
            idx = rng.choice(len(f), size=len(f), replace=True, p=w)
            return np.asarray(self.P_t[idx], float)
        elif alg == "tournament":
        # TODO: tournament selection
            pass
    
    def recombine(self, mating_pool, p_c=0.8, alg="average", alpha=0.5):
        # crossover probability p is used in average and uniform crossover
        if alg == "convex":
            P2 = np.zeros_like(mating_pool)   
            u = mating_pool[-1]         
            for i, c in enumerate(mating_pool):
                P2[i] = [alpha*u[g]+(1-alpha)*mating_pool[i][g] for g in range(4)]
                u = mating_pool[i]   
            return P2             
        elif alg == "average":
            P = np.asarray(mating_pool, float).copy()
            rng.shuffle(P)
            C = P.copy()
        
            for i in range(0, len(P) - 1, 2):
                if rng.random() < p_c:
                    a = rng.random(4)  # per-gene mix
                    p1, p2 = P[i], P[i+1]
                    C[i]   = a*p1 + (1-a)*p2
                    C[i+1] = a*p2 + (1-a)*p1
            return clamp(C)
        elif alg == "discrete":
        # TODO: Uniform Crossover
            pass          
        
    def mutate(self, P2, p_m=0.1, sigma_d=5.0, sigma_n=0.05, p_reset=0.05, alg="flip"):
        if alg == "uniform":
            P2 = np.asarray(P2, float)
            idx = np.where(rng.random(len(P2)) < p_m)[0]
            if idx.size:
                d_m = rng.uniform(*d_bounds, size=(idx.size, 2))
                n_m = rng.uniform(*n_bounds, size=(idx.size, 2))
                P2[idx] = np.column_stack([d_m, n_m])
            return P2
        elif alg == "flip":
            P = np.asarray(P2, float).copy()
            m = rng.random(P.shape) < p_m
            noise = np.zeros_like(P)
            noise[:, :2] = rng.normal(0.0, sigma_d, size=P[:, :2].shape)  # nm
            noise[:, 2:] = rng.normal(0.0, sigma_n, size=P[:, 2:].shape)  # refr. idx
            P = np.where(m, P + noise, P)
        
            idx = np.where(rng.random(len(P)) < p_reset)[0]
            if idx.size:
                P[idx, :2] = rng.uniform(*d_bounds, size=(idx.size, 2))
                P[idx, 2:] = rng.uniform(*n_bounds, size=(idx.size, 2))
            return clamp(P)      
    
    
    # max pressure
    def m_select(self, P3, elite_frac=0.10, pressure=0.80):
        # Do I even select P3? even then do I take a portion of P3 to survive regardless of f or not??
        P = np.vstack([self.P_t, np.asarray(P3, float)])
        f = np.asarray([fitness(c) for c in P], float)
    
        idx = np.argsort(f)[:len(self.P_t)]
        self.P_t = P[idx]
        self.fitness_t = f[idx]
    
        best = float(self.fitness_t[0])
        self.curr_patience = 1 if (self.y_star - best) >self.tol else (self.curr_patience + 1)
        self.y_star = best

    #tuned pressure
    def select(self, P3, elite_frac=0.10, pressure=0.80):
        mu = len(self.P_t)
    
        P3 = clamp(P3)
        f3 = np.array([fitness(c) for c in P3], float)
    
        P = np.vstack([self.P_t, P3])
        f = np.concatenate([self.fitness_t, f3])
    
        order = np.argsort(f)  # ascending (minimization)
    
        k_elite = max(1, int(round(elite_frac * mu)))
        k_det   = int(round(pressure * (mu - k_elite)))
    
        elite = order[:k_elite]
        det   = order[k_elite:k_elite + k_det]
        rest  = order[k_elite + k_det:]
    
        k_fill = mu - (k_elite + k_det)
        if k_fill > 0 and rest.size:
            # low-pressure stochastic fill (rank-based)
            ranks = np.arange(1, rest.size + 1, dtype=float)
            w = 1.0 / ranks
            w /= w.sum()
            fill = rng.choice(rest, size=k_fill, replace=False if rest.size >= k_fill else True, p=w)
            new_idx = np.concatenate([elite, det, fill])
        else:
            new_idx = np.concatenate([elite, det])
    
        Q = P[new_idx]
        Q = np.unique(np.round(Q, 6), axis=0)   # de-dup with tolerance
        
        if len(Q) < mu:
            k = mu - len(Q)
            d_new = rng.uniform(*d_bounds, size=(k, 2))
            n_new = rng.uniform(*n_bounds, size=(k, 2))
            Q = np.vstack([Q, np.column_stack([d_new, n_new])])
        
        Q = Q[:mu]
        self.P_t = Q
        self.fitness_t = np.array([fitness(c) for c in self.P_t], float)
        
        best = float(self.fitness_t.min())
        self.curr_patience = 0 if (self.y_star - best) > self.tol else (self.curr_patience + 1)
        self.y_star = best

        
            
    def run_GA(self, n_gen, patience, p_c=0.8, p_m=0.10, elite_frac=0.10, pressure=0.80):
        while self.t < n_gen and self.curr_patience < patience:
            P1 = self.MP()
            P2 = self.recombine(P1, p_c=p_c)
            P3 = self.mutate(P2, p_m=p_m)
            self.select(P3, elite_frac=elite_frac, pressure=pressure)
            self.t += 1
            print(self.y_star)


g = GA(P_size, tol=tol)
g.run_GA(generations, patience, elite_frac=0.05, pressure=0.60, p_m=0.15)
print("the sol: ",g.y_star)
print(g.P_t)
print(g.fitness_t)
idx = np.argsort(g.fitness_t)
P_t = g.P_t[idx]
fitness_t = g.fitness_t[idx]
print("the params: ", P_t)
print("lowest R: ", fitness_t)



def spectrum(lams, theta, d1, d2, n1, n2):
    set_layers(d1, d2, n1, n2)
    return np.array([refl(lam, theta) for lam in lams], float)

def map_R(lams, thetas, d1, d2, n1, n2):
    set_layers(d1, d2, n1, n2)
    Z = np.empty((thetas.size, lams.size), float)
    for i, th in enumerate(thetas):
        for j, lam in enumerate(lams):
            Z[i, j] = refl(lam, th)
    return Z

def heatmap(lams, thetas, Z, title, vmax):
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(thetas, lams, Z.T, shading="auto", vmin=0, vmax=vmax)
    plt.xlabel("angle [deg]")
    plt.ylabel("wavelength [nm]")
    plt.title(title)
    plt.colorbar(label="reflectance [%]")
    plt.tight_layout()

w_lam_plot = np.ones_like(lams_plot, float)
w_th_plot  = np.cos(np.deg2rad(angles_plot)).clip(min=0.0)

# baseline
J_bare = mean_R(eps, eps, 1.0, 1.0, lams=lams_plot, thetas=angles_plot)
print_log("result", f"bare meanR (plot grid) = {J_bare:.6f}%")

d1_star, d2_star, n1_star, n2_star = g.P_t[0]
J_opt = mean_R(d1_star, d2_star, n1_star, n2_star, lams=lams_plot, thetas=angles_plot)
print_log("result", f"opt  meanR (plot grid) = {J_opt:.6f}%")

J_bare = mean_R(eps, eps, 1.0, 1.0, lams=lams_plot, thetas=angles_plot, w_lam=w_lam_plot, w_th=w_th_plot)
print_log("result", f"bare meanR (weighted, plot grid) = {J_bare:.6f}%")

d1_star, d2_star, n1_star, n2_star = g.P_t[0]
J_opt = mean_R(d1_star, d2_star, n1_star, n2_star, lams=lams_plot, thetas=angles_plot, w_lam=w_lam_plot, w_th=w_th_plot)
print_log("result", f"opt  meanR (weighted, plot grid) = {J_opt:.6f}%")



# heatmaps
Rb_map = map_R(lams_plot, angles_plot, eps, eps, 1.0, 1.0)
heatmap(lams_plot, angles_plot, Rb_map, "Bare Si: unpolarized R(λ,θ)", vmax=100)

R_opt_map = map_R(lams_plot, angles_plot, d1_star, d2_star, n1_star, n2_star)
heatmap(lams_plot, angles_plot, R_opt_map, "Optimized unpolarized R(λ,θ)", vmax=100)

plt.show()


]]></script>

</plask>
