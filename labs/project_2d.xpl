<plask loglevel="result">

<defines>
  <define name="nref" value="4.0"/>
  <define name="lam_min" value="550"/>
  <define name="n_layer" value="nref**0.5"/>
</defines>

<materials>
  <material name="mySi" base="Si">
    <absp>0</absp>
    <nr>nref</nr>
  </material>
  <material name="mySi02" base="SiO2">
    <Nr>n_layer</Nr>
  </material>
  <module name="si_nk"/>
</materials>

<geometry>
  <cartesian2d name="main" axes="x,y" left="periodic" right="periodic" bottom="Si_Shinke">
    <stack>
      <rectangle name="layer1" material="mySi02" dx="1" dy="1"/>
      <rectangle name="layer2" material="mySi02" dx="1" dy="1"/>
    </stack>
  </cartesian2d>
</geometry>

<solvers>
  <optical name="OPTICAL" solver="Fourier2D" lib="modal">
    <geometry ref="main"/>
    <expansion size="0"/>
  </optical>
</solvers>

<script><![CDATA[
import numpy as np
import matplotlib.pyplot as plt
import os

N_LAYERS = 2

# opt
lams = np.arange(400.0, 700.1, 1.0)
eps = 5e-4
side = "top"

# GA
d_bounds = (20.0, 200.0)     # nm
n_bounds = (1.16, 3.00)
P_size= 8
generations = 2
patience = 8
tol = 0.1
seed = 8
rng = np.random.default_rng(seed)

lo = np.concatenate([np.full(N_LAYERS, d_bounds[0]), np.full(N_LAYERS, n_bounds[0])])
hi = np.concatenate([np.full(N_LAYERS, d_bounds[1]), np.full(N_LAYERS, n_bounds[1])])


# plots
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# import ASTM G173 data
path = os.path.join("ASTMG173.csv")
sun_data = np.genfromtxt(path, delimiter=',', skip_header=2)
SUN_WLS, SUN_INTENS = sun_data[:, 0], sun_data[:, 2]

def get_solar_weight(lams):
    return np.interp(lams, SUN_WLS, SUN_INTENS)
weights = get_solar_weight(lams)



class LayerMaterial(material.Material):
    def __init__(self, n):
        super().__init__()
        self._nr = n

    def absp(self, lam, T=300.):
        return 0
    
    def nr(self, lam, T=300., n=0.):
        return float(self._nr)

# def set_layers(d1, d2, n1, n2):
#     GEO.layer1.height = d1 / 1000
#     GEO.layer2.height = d2 / 1000
#     GEO.layer1.material = LayerMaterial(n1)
#     GEO.layer2.material = LayerMaterial(n2)
    
def set_stack(d_list, n_list):
    for i, (d, n) in enumerate(zip(d_list, n_list)):
        layer_name = f"layer{i+1}"
        layer = getattr(GEO, layer_name)
        layer.height = d / 1000.0
        layer.material = LayerMaterial(n)

def refls(lams):
    r_TE = OPTICAL.compute_reflectivity(lams, side, "TE")
    return r_TE 

def mean_R(d_arr, n_arr, lams=lams):
    set_stack(d_arr, n_arr)
    R = np.asarray(refls(lams), float)
    return float(np.sum(R * weights) / np.sum(weights))

def fitness(x):
    mid = len(x) // 2
    d_arr = x[:mid]
    n_arr = x[mid:]
    if np.any(d_arr < d_bounds[0]) or np.any(d_arr > d_bounds[1]): return np.inf
    if np.any(n_arr < n_bounds[0]) or np.any(n_arr > n_bounds[1]): return np.inf
    return mean_R(d_arr, n_arr)


def clamp(P):
    return np.clip(P, lo, hi)

class GA:
    def __init__(self, P_size, tol=0.1):
        d_0 = rng.uniform(d_bounds[0], d_bounds[1], [P_size, N_LAYERS])
        n_0 = rng.uniform(n_bounds[0], n_bounds[1], [P_size, N_LAYERS])
        self.t = 0
        self.curr_patience = 1
        self.y_star = 50
        self.tol = float(tol)
        self.P_t = np.column_stack([d_0, n_0]).astype(float)
        self.fitness_t = np.array([fitness(c) for c in self.P_t], float)
        self.y_star = float(self.fitness_t.min())
    
    def MP(self, alg="roulette"):
        if alg == "roulette":   # roulette wheel selection
            f = np.asarray(self.fitness_t, float)
            w = 1.0 / (f + 1e-12)
            w[~np.isfinite(w)] = 0.0
            w /= w.sum()
            idx = rng.choice(len(f), size=len(f), replace=True, p=w)
            return np.asarray(self.P_t[idx], float)
    
    def recombine(self, mating_pool, p_c=0.8, alg="average", alpha=0.5):
        # crossover probability p_c is used in average and uniform crossover
        if alg == "convex":
            P2 = np.zeros_like(mating_pool)   
            u = mating_pool[-1]
            n_genes = mating_pool.shape[1]  # Get dynamic length
            for i, c in enumerate(mating_pool):
                P2[i] = [alpha*u[g]+(1-alpha)*mating_pool[i][g] for g in range(n_genes)]
                u = mating_pool[i]   
            return clamp(P2)         
        elif alg == "average":
            P = np.asarray(mating_pool, float).copy()
            rng.shuffle(P)
            C = P.copy()
            n_genes = P.shape[1]
        
            for i in range(0, len(P) - 1, 2):
                if rng.random() < p_c:
                    a = rng.random(n_genes)
                    p1, p2 = P[i], P[i+1]
                    C[i]   = a*p1 + (1-a)*p2
                    C[i+1] = a*p2 + (1-a)*p1
            return clamp(C)  
        
    def mutate(self, P2, p_m=0.1, sigma_d=5.0, sigma_n=0.05, alg="flip"):
        if alg == "uniform":
            P2 = np.asarray(P2, float)
            idx = np.where(rng.random(len(P2)) < p_m)[0]
            if idx.size:
                N = P2.shape[1] // 2 
                d_m = rng.uniform(*d_bounds, size=(idx.size, N))
                n_m = rng.uniform(*n_bounds, size=(idx.size, N))
                P2[idx] = np.column_stack([d_m, n_m])
            return P2
        elif alg == "flip":
            P = np.asarray(P2, float).copy()
            m = rng.random(P.shape) < p_m
            noise = np.zeros_like(P)
            mid = P.shape[1] // 2 
            noise[:, :mid] = rng.normal(0.0, sigma_d, size=P[:, :mid].shape)  # nm
            noise[:, mid:] = rng.normal(0.0, sigma_n, size=P[:, mid:].shape)  # refractive idx
            P = np.where(m, P + noise, P)
            return clamp(P)



    def select(self, P3, elite_frac=0.10, pressure=0.80):
        mu = len(self.P_t)
        P3 = clamp(P3)
        f3 = np.array([fitness(c) for c in P3], float)
        P = np.vstack([self.P_t, P3])
        f = np.concatenate([self.fitness_t, f3])
        order = np.argsort(f)    
        k_elite = max(1, int(round(elite_frac * mu)))
        k_det   = int(round(pressure * (mu - k_elite)))
        
        new_idx = order[:k_elite + k_det]

    
        if len(new_idx) < mu:
            rest = order[k_elite + k_det:]
            fill = rng.choice(rest, size=mu - len(new_idx), replace=True)
            new_idx = np.concatenate([new_idx, fill])
            
        Q = P[new_idx]

        self.P_t = Q[:mu]
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
g.run_GA(generations, patience, elite_frac=0.2, p_m=0.05)
print("Best Solution R:", g.y_star)
best_chrome = g.P_t[np.argmin(g.fitness_t)]
n_layers = len(best_chrome) // 2
d_star = best_chrome[:n_layers]
n_star = best_chrome[n_layers:]
print("Optimal d:", d_star)
print("Optimal n:", n_star)



def save_plot(d_opt, n_opt, lams, filename_prefix="2_layers"):
    set_stack(d_opt, n_opt)
    R_curve = refls(lams)
    R_mean = float(np.sum(R_curve * weights) / np.sum(weights))

    d_str = "[" + ", ".join([f"{x:.1f}" for x in d_opt]) + "]"
    n_str = "[" + ", ".join([f"{x:.2f}" for x in n_opt]) + "]"
    
    n_layers = len(d_opt)
    
    title_main = f"{n_layers} Layers Optimization $R_{{w}}$ = {R_mean:.3f}%"
    param_str = f"d={d_str} nm\nn={n_str}"
    
    plt.figure(figsize=(10, 6))
    plt.plot(lams, R_curve, linewidth=2, label="Optimized R")
    plt.axhline(0, color='k', linewidth=0.5)
    
    plt.xlabel(r"$\lambda$ [nm]", fontsize=12)
    plt.ylabel("Reflectivity [%]", fontsize=12)
    plt.grid(True, alpha=0.3)

    plt.title(title_main, fontsize=13)
    plt.plot([], [], ' ', label=param_str)
    plt.legend(loc='upper right', framealpha=0.95)

    filename = f"{output_dir}/{filename_prefix}_optimized.png"
    plt.savefig(filename, dpi=300)
    print(f"Plot saved to {filename}")
    plt.close()

save_plot(d_star, n_star, lams, filename_prefix=f"{n_layers}_layers")
]]></script>

</plask>
