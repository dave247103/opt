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
      <rectangle name="layer2" material="Si_Shinke" dx="1" dy="1"/>
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


# opt
lams = np.arange(400.0, 700.1, 1.0)
eps = 5e-4
side = "top"

# GA
d_bounds = (10.0, 200.0)     # nm
n_bounds = (1.10, 3.40)     # (-)
P_size= 10
generations = 3
patience = 8
tol = 0.1
seed = 8
rng = np.random.default_rng(seed)
lo = np.array([d_bounds[0], d_bounds[0], n_bounds[0], n_bounds[0]], float)
hi = np.array([d_bounds[1], d_bounds[1], n_bounds[1], n_bounds[1]], float)


path = os.path.join("ASTMG173.csv")
sun_data = np.genfromtxt(path, delimiter=',', skip_header=2)
SUN_WLS, SUN_INTENS = sun_data[:, 0], sun_data[:, 2]

def get_solar_weight(lams):
    return np.interp(lams, SUN_WLS, SUN_INTENS)
weights = get_solar_weight(lams)

def clamp(P):
    return np.clip(P, lo, hi)


class LayerMaterial(material.Material):
    
    def __init__(self, nr):
        super().__init__()
        self._nr = nr

    def abs(self, lam, T=300.):
        return 0
    
    def nr(self, lam, T=300., n=0.):
        return self._nr
    

def set_layers(d1, d2, n1, n2):
    GEO.layer1.height = d1
    GEO.layer2.height = d2
    # GEO.layer1.nr = n1
    # GEO.layer2.nr = n2
    GEO.layer1.material = LayerMaterial(n1)
    GEO.layer2.material = LayerMaterial(n2)


def set_bare():
    set_layers(eps, eps, 1.0, 1.0)


def refls(lams):
    r_TE = OPTICAL.compute_reflectivity(lams, side, "TE")
    return r_TE
    
    
def mean_R(d1, d2, n1, n2, lams=lams):
    set_layers(d1, d2, n1, n2)
    R = np.asarray(refls(lams), float)
    return float(np.sum(R * weights) / np.sum(weights))


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


# baseline
J_bare = mean_R(eps, eps, 1.0, 1.0, lams=lams)
print_log("result", f"bare meanR (plot grid) = {J_bare:.6f}%")

d1_star, d2_star, n1_star, n2_star = g.P_t[0]
J_opt = mean_R(d1_star, d2_star, n1_star, n2_star, lams=lams)
print_log("result", f"opt  meanR (plot grid) = {
J_opt:.6f}%")
    

def plot2d(d1, d2, n1, n2, lams, title):
    
    set_layers(d1, d2, n1, n2)
    plt.figure(figsize=(8, 6))
    plt.plot(lams, refls(lams), label=f"n1={n1:.3f}, n2={n2:.3f}")
    plt.legend()
    plt.xlabel("wavelength [nm]")
    plt.ylabel("Reflectivity (%)")
    plt.title(title)
    plt.tight_layout()

plot2d(eps, eps, 1.0, 1.0, lams, "Bare Si: TE polarized R(λ,θ)")
plot2d(d1_star, d2_star, n1_star, n2_star, lams, "Optimized TE polarized R(λ,θ)")

plt.show()
]]></script>

</plask>
