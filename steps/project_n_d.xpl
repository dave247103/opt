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
    <expansion size="8"/>
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
lams   = np.arange(400.0, 700.1, 2.0)  # nm
angles = np.arange(0.0, 90.1, 3.0)      # deg

lams_plot   = np.arange(400.0, 700.1, 1.0)  # nm
angles_plot = np.arange(0.0, 90.1, 1.0)      # deg


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
        d_0 = rng.uniform(d_bounds[0], d_bounds[1], [P_size, 2])
        n_0 = rng.uniform(n_bounds[0], n_bounds[1], [P_size, 2])
        self.t = 0
        self.curr_patience = 1
        self.y_star = 50
        self.P_t = np.zeros((P_size, 4))
        for i in range(P_size):
            self.P_t[i] = d_0[i][0], d_0[i][1], n_0[i][0], n_0[i][1]
        self.P_t = np.asarray(self.P_t)
        self.fitness_t = [fitness(c) for c in self.P_t]
        self.fitness_t = np.asarray(self.fitness_t)
        self.y_star = sorted(self.fitness_t)[0]
    
    def MP(self, alg="roulette"):
        # roulette wheel selection
        if alg == "roulette":
            fitness_sum = np.sum(self.fitness_t)
            p = [f/fitness_sum for f in self.fitness_t]
            id_P = np.random.choice(len(self.fitness_t), size=len(p), p=p)
            return [self.P_t[i] for i in id_P]
        elif alg == "tournament":
        # TODO: tournament selection
            pass
    
    def recombine(self, mating_pool, p=0.5, alg="convex", alpha=0.5):
        # crossover probability p is used in average and uniform crossover
        if alg == "convex":
            P2 = np.zeros_like(mating_pool)   
            u = mating_pool[-1]         
            for i, c in enumerate(mating_pool):
                P2[i] = [alpha*u[g]+(1-alpha)*mating_pool[i][g] for g in range(4)]
                u = mating_pool[i]   
            return P2             
        elif alg == "average":
        # TODO: Average Crossover
            pass
        elif alg == "discrete":
        # TODO: Uniform Crossover
            pass          
        
    def mutate(self, P2, p=0.2, alg="uniform"):
        if alg == "uniform":
           ps = rng.random(len(P2))
           idx = np.where(ps > p)[0]
           d_m = rng.uniform(d_bounds[0], d_bounds[1], [len(idx), 2])
           n_m = rng.uniform(n_bounds[0], n_bounds[1], [len(idx), 2])
           for j, i in enumerate(idx):
               P2[i] = d_m[j][0], d_m[j][1], n_m[j][0], n_m[j][1]
           return P2
        elif alg == "flip":
        #TODO: flip-flop mutation
            pass           
    
    def select(self, P3, pressure=0.5):
        self.P_t = np.asarray(self.P_t)
        self.fitness_t = np.asarray(self.fitness_t)
        idx = np.argsort(self.fitness_t)
        self.P_t = self.P_t[idx]
        self.fitness_t = self.fitness_t[idx]
        l = int(len(self.P_t)*pressure)
        if l%2 != 0:
            l+=1            
        self.P_t = (self.P_t)[:l]
        self.P_t = np.concatenate((self.P_t, P3[:l]))
        rng.shuffle(self.P_t)
        self.fitness_t = [fitness(c) for c in self.P_t]
        if self.y_star - sorted(self.fitness_t)[0] > tol:
           self.curr_patience = 1
        else:
            self.curr_patience += 1
        self.y_star = sorted(self.fitness_t)[0]
        
            
    def run_GA(self, n_gen, patience, tol):
        while (self.t<n_gen and self.curr_patience < patience):
            print(self.y_star)
            P1 = GA.MP(self)
            P2 = GA.recombine(self, P1)
            P3 = GA.mutate(self, P2, 0.5)
            GA.select(self, P3, 0.5)
            self.t+=1
        idx = np.argsort(self.fitness_t)
        self.P_t = self.P_t[idx]

g = GA(P_size)
print(g.y_star)
g.run_GA(generations, patience, tol)
print(g.y_star)




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

def heatmap(Z, title, vmax):
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(angles, lams, Z.T, shading="auto", vmin=0, vmax=vmax)
    plt.xlabel("angle [deg]")
    plt.ylabel("wavelength [nm]")
    plt.title(title)
    plt.colorbar(label="reflectance [%]")
    plt.tight_layout()

# baseline
J_bare = mean_R(eps, eps, 1.0, 1.0)
print_log("result", f"bare meanR (opt grid) = {J_bare:.6f}%")
d1_star, d2_star, n1_star, n2_star = g.P_t[0]
J_opt = mean_R(d1_star, d2_star, n1_star, n2_star)
print_log("result", f"bare meanR (opt grid) = {J_opt:.6f}%")


# heatmaps
# Rb_map = map_R(lams, angles, eps, eps, 1.0, 1.0)
# heatmap(Rb_map, "Bare Si: unpolarized R(λ,θ)", vmax=100)
# R_opt_map = map_R(lams_plot, angles_plot, d1_star, d2_star, n1_star, n2_star)
# heatmap(R_opt_map, "Optimized unpolarized R(λ,θ)", vmax=50)
# 
# plt.show()

]]></script>

</plask>
