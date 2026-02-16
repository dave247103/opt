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
# this file plots mean refl with respect to angle for the optimized d1,d2,n1,n2

import numpy as np
import matplotlib.pyplot as plt
import os

# optimal params
d_star = [100.9, 51.3]  # nm
n_star = [1.41, 2.64] 

path = "ASTMG173.csv"
sun_data = np.genfromtxt(path, delimiter=',', skip_header=2)
SUN_WLS, SUN_INTENS = sun_data[:, 0], sun_data[:, 2]

class LayerMaterial(material.Material):
    def __init__(self, n):
        super().__init__()
        self._nr = n

    def absp(self, lam, T=300.):
        return 0
    
    def nr(self, lam, T=300., n=0.):
        return float(self._nr)

def set_stack(d_list, n_list):
    for i, (d, n) in enumerate(zip(d_list, n_list)):
        layer_name = f"layer{i+1}"
        layer = getattr(GEO, layer_name)
        layer.height = d / 1000.0
        layer.material = LayerMaterial(n)

lams_nm = np.arange(400.0, 700.1, 1.0)
angles = np.arange(0.0, 45.1, 1.0)
X, Y = np.meshgrid(lams_nm, angles)


def refl(polarization):
    R = np.empty((angles.size, lams_nm.size), float)
    k0_um = 2e3 * np.pi / lams_nm
    for i, th in enumerate(angles):
        s = np.sin(np.deg2rad(th))
        for j, lam in enumerate(lams_nm):
            OPTICAL.ktran = float(2e3 * np.pi / lam * s)
            R[i, j] = OPTICAL.compute_reflectivity(float(lam), "top", polarization)
    return R


set_stack(d_star, n_star)



R_TE = refl("TE")
R_TM = refl("TM")

w_lam = np.interp(lams_nm, SUN_WLS, SUN_INTENS)
w_ang = np.cos(np.radians(angles))
W_2D = w_ang[:, None] * w_lam[None, :]

def calc_weighted_mean(R):
    return np.sum(R * W_2D) / np.sum(W_2D)

avg_te = calc_weighted_mean(R_TE)
avg_tm = calc_weighted_mean(R_TM)

os.makedirs("plots", exist_ok=True)

def plot_map(R, pol, avg):
    plt.figure(figsize=(8, 5))
    cp = plt.pcolormesh(X, Y, R, shading='auto', cmap='magma', vmin=0, vmax=3)
    cb = plt.colorbar(cp)
    cb.set_label('Reflectivity [%]')
    
    plt.title(f"2-Layer AR Coating - {pol} Polarization\nSolar-Weighted Mean $R_w$ = {avg:.3f}%")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Angle of Incidence [°]")
    plt.tight_layout()
    plt.savefig(f"plots/2D_Map_{pol}.png", dpi=300)
    plt.show()

print(f"TE Weighted Mean: {avg_te:.3f}%")
print(f"TM Weighted Mean: {avg_tm:.3f}%")

plot_map(R_TE, "TE", avg_te)
plot_map(R_TM, "TM", avg_tm)

]]></script>

</plask>
