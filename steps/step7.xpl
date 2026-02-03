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

LAMS   = np.arange(400.0, 700.1, 2.0)   # nm
THETAS = np.arange(0.0, 70.1, 2.0)      # deg
SIDE   = "top"                          # incidence from air

D1_NM = 27.568   # SiO2 (top)
D2_NM = 53.222   # Si3N4 (below)
EPS_NM = 5e-4                           # nonzero to keep mesh valid


L1 = GEO.layer_SiO2
L2 = GEO.layer_Si3N4

def set_stack(d1_nm, d2_nm):
    L1.height = float(d1_nm) * 1e-3 
    L2.height = float(d2_nm) * 1e-3

def set_bare():
    set_stack(EPS_NM, EPS_NM)

def maps_TE_TM():
    Rte = np.empty((THETAS.size, LAMS.size), float)
    Rtm = np.empty((THETAS.size, LAMS.size), float)

    k0_um = 2e3 * np.pi / LAMS          # 2π/λ(um), with λ in nm

    for i, th in enumerate(THETAS):
        s = np.sin(np.deg2rad(th))
        for j, lam in enumerate(LAMS):
            OPTICAL.ktran = float(2e3 * np.pi / lam * s)   # 1/um, n0≈1
            Rte[i, j] = OPTICAL.compute_reflectivity(float(lam), SIDE, "TE")
            Rtm[i, j] = OPTICAL.compute_reflectivity(float(lam), SIDE, "TM")

    return Rte, Rtm, 0.5 * (Rte + Rtm)

def heatmap(Z, title, vmax=45):
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(THETAS, LAMS, Z.T, shading="auto", vmin=0, vmax=vmax)
    plt.xlabel("Angle (deg)")
    plt.ylabel("Wavelength (nm)")
    plt.title(title)
    plt.colorbar(label="Reflectivity (%)")
    plt.tight_layout()

# Si
set_bare()
_, _, Ru_bare = maps_TE_TM()

# 2 layers
set_stack(D1_NM, D2_NM)
_, _, Ru_stack = maps_TE_TM()

heatmap(Ru_bare,  "Bare Si: unpolarized R(λ,θ)", vmax=45)
heatmap(Ru_stack, f"2 layer AR: unpolarized R(λ,θ)  d1={D1_NM:.1f} nm, d2={D2_NM:.1f} nm", vmax=10)
plt.show()
]]></script>

</plask>
