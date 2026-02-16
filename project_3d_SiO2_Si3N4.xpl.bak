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
      <rectangle name="layer1" material="SiO2" dx="1" dy="0.1009"/>
      <rectangle name="layer2" material="Si3N4" dx="1" dy="0.0513"/>
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
# this file plots standard mean refl with respect to angle for the SiO2 Si3N4 with the optimized d1, d2


import numpy as np
import matplotlib.pyplot as plt
import os

# optimal params
d_star = [100.9, 51.3]  # nm
n_star = [1.41, 2.64] 


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


R_TE = refl("TE")
R_TM = refl("TM")

def calc_mean(R):
    return np.average(R)

avg_te = calc_mean(R_TE)
avg_tm = calc_mean(R_TM)


print(f"TE Mean: {avg_te:.3f}%")
print(f"TM Mean: {avg_tm:.3f}%")


def plot_map_mean(R, pol, avg):
    plt.figure(figsize=(8, 5))
    cp = plt.pcolormesh(X, Y, R, shading='auto', cmap='magma', vmin=0, vmax=30)
    cb = plt.colorbar(cp)
    cb.set_label('Reflectivity [%]')
    
    plt.title(f"2-Layer AR Coating - {pol} Polarization\n Mean $R$ = {avg:.3f}%")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Angle of Incidence [°]")
    plt.tight_layout()
    plt.savefig(f"plots/2D_mean_sio2_si3n4_{pol}.png", dpi=300)
    plt.show()

plot_map_mean(R_TE, "TE", avg_te)
plot_map_mean(R_TM, "TM", avg_tm)
]]></script>

</plask>
