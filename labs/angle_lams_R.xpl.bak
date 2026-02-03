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
  <module name="si_asp"/>
</materials>

<geometry>
  <cartesian2d name="main" axes="x,y" left="periodic" right="periodic" bottom="Si_Shinke">
    <stack>
      <rectangle material="Si_Shinke" dx="1" dy="1"/>
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


angles = np.arange(0., 90., 2.)        # deg
lams = np.arange(400., 700., 2.)        # nm


def refl(angle, lam, polarization, side='top'):
        OPTICAL.ktran = 2e3 * np.pi / lam * np.sin(np.pi / 180. * angle)
        if side == 'bottom':
            OPTICAL.ktran *= material.get('Si_Shinke').nr(lam)
        return OPTICAL.compute_reflectivity(lam, side, polarization)
            
x, y = np.meshgrid(angles, lams, indexing='ij')
refls_TE = np.zeros((size(angles), size(lams)))
refls_TM = np.zeros((size(angles), size(lams))) 

for i, angle in enumerate(angles):
    for j, lam in enumerate(lams):
        refls_TE[i, j] = refl(angle, lam , 'TE', 'top')
        refls_TM[i, j] = refl(angle, lam , 'TM', 'top')
        
cmap = plt.cm.viridis
cmap = cmap.copy()
cmap.set_over('white')

plt.figure(figsize=(8, 6))
heatmap = plt.pcolormesh(x, y, refls_TE, shading='auto', cmap=cmap, vmin = 0, vmax = 100)
plt.colorbar(heatmap, label='Reflectivity (%)',  extend='max')
plt.xlabel("Angle (deg)")
plt.ylabel("Wavelength (nm)")
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()

]]></script>

</plask>
