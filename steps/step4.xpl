<plask loglevel="detail">

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
</materials>

<geometry>
  <cartesian2d name="simple" axes="x,y " left="periodic" right="periodic" bottom="mySi">
    <rectangle material="mySiO2" dx="1" dy="{lam_min * 1e-3 / (4 * n_layer)}"/>
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

# PLaSK result (percent; compute_reflectivity returns [%])
refl_TE = OPTICAL.compute_reflectivity(550.0, 'top', 'TE')
refl_TM = OPTICAL.compute_reflectivity(550.0, 'top', 'TM')

# Analytic Fresnel
air = plask.material.Air()
mysi  = plask.material.mySi()
mysio2 = plask.material.mySiO2()

n0 = air.Nr(550.0)   # complex n+ik
nsi = mysi.Nr(550.0)
nsio2 = mysio2.Nr(550)

R_fresnel = abs((nsi - n0) / (nsi + n0))**2 * 100.0
print_log('result', f"n_SiO2(550)={nsio2}")
print_log('result', f"Air Nr(550)={n0}, my_Si Nr(550)={nsi}")
print_log('result', f"Fresnel R(550nm)={R_fresnel:.3f}%")
print_log('result', f"PLaSK R_TE={refl_TE:.33f}%, R_TM={refl_TM:.33f}%")

lams = np.arange(400.0, 700.1, 1.0)
refls_TE = OPTICAL.compute_reflectivity(lams, 'top', 'TE')
refls_TM = OPTICAL.compute_reflectivity(lams, 'top', 'TM')



plt.figure()
plt.plot(lams, refls_TE, label="TE")
plt.plot(lams, refls_TM, label="TM", linestyle="--")
plt.axhline(0, color='k', lw=0.7)
plt.axvline(550., color='b', lw=0.5)
plt.xlabel("Walength (nm)")
plt.ylabel("Reflectivity (%)")
plt.legend()
plt.show()
]]></script>

</plask>
