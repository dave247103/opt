<plask loglevel="detail">

<geometry>
  <cartesian2d name="simple" axes="x,y " left="periodic" right="periodic" bottom="Si">
    <rectangle material="Si" dx="1" dy="1"/>
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

# PLaSK result (percent; compute_reflectivity returns [%])
refl_TE = OPTICAL.compute_reflectivity(500.0, 'top', 'TE')
refl_TM = OPTICAL.compute_reflectivity(500.0, 'top', 'TM')

# Analytic Fresnel
air = plask.material.Air()
si  = plask.material.Si()

n0 = air.Nr(500.0)   # complex n+ik
ns = si.Nr(500.0)

R_fresnel = abs((ns - n0) / (ns + n0))**2 * 100.0

print_log('result', f"Air Nr(500)={n0}, Si Nr(500)={ns}")
print_log('result', f"Fresnel R(500nm)={R_fresnel:.3f}%")
print_log('result', f"PLaSK R_TE={refl_TE:.3f}%, R_TM={refl_TM:.3f}%")
]]></script>

</plask>
