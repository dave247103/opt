<plask loglevel="detail">

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
      <rectangle material="Si_Shinke" dx="1" dy="1"/>
    </stack>
  </cartesian2d>
</geometry>

<solvers>
  <optical name="OPT" solver="Fourier2D" lib="modal">
    <geometry ref="main"/>
    <expansion size="0"/>
  </optical>
</solvers>

<script><![CDATA[
import numpy as np
import matplotlib.pyplot as plt

#lams = np.arange(400.0, 700.1, 1.0)

#lam = 550. # nm

angles = np.arange(0., 90., 0.1)        # deg


#refls = OPT.compute_reflectivity(lams, 'top', 'TE')

def refl(angle, lam, polarization, side='top'):
        # OPT.ktran = 2e3 * np.pi / lam * np.sin(np.pi / 180. * angle) * material.get('mySi').nr(lam)
        OPT.ktran = 2e3 * np.pi / lam * np.sin(np.pi / 180. * angle)
        if side == 'bottom': OPT.ktran *= material.get('mySi').nr(lam)
        return OPT.compute_reflectivity(lam, 'bottom', polarization)
        #return OPT.compute_reflectivity(lam, 'top', polarization)

refls_TE = [refl(angle, 550., 'TE') for angle in angles]
refls_TM = [refl(angle, 550., 'TM') for angle in angles]

plt.axhline(0, color='k', lw=0.7)
plt.plot(angles, refls_TE, label='TE')
plt.plot(angles, refls_TM, label='TM')
plt.xlabel("Angle (deg)")
plt.ylabel("Reflectivity (%)")
plt.legend()

plt.show()




# print_log('result', f"Reflevtivity = {refl:.2f}%")
# plt.axhline(0, color='k', lw=0.7)
# plt.axvline(550., color='b', lw=0.5)
# plt.plot(lams, refls)
# plt.xlabel("Walength (nm)")
# plt.ylabel("Reflectivity (%)")
# plt.show()
]]></script>

</plask>
