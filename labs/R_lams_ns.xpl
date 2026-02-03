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
  <cartesian2d name="main" axes="x,y" left="periodic" right="periodic" bottom="mySi">
    <stack>
      <rectangle name="layer1" material="mySi02" dx="1" dy="0.060"/>
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

plt.axhline(0, color='k', lw=0.7)
plt.axvline(lam_min, color='b', lw=0.5)

lams = np.arange(400.0, 700.1, 2.0)

class LayerMaterial(material.Material):
    
    def __init__(self, nr):
        super().__init__()
        self._nr = nr

    def abs(self, lam, T=300.):
        return 0
    
    def nr(self, lam, T=300., n=0.):
        return self._nr
        
def plot_refls(thickness, nr):
    GEO.layer1.height = thickness
    GEO.layer1.material = LayerMaterial(nr)
    refls = OPT.compute_reflectivity(lams, 'top', 'TE')
    # plt.plot(lams, refls, label= f"$h$ = {thickness} um")
    plt.plot(lams, refls, label= f"$nr$ = {nr}")

    
# for h in (0.060, 0.070, 0.080):
#     plot_refls(h)

for nr in (2.1, 2.2, 2.3):
    plot_refls(0.070, nr)

# print_log('result', f"Reflevtivity = {refl:.2f}%")
plt.legend()
plt.xlabel("Walength (nm)")
plt.ylabel("Reflectivity (%)")


plt.show()
]]></script>

</plask>
