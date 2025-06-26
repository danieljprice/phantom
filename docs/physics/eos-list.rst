.. table:: Equations of state implemented in phantom
   :widths: auto

   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | ieos      | Description                                                                                                              | 
   +===========+==========================================================================================================================+
   | 1         | **Isothermal eos**                                                                                                       |
   |           |                                                                                                                          |
   |           | :math:`P = c_s^2 \rho`                                                                                                   |
   |           |                                                                                                                          |
   |           | where :math:`c_s^2 \equiv K` is a constant stored in the dump file header                                                |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 2         | **Adiabatic equation of state (code default)**                                                                           |
   |           |                                                                                                                          |
   |           | :math:`P = (\gamma - 1) \rho u`                                                                                          |
   |           |                                                                                                                          |
   |           | if the code is compiled with ISOTHERMAL=yes, ieos=2 gives a polytropic eos:                                              |
   |           |                                                                                                                          |
   |           | :math:`P = K \rho^\gamma`                                                                                                |
   |           |                                                                                                                          |
   |           | where K is a global constant specified in the dump header                                                                |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 3         | **Locally isothermal disc as in Lodato & Pringle (2007) where**                                                          |
   |           |                                                                                                                          |
   |           | :math:`P = c_s^2 (r) \rho`                                                                                               |
   |           |                                                                                                                          |
   |           | sound speed (temperature) is prescribed as a function of radius using:                                                   |
   |           |                                                                                                                          |
   |           | :math:`c_s = c_{s,0} r^{-q}` where :math:`r = \sqrt{x^2 + y^2 + z^2}`                                                    |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 4         | **Isothermal equation of state for GR, enforcing cs = constant**                                                         |
   |           |                                                                                                                          |
   |           | .. WARNING:: this is experimental: use with caution                                                                      |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 6         | **Locally isothermal disc centred on sink particle**                                                                     |
   |           |                                                                                                                          |
   |           | As in ieos=3 but in this version radius is taken with respect to a designated                                            |
   |           | sink particle (by default the first sink particle in the simulation)                                                     |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 7         | **Vertically stratified equation of state**                                                                              |
   |           |                                                                                                                          |
   |           | sound speed is prescribed as a function of (cylindrical) radius R and                                                    |
   |           | height z above the x-y plane                                                                                             |
   |           |                                                                                                                          |
   |           | .. WARNING:: should not be used for misaligned discs                                                                     |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 8         | **Barotropic equation of state**                                                                                         |
   |           |                                                                                                                          |
   |           | :math:`P = K \rho^\gamma`                                                                                                |
   |           |                                                                                                                          |
   |           | where the value of gamma (and K) are prescribed functions of density                                                     |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 9         | **Piecewise Polytropic equation of state**                                                                               |
   |           |                                                                                                                          |
   |           | :math:`P = K \rho^\gamma`                                                                                                |
   |           |                                                                                                                          |
   |           | where the value of gamma (and K) are a prescribed function of density.                                                   |
   |           | Similar to ieos=8 but with different defaults and slightly different                                                     |
   |           | functional form                                                                                                          |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 10        | **MESA equation of state**                                                                                               |
   |           |                                                                                                                          |
   |           | a tabulated equation of state including gas, radiation pressure                                                          |
   |           | and ionisation/dissociation. MESA is a stellar evolution code, so                                                        |
   |           | this equation of state is designed for matter inside stars                                                               |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 11        | **Equation of state with pressure and temperature equal to zero**                                                        |
   |           |                                                                                                                          |
   |           | :math:`P = 0`                                                                                                            |
   |           |                                                                                                                          |
   |           | useful for simulating test particle dynamics using SPH particles                                                         |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 12        | **Ideal gas plus radiation pressure**                                                                                    |
   |           |                                                                                                                          |
   |           | :math:`P = (\gamma - 1) \rho u`                                                                                          |
   |           |                                                                                                                          |
   |           | but solved by first solving the quartic equation:                                                                        |
   |           |                                                                                                                          |
   |           | :math:`u = \frac32 \frac{k_b T}{\mu m_H} + \frac{a T^4}{\rho}`                                                           |
   |           |                                                                                                                          |
   |           | for temperature (given u), then solving for pressure using                                                               |
   |           |                                                                                                                          |
   |           | :math:`P = \frac{k_b T}{\mu m_H} + \frac13 a T^4`                                                                        |
   |           |                                                                                                                          |
   |           | hence in this equation of state gamma (and temperature) are an output                                                    |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 13        | **Locally isothermal eos for generic hierarchical system**                                                               |
   |           |                                                                                                                          |
   |           | Assuming all sink particles are stars.                                                                                   |
   |           | Generalisation of Farris et al. (2014; for binaries) to N stars.                                                         |
   |           | For two sink particles this is identical to ieos=14                                                                      |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 14        | **Locally isothermal eos from Farris et al. (2014) for binary system**                                                   |
   |           |                                                                                                                          |
   |           | uses the locations of the first two sink particles                                                                       |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 15        | **Helmholtz equation of state for fully ionized gas**                                                                    |
   |           |                                                                                                                          |
   |           | .. WARNING:: not widely tested in phantom, better to use ieos=10                                                         |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 16        | **Shen (2012) equation of state for neutron stars**                                                                      |
   |           |                                                                                                                          |
   |           | this equation of state requires evolving temperature as the energy variable                                              |
   |           |                                                                                                                          |
   |           | .. WARNING:: not tested: use with caution                                                                                |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 20        | **Gas + radiation + various forms of recombination**                                                                     |
   |           |                                                                                                                          |
   |           | from HORMONE, Hirai+2020, as used in Lau+2022b                                                                           |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 21        | **HII region two temperature "equation of state"**                                                                       |
   |           |                                                                                                                          |
   |           | flips the temperature depending on whether a particle is ionised or not,                                                 |
   |           | use with ISOTHERMAL=yes                                                                                                  |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 22        | **Same as ieos=21 but sets the thermal energy**                                                                          |
   |           |                                                                                                                          |
   |           | for use when u is stored (ISOTHERMAL=no)                                                                                 |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 23        | **Tillotson (1962) equation of state for solids (basalt, granite, ice, etc.)**                                           |
   |           |                                                                                                                          |
   |           | Implementation from Benz et al. (1986), Asphaug & Melosh (1993) and Kegerreis et al. (2019)                              |
   |           |                                                                                                                          |
   |           | In the compressed (:math:`\rho > \rho_0`) or cold (:math:`u < u_{\rm iv}`) state gives                                   |
   |           |                                                                                                                          |
   |           | :math:`P_c = \left[a + \frac{b}{(u/(u_0 \eta^2) + 1}\right]\rho u + A \mu + B\mu^2`                                      |
   |           |                                                                                                                          |
   |           | where :math:`\eta = rho/rho_0`, :math:`\mu = \eta - 1`, u is the specific internal energy                                |
   |           | and a,b,A,B and :math:`u_0` are input parameters chosen for a particular material                                        |
   |           |                                                                                                                          |
   |           | In the hot, expanded state (:math:`\rho < \rho_0` and :math:`u > u_{\rm iv}`) gives                                      |
   |           |                                                                                                                          |
   |           | :math:`P_e = a\rho u + \left[\frac{b\rho u}{u/(u_0 \eta^2) + 1} + A\mu \exp{-\beta \nu} \right] \exp(-\alpha \nu^2)`     |
   |           |                                                                                                                          |
   |           | where :math:`\nu = \rho/\rho_0 - 1`. In the intermediate state pressure is interpolated using                            |
   |           |                                                                                                                          |
   |           | :math:`P = \frac{(u - u_{\rm iv}) P_e + (u_{\rm cv} - u) P_c}{u_{\rm cv} - u_{\rm iv}}.`                                 |
   |           |                                                                                                                          |
   |           | When using this equation of state bodies should be set up with uniform density equal, or close to the                    |
   |           | reference density :math:`\rho_0`, e.g. 2.7 g/cm^3 for basalt                                                             |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
   | 24        | **Tabulated EoS of Stamatellos et al. 2007 (includes opacities)**                                                        |
   |           |                                                                                                                          |
   |           | Tabulated equation of state with opacities from Lombardi et al. 2015. For use                                            |
   |           | with icooling = 9, the radiative cooling approximation (Young et al. 2024).                                              |
   |           |                                                                                                                          |
   +-----------+--------------------------------------------------------------------------------------------------------------------------+
