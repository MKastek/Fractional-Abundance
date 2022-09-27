# Fractional Abundance

As an example, we will determine the FA for Li, with the following assumptions for the plasma:
- The plasma is rare, so collisional recombination involving two electrons is
unlikely (negligible),
- The plasma is optically thin, i.e., the plasma does not absorb radiation, so ionization
due to absorption of energy hν is negligible,
- We assume a stationary state i.e. that changes occur slowly.

All kinetic equations describing the incomplete thermodynamic equilibrium of ionization and
recombination are as follows:

### Collision ionization
Associated with collisional ionization is the $k_{SCD}$ $\frac{cm^{3}}{s}$ (Effective Ionization Coefficient):

$$Li^{0} + e^{-} \rightarrow Li^{1+} + 2e^{-} \quad (k^{01}_{SCD})$$

$$Li^{1+} + e^{-} \rightarrow Li^{2+} + 2e^{-}  \quad (k^{12}_{SCD})$$

$$Li^{2+} + e^{-} \rightarrow Li^{3+} + 2e^{-}  \quad  (k^{23}_{SCD})$$


### Collision recombination
Associated with collisional recombination is the $k_{ACD}$ $\frac{cm^{3}}{s}$ (Effective Recombination Coefficient):

$$Li^{1+} + e^{-} \rightarrow Li^{0} + h\nu \quad (k^{10}_{ACD})$$

$$Li^{2+} + e^{-} \rightarrow Li^{1+} + h\nu  \quad (k^{21}_{ACD})$$

$$Li^{3+} + e^{-} \rightarrow Li^{2+} + h\nu  \quad  (k^{32}_{ACD})$$


### Rate equations

A set of kinetic equations describing the change in the value of a given ion over time is represented by the following equations, which for coronal equilibrium should be equal to 0.


$$\large  \frac{dn_{0}}{dt} = k_{ACD}^{10} n_1 n_e - k_{SCD}^{01} n_{0} n_{e} = 0 $$

$$\large  \frac{dn_{1}}{dt} = k_{SCD}^{01} n_{0}n_{e} - k_{SCD}^{12} n_{1}n_{e} + k_{ACD}^{21}n_{2}n_{e} - k_{ACD}^{10}n_{1}n_{e} = 0 $$

$$\large  \frac{dn_{2}}{dt} = k_{SCD}^{12} n_{1}n_{e} - k_{SCD}^{23} n_{2}n_{e} + k_{ACD}^{32}n_{3}n_{e} - k_{ACD}^{21}n_{2}n_{e} = 0 $$

$$\large \frac{dn_{3}}{dt} = k_{SCD}^{23} n_{2}n_{e} - k_{ACD}^{32} n_{3}n_{e} = 0$$

Assuming $K_{0}=1$

$$k_{SCD}^{01} n_{0} = k_{ACD}^{10} n_{1}, \quad K_{1} = \frac{k_{SCD}^{01}}{k_{ACD}^{10}}=\frac{n_{1}}{n_{0}}, \quad K_{0}\cdot K_{1}=\frac{n_{1}}{n_{0}}$$

$$k_{SCD}^{12} n_{1} = k_{ACD}^{21} n_{2}, \quad K_{2} = \frac{k_{SCD}^{12}}{k_{ACD}^{21}}=\frac{n_{2}}{n_{1}}, \quad K_{1}\cdot K_{2}=\frac{n_{2}}{n_{0}}$$

$$k_{SCD}^{23} n_{2} = k_{ACD}^{32} n_{3}, \quad K_{3} = \frac{k_{SCD}^{23}}{k_{ACD}^{32}}=\frac{n_{3}}{n_{2}}, \quad K_{1}\cdot K_{2}\cdot K_{3}=\frac{n_{3}}{n_{0}}$$

Each of the coefficients $k_{SCD}$, $k_{ACD}$ and thus $K$ depends on the temperature $T_{e}$ and density $N_{e}$.

### Fractional Abundance

$$\large FA(Li^{0})=\frac{n_{0}}{n_{0}+n_{1}+n_{2}+n_{3}}=\frac{K_{0}}{K_{0}+K_{0}\cdot K_{1}+K_{0}\cdot K_{1}\cdot K_{2}+K_{0} \cdot K_{1}\cdot K_{2} \cdot K_{3}}$$

$$\large FA(Li^{1+})=\frac{n_{1}}{n_{0}+n_{1}+n_{2}+n_{3}}=\frac{K_{0}\cdot K_{1}}{K_{0}+K_{0}\cdot K_{1}+K_{0}\cdot K_{1}\cdot K_{2}+K_{0} \cdot K_{1}\cdot K_{2} \cdot K_{3}}$$

$$\large FA(Li^{2+})=\frac{n_{2}}{n_{0}+n_{1}+n_{2}+n_{3}}=\frac{K_{0}\cdot K_{1} \cdot K_{2}}{K_{0}+K_{0}\cdot K_{1}+K_{0}\cdot K_{1}\cdot K_{2}+K_{0} \cdot K_{1}\cdot K_{2} \cdot K_{3}}$$

$$\large FA(Li^{3+})=\frac{n_{3}}{n_{0}+n_{1}+n_{2}+n_{3}}=\frac{K_{0}\cdot K_{1} \cdot K_{2} \cdot K_{3}}{K_{0}+K_{0}\cdot K_{1}+K_{0}\cdot K_{1}\cdot K_{2}+K_{0} \cdot K_{1}\cdot K_{2} \cdot K_{3}}$$

The above formulas can be written in abbreviated notation:

$$\Huge FA(Li^{ i+})=\prod_{j=0} ^{i} K_{j}\left(\sum_{k=0} ^{Z} \prod_{j=0} ^{k} K_{j} \right)^{-1} $$



### Implementation
Effective Recombination Coefficient and Effective Ionization Coefficient can be found in [OPEN-ADAS](https://open.adas.ac.uk/adf11?element=&acd=1&scd=1&year=&searching=1#searchbutton) database. Fractional Abundace is calculated with above equation in temperature [eV] and density [1/m^3] range available
in OPEN-ADAS files.

Calculation are performed with use of Numba in multithreaded approach, FA for every ion is calculated in another thread with use of ThreadPoolExecutor.

### GUI - tkinter
Simple GUI application allows to study FA for different atoms.

![](images/plot.PNG)
