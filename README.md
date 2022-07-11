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

$Li^{0} + e^{-} \rightarrow Li^{1+} + 2e^{-} \quad (k^{01}_{SCD})$  

$Li^{1+} + e^{-} \rightarrow Li^{2+} + 2e^{-}  \quad (k^{12}_{SCD})$

$Li^{2+} + e^{-} \rightarrow Li^{3+} + 2e^{-}  \quad  (k^{23}_{SCD})$


### Collision recombination  

$Li^{1+} + e^{-} \rightarrow Li^{0} + h\nu \quad (k^{10}_{ACD})$  

$Li^{2+} + e^{-} \rightarrow Li^{1+} + h\nu  \quad (k^{21}_{ACD})$  

$Li^{3+} + e^{-} \rightarrow Li^{2+} + h\nu  \quad  (k^{32}_{ACD})$


### Rate equations


$\frac{dn_{0}}{dt} = k^{10}_{ACD} n_{1} n_{e} - k^{01}_{SCD} n_{0} n_{e} = 0$

$\frac{dn_{1}}{dt} = k^{01}_{SCD}n_{0}n_{e} - k^{12}_{SCD}n_{1}n_{e} + k^{21}_{ACD}n_{2}n_{e} - k^{10}_{ACD}n_{1}n_{e} = 0$

$\frac{dn_{2}}{dt} = k^{12}_{SCD}n_{1}n_{e} - k^{23}_{SCD}n_{2}n_{e} + k^{32}_{ACD}n_{3}n_{e} - k^{21}_{ACD}n_{2}n_{e} = 0$

$\frac{dn_{3}}{dt} = k^{23}_{SCD}n_{2}n_{e} - k^{32}_{ACD}n_{3}n_{e} = 0 $

### GUI - tkinter
![](images/plot.PNG)