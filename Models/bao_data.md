## Data
| $z_\mathrm{eff}$ | $d_z$   | error  | Survey                 |
|------------------|---------|--------|------------------------|
| 0.106            | 0.3360  | 0.0150 | [6dFGS (2011)][1]      |
| 0.15             | 0.2239  | 0.0084 | [MGS (2015)][2]        |
| 0.32             | 0.1181  | 0.0024 | [BOSS LOWZ (2014)][3]  |
| 0.35             | 0.1126  | 0.0022 | [SDSS(R) (2012)][4]    |
| 0.57             | 0.07261 | 0.00071| [BOSS CMASS (2014)][3] |

Note: LOWZ = low-$z$ galaxies, CMASS = constant stellar mass  galaxies

The WiggleZ dataset comprises 3 measurements of the distilled parameter $d_z(z)$ at three
effective redshifts.

|Redshift slice | $z_\mathrm{eff}$ | $d_z$|
|---------------|------------------|------|
|0.2<$z$<0.6    |             0.44 |0.073 |
|0.4<$z$<0.8    |             0.6  |0.0726|
|0.6<$z$<1.0    |             0.73 |0.0592|

As the data points are not independent measurements (the three redshifts slices are overlapping)
there is some covariance between points.
WiggleZ inverse covariance matrix is
$$
(C_\mathrm{WiggleZ})^{-1}
=\left(
\begin{array}{ccc}
 1040.3 & -807.5  &  336.8 \\
 -807.5  &  3720.3 & -1551.9 \\
  336.8  & -1551.9 &  2914.9
  \end{array}\right)
$$
See table 2 in C. Blake et al., Mon. Not. R. Astron. Soc. __418__, 1707 (2011).

[1]: https://arxiv.org/abs/1106.3366 "6dF Galaxy Survey"
[2]: https://arxiv.org/abs/1409.3242 "Main Galaxy Sample"
[3]: https://arxiv.org/abs/1312.4877 "BOSS LOWZ and CMASS (DR10+DR11)"
[4]: https://arxiv.org/abs/1202.0090 "SDSS: reconstruction method"
[5]: https://arxiv.org/abs/1501.00963 "BOSS LOWZ and CMASS (DR12)"
