CDPM2 user material for LS-DYNA
===============================

This is a description of the UMAT implementation of CDPM2 for LS-DYNA
provided by Peter Grassl on https://github.com/githubgrasp/cdpm2.

CDPM2 is a damage-plastic concrete model based on work published in Grassl et
al. (2011, 2013). It is an extension of CDPM proposed in Grassl and Jirásek
(2006). The model has been developed to simulate the failure process of
concrete subjected to multi-axial and rate-dependent loading. It is based on
effective stress plasticity combined with a damage model which is based on
both plastic and elastic strain measures. This material model is available
only for solids.

The LS-DYNA user material was originally implemented by Dimitrios Xenos and
Peter Grassl in 2014. Since then, updates and corrections of the model have
been made, and the implementation is kept consistent with the CDPM2
implementation in OOFEM (``concretedpm2.C``). Because of these updates and
corrections, the user material is not exactly the same as the version
available in LS-DYNA (``*MAT_CDPM``).

Input cards
-----------

The material is defined with ``*MAT_USER_DEFINED_MATERIAL_MODELS`` using six
cards. The first two are related to the use of user materials in LS-DYNA. The
last four contain the CDPM2 specific parameters. An example is given in
``materialUser.k``.

Note that LS-DYNA reads a blank field as 0.0. The defaults below are the
values used for a blank (zero) field.

======== ==== ==== == === === ====== ===== ==
Card 1   1    2    3  4   5   6      7     8
======== ==== ==== == === === ====== ===== ==
Variable MID  RO   MT LMC NHV IORTHO IBULK IG
Type     A8   F    F  F   F   F      F     F
Input    none none 50 26  28  0      25    26
======== ==== ==== == === === ====== ===== ==

======== ===== ===== ====== ====== ==== ==== = =
Card 2   1     2     3      4      5    6    7 8
======== ===== ===== ====== ====== ==== ==== = =
Variable IVECT IFAIL ITHERM IHYPER IEOS LMCA
Type     F     F     F      F      F    F
Input    1     0     0      0      0    0
======== ===== ===== ====== ====== ==== ==== = =

======== ==== ==== ===== === ==== ==== === ====
Card 3   1    2    3     4   5    6    7   8
======== ==== ==== ===== === ==== ==== === ====
Variable E    PR   ECC   QH0 FT   FC   HP  AH
Type     F    F    F     F   F    F    F   F
Default  none 0.0  0.525 0.3 none none 0.0 0.08
======== ==== ==== ===== === ==== ==== === ====

======== ===== == ====== == ==== === ==== ==
Card 4   1     2  3      4  5    6   7    8
======== ===== == ====== == ==== === ==== ==
Variable BH    CH DH     AS DF   ERT TYPE BS
Type     F     F  F      F  F    F   F    F
Default  0.003 2  1.0E-6 15 0.85 0   0    1
======== ===== == ====== == ==== === ==== ==

======== ==== ======= ====== ====== ======= ====== ====== ====
Card 5   1    2       3      4      5       6      7      8
======== ==== ======= ====== ====== ======= ====== ====== ====
Variable WF   WF1     FT1    STRFLG FAILFLG EFC    DAMFLG PFLG
Type     F    F       F      F      F       F      F      F
Default  none 0.15*WF 0.3*FT 0      0       1.0E-4 0      0
======== ==== ======= ====== ====== ======= ====== ====== ====

======== ==== ==== = = = = = =
Card 6   1    2    3 4 5 6 7 8
======== ==== ==== = = = = = =
Variable K    G
Type     F    F
Default  none none
======== ==== ==== = = = = = =

Parameters
----------

Cards 1 and 2
~~~~~~~~~~~~~

MID
   Material identification. A unique number or label not exceeding 8
   characters must be specified.

RO
   Mass density.

MT
   Type of user material routine used. CDPM2 is implemented in ``umat50v``,
   so MT = 50.

LMC
   Number of model input parameters. There are 26 including K and G.

NHV
   Number of history variables. The model has 28 history variables, which are
   listed below.

IORTHO
   Flag which allows for a local coordinate system. Not needed here, so set
   to zero.

IBULK
   Input parameter number where K is specified. Here it is 25.

IG
   Input parameter number where G is specified. Here it is 26.

IVECT
   Flag specifying if the vectorised format is used. It is used here, so set
   to 1.

IFAIL
   Flag which allows for element deletion. It must be set to a nonzero value
   if integration points flagged as failed (see FAILFLG) should lead to
   element deletion. Otherwise, set to zero.

ITHERM
   Flag which allows for temperatures in the material model. Not used here,
   so set to zero.

IHYPER
   Flag which allows for the use of the deformation gradient. Not used here,
   so set to zero.

IEOS
   Flag which allows for the use of an equation of state. Not used here, so
   set to zero.

LMCA
   Additional memory for material data. Not used here, so set to zero.

Cards 3 to 6
~~~~~~~~~~~~

E
   Young's modulus. If the value of the Young's modulus is not available, it
   can be estimated from the uniaxial compressive strength according to, for
   instance, the fib Model Code 2010.

PR
   Poisson's ratio.

ECC
   Eccentricity parameter. If not provided (or set to zero), the default
   value 0.525 is used.

QH0
   Onset of hardening. Default is 0.3.

FT
   Uniaxial tensile strength. If the value of the uniaxial tensile strength
   is not available, it can be estimated from the uniaxial compressive
   strength according to, for instance, the fib Model Code 2010.

FC
   Uniaxial compressive strength.

HP
   Hardening parameter. Default is 0.0. HP must lie between 0 and 1-QH0.
   Note the following:

   1. For applications without strain rate effects (STRFLG = 0), HP = 0.01
      is recommended.
   2. For strain rate dependent material response (STRFLG > 0), HP = 0.5 is
      recommended.

AH, BH, CH, DH
   Hardening ductility parameters 1 to 4. Defaults are 0.08, 0.003, 2 and
   1.0E-6.

AS
   Ductility parameter during damage. Default is 15.

DF
   Flow rule parameter. Default is 0.85.

ERT
   Flag for the way the fracture energy is related to the dynamic increase
   factor (DIF). Only used if STRFLG > 0.

   - EQ.0.0: fracture energy is not increased.
   - EQ.1.0: fracture energy is increased by DIF.
   - EQ.2.0: fracture energy is increased by the square of DIF.

TYPE
   Flag for the tensile damage type.

   - EQ.0.0: linear damage formulation (default).
   - EQ.1.0: bilinear damage formulation.
   - EQ.2.0: exponential damage formulation.

   The best results are obtained with the bilinear formulation.

BS
   Damage ductility exponent during damage. Default is 1.

WF
   Parameter controlling the tensile softening branch. For TYPE = 0 (linear
   softening), WF = 2 GF/FT, where GF is the total fracture energy. For
   TYPE = 1 (bilinear softening) with the default values of WF1 and FT1,
   WF = 4.444 GF/FT. For TYPE = 2 (exponential softening), WF = GF/FT. For
   tetrahedral elements, the value of WF should be multiplied by 0.56.

   If the value of the fracture energy GF is not available, it can be
   estimated from the uniaxial compressive strength according to, for
   instance, the fib Model Code 2010.

   In some cases, softening down to zero stress might result in numerical
   problems, which are not directly related to the constitutive model. To
   overcome these problems, it might be necessary to use an exponential
   softening law.

WF1
   Crack opening threshold for the second part of the bilinear damage
   formulation. Default is 0.15 WF.

FT1
   Stress threshold for the bilinear damage formulation. Default is 0.3 FT.

STRFLG
   Strain rate flag for strength.

   - EQ.0.0: no strain rate dependence (default).
   - EQ.1.0: strain rate dependence using the first branch of the MC2010 DIF.
   - EQ.2.0: strain rate dependence using the first and second branch of the
     MC2010 DIF.

FAILFLG
   Failure flag.

   - EQ.0.0: not active (default).
   - GT.0.0: active. An integration point is flagged as failed if wt and wc
     both exceed 0.9995. Whether and when the element is eroded is controlled
     by LS-DYNA (see IFAIL).

   Independently of FAILFLG, an integration point is flagged as failed if the
   plasticity return does not converge, and its stress is set to zero.

EFC
   Parameter controlling the compressive softening branch in the exponential
   compressive damage formulation. Default is 1.0E-4.

DAMFLG
   Flag to determine the type of damage formulation that is used.

   - EQ.0.0: no damage.
   - EQ.1.0: two damage variables with a split of the effective stress into a
     tensile and a compressive part (recommended).
   - EQ.2.0: one damage variable which still takes into account tensile and
     compressive damage.
   - EQ.3.0: one damage variable (tensile damage).

   Note that a blank field is read as 0.0, i.e. no damage, so DAMFLG has to
   be set explicitly.

PFLG
   Print flag. If zero, the input parameters are written once to the output.

K
   Bulk modulus, calculated as E/(3(1-2 PR)).

G
   Shear modulus, calculated as E/(2(1+PR)).

History variables
-----------------

== =================================================
1  kappa
2  equivalent strain - tension
3  plastic strain xx
4  plastic strain yy
5  plastic strain zz
6  plastic strain xy
7  plastic strain yz
8  plastic strain xz
9  kappa tension kdt
10 kappa tension kdt1
11 kappa tension kdt2
12 kappa compression kdc
13 kappa compression kdc1
14 kappa compression kdc2
15 damage variable tension omegaT
16 damage variable compression omegaC
17 strain rate factor
18 compression factor alphac
19 equivalent strain - tension
20 equivalent strain - compression
21 total strain xx
22 total strain yy
23 total strain zz
24 total strain xy
25 total strain yz
26 total strain xz
27 equivalent strain (without rate factor)
28 principal strain used for the strain rate
== =================================================

References
----------

P. Grassl and M. Jirásek. "Damage-plastic model for concrete failure".
International Journal of Solids and Structures, Vol. 43, pp. 7166-7196, 2006.

P. Grassl, U. Nyström, R. Rempling and K. Gylltoft. "A damage-plasticity
model for the dynamic failure of concrete". 8th International Conference on
Structural Dynamics, Leuven, Belgium, 2011.

P. Grassl, D. Xenos, U. Nyström, R. Rempling and K. Gylltoft. "CDPM2: A
damage-plasticity approach to modelling the failure of concrete".
International Journal of Solids and Structures, Vol. 50, Issue 24,
pp. 3805-3816, 2013.
