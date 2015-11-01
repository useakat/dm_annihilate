# SOFTSUSY3.4.1 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.1       # version number
Block MODSEL  # Select model
     1    0   # nonUniversal
Block SMINPUTS             # Standard Model inputs
     1    1.27930117e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637900e-05   # G_Fermi
     3    1.18400000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.18000000e+00   # mb(mb)
     6    1.73100000e+02   # Mtop(pole)
     7    1.77682000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    5.00000000e+00   # tanb, DRbar, Feynman gauge
Block EXTPAR               # non-universal SUSY breaking parameters
     0    -1.00000000e+00  # Set MX=MSUSY
     1     3.00000000e+02  # M_1(MX)
     2     5.00000000e+02  # M_2(MX)
     3     1.50000000e+03  # M_3(MX)
     11    1.00000000e+03  # At(MX)
     12    1.00000000e+02  # Ab(MX)
     13    1.00000000e+02  # Atau(MX)
     23    2.00000000e+05  # mu(MX)
     26    2.00000000e+05  # mA(pole)
     31    2.00000000e+05  # meL(MX)
     32    2.00000000e+05  # mmuL(MX)
     33    2.00000000e+05  # mtauL(MX)
     34    2.00000000e+05  # meR(MX)
     35    2.00000000e+05  # mmuR(MX)
     36    2.00000000e+05  # mtauR(MX)
     41    2.00000000e+05  # mqL1(MX)
     42    2.00000000e+05  # mqL2(MX)
     43    2.00000000e+05  # mqL3(MX)
     44    2.00000000e+05  # muR(MX)
     45    2.00000000e+05  # mcR(MX)
     46    2.00000000e+05  # mtR(MX)
     47    2.00000000e+05  # mdR(MX)
     48    2.00000000e+05  # msR(MX)
     49    2.00000000e+05  # mbR(MX)
# SOFTSUSY-specific non SLHA information:
# MIXING=1 Desired accuracy=1.00000000e-01 Achieved accuracy=4.99958767e-02
# MX=2.00000030e+05 GeV
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.09446060e+01   # MW
        25     1.24720735e+02   # h0
        35     1.99978512e+05   # H0
        36     1.99977800e+05   # A0
        37     1.99977761e+05   # H+
   1000021     2.27625840e+03   # ~g
   1000022     2.27452403e+02   # ~neutralino(1)
   1000023     3.78262939e+02   # ~neutralino(2)
   1000024     3.78437151e+02   # ~chargino(1)
   1000025     2.03003308e+05   # ~neutralino(3)
   1000035    -2.03003326e+05   # ~neutralino(4)
   1000037     2.03003702e+05   # ~chargino(2)
   1000001     2.01797466e+05   # ~d_L
   1000002     2.01797574e+05   # ~u_L
   1000003     2.01797466e+05   # ~s_L
   1000004     2.01797574e+05   # ~c_L
   1000005     2.01408117e+05   # ~b_1
   1000006     2.01221249e+05   # ~t_1
   1000011     2.00412274e+05   # ~e_L
   1000012     2.00411884e+05   # ~nue_L
   1000013     2.00412274e+05   # ~mu_L
   1000014     2.00411884e+05   # ~numu_L
   1000015     2.00241125e+05   # ~stau_1
   1000016     2.00407724e+05   # ~nu_tau_L
   2000001     2.01419252e+05   # ~d_R
   2000002     2.01398507e+05   # ~u_R
   2000003     2.01419252e+05   # ~s_R
   2000004     2.01398507e+05   # ~c_R
   2000005     2.01707182e+05   # ~b_2
   2000006     2.01709832e+05   # ~t_2
   2000011     2.00249679e+05   # ~e_R
   2000013     2.00249679e+05   # ~mu_R
   2000015     2.00408509e+05   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.14505605e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.99999974e-01   # N_{1,1}
  1  2    -4.85179167e-05   # N_{1,2}
  1  3     2.17214246e-04   # N_{1,3}
  1  4    -4.83667857e-05   # N_{1,4}
  2  1     4.86049506e-05   # N_{2,1}
  2  2     9.99999922e-01   # N_{2,2}
  2  3    -3.81898003e-04   # N_{2,3}
  2  4     8.44105973e-05   # N_{2,4}
  3  1    -1.87778137e-04   # N_{3,1}
  3  2     3.29739095e-04   # N_{3,2}
  3  3     7.07106706e-01   # N_{3,3}
  3  4    -7.07106754e-01   # N_{3,4}
  4  1    -1.19382966e-04   # N_{4,1}
  4  2     2.10361155e-04   # N_{4,2}
  4  3     7.07106720e-01   # N_{4,3}
  4  4     7.07106801e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.99999854e-01   # U_{1,1}
  1  2    -5.39643452e-04   # U_{1,2}
  2  1     5.39643452e-04   # U_{2,1}
  2  2     9.99999854e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.99999993e-01   # V_{1,1}
  1  2    -1.19339851e-04   # V_{1,2}
  2  1     1.19339851e-04   # V_{2,1}
  2  2     9.99999993e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.03868149e-02   # F_{11}
  1  2     9.99792167e-01   # F_{12}
  2  1     9.99792167e-01   # F_{21}
  2  2    -2.03868149e-02   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     1.11849922e-02   # F_{11}
  1  2     9.99937446e-01   # F_{12}
  2  1     9.99937446e-01   # F_{21}
  2  2    -1.11849922e-02   # F_{22}
Block staumix               # stau mixing matrix
  1  1     2.42045994e-02   # F_{11}
  1  2     9.99707026e-01   # F_{12}
  2  1     9.99707026e-01   # F_{21}
  2  2    -2.42045994e-02   # F_{22}
Block gauge Q= 2.00000030e+05  # SM gauge couplings
     1     3.68332995e-01   # g'(Q)MSSM DRbar
     2     6.56682976e-01   # g(Q)MSSM DRbar
     3     9.03675451e-01   # g3(Q)MSSM DRbar
Block yu Q= 2.00000030e+05  
  3  3     7.06162700e-01   # Yt(Q)MSSM DRbar
Block yd Q= 2.00000030e+05  
  3  3     5.23360630e-02   # Yb(Q)MSSM DRbar
Block ye Q= 2.00000030e+05  
  3  3     4.78773414e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 2.00000030e+05 # Higgs mixing parameters
     1     2.02431061e+05    # mu(Q)MSSM DRbar
     2     4.59016350e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43738820e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     3.95464641e+10    # mA^2(Q)MSSM DRbar
Block msoft Q= 2.00000030e+05  # MSSM DRbar SUSY breaking parameters
     1     3.00000000e+02      # M_1(Q)
     2     5.00000001e+02      # M_2(Q)
     3     1.50000001e+03      # M_3(Q)
    21    -3.41346248e+09      # mH1^2(Q)
    22    -3.87608555e+10      # mH2^2(Q)
    31     2.00000000e+05      # meL(Q)
    32     2.00000000e+05      # mmuL(Q)
    33     2.00000000e+05      # mtauL(Q)
    34     2.00000000e+05      # meR(Q)
    35     2.00000000e+05      # mmuR(Q)
    36     2.00000000e+05      # mtauR(Q)
    41     2.00000000e+05      # mqL1(Q)
    42     2.00000000e+05      # mqL2(Q)
    43     2.00000000e+05      # mqL3(Q)
    44     2.00000000e+05      # muR(Q)
    45     2.00000000e+05      # mcR(Q)
    46     1.99999999e+05      # mtR(Q)
    47     2.00000000e+05      # mdR(Q)
    48     2.00000000e+05      # msR(Q)
    49     2.00000000e+05      # mbR(Q)
Block au Q= 2.00000030e+05  
  1  1     1.00000025e+03      # Au(Q)MSSM DRbar
  2  2     1.00000025e+03      # Ac(Q)MSSM DRbar
  3  3     1.00000039e+03      # At(Q)MSSM DRbar
Block ad Q= 2.00000030e+05  
  1  1     1.00000101e+02      # Ad(Q)MSSM DRbar
  2  2     1.00000101e+02      # As(Q)MSSM DRbar
  3  3     1.00000154e+02      # Ab(Q)MSSM DRbar
Block ae Q= 2.00000030e+05  
  1  1     1.00000001e+02      # Ae(Q)MSSM DRbar
  2  2     1.00000001e+02      # Amu(Q)MSSM DRbar
  3  3     1.00000001e+02      # Atau(Q)MSSM DRbar
