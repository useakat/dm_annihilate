# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.3.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 09.11.2014,  12:19
Block SPINFO         # Program information
     1   SPheno      # spectrum calculator
     2   v3.3.3      # version number
#
Block SPhenoINFO     # SPheno specific information
    1      2         # using 2-loop RGEs
    2      1         # using running masses for boundary conditions at mZ
# Either the general MSSM or a model has been used
# which has not yet been implemented in the LesHouches standard
Block MINPAR  # Input parameters
    3    4.13322294E+00  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    1.00000000E+02  # M_1
    2    1.00000000E+02  # M_2
    3    2.00000000E+03  # M_3
   11    1.00000000E+03  # A_t
   12    1.00000000E+02  # A_b
   13    1.00000000E+02  # A_tau
   23    1.00000000E+05  # mu
   24    4.00000000E+10  # M^2_A(Q)
   25    4.00000000E+00  # tan(beta)
   31    2.00000000E+05  # M^2_L11
   32    2.00000000E+05  # M^2_L22
   33    2.00000000E+05  # M^2_L33
   34    2.00000000E+05  # M^2_E11
   35    2.00000000E+05  # M^2_E22
   36    2.00000000E+05  # M^2_E33
   41    2.00000000E+05  # M^2_Q11
   42    2.00000000E+05  # M^2_Q22
   43    2.00000000E+05  # M^2_Q33
   44    2.00000000E+05  # M^2_U11
   45    2.00000000E+05  # M^2_U22
   46    2.00000000E+05  # M^2_U33
   47    2.00000000E+05  # M^2_D11
   48    2.00000000E+05  # M^2_D22
   49    2.00000000E+05  # M^2_D33
Block SMINPUTS  # SM parameters
         1     1.27929502E+02  # alpha_em^-1(MZ)^MSbar
         2     1.16637900E-05  # G_mu [GeV^-2]
         3     1.18400000E-01  # alpha_s(MZ)^MSbar
         4     9.11876000E+01  # m_Z(pole)
         5     4.18000000E+00  # m_b(m_b), MSbar
         6     1.73100000E+02  # m_t(pole)
         7     1.77682000E+00  # m_tau(pole)
         8     0.00000000E+00  # m_nu_3
        11     5.10998930E-04  # m_e(pole)
        12     0.00000000E+00  # m_nu_1
        13     1.05658372E-01  # m_muon(pole)
        14     0.00000000E+00  # m_nu_2
        21     5.00000000E-03  # m_d(2 GeV), MSbar
        22     2.50000000E-03  # m_u(2 GeV), MSbar
        23     9.50000000E-02  # m_s(2 GeV), MSbar
        24     1.27000000E+00  # m_c(m_c), MSbar
Block gauge Q=  1.00000000E+03  # (SUSY scale)
   1    3.55353177E-01  # g'(Q)^DRbar
   2    6.24196313E-01  # g(Q)^DRbar
   3    9.83283424E-01  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     6.93648116E-06   # Y_u(Q)^DRbar
  2  2     3.52372999E-03   # Y_c(Q)^DRbar
  3  3     8.15253811E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     5.59300251E-05   # Y_d(Q)^DRbar
  2  2     1.06266884E-03   # Y_s(Q)^DRbar
  3  3     5.22410913E-02   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.20592621E-05   # Y_e(Q)^DRbar
  2  2     2.49347128E-03   # Y_mu(Q)^DRbar
  3  3     4.19240919E-02   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_u(Q)^DRbar
  2  2     0.00000000E+00   # A_c(Q)^DRbar
  3  3     1.00000000E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_d(Q)^DRbar
  2  2     0.00000000E+00   # A_s(Q)^DRbar
  3  3     1.00000000E+02   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     0.00000000E+00   # A_e(Q)^DRbar
  2  2     0.00000000E+00   # A_mu(Q)^DRbar
  3  3     1.00000000E+02   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    1.00000000E+02  # M_1
   2    1.00000000E+02  # M_2
   3    2.00000000E+03  # M_3
  21    2.45543828E+10  # M^2_(H,d)
  22   -1.56749923E+10  # M^2_(H,u)
  31    2.00000000E+05  # M_(L,11)
  32    2.00000000E+05  # M_(L,22)
  33    2.00000000E+05  # M_(L,33)
  34    2.00000000E+05  # M_(E,11)
  35    2.00000000E+05  # M_(E,22)
  36    2.00000000E+05  # M_(E,33)
  41    2.00000000E+05  # M_(Q,11)
  42    2.00000000E+05  # M_(Q,22)
  43    2.00000000E+05  # M_(Q,33)
  44    2.00000000E+05  # M_(U,11)
  45    2.00000000E+05  # M_(U,22)
  46    2.00000000E+05  # M_(U,33)
  47    2.00000000E+05  # M_(D,11)
  48    2.00000000E+05  # M_(D,22)
  49    2.00000000E+05  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73100000E+02  # m_t(pole)
        23     9.11876000E+01  # m_Z(pole)
        24     8.05291406E+01  # W+
        15     1.77682000E+00  # m_tau(pole)
        25     1.41271995E+02  # h0
        35     1.96338786E+05  # H0
        36     1.96266467E+05  # A0
        37     1.96266486E+05  # H+
   1000001     2.02247003E+05  # ~d_L
   2000001     2.01627232E+05  # ~d_R
   1000002     2.02246983E+05  # ~u_L
   2000002     2.02568207E+05  # ~u_R
   1000003     2.02247161E+05  # ~s_L
   2000003     2.01627271E+05  # ~s_R
   1000004     2.02247141E+05  # ~c_L
   2000004     2.02568475E+05  # ~c_R
   1000005     2.01721831E+05  # ~b_1
   2000005     2.09530004E+05  # ~b_2
   1000006     2.09530018E+05  # ~t_1
   2000006     2.16896542E+05  # ~t_2
   1000011     2.00835735E+05  # ~e_L-
   2000011     1.99365286E+05  # ~e_R-
   1000012     2.00835763E+05  # ~nu_eL
   1000013     2.00835834E+05  # ~mu_L-
   2000013     1.99365480E+05  # ~mu_R-
   1000014     2.00835855E+05  # ~nu_muL
   1000015     1.99420185E+05  # ~tau_1-
   2000015     2.00863655E+05  # ~tau_2-
   1000016     2.00861601E+05  # ~nu_tauL
   1000021     2.79288725E+03  # ~g
   1000022     3.11705366E+01  # ~chi_10
   1000023    -1.14730576E+02  # ~chi_20
   1000025     1.00322802E+05  # ~chi_30
   1000035    -1.00322841E+05  # ~chi_40
   1000024     1.14560278E+02  # ~chi_1+
   1000037     1.00322844E+05  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -2.44978837E-01   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    1.00000000E+05  # mu
   2    4.00000000E+00  # tan[beta](Q)
   3    2.49409344E+02  # v(Q)
   4    4.00000000E+10  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     9.99999577E-01   # Re[R_st(1,1)]
   1  2     9.19884147E-04   # Re[R_st(1,2)]
   2  1    -9.19884147E-04   # Re[R_st(2,1)]
   2  2     9.99999577E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     1.72302679E-04   # Re[R_sb(1,1)]
   1  2     9.99999985E-01   # Re[R_sb(1,2)]
   2  1    -9.99999985E-01   # Re[R_sb(2,1)]
   2  2     1.72302679E-04   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     1.17382916E-03   # Re[R_sta(1,1)]
   1  2     9.99999311E-01   # Re[R_sta(1,2)]
   2  1    -9.99999311E-01   # Re[R_sta(2,1)]
   2  2     1.17382916E-03   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1     9.99999894E-01   # Re[N(1,1)]
   1  2     1.11524413E-04   # Re[N(1,2)]
   1  3     4.31185825E-04   # Re[N(1,3)]
   1  4    -1.16070255E-04   # Re[N(1,4)]
   2  1    -1.11190381E-04   # Re[N(2,1)]
   2  2     9.99999714E-01   # Re[N(2,2)]
   2  3    -7.22576513E-04   # Re[N(2,3)]
   2  4     1.93387906E-04   # Re[N(2,4)]
   3  1     3.87040612E-04   # Re[N(3,1)]
   3  2    -6.47641565E-04   # Re[N(3,2)]
   3  3    -7.07106479E-01   # Re[N(3,3)]
   3  4     7.07106681E-01   # Re[N(3,4)]
   4  1    -2.22862016E-04   # Re[N(4,1)]
   4  2     3.74168024E-04   # Re[N(4,2)]
   4  3     7.07106583E-01   # Re[N(4,3)]
   4  4     7.07106846E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1     9.99999479E-01   # Re[U(1,1)]
   1  2    -1.02120506E-03   # Re[U(1,2)]
   2  1     1.02120506E-03   # Re[U(2,1)]
   2  2     9.99999479E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.99999963E-01   # Re[V(1,1)]
   1  2     2.73394147E-04   # Re[V(1,2)]
   2  1     2.73394147E-04   # Re[V(2,1)]
   2  2     9.99999963E-01   # Re[V(2,2)]
DECAY        23     2.49520000E+00   # Z
DECAY        24     2.08500000E+00   # W
DECAY   2000011     1.00168236E+03   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.99999876E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
DECAY   1000011     2.58736160E+03   # ~e^-_L
#    BR                NDA      ID1      ID2
     9.75381337E-02    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     3.00795195E-01    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     6.01666277E-01    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
DECAY   2000013     1.00171085E+03   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.99972408E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
DECAY   1000013     2.58737693E+03   # ~mu^-_L
#    BR                NDA      ID1      ID2
     9.75376049E-02    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     3.00793559E-01    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     6.01663004E-01    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
DECAY   1000015     1.00974094E+03   # ~tau^-_1
#    BR                NDA      ID1      ID2
     9.92291224E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     1.92636156E-03    2     1000025        15   # BR(~tau^-_1 -> chi^0_3 tau^-)
     1.92633002E-03    2     1000035        15   # BR(~tau^-_1 -> chi^0_4 tau^-)
     3.85252796E-03    2    -1000037        16   # BR(~tau^-_1 -> chi^-_2 nu_tau)
DECAY   2000015     2.59169715E+03   # ~tau^-_2
#    BR                NDA      ID1      ID2
     9.73888891E-02    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     3.00333291E-01    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     7.63353344E-04    2     1000025        15   # BR(~tau^-_2 -> chi^0_3 tau^-)
     7.63322614E-04    2     1000035        15   # BR(~tau^-_2 -> chi^0_4 tau^-)
     6.00742350E-01    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     2.58736236E+03   # ~nu_e
#    BR                NDA      ID1      ID2
     9.74617183E-02    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     3.00871320E-01    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     6.01666766E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
DECAY   1000014     2.58737760E+03   # ~nu_mu
#    BR                NDA      ID1      ID2
     9.74611886E-02    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     3.00869685E-01    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     6.01663496E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
DECAY   1000016     2.59167104E+03   # ~nu_tau
#    BR                NDA      ID1      ID2
     9.73122047E-02    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     3.00409761E-01    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     6.00743767E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
     1.52660164E-03    2     1000037        15   # BR(~nu_tau -> chi^+_2 tau^-)
DECAY   2000001     1.04506156E+04   # ~d_R
#    BR                NDA      ID1      ID2
     1.07707324E-02    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     9.89229265E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000001     1.27495859E+04   # ~d_L
#    BR                NDA      ID1      ID2
     2.21132884E-03    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     6.14817700E-02    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     1.22958293E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
     8.13348510E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000003     1.04506227E+04   # ~s_R
#    BR                NDA      ID1      ID2
     1.07707272E-02    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     9.89228780E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000003     1.27496268E+04   # ~s_L
#    BR                NDA      ID1      ID2
     2.21132348E-03    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     6.14816204E-02    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     1.22957994E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
     8.13346534E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   1000005     1.04679316E+04   # ~b_1
#    BR                NDA      ID1      ID2
     1.07579613E-02    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     2.96356421E-04    2     1000025         5   # BR(~b_1 -> chi^0_3 b)
     2.96356155E-04    2     1000035         5   # BR(~b_1 -> chi^0_4 b)
     5.92709272E-04    2    -1000037         6   # BR(~b_1 -> chi^-_2 t)
     9.88056602E-01    2     1000021         5   # BR(~b_1 -> ~g b)
DECAY   2000005     1.48617113E+04   # ~b_2
#    BR                NDA      ID1      ID2
     1.96537155E-03    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     5.46434035E-02    2     1000023         5   # BR(~b_2 -> chi^0_2 b)
     2.27386238E-04    2     1000025         5   # BR(~b_2 -> chi^0_3 b)
     2.27374450E-04    2     1000035         5   # BR(~b_2 -> chi^0_4 b)
     1.09282008E-01    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     1.10744326E-01    2    -1000037         6   # BR(~b_2 -> chi^-_2 t)
     7.22901915E-01    2     1000021         5   # BR(~b_2 -> ~g b)
DECAY   1000002     1.27495854E+04   # ~u_L
#    BR                NDA      ID1      ID2
     2.21653318E-03    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     6.14765775E-02    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     1.22958405E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
     8.13348463E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000002     1.08386829E+04   # ~u_R
#    BR                NDA      ID1      ID2
     4.17342573E-02    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     9.58265737E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000004     1.27496263E+04   # ~c_L
#    BR                NDA      ID1      ID2
     2.21652801E-03    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     6.14764278E-02    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     1.22958106E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
     8.13346488E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   2000004     1.08387542E+04   # ~c_R
#    BR                NDA      ID1      ID2
     4.17340374E-02    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     9.58260696E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000006     5.55180138E+04   # ~t_1
#    BR                NDA      ID1      ID2
     5.27361900E-04    2     1000022         6   # BR(~t_1 -> chi^0_1 t)
     1.46263318E-02    2     1000023         6   # BR(~t_1 -> chi^0_2 t)
     1.48226118E-02    2     1000025         6   # BR(~t_1 -> chi^0_3 t)
     1.48226889E-02    2     1000035         6   # BR(~t_1 -> chi^0_4 t)
     2.92539491E-02    2     1000024         5   # BR(~t_1 -> chi^+_1 b)
     1.21744025E-04    2     1000037         5   # BR(~t_1 -> chi^+_2 b)
     1.93514598E-01    2     1000021         6   # BR(~t_1 -> ~g t)
#    BR                NDA      ID1      ID2       ID3
     7.32309564E-01    3     1000022        24         5   # BR(~t_1 -> chi^0_1 W^+ b)
DECAY   2000006     1.51534147E+04   # ~t_2
#    BR                NDA      ID1      ID2
     3.19623720E-02    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     5.84707523E-02    2     1000025         6   # BR(~t_2 -> chi^0_3 t)
     5.84704277E-02    2     1000035         6   # BR(~t_2 -> chi^0_4 t)
     1.16941264E-01    2     1000037         5   # BR(~t_2 -> chi^+_2 b)
     7.33928633E-01    2     1000021         6   # BR(~t_2 -> ~g t)
     1.01362989E-04    2     2000005        24   # BR(~t_2 -> ~b_2 W^+)
DECAY   1000024     7.62108197E-09   # chi^+_1
#    BR                NDA      ID1      ID2
     9.99999999E-01    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000037     6.26471140E+02   # chi^+_2
#    BR                NDA      ID1      ID2
     9.61035905E-02    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     2.68785887E-01    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     2.77917457E-01    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     3.09985934E-01    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.17214116E-02    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     3.53561271E-02    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
     1.28777142E-04    3     1000021        -5         6   # BR(chi^+_2 -> ~g b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     6.64812124E-12   # chi^0_2
#    BR                NDA      ID1      ID2
#    BR                NDA      ID1      ID2       ID3
     1.38183286E-02    3     1000022         2        -2   # BR(chi^0_2 -> chi^0_1 u u_bar)
     1.79442352E-03    3     1000022         4        -4   # BR(chi^0_2 -> chi^0_1 c c_bar)
     4.72412883E-03    3     1000022         1        -1   # BR(chi^0_2 -> chi^0_1 d d_bar)
     7.96248822E-03    3     1000022         3        -3   # BR(chi^0_2 -> chi^0_1 s s_bar)
     1.66312550E-02    3     1000022        11       -11   # BR(chi^0_2 -> chi^0_1 e^- e^+)
     1.10885192E-02    3     1000022        13       -13   # BR(chi^0_2 -> chi^0_1 mu^- mu^+)
     1.17460203E-02    3     1000022        15       -15   # BR(chi^0_2 -> chi^0_1 tau^- tau^+)
     9.19430801E-02    3     1000022        12       -12   # BR(chi^0_2 -> chi^0_1 nu_e nu_bar_e)
     2.81310470E-01    3     1000024         1        -2   # BR(chi^0_2 -> chi^+_1 d u_bar)
     2.81310470E-01    3    -1000024        -1         2   # BR(chi^0_2 -> chi^-_1 d_bar u)
     1.38835377E-01    3     1000024        13       -14   # BR(chi^0_2 -> chi^+_1 mu^- nu_bar_mu)
     1.38835377E-01    3    -1000024       -13        14   # BR(chi^0_2 -> chi^-_1 mu^+ nu_mu)
DECAY   1000025     6.30458427E+02   # chi^0_3
#    BR                NDA      ID1      ID2
     2.67850909E-01    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     2.67850909E-01    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     2.45365312E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     6.91978870E-02    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     7.35008286E-02    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     2.26242912E-01    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     3.53065537E-02    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     3.53065537E-02    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
     1.90573173E-04    3     1000021         6        -6   # BR(chi^0_3 -> ~g t t_bar)
DECAY   1000035     6.27696024E+02   # chi^0_4
#    BR                NDA      ID1      ID2
     2.68793120E-01    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     2.68793120E-01    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     7.43297204E-02    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     2.08226951E-01    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     2.65437886E-02    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     8.21811167E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     3.54621663E-02    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     3.54621663E-02    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
     1.91411696E-04    3     1000021         6        -6   # BR(chi^0_4 -> ~g t t_bar)
DECAY   1000021     3.56186327E-09   # ~g
#    BR                NDA      ID1      ID2
#    BR                NDA      ID1      ID2       ID3
     3.50475592E-02    3     1000022         2        -2   # BR(~g -> chi^0_1 u u_bar)
     3.50473780E-02    3     1000022         4        -4   # BR(~g -> chi^0_1 c c_bar)
     2.51084387E-02    3     1000022         6        -6   # BR(~g -> chi^0_1 t t_bar)
     1.04691812E-02    3     1000022         1        -1   # BR(~g -> chi^0_1 d d_bar)
     1.04691682E-02    3     1000022         3        -3   # BR(~g -> chi^0_1 s s_bar)
     1.01797742E-02    3     1000022         5        -5   # BR(~g -> chi^0_1 b b_bar)
     5.12388293E-02    3     1000023         2        -2   # BR(~g -> chi^0_2 u u_bar)
     5.12386696E-02    3     1000023         4        -4   # BR(~g -> chi^0_2 c c_bar)
     4.17628455E-02    3     1000023         6        -6   # BR(~g -> chi^0_2 t t_bar)
     5.12431334E-02    3     1000023         1        -1   # BR(~g -> chi^0_2 d d_bar)
     5.12429736E-02    3     1000023         3        -3   # BR(~g -> chi^0_2 s s_bar)
     4.44792776E-02    3     1000023         5        -5   # BR(~g -> chi^0_2 b b_bar)
     1.02498331E-01    3     1000024         1        -2   # BR(~g -> chi^+_1 d u_bar)
     1.02498331E-01    3    -1000024        -1         2   # BR(~g -> chi^-_1 d_bar u)
     1.02498011E-01    3     1000024         3        -4   # BR(~g -> chi^+_1 s c_bar)
     1.02498011E-01    3    -1000024        -3         4   # BR(~g -> chi^-_1 s_bar c)
     8.62391074E-02    3     1000024         5        -6   # BR(~g -> chi^+_1 b t_bar)
     8.62391074E-02    3    -1000024        -5         6   # BR(~g -> chi^-_1 b_bar t)
DECAY        25     6.50770026E-03   # h^0
#    BR                NDA      ID1      ID2
     1.63980957E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     4.63066346E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     1.11673060E-04    2           3        -3   # BR(h^0 -> s s_bar)
     2.64684226E-01    2           5        -5   # BR(h^0 -> b b_bar)
     1.86877827E-02    2           4        -4   # BR(h^0 -> c c_bar)
     6.68342189E-02    2          21        21   # BR(h^0 -> g g)
     2.96092081E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     2.63555094E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     2.63555094E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     2.63555094E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     9.22442830E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     9.22442830E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     2.63555094E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     2.63555094E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     2.63555094E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     9.22442830E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     9.22442830E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     2.45744120E-03    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     2.45744120E-03    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     2.47206883E-03    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     1.46276262E-02    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     1.14095484E-02    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     1.14095484E-02    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     1.10438578E-02    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     8.48402319E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     8.77657571E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35                NaN   # H^0
DECAY        36                NaN   # A^0
DECAY        37     1.87391369E+03   # H^+
#    BR                NDA      ID1      ID2
     3.44687154E-03    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     2.60446540E-01    2          -5         6   # BR(H^+ -> b_bar t)
     2.21427799E-01    2     1000024   1000025   # BR(H^+ -> chi^+_1 chi^0_3)
     2.21427663E-01    2     1000024   1000035   # BR(H^+ -> chi^+_1 chi^0_4)
     7.18273058E-02    2     1000037   1000022   # BR(H^+ -> chi^+_2 chi^0_1)
     2.21399819E-01    2     1000037   1000023   # BR(H^+ -> chi^+_2 chi^0_2)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.00000139E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    1.59999986E+01    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    1.60000000E+01        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99999913E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    6.25000868E-02    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    6.25000000E-02        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.00000139E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    1.59999986E+01    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    1.60000000E+01        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    1.00000000E+00        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    3.01115290E-14        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    1.00000000E+00        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    3.01115290E-14        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.13774713E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    5.92837637E-02        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    2.54617431E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    3.01115290E-14        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    1.00000000E+00        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.20002583E-04   # BR(b -> s gamma)
    2    1.59099811E-06   # BR(b -> s mu+ mu-)
    3    3.51716525E-05   # BR(b -> s nu nu)
    4    2.05658423E-15   # BR(Bd -> e+ e-)
    5    8.78549428E-11   # BR(Bd -> mu+ mu-)
    6    1.83893972E-08   # BR(Bd -> tau+ tau-)
    7    6.93298940E-14   # BR(Bs -> e+ e-)
    8    2.96177084E-09   # BR(Bs -> mu+ mu-)
    9    6.28140125E-07   # BR(Bs -> tau+ tau-)
   10    9.67944567E-05   # BR(B_u -> tau nu)
   11    9.99999986E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
   12    5.41634197E-01   # |Delta(M_Bd)| [ps^-1] 
   13    1.93532378E+01   # |Delta(M_Bs)| [ps^-1] 
   16    2.15529998E-03   # epsilon_K
   17    2.28163865E-15   # Delta(M_K)
   18    2.47418372E-11   # BR(K^0 -> pi^0 nu nu)
   19    8.27669108E-11   # BR(K^+ -> pi^+ nu nu)
   20   -8.27315838E-21   # Delta(g-2)_electron/2
   21   -3.53702914E-16   # Delta(g-2)_muon/2
   22   -1.00076376E-13   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39   -1.82452077E-03   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.88674011E-01   # C7
     0305 4422   00   2    -1.88676506E-01   # C7
     0305 4322   00   2    -5.16821342E-08   # C7'
     0305 6421   00   0    -9.51816942E-02   # C8
     0305 6421   00   2    -9.51852808E-02   # C8
     0305 6321   00   2    -7.39280712E-08   # C8'
 03051111 4133   00   0     1.17242874E+00   # C9 e+e-
 03051111 4133   00   2     1.17242877E+00   # C9 e+e-
 03051111 4233   00   2     1.04304069E-07   # C9' e+e-
 03051111 4137   00   0    -3.98780252E+00   # C10 e+e-
 03051111 4137   00   2    -3.98780337E+00   # C10 e+e-
 03051111 4237   00   2    -2.83602457E-06   # C10' e+e-
 03051313 4133   00   0     1.17242874E+00   # C9 mu+mu-
 03051313 4133   00   2     1.17242877E+00   # C9 mu+mu-
 03051313 4233   00   2     1.04304069E-07   # C9' mu+mu-
 03051313 4137   00   0    -3.98780252E+00   # C10 mu+mu-
 03051313 4137   00   2    -3.98780337E+00   # C10 mu+mu-
 03051313 4237   00   2    -2.83602457E-06   # C10' mu+mu-
 03051212 4137   00   0     1.50325945E+00   # C11 nu_1 nu_1
 03051212 4137   00   2     1.50325965E+00   # C11 nu_1 nu_1
 03051212 4237   00   2     6.82930139E-07   # C11' nu_1 nu_1
 03051414 4137   00   0     1.50325945E+00   # C11 nu_2 nu_2
 03051414 4137   00   2     1.50325965E+00   # C11 nu_2 nu_2
 03051414 4237   00   2     6.82930139E-07   # C11' nu_2 nu_2
 03051616 4137   00   0     1.50325945E+00   # C11 nu_3 nu_3
 03051616 4137   00   2     1.50325965E+00   # C11 nu_3 nu_3
 03051616 4237   00   2     6.82930139E-07   # C11' nu_3 nu_3
Block IMFWCOEF Q=  1.60000000E+02  # Im(Wilson coefficients) at scale Q
#    id        order  M        value         comment
     0305 4422   00   0     3.84876572E-07   # C7
     0305 4422   00   2     3.84869599E-07   # C7
     0305 4322   00   2    -2.33100003E-13   # C7'
     0305 6421   00   0     3.29672707E-07   # C8
     0305 6421   00   2     3.29678986E-07   # C8
     0305 6321   00   2    -2.10779586E-12   # C8'
 03051111 4133   00   2     5.21309224E-10   # C9 e+e-
 03051111 4233   00   2     1.86254040E-09   # C9' e+e-
 03051111 4137   00   2     5.12780814E-11   # C10 e+e-
 03051111 4237   00   2    -5.06424184E-08   # C10' e+e-
 03051313 4133   00   2     5.21309215E-10   # C9 mu+mu-
 03051313 4233   00   2     1.86254040E-09   # C9' mu+mu-
 03051313 4137   00   2     5.12780922E-11   # C10 mu+mu-
 03051313 4237   00   2    -5.06424184E-08   # C10' mu+mu-
 03051212 4137   00   2    -4.36489906E-13   # C11 nu_1 nu_1
 03051212 4237   00   2     1.21949698E-08   # C11' nu_1 nu_1
 03051414 4137   00   2    -4.36489288E-13   # C11 nu_2 nu_2
 03051414 4237   00   2     1.21949698E-08   # C11' nu_2 nu_2
 03051616 4137   00   2    -4.36315302E-13   # C11 nu_3 nu_3
 03051616 4237   00   2     1.21949698E-08   # C11' nu_3 nu_3
