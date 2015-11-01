# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.3.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 09.11.2014,  02:35
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
    3    4.12954225E+00  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    1.00000000E+02  # M_1
    2    1.00000000E+02  # M_2
    3    1.00000000E+02  # M_3
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
         1     1.27929506E+02  # alpha_em^-1(MZ)^MSbar
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
   1    3.55367667E-01  # g'(Q)^DRbar
   2    6.24157815E-01  # g(Q)^DRbar
   3    1.02158855E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     6.80580157E-06   # Y_u(Q)^DRbar
  2  2     3.45734489E-03   # Y_c(Q)^DRbar
  3  3     7.99232745E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     5.48924132E-05   # Y_d(Q)^DRbar
  2  2     1.04295434E-03   # Y_s(Q)^DRbar
  3  3     5.13676636E-02   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.20490406E-05   # Y_e(Q)^DRbar
  2  2     2.49135779E-03   # Y_mu(Q)^DRbar
  3  3     4.18885730E-02   # Y_tau(Q)^DRbar
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
   3    1.00000000E+02  # M_3
  21    2.47161819E+10  # M^2_(H,d)
  22   -1.52222884E+10  # M^2_(H,u)
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
        24     8.05288150E+01  # W+
        15     1.77682000E+00  # m_tau(pole)
        25     1.37203125E+02  # h0
        35     1.96280106E+05  # H0
        36     1.96213666E+05  # A0
        37     1.96213686E+05  # H+
   1000001     2.02428640E+05  # ~d_L
   2000001     2.01802083E+05  # ~d_R
   1000002     2.02428620E+05  # ~u_L
   2000002     2.02741926E+05  # ~u_R
   1000003     2.02428793E+05  # ~s_L
   2000003     2.01802121E+05  # ~s_R
   1000004     2.02428773E+05  # ~c_L
   2000004     2.02742188E+05  # ~c_R
   1000005     2.01894187E+05  # ~b_1
   2000005     2.09500497E+05  # ~b_2
   1000006     2.09500508E+05  # ~t_1
   2000006     2.16633676E+05  # ~t_2
   1000011     2.00832503E+05  # ~e_L-
   2000011     1.99371583E+05  # ~e_R-
   1000012     2.00832531E+05  # ~nu_eL
   1000013     2.00832602E+05  # ~mu_L-
   2000013     1.99371777E+05  # ~mu_R-
   1000014     2.00832623E+05  # ~nu_muL
   1000015     1.99426444E+05  # ~tau_1-
   2000015     2.00860405E+05  # ~tau_2-
   1000016     2.00858472E+05  # ~nu_tauL
   1000021     1.76020347E+02  # ~g
   1000022     3.15544455E+01  # ~chi_10
   1000023    -1.13516349E+02  # ~chi_20
   1000025     1.00229822E+05  # ~chi_30
   1000035    -1.00229861E+05  # ~chi_40
   1000024     1.13344647E+02  # ~chi_1+
   1000037     1.00229863E+05  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -2.44978832E-01   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    1.00000000E+05  # mu
   2    4.00000000E+00  # tan[beta](Q)
   3    2.49588917E+02  # v(Q)
   4    4.00000000E+10  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     9.99999567E-01   # Re[R_st(1,1)]
   1  2     9.30366352E-04   # Re[R_st(1,2)]
   2  1    -9.30366352E-04   # Re[R_st(2,1)]
   2  2     9.99999567E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     1.69589407E-04   # Re[R_sb(1,1)]
   1  2     9.99999986E-01   # Re[R_sb(1,2)]
   2  1    -9.99999986E-01   # Re[R_sb(2,1)]
   2  2     1.69589407E-04   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     1.18300371E-03   # Re[R_sta(1,1)]
   1  2     9.99999300E-01   # Re[R_sta(1,2)]
   2  1    -9.99999300E-01   # Re[R_sta(2,1)]
   2  2     1.18300371E-03   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1     9.99999893E-01   # Re[N(1,1)]
   1  2     1.12710440E-04   # Re[N(1,2)]
   1  3     4.32903366E-04   # Re[N(1,3)]
   1  4    -1.16273902E-04   # Re[N(1,4)]
   2  1    -1.12373712E-04   # Re[N(2,1)]
   2  2     9.99999712E-01   # Re[N(2,2)]
   2  3    -7.25743367E-04   # Re[N(2,3)]
   2  4     1.93767294E-04   # Re[N(2,4)]
   3  1     3.88400146E-04   # Re[N(3,1)]
   3  2    -6.50148525E-04   # Re[N(3,2)]
   3  3    -7.07106477E-01   # Re[N(3,3)]
   3  4     7.07106680E-01   # Re[N(3,4)]
   4  1    -2.23933166E-04   # Re[N(4,1)]
   4  2     3.76138675E-04   # Re[N(4,2)]
   4  3     7.07106581E-01   # Re[N(4,3)]
   4  4     7.07106846E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1     9.99999474E-01   # Re[U(1,1)]
   1  2    -1.02570721E-03   # Re[U(1,2)]
   2  1     1.02570721E-03   # Re[U(2,1)]
   2  2     9.99999474E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.99999962E-01   # Re[V(1,1)]
   1  2     2.73934261E-04   # Re[V(1,2)]
   2  1     2.73934261E-04   # Re[V(2,1)]
   2  2     9.99999962E-01   # Re[V(2,2)]
DECAY        23     2.49520000E+00   # Z
DECAY        24     2.08500000E+00   # W
DECAY   2000011     1.00179569E+03   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.99999874E-01    2     1000022        11   # BR(~e^-_R -> chi^0_1 e^-)
DECAY   1000011     2.58705253E+03   # ~e^-_L
#    BR                NDA      ID1      ID2
     9.75565740E-02    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     3.00788775E-01    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     6.01654252E-01    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
DECAY   2000013     1.00182417E+03   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.99972419E-01    2     1000022        13   # BR(~mu^-_R -> chi^0_1 mu^-)
DECAY   1000013     2.58706786E+03   # ~mu^-_L
#    BR                NDA      ID1      ID2
     9.75560454E-02    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     3.00787140E-01    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     6.01650980E-01    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
DECAY   1000015     1.00985129E+03   # ~tau^-_1
#    BR                NDA      ID1      ID2
     9.92294835E-01    2     1000022        15   # BR(~tau^-_1 -> chi^0_1 tau^-)
     1.92544589E-03    2     1000025        15   # BR(~tau^-_1 -> chi^0_3 tau^-)
     1.92541420E-03    2     1000035        15   # BR(~tau^-_1 -> chi^0_4 tau^-)
     3.85069478E-03    2    -1000037        16   # BR(~tau^-_1 -> chi^-_2 nu_tau)
DECAY   2000015     2.59138566E+03   # ~tau^-_2
#    BR                NDA      ID1      ID2
     9.74073705E-02    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     3.00327069E-01    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     7.63062857E-04    2     1000025        15   # BR(~tau^-_2 -> chi^0_3 tau^-)
     7.63031868E-04    2     1000035        15   # BR(~tau^-_2 -> chi^0_4 tau^-)
     6.00730719E-01    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     2.58705330E+03   # ~nu_e
#    BR                NDA      ID1      ID2
     9.74793396E-02    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     3.00865717E-01    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     6.01654746E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
DECAY   1000014     2.58706854E+03   # ~nu_mu
#    BR                NDA      ID1      ID2
     9.74788101E-02    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     3.00864083E-01    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     6.01651478E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
DECAY   1000016     2.59136121E+03   # ~nu_tau
#    BR                NDA      ID1      ID2
     9.73298586E-02    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     3.00404350E-01    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     6.00732133E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
     1.52602087E-03    2     1000037        15   # BR(~nu_tau -> chi^+_2 tau^-)
DECAY   2000001     1.12858245E+04   # ~d_R
#    BR                NDA      ID1      ID2
     9.98310624E-03    2     1000022         1   # BR(~d_R -> chi^0_1 d)
     9.90016891E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000001     1.35894259E+04   # ~d_L
#    BR                NDA      ID1      ID2
     2.07667340E-03    2     1000022         1   # BR(~d_L -> chi^0_1 d)
     5.77268450E-02    2     1000023         1   # BR(~d_L -> chi^0_2 d)
     1.15448711E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
     8.24747678E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000003     1.12858315E+04   # ~s_R
#    BR                NDA      ID1      ID2
     9.98310185E-03    2     1000022         3   # BR(~s_R -> chi^0_1 s)
     9.90016457E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000003     1.35894661E+04   # ~s_L
#    BR                NDA      ID1      ID2
     2.07666884E-03    2     1000022         3   # BR(~s_L -> chi^0_1 s)
     5.77267178E-02    2     1000023         3   # BR(~s_L -> chi^0_2 s)
     1.15448457E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
     8.24745862E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   1000005     1.13030113E+04   # ~b_1
#    BR                NDA      ID1      ID2
     9.97247643E-03    2     1000022         5   # BR(~b_1 -> chi^0_1 b)
     2.66208334E-04    2     1000025         5   # BR(~b_1 -> chi^0_3 b)
     2.66208074E-04    2     1000035         5   # BR(~b_1 -> chi^0_4 b)
     5.32413323E-04    2    -1000037         6   # BR(~b_1 -> chi^-_2 t)
     9.88962680E-01    2     1000021         5   # BR(~b_1 -> ~g b)
DECAY   2000005     1.56538792E+04   # ~b_2
#    BR                NDA      ID1      ID2
     1.86577979E-03    2     1000022         5   # BR(~b_2 -> chi^0_1 b)
     5.18644828E-02    2     1000023         5   # BR(~b_2 -> chi^0_2 b)
     2.08887422E-04    2     1000025         5   # BR(~b_2 -> chi^0_3 b)
     2.08876170E-04    2     1000035         5   # BR(~b_2 -> chi^0_4 b)
     1.03724363E-01    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     1.01128596E-01    2    -1000037         6   # BR(~b_2 -> chi^-_2 t)
     7.40991760E-01    2     1000021         5   # BR(~b_2 -> ~g b)
DECAY   1000002     1.35894253E+04   # ~u_L
#    BR                NDA      ID1      ID2
     2.08161237E-03    2     1000022         2   # BR(~u_L -> chi^0_1 u)
     5.77219175E-02    2     1000023         2   # BR(~u_L -> chi^0_2 u)
     1.15448817E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
     8.24747633E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000002     1.16779625E+04   # ~u_R
#    BR                NDA      ID1      ID2
     3.87712517E-02    2     1000022         2   # BR(~u_R -> chi^0_1 u)
     9.61228743E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000004     1.35894655E+04   # ~c_L
#    BR                NDA      ID1      ID2
     2.08160799E-03    2     1000022         4   # BR(~c_L -> chi^0_1 c)
     5.77217901E-02    2     1000023         4   # BR(~c_L -> chi^0_2 c)
     1.15448563E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
     8.24745818E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   2000004     1.16780326E+04   # ~c_R
#    BR                NDA      ID1      ID2
     3.87710687E-02    2     1000022         4   # BR(~c_R -> chi^0_1 c)
     9.61224211E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000006     5.62911250E+04   # ~t_1
#    BR                NDA      ID1      ID2
     5.20094726E-04    2     1000022         6   # BR(~t_1 -> chi^0_1 t)
     1.44216339E-02    2     1000023         6   # BR(~t_1 -> chi^0_2 t)
     1.40612626E-02    2     1000025         6   # BR(~t_1 -> chi^0_3 t)
     1.40613371E-02    2     1000035         6   # BR(~t_1 -> chi^0_4 t)
     2.88445486E-02    2     1000024         5   # BR(~t_1 -> chi^+_1 b)
     1.16182866E-04    2     1000037         5   # BR(~t_1 -> chi^+_2 b)
     2.06060552E-01    2     1000021         6   # BR(~t_1 -> ~g t)
#    BR                NDA      ID1      ID2       ID3
     7.21913371E-01    3     1000022        24         5   # BR(~t_1 -> chi^0_1 W^+ b)
DECAY   2000006     1.58823172E+04   # ~t_2
#    BR                NDA      ID1      ID2
     3.04610168E-02    2     1000022         6   # BR(~t_2 -> chi^0_1 t)
     5.35345997E-02    2     1000025         6   # BR(~t_2 -> chi^0_3 t)
     5.35342960E-02    2     1000035         6   # BR(~t_2 -> chi^0_4 t)
     1.07068970E-01    2     1000037         5   # BR(~t_2 -> chi^+_2 b)
     7.55199836E-01    2     1000021         6   # BR(~t_2 -> ~g t)
DECAY   1000024     5.11167055E-09   # chi^+_1
#    BR                NDA      ID1      ID2
     9.99999999E-01    2     1000022        24   # BR(chi^+_1 -> chi^0_1 W^+)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000037     6.27301033E+02   # chi^+_2
#    BR                NDA      ID1      ID2
     9.64346098E-02    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     2.69945090E-01    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     2.79113847E-01    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     3.09254057E-01    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
     1.12348771E-02    3     1000022        -5         6   # BR(chi^+_2 -> chi^0_1 b_bar t)
     3.38828401E-02    3     1000023        -5         6   # BR(chi^+_2 -> chi^0_2 b_bar t)
     1.33887492E-04    3     1000021        -5         6   # BR(chi^+_2 -> ~g b_bar t)
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     2.23934658E-12   # chi^0_2
#    BR                NDA      ID1      ID2
#    BR                NDA      ID1      ID2       ID3
     1.15726310E-01    3     1000022         3        -3   # BR(chi^0_2 -> chi^0_1 s s_bar)
     1.97706742E-02    3     1000022        11       -11   # BR(chi^0_2 -> chi^0_1 e^- e^+)
     2.06042719E-02    3     1000022        13       -13   # BR(chi^0_2 -> chi^0_1 mu^- mu^+)
     2.06413480E-04    3     1000022        15       -15   # BR(chi^0_2 -> chi^0_1 tau^- tau^+)
     1.52326499E-01    3     1000022        12       -12   # BR(chi^0_2 -> chi^0_1 nu_e nu_bar_e)
     3.45682828E-01    3     1000024         1        -2   # BR(chi^0_2 -> chi^+_1 d u_bar)
     3.45682828E-01    3    -1000024        -1         2   # BR(chi^0_2 -> chi^-_1 d_bar u)
DECAY   1000025     6.30782925E+02   # chi^0_3
#    BR                NDA      ID1      ID2
     2.69190880E-01    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     2.69190880E-01    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     2.46895250E-02    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     6.96935738E-02    2     1000023        23   # BR(chi^0_3 -> chi^0_2 Z)
     7.34012325E-02    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
     2.25894439E-01    2     1000023        25   # BR(chi^0_3 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     3.38626937E-02    3     1000024         5        -6   # BR(chi^0_3 -> chi^+_1 b t_bar)
     3.38626937E-02    3    -1000024        -5         6   # BR(chi^0_3 -> chi^-_1 b_bar t)
     1.98306931E-04    3     1000021         6        -6   # BR(chi^0_3 -> ~g t t_bar)
DECAY   1000035     6.28542344E+02   # chi^0_4
#    BR                NDA      ID1      ID2
     2.69924440E-01    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     2.69924440E-01    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     7.45386236E-02    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     2.08962456E-01    2     1000023        23   # BR(chi^0_4 -> chi^0_2 Z)
     2.64853187E-02    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
     8.19826065E-02    2     1000023        25   # BR(chi^0_4 -> chi^0_2 h^0)
#    BR                NDA      ID1      ID2       ID3
     3.39836219E-02    3     1000024         5        -6   # BR(chi^0_4 -> chi^+_1 b t_bar)
     3.39836219E-02    3    -1000024        -5         6   # BR(chi^0_4 -> chi^-_1 b_bar t)
     1.99014337E-04    3     1000021         6        -6   # BR(chi^0_4 -> ~g t t_bar)
DECAY   1000021     4.33573545E-16   # ~g
#    BR                NDA      ID1      ID2
#    BR                NDA      ID1      ID2       ID3
     3.01320520E-01    3     1000022         2        -2   # BR(~g -> chi^0_1 u u_bar)
     3.01055385E-01    3     1000022         4        -4   # BR(~g -> chi^0_1 c c_bar)
     9.00024392E-02    3     1000022         1        -1   # BR(~g -> chi^0_1 d d_bar)
     9.00026135E-02    3     1000022         3        -3   # BR(~g -> chi^0_1 s s_bar)
     8.69735942E-02    3     1000022         5        -5   # BR(~g -> chi^0_1 b b_bar)
     1.02687388E-02    3     1000023         2        -2   # BR(~g -> chi^0_2 u u_bar)
     1.00320572E-02    3     1000023         4        -4   # BR(~g -> chi^0_2 c c_bar)
     1.02105955E-02    3     1000023         1        -1   # BR(~g -> chi^0_2 d d_bar)
     1.02795271E-02    3     1000023         3        -3   # BR(~g -> chi^0_2 s s_bar)
     8.46145766E-03    3     1000023         5        -5   # BR(~g -> chi^0_2 b b_bar)
     2.03236594E-02    3     1000024         1        -2   # BR(~g -> chi^+_1 d u_bar)
     2.03236594E-02    3    -1000024        -1         2   # BR(~g -> chi^-_1 d_bar u)
     2.03636109E-02    3     1000024         3        -4   # BR(~g -> chi^+_1 s c_bar)
     2.03636109E-02    3    -1000024        -3         4   # BR(~g -> chi^-_1 s_bar c)
DECAY        25     5.06332448E-03   # h^0
#    BR                NDA      ID1      ID2
     2.04454981E-04    2          13       -13   # BR(h^0 -> mu^- mu^+)
     5.77327597E-02    2          15       -15   # BR(h^0 -> tau^- tau^+)
     1.37116840E-04    2           3        -3   # BR(h^0 -> s s_bar)
     3.26303069E-01    2           5        -5   # BR(h^0 -> b b_bar)
     2.29496568E-02    2           4        -4   # BR(h^0 -> c c_bar)
     8.10913973E-02    2          21        21   # BR(h^0 -> g g)
     3.28886615E-03    2          22        22   # BR(h^0 -> photon photon)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     2.23340881E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     2.23340881E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     2.23340881E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     7.81693084E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     7.81693084E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     2.23340881E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     2.23340881E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     2.23340881E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     7.81693084E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     7.81693084E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     2.07003498E-03    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     2.07003498E-03    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     2.08235662E-03    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     1.23216368E-02    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     9.61087670E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     9.61087670E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     9.30283578E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     7.14654934E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     7.39298208E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        35                NaN   # H^0
DECAY        36                NaN   # A^0
DECAY        37     1.85572895E+03   # H^+
#    BR                NDA      ID1      ID2
     3.47381823E-03    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     2.52788916E-01    2          -5         6   # BR(H^+ -> b_bar t)
     2.23718914E-01    2     1000024   1000025   # BR(H^+ -> chi^+_1 chi^0_3)
     2.23718777E-01    2     1000024   1000035   # BR(H^+ -> chi^+_1 chi^0_4)
     7.25854408E-02    2     1000037   1000022   # BR(H^+ -> chi^+_2 chi^0_1)
     2.23690341E-01    2     1000037   1000023   # BR(H^+ -> chi^+_2 chi^0_2)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.00000135E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    1.59999987E+01    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    1.60000000E+01        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99999916E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    6.25000843E-02    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    6.25000000E-02        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.00000135E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    1.59999987E+01    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    1.60000000E+01        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    1.00000000E+00        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    2.84549436E-14        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    1.00000000E+00        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    2.84549436E-14        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.15085989E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    5.92968391E-02        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    2.54625280E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    2.84549436E-14        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    1.00000000E+00        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.20002522E-04   # BR(b -> s gamma)
    2    1.59099812E-06   # BR(b -> s mu+ mu-)
    3    3.51719734E-05   # BR(b -> s nu nu)
    4    2.05592758E-15   # BR(Bd -> e+ e-)
    5    8.78268917E-11   # BR(Bd -> mu+ mu-)
    6    1.83835265E-08   # BR(Bd -> tau+ tau-)
    7    6.93091516E-14   # BR(Bs -> e+ e-)
    8    2.96088473E-09   # BR(Bs -> mu+ mu-)
    9    6.27952230E-07   # BR(Bs -> tau+ tau-)
   10    9.67944568E-05   # BR(B_u -> tau nu)
   11    9.99999986E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
   12    5.41627656E-01   # |Delta(M_Bd)| [ps^-1] 
   13    1.93530089E+01   # |Delta(M_Bs)| [ps^-1] 
   16    2.15529636E-03   # epsilon_K
   17    2.28163860E-15   # Delta(M_K)
   18    2.47420629E-11   # BR(K^0 -> pi^0 nu nu)
   19    8.27674480E-11   # BR(K^+ -> pi^+ nu nu)
   20   -8.44682256E-21   # Delta(g-2)_electron/2
   21   -3.61128067E-16   # Delta(g-2)_muon/2
   22   -1.02175736E-13   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39   -1.82684049E-03   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.88674487E-01   # C7
     0305 4422   00   2    -1.88676926E-01   # C7
     0305 4322   00   2    -5.03968591E-08   # C7'
     0305 6421   00   0    -9.51818457E-02   # C8
     0305 6421   00   2    -9.51853447E-02   # C8
     0305 6321   00   2    -7.19813594E-08   # C8'
 03051111 4133   00   0     1.17160182E+00   # C9 e+e-
 03051111 4133   00   2     1.17160185E+00   # C9 e+e-
 03051111 4233   00   2     9.92756282E-08   # C9' e+e-
 03051111 4137   00   0    -3.98699959E+00   # C10 e+e-
 03051111 4137   00   2    -3.98700044E+00   # C10 e+e-
 03051111 4237   00   2    -2.71406265E-06   # C10' e+e-
 03051313 4133   00   0     1.17160182E+00   # C9 mu+mu-
 03051313 4133   00   2     1.17160185E+00   # C9 mu+mu-
 03051313 4233   00   2     9.92756282E-08   # C9' mu+mu-
 03051313 4137   00   0    -3.98699959E+00   # C10 mu+mu-
 03051313 4137   00   2    -3.98700044E+00   # C10 mu+mu-
 03051313 4237   00   2    -2.71406265E-06   # C10' mu+mu-
 03051212 4137   00   0     1.50326630E+00   # C11 nu_1 nu_1
 03051212 4137   00   2     1.50326651E+00   # C11 nu_1 nu_1
 03051212 4237   00   2     6.53696770E-07   # C11' nu_1 nu_1
 03051414 4137   00   0     1.50326630E+00   # C11 nu_2 nu_2
 03051414 4137   00   2     1.50326651E+00   # C11 nu_2 nu_2
 03051414 4237   00   2     6.53696770E-07   # C11' nu_2 nu_2
 03051616 4137   00   0     1.50326630E+00   # C11 nu_3 nu_3
 03051616 4137   00   2     1.50326651E+00   # C11 nu_3 nu_3
 03051616 4237   00   2     6.53696770E-07   # C11' nu_3 nu_3
Block IMFWCOEF Q=  1.60000000E+02  # Im(Wilson coefficients) at scale Q
#    id        order  M        value         comment
     0305 4422   00   0     3.84879756E-07   # C7
     0305 4422   00   2     3.84869586E-07   # C7
     0305 4322   00   2    -3.58287122E-13   # C7'
     0305 6421   00   0     3.29675432E-07   # C8
     0305 6421   00   2     3.29771173E-07   # C8
     0305 6321   00   2    -4.48077073E-13   # C8'
 03051111 4133   00   2     5.08727257E-10   # C9 e+e-
 03051111 4233   00   2     1.77279697E-09   # C9' e+e-
 03051111 4137   00   2     4.96657090E-11   # C10 e+e-
 03051111 4237   00   2    -4.84658936E-08   # C10' e+e-
 03051313 4133   00   2     5.08727248E-10   # C9 mu+mu-
 03051313 4233   00   2     1.77279697E-09   # C9' mu+mu-
 03051313 4137   00   2     4.96657195E-11   # C10 mu+mu-
 03051313 4237   00   2    -4.84658936E-08   # C10' mu+mu-
 03051212 4137   00   2    -3.23937339E-13   # C11 nu_1 nu_1
 03051212 4237   00   2     1.16732744E-08   # C11' nu_1 nu_1
 03051414 4137   00   2    -3.23936735E-13   # C11 nu_2 nu_2
 03051414 4237   00   2     1.16732744E-08   # C11' nu_2 nu_2
 03051616 4137   00   2    -3.23766742E-13   # C11 nu_3 nu_3
 03051616 4237   00   2     1.16732744E-08   # C11' nu_3 nu_3
