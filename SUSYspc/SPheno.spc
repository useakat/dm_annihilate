# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.3.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 11.11.2014,  12:54
Block SPINFO         # Program information
     1   SPheno      # spectrum calculator
     2   v3.3.3      # version number
#
Block SPhenoINFO     # SPheno specific information
    1      1         # using 1-loop RGEs
    2      1         # using running masses for boundary conditions at mZ
# Either the general MSSM or a model has been used
# which has not yet been implemented in the LesHouches standard
Block MINPAR  # Input parameters
    3    1.08883557E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    1.00000000E+02  # M_1
    2    1.00000000E+02  # M_2
    3    1.50000000E+03  # M_3
   11    1.00000000E+03  # A_t
   12    1.00000000E+02  # A_b
   13    1.00000000E+02  # A_tau
   23    2.00000000E+05  # mu
   25    1.00000000E+01  # tan(beta)
   26    2.00000000E+05  # m_A, pole mass
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
         1     1.27929573E+02  # alpha_em^-1(MZ)^MSbar
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
Block gauge Q=  2.00000037E+05  # (SUSY scale)
   1    3.72868285E-01  # g'(Q)^DRbar
   2    6.30898769E-01  # g(Q)^DRbar
   3    8.99517252E-01  # g3(Q)^DRbar
Block Yu Q=  2.00000037E+05  # (SUSY scale)
  1  1     5.82904019E-06   # Y_u(Q)^DRbar
  2  2     2.96115349E-03   # Y_c(Q)^DRbar
  3  3     7.27682611E-01   # Y_t(Q)^DRbar
Block Yd Q=  2.00000037E+05  # (SUSY scale)
  1  1     1.17578480E-04   # Y_d(Q)^DRbar
  2  2     2.23398846E-03   # Y_s(Q)^DRbar
  3  3     1.12261103E-01   # Y_b(Q)^DRbar
Block Ye Q=  2.00000037E+05  # (SUSY scale)
  1  1     2.94744358E-05   # Y_e(Q)^DRbar
  2  2     6.09437494E-03   # Y_mu(Q)^DRbar
  3  3     1.02470792E-01   # Y_tau(Q)^DRbar
Block Au Q=  2.00000037E+05  # (SUSY scale)
  1  1     0.00000000E+00   # A_u(Q)^DRbar
  2  2     0.00000000E+00   # A_c(Q)^DRbar
  3  3     1.00000000E+03   # A_t(Q)^DRbar
Block Ad Q=  2.00000037E+05  # (SUSY scale)
  1  1     0.00000000E+00   # A_d(Q)^DRbar
  2  2     0.00000000E+00   # A_s(Q)^DRbar
  3  3     1.00000000E+02   # A_b(Q)^DRbar
Block Ae Q=  2.00000037E+05  # (SUSY scale)
  1  1     0.00000000E+00   # A_e(Q)^DRbar
  2  2     0.00000000E+00   # A_mu(Q)^DRbar
  3  3     1.00000000E+02   # A_tau(Q)^DRbar
Block MSOFT Q=  2.00000037E+05  # soft SUSY breaking masses at Q
   1    1.00000000E+02  # M_1
   2    1.00000000E+02  # M_2
   3    1.50000000E+03  # M_3
  21   -9.63153348E+08  # M^2_(H,d)
  22   -3.91194283E+10  # M^2_(H,u)
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
        24     8.05239790E+01  # W+
        15     1.77682000E+00  # m_tau(pole)
        25     1.35270522E+02  # h0
        35     2.00000314E+05  # H0
        36     2.00000000E+05  # A0
        37     1.99999989E+05  # H+
   1000001     2.02068914E+05  # ~d_L
   2000001     2.01614736E+05  # ~d_R
   1000002     2.02068891E+05  # ~u_L
   2000002     2.01593897E+05  # ~u_R
   1000003     2.02068901E+05  # ~s_L
   2000003     2.01614714E+05  # ~s_R
   1000004     2.02068878E+05  # ~c_L
   2000004     2.01593893E+05  # ~c_R
   1000005     2.01560496E+05  # ~b_1
   2000005     2.01929637E+05  # ~b_2
   1000006     2.01382327E+05  # ~t_1
   2000006     2.01929458E+05  # ~t_2
   1000011     2.00398983E+05  # ~e_L-
   2000011     2.00268273E+05  # ~e_R-
   1000012     2.00399041E+05  # ~nu_eL
   1000013     2.00398911E+05  # ~mu_L-
   2000013     2.00268128E+05  # ~mu_R-
   1000014     2.00398967E+05  # ~nu_muL
   1000015     2.00227212E+05  # ~tau_1-
   2000015     2.00378727E+05  # ~tau_2-
   1000016     2.00378072E+05  # ~nu_tauL
   1000021     2.27346079E+03  # ~g
   1000022     2.38563630E+01  # ~chi_10
   1000023     6.50806412E+01  # ~chi_20
   1000025     2.00484261E+05  # ~chi_30
   1000035    -2.00484268E+05  # ~chi_40
   1000024     2.39161474E+01  # ~chi_1+
   1000037     2.00484247E+05  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -9.96687178E-02   # alpha
Block Hmix Q=  2.00000037E+05  # Higgs mixing parameters
   1    2.00000000E+05  # mu
   2    1.00000000E+01  # tan[beta](Q)
   3    2.44365238E+02  # v(Q)
   4    3.94505188E+10  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     1.11284390E-02   # Re[R_st(1,1)]
   1  2     9.99938077E-01   # Re[R_st(1,2)]
   2  1    -9.99938077E-01   # Re[R_st(2,1)]
   2  2     1.11284390E-02   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     2.68617685E-02   # Re[R_sb(1,1)]
   1  2     9.99639158E-01   # Re[R_sb(1,2)]
   2  1    -9.99639158E-01   # Re[R_sb(2,1)]
   2  2     2.68617685E-02   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     6.07935639E-02   # Re[R_sta(1,1)]
   1  2     9.98150361E-01   # Re[R_sta(1,2)]
   2  1    -9.98150361E-01   # Re[R_sta(2,1)]
   2  2     6.07935639E-02   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1    -8.15475205E-05   # Re[N(1,1)]
   1  2     9.99999924E-01   # Re[N(1,2)]
   1  3    -3.78858787E-04   # Re[N(1,3)]
   1  4     3.82370656E-05   # Re[N(1,4)]
   2  1     9.99999971E-01   # Re[N(2,1)]
   2  2     8.16344621E-05   # Re[N(2,2)]
   2  3     2.27125864E-04   # Re[N(2,3)]
   2  4    -2.32543270E-05   # Re[N(2,4)]
   3  1    -1.77069595E-04   # Re[N(3,1)]
   3  2     2.94916860E-04   # Re[N(3,2)]
   3  3     7.07106707E-01   # Re[N(3,3)]
   3  4    -7.07106772E-01   # Re[N(3,4)]
   4  1    -1.44178597E-04   # Re[N(4,1)]
   4  2     2.40844165E-04   # Re[N(4,2)]
   4  3     7.07106718E-01   # Re[N(4,3)]
   4  4     7.07106789E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.99999856E-01   # Re[U(1,1)]
   1  2     5.35825740E-04   # Re[U(1,2)]
   2  1     5.35825740E-04   # Re[U(2,1)]
   2  2     9.99999856E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.99999999E-01   # Re[V(1,1)]
   1  2     5.40806215E-05   # Re[V(1,2)]
   2  1     5.40806215E-05   # Re[V(2,1)]
   2  2     9.99999999E-01   # Re[V(2,2)]
DECAY        23     2.49520000E+00   # Z
DECAY        24     2.08500000E+00   # W
DECAY   2000011     1.10785538E+03   # ~e^-_R
#    BR                NDA      ID1      ID2
     9.99999992E-01    2     1000023        11   # BR(~e^-_R -> chi^0_2 e^-)
DECAY   1000011     2.65746915E+03   # ~e^-_L
#    BR                NDA      ID1      ID2
     2.98541594E-01    2     1000022        11   # BR(~e^-_L -> chi^0_1 e^-)
     1.04317740E-01    2     1000023        11   # BR(~e^-_L -> chi^0_2 e^-)
     5.97140665E-01    2    -1000024        12   # BR(~e^-_L -> chi^-_1 nu_e)
DECAY   2000013     1.10788174E+03   # ~mu^-_R
#    BR                NDA      ID1      ID2
     9.99962373E-01    2     1000023        13   # BR(~mu^-_R -> chi^0_2 mu^-)
DECAY   1000013     2.65744105E+03   # ~mu^-_L
#    BR                NDA      ID1      ID2
     2.98539413E-01    2     1000022        13   # BR(~mu^-_L -> chi^0_1 mu^-)
     1.04324269E-01    2     1000023        13   # BR(~mu^-_L -> chi^0_2 mu^-)
     5.97136303E-01    2    -1000024        14   # BR(~mu^-_L -> chi^-_1 nu_mu)
DECAY   1000015     1.11337591E+03   # ~tau^-_1
#    BR                NDA      ID1      ID2
     9.92082942E-01    2     1000023        15   # BR(~tau^-_1 -> chi^0_2 tau^-)
DECAY   2000015     2.65147914E+03   # ~tau^-_2
#    BR                NDA      ID1      ID2
     2.98076886E-01    2     1000022        15   # BR(~tau^-_2 -> chi^0_1 tau^-)
     1.05700419E-01    2     1000023        15   # BR(~tau^-_2 -> chi^0_2 tau^-)
     5.96211155E-01    2    -1000024        16   # BR(~tau^-_2 -> chi^-_1 nu_tau)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000012     2.65747021E+03   # ~nu_e
#    BR                NDA      ID1      ID2
     2.98599121E-01    2     1000022        12   # BR(~nu_e -> chi^0_1 nu_e)
     1.04260109E-01    2     1000023        12   # BR(~nu_e -> chi^0_2 nu_e)
     5.97140770E-01    2     1000024        11   # BR(~nu_e -> chi^+_1 e^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000014     2.65746932E+03   # ~nu_mu
#    BR                NDA      ID1      ID2
     2.98599111E-01    2     1000022        14   # BR(~nu_mu -> chi^0_1 nu_mu)
     1.04260105E-01    2     1000023        14   # BR(~nu_mu -> chi^0_2 nu_mu)
     5.97140749E-01    2     1000024        13   # BR(~nu_mu -> chi^+_1 mu^-)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000016     2.65722959E+03   # ~nu_tau
#    BR                NDA      ID1      ID2
     2.98594912E-01    2     1000022        16   # BR(~nu_tau -> chi^0_1 nu_tau)
     1.04258639E-01    2     1000023        16   # BR(~nu_tau -> chi^0_2 nu_tau)
     5.97132358E-01    2     1000024        15   # BR(~nu_tau -> chi^+_1 tau^-)
DECAY   2000001     8.77618521E+03   # ~d_R
#    BR                NDA      ID1      ID2
     9.85879668E-01    2     1000021         1   # BR(~d_R -> ~g d)
DECAY   1000001     1.11029737E+04   # ~d_L
#    BR                NDA      ID1      ID2
     1.44115126E-01    2    -1000024         2   # BR(~d_L -> chi^-_1 u)
     7.81030715E-01    2     1000021         1   # BR(~d_L -> ~g d)
DECAY   2000003     8.77618474E+03   # ~s_R
#    BR                NDA      ID1      ID2
     9.85879617E-01    2     1000021         3   # BR(~s_R -> ~g s)
DECAY   1000003     1.11029726E+04   # ~s_L
#    BR                NDA      ID1      ID2
     1.44115102E-01    2    -1000024         4   # BR(~s_L -> chi^-_1 c)
     7.81030739E-01    2     1000021         3   # BR(~s_L -> ~g s)
DECAY   1000005     8.77550863E+03   # ~b_1
#    BR                NDA      ID1      ID2
     9.85690320E-01    2     1000021         5   # BR(~b_1 -> ~g b)
DECAY   2000005     1.10943416E+04   # ~b_2
#    BR                NDA      ID1      ID2
     1.44022901E-01    2    -1000024         6   # BR(~b_2 -> chi^-_1 t)
     7.81099358E-01    2     1000021         5   # BR(~b_2 -> ~g b)
DECAY   2000002     9.14700714E+03   # ~u_R
#    BR                NDA      ID1      ID2
     9.45814041E-01    2     1000021         2   # BR(~u_R -> ~g u)
DECAY   1000002     1.11029729E+04   # ~u_L
#    BR                NDA      ID1      ID2
     1.44115160E-01    2     1000024         1   # BR(~u_L -> chi^+_1 d)
     7.81030679E-01    2     1000021         2   # BR(~u_L -> ~g u)
DECAY   2000004     9.14700700E+03   # ~c_R
#    BR                NDA      ID1      ID2
     9.45814040E-01    2     1000021         4   # BR(~c_R -> ~g c)
DECAY   1000004     1.11029723E+04   # ~c_L
#    BR                NDA      ID1      ID2
     1.44115159E-01    2     1000024         3   # BR(~c_L -> chi^+_1 s)
     7.81030672E-01    2     1000021         4   # BR(~c_L -> ~g c)
DECAY   1000006     9.26373184E+03   # ~t_1
#    BR                NDA      ID1      ID2
     9.32915030E-01    2     1000021         6   # BR(~t_1 -> ~g t)
#    BR                NDA      ID1      ID2       ID3
DECAY   2000006     1.10957504E+04   # ~t_2
#    BR                NDA      ID1      ID2
     1.44091415E-01    2     1000024         5   # BR(~t_2 -> chi^+_1 b)
     7.80998032E-01    2     1000021         6   # BR(~t_2 -> ~g t)
DECAY   1000024     1.75116912E-11   # chi^+_1
#    BR                NDA      ID1      ID2       ID3
     1.00000000E+00    3     1000022        -1         2   # BR(chi^+_1 -> chi^0_1 d_bar u)
DECAY   1000037     1.27343926E+03   # chi^+_2
#    BR                NDA      ID1      ID2
     2.80209606E-01    2     1000022        24   # BR(chi^+_2 -> chi^0_1 W^+)
     1.00762726E-01    2     1000023        24   # BR(chi^+_2 -> chi^0_2 W^+)
     2.94785314E-01    2     1000024        23   # BR(chi^+_2 -> chi^+_1 Z)
     3.11681857E-01    2     1000024        25   # BR(chi^+_2 -> chi^+_1 h^0)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000022     0.00000000E+00   # chi^0_1
DECAY   1000023     3.48790615E-10   # chi^0_2
#    BR                NDA      ID1      ID2
#    BR                NDA      ID1      ID2       ID3
     2.23150977E-01    3     1000024         1        -2   # BR(chi^0_2 -> chi^+_1 d u_bar)
     2.23150977E-01    3    -1000024        -1         2   # BR(chi^0_2 -> chi^-_1 d_bar u)
     1.48223375E-01    3     1000024        13       -14   # BR(chi^0_2 -> chi^+_1 mu^- nu_bar_mu)
     1.48223375E-01    3    -1000024       -13        14   # BR(chi^0_2 -> chi^-_1 mu^+ nu_mu)
     1.28429028E-01    3     1000024        15       -16   # BR(chi^0_2 -> chi^+_1 tau^- nu_bar_tau)
     1.28429028E-01    3    -1000024       -15        16   # BR(chi^0_2 -> chi^-_1 tau^+ nu_tau)
DECAY   1000025     1.25723778E+03   # chi^0_3
#    BR                NDA      ID1      ID2
     2.83695828E-01    2     1000024       -24   # BR(chi^0_3 -> chi^+_1 W^-)
     2.83695828E-01    2    -1000024        24   # BR(chi^0_3 -> chi^-_1 W^+)
     1.19443529E-01    2     1000022        23   # BR(chi^0_3 -> chi^0_1 Z)
     1.89160787E-01    2     1000022        25   # BR(chi^0_3 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000035     1.25342964E+03   # chi^0_4
#    BR                NDA      ID1      ID2
     2.84588192E-01    2     1000024       -24   # BR(chi^0_4 -> chi^+_1 W^-)
     2.84588192E-01    2    -1000024        24   # BR(chi^0_4 -> chi^-_1 W^+)
     1.79641731E-01    2     1000022        23   # BR(chi^0_4 -> chi^0_1 Z)
     1.26952575E-01    2     1000022        25   # BR(chi^0_4 -> chi^0_1 h^0)
#    BR                NDA      ID1      ID2       ID3
DECAY   1000021     1.27818479E-09   # ~g
#    BR                NDA      ID1      ID2
#    BR                NDA      ID1      ID2       ID3
DECAY        25     4.76289723E-03   # h^0
#    BR                NDA      ID1      ID2
     3.85029822E-01    2           5        -5   # BR(h^0 -> b b_bar)
# writing decays into V V* as 3-body decays
#    BR                NDA      ID1      ID2       ID3
     1.93890006E-02    3          24        11        12   # BR(h^0 -> W^+ e^- nu_e)
     1.93890006E-02    3          24        13        14   # BR(h^0 -> W^+ mu^- nu_mu)
     1.93890006E-02    3          24        15        16   # BR(h^0 -> W^+ tau^- nu_tau)
     6.78615022E-02    3          24         1        -2   # BR(h^0 -> W^+ d u_bar)
     6.78615022E-02    3          24         3        -4   # BR(h^0 -> W^+ s c_bar)
     1.93890006E-02    3         -24       -11       -12   # BR(h^0 -> W^- e^+ nu_bar_e)
     1.93890006E-02    3         -24       -13       -14   # BR(h^0 -> W^- mu^+ nu_bar_mu)
     1.93890006E-02    3         -24       -15       -16   # BR(h^0 -> W^- tau^+ nu_bar_tau)
     6.78615022E-02    3         -24        -1         2   # BR(h^0 -> W^- d_bar u)
     6.78615022E-02    3         -24        -3         4   # BR(h^0 -> W^- s_bar c)
     1.77376028E-03    3          23        11       -11   # BR(h^0 -> Z e^- e^+)
     1.77376028E-03    3          23        13       -13   # BR(h^0 -> Z mu^- mu^+)
     1.78431838E-03    3          23        15       -15   # BR(h^0 -> Z tau^- tau^+)
     1.05580969E-02    3          23        12       -12   # BR(h^0 -> Z nu_e nu_bar_e)
     8.23531559E-03    3          23         1        -1   # BR(h^0 -> Z d d_bar)
     8.23531559E-03    3          23         3        -3   # BR(h^0 -> Z s s_bar)
     7.97136316E-03    3          23         5        -5   # BR(h^0 -> Z b b_bar)
     6.12369621E-03    3          23         2        -2   # BR(h^0 -> Z u u_bar)
     6.33485814E-03    3          23         4        -4   # BR(h^0 -> Z c c_bar)
DECAY        37     2.53095774E+02   # H^+
#    BR                NDA      ID1      ID2
     1.63438375E-01    2         -15        12   # BR(H^+ -> tau^+ nu_e)
     8.35744414E-01    2          -5         6   # BR(H^+ -> b_bar t)
DECAY         6     2.43000000E+00   # top
#    BR                NDA      ID1      ID2
     1.00000000E+00    2           5        24   # BR(t -> b W)
Block HiggsBoundsInputHiggsCouplingsFermions
# ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP2
    1.00000131E+00    0.00000000E+00        3  25   5   5  # h0-b-b eff. coupling^2, normalised to SM
    9.99999987E+01    0.00000000E+00        3  35   5   5  # H0-b-b eff. coupling^2, normalised to SM
    0.00000000E+00    1.00000000E+02        3  36   5   5  # A0-b-b eff. coupling^2, normalised to SM
#
    9.99999987E-01    0.00000000E+00        3  25   6   6  # h0-t-t eff. coupling^2, normalised to SM
    1.00000131E-02    0.00000000E+00        3  35   6   6  # H0-t-t eff. coupling^2, normalised to SM
    0.00000000E+00    1.00000000E-02        3  36   6   6  # A0-t-t eff. coupling^2, normalised to SM
#
    1.00000131E+00    0.00000000E+00        3  25  15  15  # h0-tau-tau eff. coupling^2, normalised to SM
    9.99999987E+01    0.00000000E+00        3  35  15  15  # H0-tau-tau eff. coupling^2, normalised to SM
    0.00000000E+00    1.00000000E+02        3  36  15  15  # A0-tau-tau eff. coupling^2, normalised to SM
#
Block HiggsBoundsInputHiggsCouplingsBosons
    1.00000000E+00        3  25  24  24  # h0-W-W eff. coupling^2, normalised to SM
    4.27004274E-15        3  35  24  24  # H0-W-W eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  24  24  # A0-W-W eff. coupling^2, normalised to SM
#
    1.00000000E+00        3  25  23  23  # h0-Z-Z eff. coupling^2, normalised to SM
    4.27004274E-15        3  35  23  23  # H0-Z-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  23  23  # A0-Z-Z eff. coupling^2, normalised to SM
#
    1.13492380E+00        3  25  21  21  # h0-g-g eff. coupling^2, normalised to SM
    7.27813046E-03        3  35  21  21  # H0-g-g eff. coupling^2, normalised to SM
    1.14057109E-01        3  36  21  21  # A0-g-g eff. coupling^2, normalised to SM
#
    0.00000000E+00        3  25  25  23  # h0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  25  23  # H0-h0-Z eff. coupling^2, normalised to SM
    4.27004274E-15        3  36  25  23  # A0-h0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  35  35  23  # H0-H0-Z eff. coupling^2, normalised to SM
    1.00000000E+00        3  36  35  23  # A0-H0-Z eff. coupling^2, normalised to SM
    0.00000000E+00        3  36  36  23  # A0-A0-Z eff. coupling^2, normalised to SM
#
    0.00000000E+00        4  25  21  21  23  # h0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  35  21  21  23  # H0-g-g-Z eff. coupling^2, normalised to SM
    0.00000000E+00        4  36  21  21  23  # A0-g-g-Z eff. coupling^2, normalised to SM
Block SPhenoLowEnergy  # low energy observables
    1    3.20002504E-04   # BR(b -> s gamma)
    2    1.59099767E-06   # BR(b -> s mu+ mu-)
    3    3.51767155E-05   # BR(b -> s nu nu)
    4    2.04895882E-15   # BR(Bd -> e+ e-)
    5    8.75291929E-11   # BR(Bd -> mu+ mu-)
    6    1.83211533E-08   # BR(Bd -> tau+ tau-)
    7    6.90724324E-14   # BR(Bs -> e+ e-)
    8    2.95077207E-09   # BR(Bs -> mu+ mu-)
    9    6.25806221E-07   # BR(Bs -> tau+ tau-)
   10    9.67944490E-05   # BR(B_u -> tau nu)
   11    9.99999906E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
   12    5.41617311E-01   # |Delta(M_Bd)| [ps^-1] 
   13    1.93526358E+01   # |Delta(M_Bs)| [ps^-1] 
   16    2.15524495E-03   # epsilon_K
   17    2.28163801E-15   # Delta(M_K)
   18    2.47454081E-11   # BR(K^0 -> pi^0 nu nu)
   19    8.27754101E-11   # BR(K^+ -> pi^+ nu nu)
   20   -2.20818235E-21   # Delta(g-2)_electron/2
   21   -9.44101997E-17   # Delta(g-2)_muon/2
   22   -2.70376159E-14   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39   -2.34622199E-04   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.88681561E-01   # C7
     0305 4422   00   2    -1.88683977E-01   # C7
     0305 4322   00   2    -4.98902010E-08   # C7'
     0305 6421   00   0    -9.51840960E-02   # C8
     0305 6421   00   2    -9.51875875E-02   # C8
     0305 6321   00   2    -7.20248545E-08   # C8'
 03051111 4133   00   0     1.16528672E+00   # C9 e+e-
 03051111 4133   00   2     1.16528681E+00   # C9 e+e-
 03051111 4233   00   2     1.12223287E-06   # C9' e+e-
 03051111 4137   00   0    -3.98104090E+00   # C10 e+e-
 03051111 4137   00   2    -3.98104084E+00   # C10 e+e-
 03051111 4237   00   2    -3.20215190E-05   # C10' e+e-
 03051313 4133   00   0     1.16528672E+00   # C9 mu+mu-
 03051313 4133   00   2     1.16528681E+00   # C9 mu+mu-
 03051313 4233   00   2     1.12223287E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.98104090E+00   # C10 mu+mu-
 03051313 4137   00   2    -3.98104084E+00   # C10 mu+mu-
 03051313 4237   00   2    -3.20215190E-05   # C10' mu+mu-
 03051212 4137   00   0     1.50336814E+00   # C11 nu_1 nu_1
 03051212 4137   00   2     1.50336813E+00   # C11 nu_1 nu_1
 03051212 4237   00   2     7.72482163E-06   # C11' nu_1 nu_1
 03051414 4137   00   0     1.50336814E+00   # C11 nu_2 nu_2
 03051414 4137   00   2     1.50336813E+00   # C11 nu_2 nu_2
 03051414 4237   00   2     7.72482163E-06   # C11' nu_2 nu_2
 03051616 4137   00   0     1.50336814E+00   # C11 nu_3 nu_3
 03051616 4137   00   2     1.50336813E+00   # C11 nu_3 nu_3
 03051616 4237   00   2     7.72482163E-06   # C11' nu_3 nu_3
Block IMFWCOEF Q=  1.60000000E+02  # Im(Wilson coefficients) at scale Q
#    id        order  M        value         comment
     0305 4422   00   0     3.84926274E-07   # C7
     0305 4422   00   2     3.84939536E-07   # C7
     0305 4322   00   2     4.90928765E-13   # C7'
     0305 6421   00   0     3.29715255E-07   # C8
     0305 6421   00   2     3.29252071E-07   # C8
     0305 6321   00   2    -1.37955468E-11   # C8'
 03051111 4133   00   2     1.06421054E-09   # C9 e+e-
 03051111 4233   00   2     2.00344832E-08   # C9' e+e-
 03051111 4137   00   2     9.58902251E-11   # C10 e+e-
 03051111 4237   00   2    -5.71659056E-07   # C10' e+e-
 03051313 4133   00   2     1.06421017E-09   # C9 mu+mu-
 03051313 4233   00   2     2.00344832E-08   # C9' mu+mu-
 03051313 4137   00   2     9.58906147E-11   # C10 mu+mu-
 03051313 4237   00   2    -5.71659056E-07   # C10' mu+mu-
 03051212 4137   00   2     1.14218270E-12   # C11 nu_1 nu_1
 03051212 4237   00   2     1.37906145E-07   # C11' nu_1 nu_1
 03051414 4137   00   2     1.14220602E-12   # C11 nu_2 nu_2
 03051414 4237   00   2     1.37906145E-07   # C11' nu_2 nu_2
 03051616 4137   00   2     1.14876491E-12   # C11 nu_3 nu_3
 03051616 4237   00   2     1.37906145E-07   # C11' nu_3 nu_3
