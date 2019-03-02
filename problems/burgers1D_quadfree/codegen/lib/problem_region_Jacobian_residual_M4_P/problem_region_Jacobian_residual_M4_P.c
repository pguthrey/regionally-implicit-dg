/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * problem_region_Jacobian_residual_M4_P.c
 *
 * Code generation for function 'problem_region_Jacobian_residual_M4_P'
 *
 */

/* Include files */
#include <math.h>
#include <string.h>
#include "problem_region_Jacobian_residual_M4_P.h"

/* Function Definitions */
void problem_region_Jacobian_residual_M4_P(const double q[30], const double
  q_past[30], double nuv1, double Jacobian[900], double residual[30])
{
  double Uphi1_east_3;
  double Uphi2_east_3;
  double Uphi3_east_3;
  double Uphi1_east_2;
  double Uphi2_east_2;
  double Uphi3_east_2;
  double Uphi1_east_1;
  double Uphi2_east_1;
  double Uphi3_east_1;
  double Uphi1_west_3;
  double Uphi2_west_3;
  double Uphi3_west_3;
  double Uphi1_west_2;
  double Uphi2_west_2;
  double Uphi3_west_2;
  double Uphi1_west_1;
  double Uphi2_west_1;
  double Uphi3_west_1;
  double lambda_sample_12_1;
  double lambda_sample_12_2;
  double lambda_sample_12_3;
  double lambda_sample_12_4;
  double lambda_sample_23_1;
  double lambda_sample_23_2;
  double lambda_sample_23_3;
  double lambda_sample_23_4;
  double lambda_12_1;
  double lambda_12_2;
  double lambda_12_3;
  double lambda_12_4;
  double lambda_23_1;
  double lambda_23_2;
  double lambda_23_3;
  double lambda_23_4;
  double residual_tmp;
  double b_residual_tmp;
  double c_residual_tmp;
  double d_residual_tmp;
  double e_residual_tmp;
  double f_residual_tmp;
  double g_residual_tmp;
  double h_residual_tmp;
  double i_residual_tmp;
  double j_residual_tmp;
  double k_residual_tmp;
  double l_residual_tmp;
  double m_residual_tmp;
  double n_residual_tmp;
  double o_residual_tmp;
  double p_residual_tmp;
  double q_residual_tmp;
  double r_residual_tmp;
  double s_residual_tmp;
  double t_residual_tmp;
  double u_residual_tmp;
  double v_residual_tmp;
  double w_residual_tmp;
  double x_residual_tmp;
  double y_residual_tmp;
  double ab_residual_tmp;
  double bb_residual_tmp;
  double cb_residual_tmp;
  double db_residual_tmp;
  double eb_residual_tmp;
  double fb_residual_tmp;
  double gb_residual_tmp;
  double hb_residual_tmp;
  double ib_residual_tmp;
  double jb_residual_tmp;
  double kb_residual_tmp;
  double lb_residual_tmp;
  double mb_residual_tmp;
  double nb_residual_tmp;
  double ob_residual_tmp;
  double pb_residual_tmp;
  double qb_residual_tmp;
  double rb_residual_tmp;
  double sb_residual_tmp;
  double tb_residual_tmp;
  double ub_residual_tmp;
  double vb_residual_tmp;
  double wb_residual_tmp;
  double xb_residual_tmp;
  double yb_residual_tmp;
  double ac_residual_tmp;
  double bc_residual_tmp;
  double cc_residual_tmp;
  double dc_residual_tmp;
  double ec_residual_tmp;
  double fc_residual_tmp;
  double gc_residual_tmp;
  double hc_residual_tmp;
  double ic_residual_tmp;
  double jc_residual_tmp;
  double kc_residual_tmp;
  double lc_residual_tmp;
  double mc_residual_tmp;
  double nc_residual_tmp;
  double oc_residual_tmp;
  double pc_residual_tmp;
  double qc_residual_tmp;
  double rc_residual_tmp;
  double sc_residual_tmp;
  double tc_residual_tmp;
  double uc_residual_tmp;
  double vc_residual_tmp;
  double wc_residual_tmp;
  double xc_residual_tmp;
  double yc_residual_tmp;
  double ad_residual_tmp;
  double bd_residual_tmp;
  double cd_residual_tmp;
  double dd_residual_tmp;
  double ed_residual_tmp;
  double fd_residual_tmp;
  double gd_residual_tmp;
  double hd_residual_tmp;
  double id_residual_tmp;
  double jd_residual_tmp;
  double kd_residual_tmp;
  double ld_residual_tmp;
  double md_residual_tmp;
  double nd_residual_tmp;
  double od_residual_tmp;
  double pd_residual_tmp;
  double qd_residual_tmp;
  double rd_residual_tmp;
  double sd_residual_tmp;
  double td_residual_tmp;
  double ud_residual_tmp;
  double vd_residual_tmp;
  double wd_residual_tmp;
  double xd_residual_tmp;
  double yd_residual_tmp;
  double ae_residual_tmp;
  double be_residual_tmp;
  double ce_residual_tmp;
  double de_residual_tmp;
  double ee_residual_tmp;
  double fe_residual_tmp;
  double ge_residual_tmp;
  double he_residual_tmp;
  double ie_residual_tmp;
  double je_residual_tmp;
  double ke_residual_tmp;
  double le_residual_tmp;
  double me_residual_tmp;
  double ne_residual_tmp;
  double oe_residual_tmp;
  double pe_residual_tmp;
  double qe_residual_tmp;
  double re_residual_tmp;
  double se_residual_tmp;
  double te_residual_tmp;
  double ue_residual_tmp;
  double ve_residual_tmp;
  double we_residual_tmp;
  double xe_residual_tmp;
  double ye_residual_tmp;
  double af_residual_tmp;
  double bf_residual_tmp;
  double cf_residual_tmp;
  double df_residual_tmp;
  double ef_residual_tmp;
  double ff_residual_tmp;
  double gf_residual_tmp;
  double hf_residual_tmp;
  double if_residual_tmp;
  double jf_residual_tmp;
  double kf_residual_tmp;
  double lf_residual_tmp;
  double mf_residual_tmp;
  double nf_residual_tmp;
  double of_residual_tmp;
  double pf_residual_tmp;
  double qf_residual_tmp;
  double rf_residual_tmp;
  double sf_residual_tmp;
  double tf_residual_tmp;
  double uf_residual_tmp;
  double vf_residual_tmp;
  double wf_residual_tmp;
  double xf_residual_tmp;
  double yf_residual_tmp;
  double ag_residual_tmp;
  double bg_residual_tmp;
  double cg_residual_tmp;
  double dg_residual_tmp;
  double eg_residual_tmp;
  double fg_residual_tmp;
  double gg_residual_tmp;
  double hg_residual_tmp;
  double ig_residual_tmp;
  double jg_residual_tmp;
  double kg_residual_tmp;
  double lg_residual_tmp;
  double mg_residual_tmp;
  double ng_residual_tmp;
  double og_residual_tmp;
  double pg_residual_tmp;
  double qg_residual_tmp;
  double rg_residual_tmp;
  double sg_residual_tmp;
  double tg_residual_tmp;
  double ug_residual_tmp;
  double vg_residual_tmp;
  double wg_residual_tmp;
  double xg_residual_tmp;
  double yg_residual_tmp;
  double ah_residual_tmp;
  double bh_residual_tmp;
  double ch_residual_tmp;
  double dh_residual_tmp;
  double eh_residual_tmp;
  double fh_residual_tmp;
  double gh_residual_tmp;
  double hh_residual_tmp;
  double ih_residual_tmp;
  double jh_residual_tmp;
  double kh_residual_tmp;
  double lh_residual_tmp;
  double mh_residual_tmp;
  double nh_residual_tmp;
  double oh_residual_tmp;
  double ph_residual_tmp;
  double qh_residual_tmp;
  double rh_residual_tmp;
  double d0;
  double d1;
  double d2;
  double d3;
  double sh_residual_tmp;
  double th_residual_tmp;
  double uh_residual_tmp;
  double vh_residual_tmp;
  double wh_residual_tmp;
  double xh_residual_tmp;
  double yh_residual_tmp;
  double d4;
  double d5;
  Uphi1_east_3 = q[7] + 1.7320508075688772 * q[8];
  Uphi2_east_3 = q[17] + 1.7320508075688772 * q[18];
  Uphi3_east_3 = q[27] + 1.7320508075688772 * q[28];
  Uphi1_east_2 = (q[4] + 1.7320508075688772 * q[5]) + 2.23606797749979 * q[6];
  Uphi2_east_2 = (q[14] + 1.7320508075688772 * q[15]) + 2.23606797749979 * q[16];
  Uphi3_east_2 = (q[24] + 1.7320508075688772 * q[25]) + 2.23606797749979 * q[26];
  Uphi1_east_1 = ((q[0] + 1.7320508075688772 * q[1]) + 2.23606797749979 * q[2])
    + 2.6457513110645907 * q[3];
  Uphi2_east_1 = ((q[10] + 1.7320508075688772 * q[11]) + 2.23606797749979 * q[12])
    + 2.6457513110645907 * q[13];
  Uphi3_east_1 = ((q[20] + 1.7320508075688772 * q[21]) + 2.23606797749979 * q[22])
    + 2.6457513110645907 * q[23];
  Uphi1_west_3 = q[7] - 1.7320508075688772 * q[8];
  Uphi2_west_3 = q[17] - 1.7320508075688772 * q[18];
  Uphi3_west_3 = q[27] - 1.7320508075688772 * q[28];
  Uphi1_west_2 = (q[4] - 1.7320508075688772 * q[5]) + 2.23606797749979 * q[6];
  Uphi2_west_2 = (q[14] - 1.7320508075688772 * q[15]) + 2.23606797749979 * q[16];
  Uphi3_west_2 = (q[24] - 1.7320508075688772 * q[25]) + 2.23606797749979 * q[26];
  Uphi1_west_1 = ((q[0] - 1.7320508075688772 * q[1]) + 2.23606797749979 * q[2])
    - 2.6457513110645907 * q[3];
  Uphi2_west_1 = ((q[10] - 1.7320508075688772 * q[11]) + 2.23606797749979 * q[12])
    - 2.6457513110645907 * q[13];
  Uphi3_west_1 = ((q[20] - 1.7320508075688772 * q[21]) + 2.23606797749979 * q[22])
    - 2.6457513110645907 * q[23];
  lambda_sample_12_1 = fmax(fabs(((Uphi1_east_1 - 1.4915318439233631 *
    Uphi1_east_2) + 1.3692196008356177 * Uphi1_east_3) - 0.80628473498821784 *
    q[9]), fabs(((Uphi2_west_1 - 1.4915318439233631 * Uphi2_west_2) +
                 1.3692196008356177 * Uphi2_west_3) - 0.80628473498821784 * q[19]));
  lambda_sample_12_2 = fmax(fabs(((Uphi1_east_1 - 0.58886444109926006 *
    Uphi1_east_2) - 0.73034303583567772 * Uphi1_east_3) + 1.0893298949365755 *
    q[9]), fabs(((Uphi2_west_1 - 0.58886444109926006 * Uphi2_west_2) -
                 0.73034303583567772 * Uphi2_west_3) + 1.0893298949365755 * q[19]));
  lambda_sample_12_3 = fmax(fabs(((Uphi1_east_1 + 0.58886444109926006 *
    Uphi1_east_2) - 0.73034303583567772 * Uphi1_east_3) - 1.0893298949365755 *
    q[9]), fabs(((Uphi2_west_1 + 0.58886444109926006 * Uphi2_west_2) -
                 0.73034303583567772 * Uphi2_west_3) - 1.0893298949365755 * q[19]));
  lambda_sample_12_4 = fmax(fabs(((Uphi1_east_1 + 1.4915318439233631 *
    Uphi1_east_2) + 1.3692196008356177 * Uphi1_east_3) + 0.80628473498821784 *
    q[9]), fabs(((Uphi2_west_1 + 1.4915318439233631 * Uphi2_west_2) +
                 1.3692196008356177 * Uphi2_west_3) + 0.80628473498821784 * q[19]));
  lambda_sample_23_1 = fmax(fabs(((Uphi2_east_1 - 1.4915318439233631 *
    Uphi2_east_2) + 1.3692196008356177 * Uphi2_east_3) - 0.80628473498821784 *
    q[19]), fabs(((Uphi3_west_1 - 1.4915318439233631 * Uphi3_west_2) +
                  1.3692196008356177 * Uphi3_west_3) - 0.80628473498821784 * q
                 [29]));
  lambda_sample_23_2 = fmax(fabs(((Uphi2_east_1 - 0.58886444109926006 *
    Uphi2_east_2) - 0.73034303583567772 * Uphi2_east_3) + 1.0893298949365755 *
    q[19]), fabs(((Uphi3_west_1 - 0.58886444109926006 * Uphi3_west_2) -
                  0.73034303583567772 * Uphi3_west_3) + 1.0893298949365755 * q
                 [29]));
  lambda_sample_23_3 = fmax(fabs(((Uphi2_east_1 + 0.58886444109926006 *
    Uphi2_east_2) - 0.73034303583567772 * Uphi2_east_3) - 1.0893298949365755 *
    q[19]), fabs(((Uphi3_west_1 + 0.58886444109926006 * Uphi3_west_2) -
                  0.73034303583567772 * Uphi3_west_3) - 1.0893298949365755 * q
                 [29]));
  lambda_sample_23_4 = fmax(fabs(((Uphi2_east_1 + 1.4915318439233631 *
    Uphi2_east_2) + 1.3692196008356177 * Uphi2_east_3) + 0.80628473498821784 *
    q[19]), fabs(((Uphi3_west_1 + 1.4915318439233631 * Uphi3_west_2) +
                  1.3692196008356177 * Uphi3_west_3) + 0.80628473498821784 * q
                 [29]));
  lambda_12_1 = ((0.1739274225687269 * lambda_sample_12_1 + 0.32607257743127305 *
                  lambda_sample_12_2) + 0.32607257743127305 * lambda_sample_12_3)
    + 0.1739274225687269 * lambda_sample_12_4;
  lambda_12_2 = ((0.1920125460668618 * lambda_sample_12_3 - 0.1920125460668618 *
                  lambda_sample_12_2) - 0.25941828929277116 * lambda_sample_12_1)
    + 0.25941828929277116 * lambda_sample_12_4;
  lambda_12_3 = ((0.23814483610392004 * lambda_sample_12_1 - 0.23814483610392004
                  * lambda_sample_12_2) - 0.23814483610392004 *
                 lambda_sample_12_3) + 0.23814483610392004 * lambda_sample_12_4;
  lambda_12_4 = ((0.35520060651490704 * lambda_sample_12_2 - 0.14023502581300973
                  * lambda_sample_12_1) - 0.35520060651490704 *
                 lambda_sample_12_3) + 0.14023502581300973 * lambda_sample_12_4;
  lambda_23_1 = ((0.1739274225687269 * lambda_sample_23_1 + 0.32607257743127305 *
                  lambda_sample_23_2) + 0.32607257743127305 * lambda_sample_23_3)
    + 0.1739274225687269 * lambda_sample_23_4;
  lambda_23_2 = ((0.1920125460668618 * lambda_sample_23_3 - 0.1920125460668618 *
                  lambda_sample_23_2) - 0.25941828929277116 * lambda_sample_23_1)
    + 0.25941828929277116 * lambda_sample_23_4;
  lambda_23_3 = ((0.23814483610392004 * lambda_sample_23_1 - 0.23814483610392004
                  * lambda_sample_23_2) - 0.23814483610392004 *
                 lambda_sample_23_3) + 0.23814483610392004 * lambda_sample_23_4;
  lambda_23_4 = ((0.35520060651490704 * lambda_sample_23_2 - 0.14023502581300973
                  * lambda_sample_23_1) - 0.35520060651490704 *
                 lambda_sample_23_3) + 0.14023502581300973 * lambda_sample_23_4;
  lambda_sample_12_1 = Uphi2_east_1 * q[19];
  lambda_sample_12_2 = 0.87831006565367986 * Uphi2_east_2 * Uphi2_east_3;
  lambda_sample_12_3 = 0.59628479399994394 * Uphi2_east_3 * q[19];
  lambda_sample_12_4 = Uphi3_west_1 * q[29];
  lambda_sample_23_1 = 0.87831006565367986 * Uphi3_west_2 * Uphi3_west_3;
  lambda_sample_23_2 = 0.59628479399994394 * Uphi3_west_3 * q[29];
  lambda_sample_23_3 = Uphi2_east_1 * lambda_23_4;
  lambda_sample_23_4 = 0.87831006565367986 * Uphi2_east_2 * lambda_23_3;
  residual_tmp = 0.87831006565367986 * Uphi2_east_3 * lambda_23_2;
  b_residual_tmp = q[19] * lambda_23_1;
  c_residual_tmp = 0.59628479399994394 * Uphi2_east_3 * lambda_23_4;
  d_residual_tmp = 0.59628479399994394 * q[19] * lambda_23_3;
  e_residual_tmp = Uphi3_west_1 * lambda_23_4;
  f_residual_tmp = 0.87831006565367986 * Uphi3_west_2 * lambda_23_3;
  g_residual_tmp = 0.87831006565367986 * Uphi3_west_3 * lambda_23_2;
  h_residual_tmp = q[29] * lambda_23_1;
  i_residual_tmp = 0.59628479399994394 * Uphi3_west_3 * lambda_23_4;
  j_residual_tmp = 0.59628479399994394 * q[29] * lambda_23_3;
  residual[29] = (((((((9.16515138991168 * q[24] - 5.2915026221291814 * q[20]) -
                       11.832159566199232 * q[27]) + 14.0 * q[29]) +
                     5.2915026221291814 * q_past[20]) + 9.16515138991168 *
                    q_past[24]) + 11.832159566199232 * q_past[27]) + 14.0 *
                  q_past[29]) + nuv1 * ((((((((((((((((((((2.0 * Uphi3_east_1 *
    q[29] - lambda_sample_12_2) - lambda_sample_12_3) - lambda_sample_12_1) +
    1.7566201313073597 * Uphi3_east_2 * Uphi3_east_3) + 1.1925695879998879 *
    Uphi3_east_3 * q[29]) - lambda_sample_12_4) - lambda_sample_23_1) -
    lambda_sample_23_2) - lambda_sample_23_3) - lambda_sample_23_4) -
    residual_tmp) - b_residual_tmp) - c_residual_tmp) - d_residual_tmp) +
    e_residual_tmp) + f_residual_tmp) + g_residual_tmp) + h_residual_tmp) +
    i_residual_tmp) + j_residual_tmp);
  k_residual_tmp = q[29] * q[29];
  l_residual_tmp = Uphi2_east_2 * Uphi2_east_2;
  m_residual_tmp = Uphi2_east_3 * Uphi2_east_3;
  n_residual_tmp = q[19] * q[19];
  o_residual_tmp = Uphi3_east_2 * Uphi3_east_2;
  p_residual_tmp = Uphi3_east_3 * Uphi3_east_3;
  q_residual_tmp = Uphi3_west_2 * Uphi3_west_2;
  r_residual_tmp = Uphi3_west_3 * Uphi3_west_3;
  s_residual_tmp = q[24] * q[24];
  t_residual_tmp = q[25] * q[25];
  u_residual_tmp = q[26] * q[26];
  v_residual_tmp = q[27] * q[27];
  w_residual_tmp = q[28] * q[28];
  x_residual_tmp = 1.7320508075688772 * Uphi2_east_1 * Uphi2_east_3;
  y_residual_tmp = 1.52127765851133 * Uphi2_east_2 * q[19];
  ab_residual_tmp = 1.7320508075688772 * Uphi3_west_1 * Uphi3_west_3;
  bb_residual_tmp = 1.52127765851133 * Uphi3_west_2 * q[29];
  cb_residual_tmp = 1.7320508075688772 * Uphi2_east_1 * lambda_23_3;
  db_residual_tmp = 1.5491933384829668 * Uphi2_east_2 * lambda_23_2;
  eb_residual_tmp = 1.7320508075688772 * Uphi2_east_3 * lambda_23_1;
  fb_residual_tmp = 1.52127765851133 * Uphi2_east_2 * lambda_23_4;
  gb_residual_tmp = 1.1065666703449764 * Uphi2_east_3 * lambda_23_3;
  hb_residual_tmp = 1.52127765851133 * q[19] * lambda_23_2;
  ib_residual_tmp = 1.0327955589886444 * q[19] * lambda_23_4;
  jb_residual_tmp = 1.7320508075688772 * Uphi3_west_1 * lambda_23_3;
  kb_residual_tmp = 1.5491933384829668 * Uphi3_west_2 * lambda_23_2;
  lb_residual_tmp = 1.7320508075688772 * Uphi3_west_3 * lambda_23_1;
  mb_residual_tmp = 1.52127765851133 * Uphi3_west_2 * lambda_23_4;
  nb_residual_tmp = 1.1065666703449764 * Uphi3_west_3 * lambda_23_3;
  ob_residual_tmp = 1.52127765851133 * q[29] * lambda_23_2;
  pb_residual_tmp = 1.0327955589886444 * q[29] * lambda_23_4;
  qb_residual_tmp = 0.7745966692414834 * l_residual_tmp;
  rb_residual_tmp = 0.55328333517248818 * m_residual_tmp;
  sb_residual_tmp = 0.5163977794943222 * n_residual_tmp;
  tb_residual_tmp = 0.7745966692414834 * q_residual_tmp;
  ub_residual_tmp = 0.55328333517248818 * r_residual_tmp;
  vb_residual_tmp = 0.5163977794943222 * k_residual_tmp;
  residual[28] = (((((4.47213595499958 * q[21] - 7.745966692414834 * q[25]) +
                     10.0 * q[28]) - 4.47213595499958 * q_past[21]) -
                   7.745966692414834 * q_past[25]) - 10.0 * q_past[28]) + nuv1 *
    (((((((((((((((((((((((((((((((((((((x_residual_tmp + y_residual_tmp) +
    3.4641016151377544 * Uphi3_east_1 * Uphi3_east_3) + 3.04255531702266 *
    Uphi3_east_2 * q[29]) + ab_residual_tmp) + bb_residual_tmp) -
    6.9282032302755088 * q[20] * q[27]) - 6.9282032302755088 * q[21] * q[28]) -
    6.08511063404532 * q[24] * q[29]) + cb_residual_tmp) + db_residual_tmp) +
    eb_residual_tmp) + fb_residual_tmp) + gb_residual_tmp) + hb_residual_tmp) +
    ib_residual_tmp) - jb_residual_tmp) - kb_residual_tmp) - lb_residual_tmp) -
                       mb_residual_tmp) - nb_residual_tmp) - ob_residual_tmp) -
                    pb_residual_tmp) + qb_residual_tmp) + rb_residual_tmp) +
                 sb_residual_tmp) + 1.5491933384829668 * o_residual_tmp) +
               1.1065666703449764 * p_residual_tmp) + 1.0327955589886444 *
              k_residual_tmp) + tb_residual_tmp) + ub_residual_tmp) +
           vb_residual_tmp) - 3.0983866769659336 * s_residual_tmp) -
         3.0983866769659336 * t_residual_tmp) - 3.0983866769659336 *
        u_residual_tmp) - 2.2131333406899527 * v_residual_tmp) -
      2.2131333406899527 * w_residual_tmp) - 2.0655911179772892 * k_residual_tmp);
  wb_residual_tmp = Uphi2_east_1 * Uphi2_east_3;
  xb_residual_tmp = Uphi3_west_1 * Uphi3_west_3;
  yb_residual_tmp = 0.87831006565367986 * Uphi3_west_2 * q[29];
  ac_residual_tmp = Uphi2_east_1 * lambda_23_3;
  bc_residual_tmp = 0.89442719099991586 * Uphi2_east_2 * lambda_23_2;
  cc_residual_tmp = Uphi2_east_3 * lambda_23_1;
  dc_residual_tmp = 0.87831006565367986 * Uphi2_east_2 * lambda_23_4;
  ec_residual_tmp = 0.63887656499993994 * Uphi2_east_3 * lambda_23_3;
  fc_residual_tmp = 0.87831006565367986 * q[19] * lambda_23_2;
  gc_residual_tmp = 0.59628479399994394 * q[19] * lambda_23_4;
  hc_residual_tmp = Uphi3_west_1 * lambda_23_3;
  ic_residual_tmp = 0.89442719099991586 * Uphi3_west_2 * lambda_23_2;
  jc_residual_tmp = Uphi3_west_3 * lambda_23_1;
  kc_residual_tmp = 0.87831006565367986 * Uphi3_west_2 * lambda_23_4;
  lc_residual_tmp = 0.63887656499993994 * Uphi3_west_3 * lambda_23_3;
  mc_residual_tmp = 0.87831006565367986 * q[29] * lambda_23_2;
  nc_residual_tmp = 0.59628479399994394 * q[29] * lambda_23_4;
  oc_residual_tmp = 0.44721359549995793 * l_residual_tmp;
  pc_residual_tmp = 0.31943828249996997 * m_residual_tmp;
  qc_residual_tmp = 0.29814239699997197 * n_residual_tmp;
  rc_residual_tmp = 0.44721359549995793 * q_residual_tmp;
  sc_residual_tmp = 0.31943828249996997 * r_residual_tmp;
  tc_residual_tmp = 0.29814239699997197 * k_residual_tmp;
  residual[27] = (((((((4.47213595499958 * q[20] - 7.745966692414834 * q[24]) +
                       10.0 * q[27]) + 11.832159566199232 * q[29]) -
                     4.47213595499958 * q_past[20]) - 7.745966692414834 *
                    q_past[24]) - 10.0 * q_past[27]) - 11.832159566199232 *
                  q_past[29]) + nuv1 * ((((((((((((((((((((((((((((2.0 *
    Uphi3_east_1 * Uphi3_east_3 - 0.87831006565367986 * Uphi2_east_2 * q[19]) -
    wb_residual_tmp) + 1.7566201313073597 * Uphi3_east_2 * q[29]) -
    xb_residual_tmp) - yb_residual_tmp) - ac_residual_tmp) - bc_residual_tmp) -
    cc_residual_tmp) - dc_residual_tmp) - ec_residual_tmp) - fc_residual_tmp) -
    gc_residual_tmp) + hc_residual_tmp) + ic_residual_tmp) + jc_residual_tmp) +
    kc_residual_tmp) + lc_residual_tmp) + mc_residual_tmp) + nc_residual_tmp) -
    oc_residual_tmp) - pc_residual_tmp) - qc_residual_tmp) + 0.89442719099991586
    * o_residual_tmp) + 0.63887656499993994 * p_residual_tmp) +
    0.59628479399994394 * k_residual_tmp) - rc_residual_tmp) - sc_residual_tmp)
    - tc_residual_tmp);
  uc_residual_tmp = 2.0 * Uphi2_east_2 * Uphi2_east_3;
  vc_residual_tmp = 1.9639610121239315 * Uphi2_east_3 * q[19];
  wc_residual_tmp = 2.23606797749979 * Uphi3_west_1 * Uphi3_west_2;
  xc_residual_tmp = 2.0 * Uphi3_west_2 * Uphi3_west_3;
  yc_residual_tmp = 1.9639610121239315 * Uphi3_west_3 * q[29];
  ad_residual_tmp = 2.23606797749979 * Uphi2_east_1 * lambda_23_2;
  bd_residual_tmp = 2.23606797749979 * Uphi2_east_2 * lambda_23_1;
  cd_residual_tmp = 2.0 * Uphi2_east_2 * lambda_23_3;
  dd_residual_tmp = 2.0 * Uphi2_east_3 * lambda_23_2;
  ed_residual_tmp = 1.9639610121239315 * Uphi2_east_3 * lambda_23_4;
  fd_residual_tmp = 1.9639610121239315 * q[19] * lambda_23_3;
  gd_residual_tmp = 2.23606797749979 * Uphi3_west_1 * lambda_23_2;
  hd_residual_tmp = 2.23606797749979 * Uphi3_west_2 * lambda_23_1;
  id_residual_tmp = 2.0 * Uphi3_west_2 * lambda_23_3;
  jd_residual_tmp = 2.0 * Uphi3_west_3 * lambda_23_2;
  kd_residual_tmp = 1.9639610121239315 * Uphi3_west_3 * lambda_23_4;
  ld_residual_tmp = 1.9639610121239315 * q[29] * lambda_23_3;
  residual[26] = (((6.0 * q[26] - 3.4641016151377544 * q[22]) +
                   3.4641016151377544 * q_past[22]) + 6.0 * q_past[26]) + nuv1 *
    (((((((((((((((((((((((((((((4.47213595499958 * Uphi3_east_1 * Uphi3_east_2
    - uc_residual_tmp) - vc_residual_tmp) - 2.23606797749979 * Uphi2_east_1 *
    Uphi2_east_2) + 4.0 * Uphi3_east_2 * Uphi3_east_3) + 3.927922024247863 *
    Uphi3_east_3 * q[29]) - wc_residual_tmp) - xc_residual_tmp) -
    yc_residual_tmp) - 15.491933384829668 * q[20] * q[25]) - 15.491933384829668 *
    q[21] * q[24]) - 13.856406460551018 * q[21] * q[26]) - 13.856406460551018 *
                      q[22] * q[25]) - 13.60672102833218 * q[23] * q[26]) -
                    13.856406460551018 * q[24] * q[28]) - 13.856406460551018 *
                   q[25] * q[27]) - 12.393546707863734 * q[26] * q[28]) -
                 13.60672102833218 * q[28] * q[29]) - ad_residual_tmp) -
               bd_residual_tmp) - cd_residual_tmp) - dd_residual_tmp) -
            ed_residual_tmp) - fd_residual_tmp) + gd_residual_tmp) +
         hd_residual_tmp) + id_residual_tmp) + jd_residual_tmp) +
      kd_residual_tmp) + ld_residual_tmp);
  md_residual_tmp = 1.7320508075688772 * Uphi2_east_1 * Uphi2_east_2;
  nd_residual_tmp = 1.5491933384829668 * Uphi2_east_2 * Uphi2_east_3;
  od_residual_tmp = 1.52127765851133 * Uphi2_east_3 * q[19];
  pd_residual_tmp = 1.7320508075688772 * Uphi3_west_1 * Uphi3_west_2;
  qd_residual_tmp = 1.5491933384829668 * Uphi3_west_2 * Uphi3_west_3;
  rd_residual_tmp = 1.52127765851133 * Uphi3_west_3 * q[29];
  sd_residual_tmp = 1.7320508075688772 * Uphi2_east_1 * lambda_23_2;
  td_residual_tmp = 1.7320508075688772 * Uphi2_east_2 * lambda_23_1;
  ud_residual_tmp = 1.5491933384829668 * Uphi2_east_2 * lambda_23_3;
  vd_residual_tmp = 1.5491933384829668 * Uphi2_east_3 * lambda_23_2;
  wd_residual_tmp = 1.52127765851133 * Uphi2_east_3 * lambda_23_4;
  xd_residual_tmp = 1.52127765851133 * q[19] * lambda_23_3;
  yd_residual_tmp = 1.7320508075688772 * Uphi3_west_1 * lambda_23_2;
  ae_residual_tmp = 1.7320508075688772 * Uphi3_west_2 * lambda_23_1;
  be_residual_tmp = 1.5491933384829668 * Uphi3_west_2 * lambda_23_3;
  ce_residual_tmp = 1.5491933384829668 * Uphi3_west_3 * lambda_23_2;
  de_residual_tmp = 1.52127765851133 * Uphi3_west_3 * lambda_23_4;
  ee_residual_tmp = 1.52127765851133 * q[29] * lambda_23_3;
  residual[25] = (((((6.0 * q[25] - 3.4641016151377544 * q[21]) +
                     7.745966692414834 * q[28]) + 3.4641016151377544 * q_past[21])
                   + 6.0 * q_past[25]) + 7.745966692414834 * q_past[28]) + nuv1 *
    ((((((((((((((((((((((((((md_residual_tmp + nd_residual_tmp) +
    od_residual_tmp) + 3.4641016151377544 * Uphi3_east_1 * Uphi3_east_2) +
    3.0983866769659336 * Uphi3_east_2 * Uphi3_east_3) + 3.04255531702266 *
    Uphi3_east_3 * q[29]) + pd_residual_tmp) + qd_residual_tmp) +
                       rd_residual_tmp) - 6.9282032302755088 * q[20] * q[24]) -
                     6.9282032302755088 * q[21] * q[25]) - 6.9282032302755088 *
                    q[22] * q[26]) - 6.1967733539318672 * q[24] * q[27]) -
                  6.1967733539318672 * q[25] * q[28]) - 6.08511063404532 * q[27]
                 * q[29]) + sd_residual_tmp) + td_residual_tmp) +
              ud_residual_tmp) + vd_residual_tmp) + wd_residual_tmp) +
           xd_residual_tmp) - yd_residual_tmp) - ae_residual_tmp) -
        be_residual_tmp) - ce_residual_tmp) - de_residual_tmp) - ee_residual_tmp);
  fe_residual_tmp = Uphi2_east_1 * Uphi2_east_2;
  ge_residual_tmp = 0.89442719099991586 * Uphi2_east_2 * Uphi2_east_3;
  he_residual_tmp = 0.87831006565367986 * Uphi2_east_3 * q[19];
  ie_residual_tmp = Uphi3_west_1 * Uphi3_west_2;
  je_residual_tmp = 0.89442719099991586 * Uphi3_west_2 * Uphi3_west_3;
  ke_residual_tmp = 0.87831006565367986 * Uphi3_west_3 * q[29];
  le_residual_tmp = Uphi2_east_1 * lambda_23_2;
  me_residual_tmp = Uphi2_east_2 * lambda_23_1;
  ne_residual_tmp = 0.89442719099991586 * Uphi2_east_2 * lambda_23_3;
  oe_residual_tmp = 0.89442719099991586 * Uphi2_east_3 * lambda_23_2;
  pe_residual_tmp = 0.87831006565367986 * Uphi2_east_3 * lambda_23_4;
  qe_residual_tmp = 0.87831006565367986 * q[19] * lambda_23_3;
  re_residual_tmp = Uphi3_west_1 * lambda_23_2;
  se_residual_tmp = Uphi3_west_2 * lambda_23_1;
  te_residual_tmp = 0.89442719099991586 * Uphi3_west_2 * lambda_23_3;
  ue_residual_tmp = 0.89442719099991586 * Uphi3_west_3 * lambda_23_2;
  ve_residual_tmp = 0.87831006565367986 * Uphi3_west_3 * lambda_23_4;
  we_residual_tmp = 0.87831006565367986 * q[29] * lambda_23_3;
  residual[24] = (((((((6.0 * q[24] - 3.4641016151377544 * q[20]) +
                       7.745966692414834 * q[27]) + 9.16515138991168 * q[29]) +
                     3.4641016151377544 * q_past[20]) + 6.0 * q_past[24]) +
                   7.745966692414834 * q_past[27]) + 9.16515138991168 * q_past
                  [29]) + nuv1 * ((((((((((((((((((((2.0 * Uphi3_east_1 *
    Uphi3_east_2 - ge_residual_tmp) - he_residual_tmp) - fe_residual_tmp) +
    1.7888543819998317 * Uphi3_east_2 * Uphi3_east_3) + 1.7566201313073597 *
    Uphi3_east_3 * q[29]) - ie_residual_tmp) - je_residual_tmp) -
    ke_residual_tmp) - le_residual_tmp) - me_residual_tmp) - ne_residual_tmp) -
    oe_residual_tmp) - pe_residual_tmp) - qe_residual_tmp) + re_residual_tmp) +
    se_residual_tmp) + te_residual_tmp) + ue_residual_tmp) + ve_residual_tmp) +
    we_residual_tmp);
  xe_residual_tmp = Uphi2_east_1 * Uphi2_east_1;
  ye_residual_tmp = Uphi3_east_1 * Uphi3_east_1;
  af_residual_tmp = Uphi3_west_1 * Uphi3_west_1;
  bf_residual_tmp = q[20] * q[20];
  cf_residual_tmp = q[21] * q[21];
  df_residual_tmp = q[22] * q[22];
  ef_residual_tmp = q[23] * q[23];
  ff_residual_tmp = 2.6457513110645907 * Uphi2_east_1 * lambda_23_1;
  gf_residual_tmp = 2.6457513110645907 * Uphi2_east_2 * lambda_23_2;
  hf_residual_tmp = 2.6457513110645907 * Uphi2_east_3 * lambda_23_3;
  if_residual_tmp = 2.6457513110645907 * q[19] * lambda_23_4;
  jf_residual_tmp = 2.6457513110645907 * Uphi3_west_1 * lambda_23_1;
  kf_residual_tmp = 2.6457513110645907 * Uphi3_west_2 * lambda_23_2;
  lf_residual_tmp = 2.6457513110645907 * Uphi3_west_3 * lambda_23_3;
  mf_residual_tmp = 2.6457513110645907 * q[29] * lambda_23_4;
  nf_residual_tmp = 1.3228756555322954 * xe_residual_tmp;
  of_residual_tmp = 1.3228756555322954 * l_residual_tmp;
  pf_residual_tmp = 1.3228756555322954 * m_residual_tmp;
  qf_residual_tmp = 1.3228756555322954 * n_residual_tmp;
  rf_residual_tmp = 1.3228756555322954 * af_residual_tmp;
  sf_residual_tmp = 1.3228756555322954 * q_residual_tmp;
  tf_residual_tmp = 1.3228756555322954 * r_residual_tmp;
  uf_residual_tmp = 1.3228756555322954 * k_residual_tmp;
  residual[23] = (2.0 * q[23] - 2.0 * q_past[23]) + nuv1 *
    ((((((((((((((((((((((((((((((((ff_residual_tmp - 20.784609690826528 * q[21]
    * q[23]) - 23.664319132398465 * q[24] * q[26]) - 23.664319132398465 * q[20] *
    q[22]) + gf_residual_tmp) + hf_residual_tmp) + if_residual_tmp) -
    jf_residual_tmp) - kf_residual_tmp) - lf_residual_tmp) - mf_residual_tmp) +
    nf_residual_tmp) + of_residual_tmp) + pf_residual_tmp) + qf_residual_tmp) +
                      2.6457513110645907 * ye_residual_tmp) + 2.6457513110645907
                     * o_residual_tmp) + 2.6457513110645907 * p_residual_tmp) +
                   2.6457513110645907 * k_residual_tmp) + rf_residual_tmp) +
                 sf_residual_tmp) + tf_residual_tmp) + uf_residual_tmp) -
              5.2915026221291814 * bf_residual_tmp) - 15.874507866387544 *
             cf_residual_tmp) - 12.850792082313726 * df_residual_tmp) -
           12.346839451634755 * ef_residual_tmp) - 5.2915026221291814 *
          s_residual_tmp) - 15.874507866387544 * t_residual_tmp) -
        12.850792082313726 * u_residual_tmp) - 5.2915026221291814 *
       v_residual_tmp) - 15.874507866387544 * w_residual_tmp) -
     5.2915026221291814 * k_residual_tmp);
  vf_residual_tmp = 2.23606797749979 * Uphi2_east_1 * lambda_23_1;
  wf_residual_tmp = 2.23606797749979 * Uphi2_east_2 * lambda_23_2;
  xf_residual_tmp = 2.23606797749979 * Uphi2_east_3 * lambda_23_3;
  yf_residual_tmp = 2.23606797749979 * q[19] * lambda_23_4;
  ag_residual_tmp = 2.23606797749979 * Uphi3_west_1 * lambda_23_1;
  bg_residual_tmp = 2.23606797749979 * Uphi3_west_2 * lambda_23_2;
  cg_residual_tmp = 2.23606797749979 * Uphi3_west_3 * lambda_23_3;
  dg_residual_tmp = 2.23606797749979 * q[29] * lambda_23_4;
  eg_residual_tmp = 1.1180339887498949 * xe_residual_tmp;
  fg_residual_tmp = 1.1180339887498949 * l_residual_tmp;
  gg_residual_tmp = 1.1180339887498949 * m_residual_tmp;
  hg_residual_tmp = 1.1180339887498949 * n_residual_tmp;
  ig_residual_tmp = 1.1180339887498949 * af_residual_tmp;
  jg_residual_tmp = 1.1180339887498949 * q_residual_tmp;
  kg_residual_tmp = 1.1180339887498949 * r_residual_tmp;
  lg_residual_tmp = 1.1180339887498949 * k_residual_tmp;
  residual[22] = (((2.0 * q[22] + 3.4641016151377544 * q[26]) - 2.0 * q_past[22])
                  - 3.4641016151377544 * q_past[26]) + nuv1 *
    (((((((((((((((((((((((((ag_residual_tmp - 13.856406460551018 * q[21] * q[22])
    - 13.60672102833218 * q[22] * q[23]) - 15.491933384829668 * q[24] * q[25]) -
    13.856406460551018 * q[25] * q[26]) - 15.491933384829668 * q[27] * q[28]) -
    vf_residual_tmp) - wf_residual_tmp) - xf_residual_tmp) - yf_residual_tmp) -
                    15.491933384829668 * q[20] * q[21]) + bg_residual_tmp) +
                  cg_residual_tmp) + dg_residual_tmp) - eg_residual_tmp) -
               fg_residual_tmp) - gg_residual_tmp) - hg_residual_tmp) +
            2.23606797749979 * ye_residual_tmp) + 2.23606797749979 *
           o_residual_tmp) + 2.23606797749979 * p_residual_tmp) +
         2.23606797749979 * k_residual_tmp) - ig_residual_tmp) - jg_residual_tmp)
      - kg_residual_tmp) - lg_residual_tmp);
  mg_residual_tmp = 1.7320508075688772 * Uphi2_east_1 * lambda_23_1;
  ng_residual_tmp = 1.7320508075688772 * Uphi2_east_2 * lambda_23_2;
  og_residual_tmp = 1.7320508075688772 * Uphi2_east_3 * lambda_23_3;
  pg_residual_tmp = 1.7320508075688772 * q[19] * lambda_23_4;
  qg_residual_tmp = 1.7320508075688772 * Uphi3_west_1 * lambda_23_1;
  rg_residual_tmp = 1.7320508075688772 * Uphi3_west_2 * lambda_23_2;
  sg_residual_tmp = 1.7320508075688772 * Uphi3_west_3 * lambda_23_3;
  tg_residual_tmp = 1.7320508075688772 * q[29] * lambda_23_4;
  ug_residual_tmp = 0.8660254037844386 * xe_residual_tmp;
  vg_residual_tmp = 0.8660254037844386 * l_residual_tmp;
  wg_residual_tmp = 0.8660254037844386 * m_residual_tmp;
  xg_residual_tmp = 0.8660254037844386 * n_residual_tmp;
  yg_residual_tmp = 0.8660254037844386 * af_residual_tmp;
  ah_residual_tmp = 0.8660254037844386 * q_residual_tmp;
  bh_residual_tmp = 0.8660254037844386 * r_residual_tmp;
  ch_residual_tmp = 0.8660254037844386 * k_residual_tmp;
  residual[21] = (((((2.0 * q[21] + 3.4641016151377544 * q[25]) +
                     4.47213595499958 * q[28]) - 2.0 * q_past[21]) -
                   3.4641016151377544 * q_past[25]) - 4.47213595499958 * q_past
                  [28]) + nuv1 * (((((((((((((((((((((((((((((mg_residual_tmp +
    ng_residual_tmp) + og_residual_tmp) + pg_residual_tmp) - qg_residual_tmp) -
    rg_residual_tmp) - sg_residual_tmp) - tg_residual_tmp) + ug_residual_tmp) +
    vg_residual_tmp) + wg_residual_tmp) + xg_residual_tmp) + 1.7320508075688772 *
    ye_residual_tmp) + 1.7320508075688772 * o_residual_tmp) + 1.7320508075688772
    * p_residual_tmp) + 1.7320508075688772 * k_residual_tmp) + yg_residual_tmp)
    + ah_residual_tmp) + bh_residual_tmp) + ch_residual_tmp) -
    3.4641016151377544 * bf_residual_tmp) - 3.4641016151377544 * cf_residual_tmp)
    - 3.4641016151377544 * df_residual_tmp) - 3.4641016151377544 *
    ef_residual_tmp) - 3.4641016151377544 * s_residual_tmp) - 3.4641016151377544
    * t_residual_tmp) - 3.4641016151377544 * u_residual_tmp) -
    3.4641016151377544 * v_residual_tmp) - 3.4641016151377544 * w_residual_tmp)
    - 3.4641016151377544 * k_residual_tmp);
  s_residual_tmp = Uphi2_east_1 * lambda_23_1;
  t_residual_tmp = Uphi2_east_2 * lambda_23_2;
  u_residual_tmp = Uphi2_east_3 * lambda_23_3;
  v_residual_tmp = q[19] * lambda_23_4;
  w_residual_tmp = Uphi3_west_1 * lambda_23_1;
  bf_residual_tmp = Uphi3_west_2 * lambda_23_2;
  cf_residual_tmp = Uphi3_west_3 * lambda_23_3;
  df_residual_tmp = q[29] * lambda_23_4;
  xe_residual_tmp *= 0.5;
  l_residual_tmp *= 0.5;
  m_residual_tmp *= 0.5;
  ef_residual_tmp = 0.5 * n_residual_tmp;
  af_residual_tmp *= 0.5;
  q_residual_tmp *= 0.5;
  r_residual_tmp *= 0.5;
  dh_residual_tmp = 0.5 * k_residual_tmp;
  residual[20] = (((((((2.0 * q[20] + 3.4641016151377544 * q[24]) +
                       4.47213595499958 * q[27]) + 5.2915026221291814 * q[29]) -
                     2.0 * q_past[20]) - 3.4641016151377544 * q_past[24]) -
                   4.47213595499958 * q_past[27]) - 5.2915026221291814 * q_past
                  [29]) + nuv1 * (((((((((((((((((((w_residual_tmp -
    t_residual_tmp) - u_residual_tmp) - v_residual_tmp) - s_residual_tmp) +
    bf_residual_tmp) + cf_residual_tmp) + df_residual_tmp) - xe_residual_tmp) -
    l_residual_tmp) - m_residual_tmp) - ef_residual_tmp) + ye_residual_tmp) +
    o_residual_tmp) + p_residual_tmp) + k_residual_tmp) - af_residual_tmp) -
    q_residual_tmp) - r_residual_tmp) - dh_residual_tmp);
  k_residual_tmp = Uphi1_east_1 * q[9];
  o_residual_tmp = 0.87831006565367986 * Uphi1_east_2 * Uphi1_east_3;
  p_residual_tmp = 0.59628479399994394 * Uphi1_east_3 * q[9];
  ye_residual_tmp = Uphi2_west_1 * q[19];
  eh_residual_tmp = 0.87831006565367986 * Uphi2_west_2 * Uphi2_west_3;
  fh_residual_tmp = 0.59628479399994394 * Uphi2_west_3 * q[19];
  gh_residual_tmp = Uphi1_east_1 * lambda_12_4;
  hh_residual_tmp = 0.87831006565367986 * Uphi1_east_2 * lambda_12_3;
  ih_residual_tmp = 0.87831006565367986 * Uphi1_east_3 * lambda_12_2;
  jh_residual_tmp = q[9] * lambda_12_1;
  kh_residual_tmp = 0.59628479399994394 * Uphi1_east_3 * lambda_12_4;
  lh_residual_tmp = 0.59628479399994394 * q[9] * lambda_12_3;
  mh_residual_tmp = Uphi2_west_1 * lambda_12_4;
  nh_residual_tmp = 0.87831006565367986 * Uphi2_west_2 * lambda_12_3;
  oh_residual_tmp = 0.87831006565367986 * Uphi2_west_3 * lambda_12_2;
  ph_residual_tmp = q[19] * lambda_12_1;
  qh_residual_tmp = 0.59628479399994394 * Uphi2_west_3 * lambda_12_4;
  rh_residual_tmp = 0.59628479399994394 * q[19] * lambda_12_3;
  residual[19] = (((((((9.16515138991168 * q[14] - 5.2915026221291814 * q[10]) -
                       11.832159566199232 * q[17]) + 14.0 * q[19]) +
                     5.2915026221291814 * q_past[10]) + 9.16515138991168 *
                    q_past[14]) + 11.832159566199232 * q_past[17]) + 14.0 *
                  q_past[19]) + nuv1 *
    (((((((((((((((((((((((((((((((((((lambda_sample_12_1 - o_residual_tmp) -
    p_residual_tmp) - k_residual_tmp) + lambda_sample_12_2) + lambda_sample_12_3)
    - ye_residual_tmp) - eh_residual_tmp) - fh_residual_tmp) +
    lambda_sample_12_4) + lambda_sample_23_1) + lambda_sample_23_2) -
    gh_residual_tmp) - hh_residual_tmp) - ih_residual_tmp) - jh_residual_tmp) -
    kh_residual_tmp) - lh_residual_tmp) + mh_residual_tmp) + nh_residual_tmp) +
                    oh_residual_tmp) + ph_residual_tmp) + qh_residual_tmp) +
                 rh_residual_tmp) + lambda_sample_23_3) + lambda_sample_23_4) +
              residual_tmp) + b_residual_tmp) + c_residual_tmp) + d_residual_tmp)
          - e_residual_tmp) - f_residual_tmp) - g_residual_tmp) - h_residual_tmp)
      - i_residual_tmp) - j_residual_tmp);
  d0 = 1.7320508075688772 * Uphi1_east_1 * Uphi1_east_3 + 1.52127765851133 *
    Uphi1_east_2 * q[9];
  d1 = 1.7320508075688772 * Uphi2_west_1 * Uphi2_west_3;
  d2 = 1.52127765851133 * Uphi2_west_2 * q[19];
  d3 = 1.7320508075688772 * Uphi1_east_1 * lambda_12_3;
  lambda_sample_12_1 = Uphi1_east_2 * Uphi1_east_2;
  lambda_sample_12_2 = Uphi1_east_3 * Uphi1_east_3;
  lambda_sample_12_3 = q[9] * q[9];
  lambda_sample_12_4 = Uphi2_west_2 * Uphi2_west_2;
  lambda_sample_23_1 = Uphi2_west_3 * Uphi2_west_3;
  lambda_sample_23_2 = q[14] * q[14];
  lambda_sample_23_3 = q[15] * q[15];
  lambda_sample_23_4 = q[16] * q[16];
  residual_tmp = q[17] * q[17];
  b_residual_tmp = q[18] * q[18];
  c_residual_tmp = 1.1065666703449764 * Uphi1_east_3 * lambda_12_3;
  d_residual_tmp = 1.52127765851133 * q[9] * lambda_12_2;
  e_residual_tmp = 1.0327955589886444 * q[9] * lambda_12_4;
  f_residual_tmp = 1.7320508075688772 * Uphi2_west_1 * lambda_12_3;
  g_residual_tmp = 1.5491933384829668 * Uphi2_west_2 * lambda_12_2;
  h_residual_tmp = 1.7320508075688772 * Uphi2_west_3 * lambda_12_1;
  i_residual_tmp = 1.52127765851133 * Uphi2_west_2 * lambda_12_4;
  j_residual_tmp = 1.1065666703449764 * Uphi2_west_3 * lambda_12_3;
  sh_residual_tmp = 1.52127765851133 * q[19] * lambda_12_2;
  th_residual_tmp = 1.0327955589886444 * q[19] * lambda_12_4;
  uh_residual_tmp = 0.7745966692414834 * lambda_sample_12_1;
  vh_residual_tmp = 0.55328333517248818 * lambda_sample_12_2;
  wh_residual_tmp = 0.5163977794943222 * lambda_sample_12_3;
  xh_residual_tmp = 0.7745966692414834 * lambda_sample_12_4;
  yh_residual_tmp = 0.55328333517248818 * lambda_sample_23_1;
  residual[18] = (((((4.47213595499958 * q[11] - 7.745966692414834 * q[15]) +
                     10.0 * q[18]) - 4.47213595499958 * q_past[11]) -
                   7.745966692414834 * q_past[15]) - 10.0 * q_past[18]) + nuv1 *
    (((((((((((((((((((((((((((((((((((((((((((((((((((((((d0 + x_residual_tmp)
    + y_residual_tmp) + d1) + d2) + ab_residual_tmp) + bb_residual_tmp) -
    6.9282032302755088 * q[10] * q[17]) - 6.9282032302755088 * q[11] * q[18]) -
    6.08511063404532 * q[14] * q[19]) + d3) + 1.5491933384829668 * Uphi1_east_2 *
    lambda_12_2) + 1.7320508075688772 * Uphi1_east_3 * lambda_12_1) +
    1.52127765851133 * Uphi1_east_2 * lambda_12_4) + c_residual_tmp) +
    d_residual_tmp) + e_residual_tmp) - f_residual_tmp) - g_residual_tmp) -
    h_residual_tmp) - i_residual_tmp) - j_residual_tmp) - sh_residual_tmp) -
    th_residual_tmp) + cb_residual_tmp) + db_residual_tmp) + eb_residual_tmp) +
    fb_residual_tmp) + gb_residual_tmp) + hb_residual_tmp) + ib_residual_tmp) -
    jb_residual_tmp) - kb_residual_tmp) - lb_residual_tmp) - mb_residual_tmp) -
    nb_residual_tmp) - ob_residual_tmp) - pb_residual_tmp) + uh_residual_tmp) +
                     vh_residual_tmp) + wh_residual_tmp) + qb_residual_tmp) +
                  rb_residual_tmp) + sb_residual_tmp) + xh_residual_tmp) +
               yh_residual_tmp) + sb_residual_tmp) + tb_residual_tmp) +
            ub_residual_tmp) + vb_residual_tmp) - 3.0983866769659336 *
          lambda_sample_23_2) - 3.0983866769659336 * lambda_sample_23_3) -
        3.0983866769659336 * lambda_sample_23_4) - 2.2131333406899527 *
       residual_tmp) - 2.2131333406899527 * b_residual_tmp) - 2.0655911179772892
     * n_residual_tmp);
  d4 = 0.87831006565367986 * Uphi1_east_2 * q[9];
  d5 = Uphi1_east_1 * Uphi1_east_3;
  x_residual_tmp = Uphi2_west_1 * Uphi2_west_3;
  y_residual_tmp = Uphi1_east_1 * lambda_12_3;
  ab_residual_tmp = 0.89442719099991586 * Uphi1_east_2 * lambda_12_2;
  bb_residual_tmp = Uphi1_east_3 * lambda_12_1;
  cb_residual_tmp = 0.87831006565367986 * Uphi1_east_2 * lambda_12_4;
  db_residual_tmp = 0.63887656499993994 * Uphi1_east_3 * lambda_12_3;
  eb_residual_tmp = 0.87831006565367986 * q[9] * lambda_12_2;
  fb_residual_tmp = 0.59628479399994394 * q[9] * lambda_12_4;
  gb_residual_tmp = Uphi2_west_1 * lambda_12_3;
  hb_residual_tmp = 0.89442719099991586 * Uphi2_west_2 * lambda_12_2;
  ib_residual_tmp = Uphi2_west_3 * lambda_12_1;
  jb_residual_tmp = 0.87831006565367986 * Uphi2_west_2 * lambda_12_4;
  kb_residual_tmp = 0.63887656499993994 * Uphi2_west_3 * lambda_12_3;
  lb_residual_tmp = 0.87831006565367986 * q[19] * lambda_12_2;
  mb_residual_tmp = 0.59628479399994394 * q[19] * lambda_12_4;
  nb_residual_tmp = 0.44721359549995793 * lambda_sample_12_1;
  ob_residual_tmp = 0.31943828249996997 * lambda_sample_12_2;
  pb_residual_tmp = 0.29814239699997197 * lambda_sample_12_3;
  qb_residual_tmp = 0.44721359549995793 * lambda_sample_12_4;
  rb_residual_tmp = 0.31943828249996997 * lambda_sample_23_1;
  residual[17] = (((((((4.47213595499958 * q[10] - 7.745966692414834 * q[14]) +
                       10.0 * q[17]) + 11.832159566199232 * q[19]) -
                     4.47213595499958 * q_past[10]) - 7.745966692414834 *
                    q_past[14]) - 10.0 * q_past[17]) - 11.832159566199232 *
                  q_past[19]) + nuv1 *
    (((((((((((((((((((((((((((((((((((((((((((((((wb_residual_tmp - d4) - d5) +
    0.87831006565367986 * Uphi2_east_2 * q[19]) - x_residual_tmp) -
    0.87831006565367986 * Uphi2_west_2 * q[19]) + xb_residual_tmp) +
    yb_residual_tmp) - y_residual_tmp) - ab_residual_tmp) - bb_residual_tmp) -
    cb_residual_tmp) - db_residual_tmp) - eb_residual_tmp) - fb_residual_tmp) +
    gb_residual_tmp) + hb_residual_tmp) + ib_residual_tmp) + jb_residual_tmp) +
    kb_residual_tmp) + lb_residual_tmp) + mb_residual_tmp) + ac_residual_tmp) +
    bc_residual_tmp) + cc_residual_tmp) + dc_residual_tmp) + ec_residual_tmp) +
    fc_residual_tmp) + gc_residual_tmp) - hc_residual_tmp) - ic_residual_tmp) -
                     jc_residual_tmp) - kc_residual_tmp) - lc_residual_tmp) -
                  mc_residual_tmp) - nc_residual_tmp) - nb_residual_tmp) -
               ob_residual_tmp) - pb_residual_tmp) + oc_residual_tmp) +
            pc_residual_tmp) + qc_residual_tmp) - qb_residual_tmp) -
         rb_residual_tmp) - qc_residual_tmp) + rc_residual_tmp) +
      sc_residual_tmp) + tc_residual_tmp);
  tb_residual_tmp = 2.23606797749979 * Uphi1_east_1 * Uphi1_east_2;
  ub_residual_tmp = 2.23606797749979 * Uphi2_west_1 * Uphi2_west_2;
  vb_residual_tmp = 2.0 * Uphi2_west_2 * Uphi2_west_3;
  wb_residual_tmp = 1.9639610121239315 * Uphi2_west_3 * q[19];
  xb_residual_tmp = 2.23606797749979 * Uphi1_east_1 * lambda_12_2;
  yb_residual_tmp = 2.23606797749979 * Uphi1_east_2 * lambda_12_1;
  ac_residual_tmp = 2.0 * Uphi1_east_2 * lambda_12_3;
  bc_residual_tmp = 2.0 * Uphi1_east_3 * lambda_12_2;
  cc_residual_tmp = 1.9639610121239315 * Uphi1_east_3 * lambda_12_4;
  dc_residual_tmp = 1.9639610121239315 * q[9] * lambda_12_3;
  ec_residual_tmp = 2.23606797749979 * Uphi2_west_1 * lambda_12_2;
  fc_residual_tmp = 2.23606797749979 * Uphi2_west_2 * lambda_12_1;
  gc_residual_tmp = 2.0 * Uphi2_west_2 * lambda_12_3;
  hc_residual_tmp = 2.0 * Uphi2_west_3 * lambda_12_2;
  ic_residual_tmp = 1.9639610121239315 * Uphi2_west_3 * lambda_12_4;
  jc_residual_tmp = 1.9639610121239315 * q[19] * lambda_12_3;
  residual[16] = (((6.0 * q[16] - 3.4641016151377544 * q[12]) +
                   3.4641016151377544 * q_past[12]) + 6.0 * q_past[16]) + nuv1 *
    ((((((((((((((((((((((((((((((((((((((((((((2.23606797749979 * Uphi2_east_1 *
    Uphi2_east_2 - 2.0 * Uphi1_east_2 * Uphi1_east_3) - 1.9639610121239315 *
    Uphi1_east_3 * q[9]) - tb_residual_tmp) + uc_residual_tmp) + vc_residual_tmp)
    - ub_residual_tmp) - vb_residual_tmp) - wb_residual_tmp) + wc_residual_tmp)
    + xc_residual_tmp) + yc_residual_tmp) - 15.491933384829668 * q[10] * q[15])
    - 15.491933384829668 * q[11] * q[14]) - 13.856406460551018 * q[11] * q[16])
    - 13.856406460551018 * q[12] * q[15]) - 13.60672102833218 * q[13] * q[16]) -
    13.856406460551018 * q[14] * q[18]) - 13.856406460551018 * q[15] * q[17]) -
    12.393546707863734 * q[16] * q[18]) - 13.60672102833218 * q[18] * q[19]) -
    xb_residual_tmp) - yb_residual_tmp) - ac_residual_tmp) - bc_residual_tmp) -
    cc_residual_tmp) - dc_residual_tmp) + ec_residual_tmp) + fc_residual_tmp) +
                    gc_residual_tmp) + hc_residual_tmp) + ic_residual_tmp) +
                 jc_residual_tmp) + ad_residual_tmp) + bd_residual_tmp) +
              cd_residual_tmp) + dd_residual_tmp) + ed_residual_tmp) +
           fd_residual_tmp) - gd_residual_tmp) - hd_residual_tmp) -
        id_residual_tmp) - jd_residual_tmp) - kd_residual_tmp) - ld_residual_tmp);
  kc_residual_tmp = (1.7320508075688772 * Uphi1_east_1 * Uphi1_east_2 +
                     1.5491933384829668 * Uphi1_east_2 * Uphi1_east_3) +
    1.52127765851133 * Uphi1_east_3 * q[9];
  lc_residual_tmp = 1.7320508075688772 * Uphi2_west_1 * Uphi2_west_2;
  mc_residual_tmp = 1.5491933384829668 * Uphi2_west_2 * Uphi2_west_3;
  nc_residual_tmp = 1.52127765851133 * Uphi2_west_3 * q[19];
  oc_residual_tmp = 1.7320508075688772 * Uphi1_east_1 * lambda_12_2;
  pc_residual_tmp = 1.7320508075688772 * Uphi1_east_2 * lambda_12_1;
  rc_residual_tmp = 1.5491933384829668 * Uphi1_east_2 * lambda_12_3;
  sc_residual_tmp = 1.5491933384829668 * Uphi1_east_3 * lambda_12_2;
  tc_residual_tmp = 1.52127765851133 * Uphi1_east_3 * lambda_12_4;
  uc_residual_tmp = 1.52127765851133 * q[9] * lambda_12_3;
  vc_residual_tmp = 1.7320508075688772 * Uphi2_west_1 * lambda_12_2;
  wc_residual_tmp = 1.7320508075688772 * Uphi2_west_2 * lambda_12_1;
  xc_residual_tmp = 1.5491933384829668 * Uphi2_west_2 * lambda_12_3;
  yc_residual_tmp = 1.5491933384829668 * Uphi2_west_3 * lambda_12_2;
  ad_residual_tmp = 1.52127765851133 * Uphi2_west_3 * lambda_12_4;
  bd_residual_tmp = 1.52127765851133 * q[19] * lambda_12_3;
  residual[15] = (((((6.0 * q[15] - 3.4641016151377544 * q[11]) +
                     7.745966692414834 * q[18]) + 3.4641016151377544 * q_past[11])
                   + 6.0 * q_past[15]) + 7.745966692414834 * q_past[18]) + nuv1 *
    (((((((((((((((((((((((((((((((((((((((kc_residual_tmp + md_residual_tmp) +
    nd_residual_tmp) + od_residual_tmp) + lc_residual_tmp) + mc_residual_tmp) +
    nc_residual_tmp) + pd_residual_tmp) + qd_residual_tmp) + rd_residual_tmp) -
    6.9282032302755088 * q[10] * q[14]) - 6.9282032302755088 * q[11] * q[15]) -
    6.9282032302755088 * q[12] * q[16]) - 6.1967733539318672 * q[14] * q[17]) -
    6.1967733539318672 * q[15] * q[18]) - 6.08511063404532 * q[17] * q[19]) +
    oc_residual_tmp) + pc_residual_tmp) + rc_residual_tmp) + sc_residual_tmp) +
    tc_residual_tmp) + uc_residual_tmp) - vc_residual_tmp) - wc_residual_tmp) -
                    xc_residual_tmp) - yc_residual_tmp) - ad_residual_tmp) -
                 bd_residual_tmp) + sd_residual_tmp) + td_residual_tmp) +
              ud_residual_tmp) + vd_residual_tmp) + wd_residual_tmp) +
           xd_residual_tmp) - yd_residual_tmp) - ae_residual_tmp) -
        be_residual_tmp) - ce_residual_tmp) - de_residual_tmp) - ee_residual_tmp);
  cd_residual_tmp = Uphi1_east_1 * Uphi1_east_2;
  dd_residual_tmp = 0.89442719099991586 * Uphi1_east_2 * Uphi1_east_3;
  ed_residual_tmp = 0.87831006565367986 * Uphi1_east_3 * q[9];
  fd_residual_tmp = Uphi2_west_1 * Uphi2_west_2;
  gd_residual_tmp = 0.89442719099991586 * Uphi2_west_2 * Uphi2_west_3;
  hd_residual_tmp = 0.87831006565367986 * Uphi2_west_3 * q[19];
  id_residual_tmp = Uphi1_east_1 * lambda_12_2;
  jd_residual_tmp = Uphi1_east_2 * lambda_12_1;
  kd_residual_tmp = 0.89442719099991586 * Uphi1_east_2 * lambda_12_3;
  ld_residual_tmp = 0.89442719099991586 * Uphi1_east_3 * lambda_12_2;
  md_residual_tmp = 0.87831006565367986 * Uphi1_east_3 * lambda_12_4;
  nd_residual_tmp = 0.87831006565367986 * q[9] * lambda_12_3;
  od_residual_tmp = Uphi2_west_1 * lambda_12_2;
  pd_residual_tmp = Uphi2_west_2 * lambda_12_1;
  qd_residual_tmp = 0.89442719099991586 * Uphi2_west_2 * lambda_12_3;
  rd_residual_tmp = 0.89442719099991586 * Uphi2_west_3 * lambda_12_2;
  sd_residual_tmp = 0.87831006565367986 * Uphi2_west_3 * lambda_12_4;
  td_residual_tmp = 0.87831006565367986 * q[19] * lambda_12_3;
  residual[14] = (((((((6.0 * q[14] - 3.4641016151377544 * q[10]) +
                       7.745966692414834 * q[17]) + 9.16515138991168 * q[19]) +
                     3.4641016151377544 * q_past[10]) + 6.0 * q_past[14]) +
                   7.745966692414834 * q_past[17]) + 9.16515138991168 * q_past
                  [19]) + nuv1 *
    (((((((((((((((((((((((((((((((((((fe_residual_tmp - dd_residual_tmp) -
    ed_residual_tmp) - cd_residual_tmp) + ge_residual_tmp) + he_residual_tmp) -
    fd_residual_tmp) - gd_residual_tmp) - hd_residual_tmp) + ie_residual_tmp) +
    je_residual_tmp) + ke_residual_tmp) - id_residual_tmp) - jd_residual_tmp) -
    kd_residual_tmp) - ld_residual_tmp) - md_residual_tmp) - nd_residual_tmp) +
                      od_residual_tmp) + pd_residual_tmp) + qd_residual_tmp) +
                   rd_residual_tmp) + sd_residual_tmp) + td_residual_tmp) +
                le_residual_tmp) + me_residual_tmp) + ne_residual_tmp) +
             oe_residual_tmp) + pe_residual_tmp) + qe_residual_tmp) -
          re_residual_tmp) - se_residual_tmp) - te_residual_tmp) -
       ue_residual_tmp) - ve_residual_tmp) - we_residual_tmp);
  ud_residual_tmp = Uphi1_east_1 * Uphi1_east_1;
  vd_residual_tmp = Uphi2_west_1 * Uphi2_west_1;
  wd_residual_tmp = q[10] * q[10];
  xd_residual_tmp = q[11] * q[11];
  yd_residual_tmp = q[12] * q[12];
  ae_residual_tmp = q[13] * q[13];
  be_residual_tmp = 2.6457513110645907 * Uphi1_east_2 * lambda_12_2;
  ce_residual_tmp = 2.6457513110645907 * Uphi1_east_3 * lambda_12_3;
  de_residual_tmp = 2.6457513110645907 * q[9] * lambda_12_4;
  ee_residual_tmp = 2.6457513110645907 * Uphi2_west_1 * lambda_12_1;
  fe_residual_tmp = 2.6457513110645907 * Uphi2_west_2 * lambda_12_2;
  ge_residual_tmp = 2.6457513110645907 * Uphi2_west_3 * lambda_12_3;
  he_residual_tmp = 2.6457513110645907 * q[19] * lambda_12_4;
  ie_residual_tmp = 1.3228756555322954 * ud_residual_tmp;
  je_residual_tmp = 1.3228756555322954 * lambda_sample_12_1;
  ke_residual_tmp = 1.3228756555322954 * lambda_sample_12_2;
  le_residual_tmp = 1.3228756555322954 * lambda_sample_12_3;
  me_residual_tmp = 1.3228756555322954 * vd_residual_tmp;
  ne_residual_tmp = 1.3228756555322954 * lambda_sample_12_4;
  oe_residual_tmp = 1.3228756555322954 * lambda_sample_23_1;
  residual[13] = (2.0 * q[13] - 2.0 * q_past[13]) + nuv1 *
    ((((((((((((((((((((((((((((((((((((((((((((2.6457513110645907 *
    Uphi1_east_1 * lambda_12_1 - 20.784609690826528 * q[11] * q[13]) -
    23.664319132398465 * q[14] * q[16]) - 23.664319132398465 * q[10] * q[12]) +
    be_residual_tmp) + ce_residual_tmp) + de_residual_tmp) - ee_residual_tmp) -
    fe_residual_tmp) - ge_residual_tmp) - he_residual_tmp) + ff_residual_tmp) +
    gf_residual_tmp) + hf_residual_tmp) + if_residual_tmp) - jf_residual_tmp) -
    kf_residual_tmp) - lf_residual_tmp) - mf_residual_tmp) + ie_residual_tmp) +
    je_residual_tmp) + ke_residual_tmp) + le_residual_tmp) + nf_residual_tmp) +
    of_residual_tmp) + pf_residual_tmp) + qf_residual_tmp) + me_residual_tmp) +
                     ne_residual_tmp) + oe_residual_tmp) + qf_residual_tmp) +
                  rf_residual_tmp) + sf_residual_tmp) + tf_residual_tmp) +
               uf_residual_tmp) - 5.2915026221291814 * wd_residual_tmp) -
             15.874507866387544 * xd_residual_tmp) - 12.850792082313726 *
            yd_residual_tmp) - 12.346839451634755 * ae_residual_tmp) -
          5.2915026221291814 * lambda_sample_23_2) - 15.874507866387544 *
         lambda_sample_23_3) - 12.850792082313726 * lambda_sample_23_4) -
       5.2915026221291814 * residual_tmp) - 15.874507866387544 * b_residual_tmp)
     - 5.2915026221291814 * n_residual_tmp);
  pe_residual_tmp = 2.23606797749979 * Uphi1_east_1 * lambda_12_1;
  qe_residual_tmp = 2.23606797749979 * Uphi1_east_2 * lambda_12_2;
  re_residual_tmp = 2.23606797749979 * Uphi1_east_3 * lambda_12_3;
  se_residual_tmp = 2.23606797749979 * q[9] * lambda_12_4;
  te_residual_tmp = 2.23606797749979 * Uphi2_west_1 * lambda_12_1;
  ue_residual_tmp = 2.23606797749979 * Uphi2_west_2 * lambda_12_2;
  ve_residual_tmp = 2.23606797749979 * Uphi2_west_3 * lambda_12_3;
  we_residual_tmp = 2.23606797749979 * q[19] * lambda_12_4;
  ff_residual_tmp = 1.1180339887498949 * ud_residual_tmp;
  gf_residual_tmp = 1.1180339887498949 * lambda_sample_12_1;
  hf_residual_tmp = 1.1180339887498949 * lambda_sample_12_2;
  if_residual_tmp = 1.1180339887498949 * lambda_sample_12_3;
  jf_residual_tmp = 1.1180339887498949 * vd_residual_tmp;
  kf_residual_tmp = 1.1180339887498949 * lambda_sample_12_4;
  lf_residual_tmp = 1.1180339887498949 * lambda_sample_23_1;
  residual[12] = (((2.0 * q[12] + 3.4641016151377544 * q[16]) - 2.0 * q_past[12])
                  - 3.4641016151377544 * q_past[16]) + nuv1 *
    (((((((((((((((((((((((((((((((((((((te_residual_tmp - 13.856406460551018 *
    q[11] * q[12]) - 13.60672102833218 * q[12] * q[13]) - 15.491933384829668 *
    q[14] * q[15]) - 13.856406460551018 * q[15] * q[16]) - 15.491933384829668 *
    q[17] * q[18]) - pe_residual_tmp) - qe_residual_tmp) - re_residual_tmp) -
    se_residual_tmp) - 15.491933384829668 * q[10] * q[11]) + ue_residual_tmp) +
    ve_residual_tmp) + we_residual_tmp) + vf_residual_tmp) + wf_residual_tmp) +
    xf_residual_tmp) + yf_residual_tmp) - ag_residual_tmp) - bg_residual_tmp) -
                      cg_residual_tmp) - dg_residual_tmp) - ff_residual_tmp) -
                   gf_residual_tmp) - hf_residual_tmp) - if_residual_tmp) +
                eg_residual_tmp) + fg_residual_tmp) + gg_residual_tmp) +
             hg_residual_tmp) - jf_residual_tmp) - kf_residual_tmp) -
          lf_residual_tmp) - hg_residual_tmp) + ig_residual_tmp) +
       jg_residual_tmp) + kg_residual_tmp) + lg_residual_tmp);
  mf_residual_tmp = ((((((1.7320508075688772 * Uphi1_east_1 * lambda_12_1 +
    1.7320508075688772 * Uphi1_east_2 * lambda_12_2) + 1.7320508075688772 *
    Uphi1_east_3 * lambda_12_3) + 1.7320508075688772 * q[9] * lambda_12_4) -
                       1.7320508075688772 * Uphi2_west_1 * lambda_12_1) -
                      1.7320508075688772 * Uphi2_west_2 * lambda_12_2) -
                     1.7320508075688772 * Uphi2_west_3 * lambda_12_3) -
    1.7320508075688772 * q[19] * lambda_12_4;
  nf_residual_tmp = 0.8660254037844386 * ud_residual_tmp;
  of_residual_tmp = 0.8660254037844386 * lambda_sample_12_1;
  pf_residual_tmp = 0.8660254037844386 * lambda_sample_12_2;
  rf_residual_tmp = 0.8660254037844386 * lambda_sample_12_3;
  sf_residual_tmp = 0.8660254037844386 * vd_residual_tmp;
  tf_residual_tmp = 0.8660254037844386 * lambda_sample_12_4;
  uf_residual_tmp = 0.8660254037844386 * lambda_sample_23_1;
  residual[11] = (((((2.0 * q[11] + 3.4641016151377544 * q[15]) +
                     4.47213595499958 * q[18]) - 2.0 * q_past[11]) -
                   3.4641016151377544 * q_past[15]) - 4.47213595499958 * q_past
                  [18]) + nuv1 *
    ((((((((((((((((((((((((((((((((((mf_residual_tmp + mg_residual_tmp) +
    ng_residual_tmp) + og_residual_tmp) + pg_residual_tmp) - qg_residual_tmp) -
    rg_residual_tmp) - sg_residual_tmp) - tg_residual_tmp) + nf_residual_tmp) +
    of_residual_tmp) + pf_residual_tmp) + rf_residual_tmp) + ug_residual_tmp) +
    vg_residual_tmp) + wg_residual_tmp) + xg_residual_tmp) + sf_residual_tmp) +
                     tf_residual_tmp) + uf_residual_tmp) + xg_residual_tmp) +
                  yg_residual_tmp) + ah_residual_tmp) + bh_residual_tmp) +
               ch_residual_tmp) - 3.4641016151377544 * wd_residual_tmp) -
             3.4641016151377544 * xd_residual_tmp) - 3.4641016151377544 *
            yd_residual_tmp) - 3.4641016151377544 * ae_residual_tmp) -
          3.4641016151377544 * lambda_sample_23_2) - 3.4641016151377544 *
         lambda_sample_23_3) - 3.4641016151377544 * lambda_sample_23_4) -
       3.4641016151377544 * residual_tmp) - 3.4641016151377544 * b_residual_tmp)
     - 3.4641016151377544 * n_residual_tmp);
  lambda_sample_23_2 = Uphi1_east_1 * lambda_12_1;
  lambda_sample_23_3 = Uphi1_east_2 * lambda_12_2;
  lambda_sample_23_4 = Uphi1_east_3 * lambda_12_3;
  residual_tmp = q[9] * lambda_12_4;
  b_residual_tmp = Uphi2_west_1 * lambda_12_1;
  n_residual_tmp = Uphi2_west_2 * lambda_12_2;
  wd_residual_tmp = Uphi2_west_3 * lambda_12_3;
  xd_residual_tmp = q[19] * lambda_12_4;
  ud_residual_tmp *= 0.5;
  lambda_sample_12_1 *= 0.5;
  lambda_sample_12_2 *= 0.5;
  yd_residual_tmp = 0.5 * lambda_sample_12_3;
  vd_residual_tmp *= 0.5;
  lambda_sample_12_4 *= 0.5;
  lambda_sample_23_1 *= 0.5;
  residual[10] = (((((((2.0 * q[10] + 3.4641016151377544 * q[14]) +
                       4.47213595499958 * q[17]) + 5.2915026221291814 * q[19]) -
                     2.0 * q_past[10]) - 3.4641016151377544 * q_past[14]) -
                   4.47213595499958 * q_past[17]) - 5.2915026221291814 * q_past
                  [19]) + nuv1 * (((((((((((((((((((((((((((((((b_residual_tmp -
    lambda_sample_23_3) - lambda_sample_23_4) - residual_tmp) -
    lambda_sample_23_2) + n_residual_tmp) + wd_residual_tmp) + xd_residual_tmp)
    + s_residual_tmp) + t_residual_tmp) + u_residual_tmp) + v_residual_tmp) -
    w_residual_tmp) - bf_residual_tmp) - cf_residual_tmp) - df_residual_tmp) -
    ud_residual_tmp) - lambda_sample_12_1) - lambda_sample_12_2) -
    yd_residual_tmp) + xe_residual_tmp) + l_residual_tmp) + m_residual_tmp) +
    ef_residual_tmp) - vd_residual_tmp) - lambda_sample_12_4) -
    lambda_sample_23_1) - ef_residual_tmp) + af_residual_tmp) + q_residual_tmp)
    + r_residual_tmp) + dh_residual_tmp);
  residual[9] = (((((((9.16515138991168 * q[4] - 5.2915026221291814 * q[0]) -
                      11.832159566199232 * q[7]) + 14.0 * q[9]) +
                    5.2915026221291814 * q_past[0]) + 9.16515138991168 * q_past
                   [4]) + 11.832159566199232 * q_past[7]) + 14.0 * q_past[9]) +
    nuv1 * ((((((((((((((((((((k_residual_tmp + o_residual_tmp) + p_residual_tmp)
    - 2.0 * Uphi1_west_1 * q[9]) - 1.7566201313073597 * Uphi1_west_2 *
    Uphi1_west_3) - 1.1925695879998879 * Uphi1_west_3 * q[9]) + ye_residual_tmp)
    + eh_residual_tmp) + fh_residual_tmp) + gh_residual_tmp) + hh_residual_tmp)
                     + ih_residual_tmp) + jh_residual_tmp) + kh_residual_tmp) +
                  lh_residual_tmp) - mh_residual_tmp) - nh_residual_tmp) -
               oh_residual_tmp) - ph_residual_tmp) - qh_residual_tmp) -
            rh_residual_tmp);
  k_residual_tmp = Uphi1_west_2 * Uphi1_west_2;
  l_residual_tmp = Uphi1_west_3 * Uphi1_west_3;
  m_residual_tmp = q[4] * q[4];
  o_residual_tmp = q[5] * q[5];
  p_residual_tmp = q[6] * q[6];
  q_residual_tmp = q[7] * q[7];
  r_residual_tmp = q[8] * q[8];
  residual[8] = (((((4.47213595499958 * q[1] - 7.745966692414834 * q[5]) + 10.0 *
                    q[8]) - 4.47213595499958 * q_past[1]) - 7.745966692414834 *
                  q_past[5]) - 10.0 * q_past[8]) + nuv1 *
    ((((((((((((((((((((((((((((((((((((d0 + 3.4641016151377544 * Uphi1_west_1 *
    Uphi1_west_3) + 3.04255531702266 * Uphi1_west_2 * q[9]) + d1) + d2) -
    6.9282032302755088 * q[0] * q[7]) - 6.9282032302755088 * q[1] * q[8]) -
    6.08511063404532 * q[4] * q[9]) + d3) + 1.5491933384829668 * Uphi1_east_2 *
    lambda_12_2) + 1.7320508075688772 * Uphi1_east_3 * lambda_12_1) +
    1.52127765851133 * Uphi1_east_2 * lambda_12_4) + c_residual_tmp) +
    d_residual_tmp) + e_residual_tmp) - f_residual_tmp) - g_residual_tmp) -
    h_residual_tmp) - i_residual_tmp) - j_residual_tmp) - sh_residual_tmp) -
                    th_residual_tmp) + uh_residual_tmp) + vh_residual_tmp) +
                 wh_residual_tmp) + 1.5491933384829668 * k_residual_tmp) +
               1.1065666703449764 * l_residual_tmp) + 1.0327955589886444 *
              lambda_sample_12_3) + xh_residual_tmp) + yh_residual_tmp) +
           sb_residual_tmp) - 3.0983866769659336 * m_residual_tmp) -
         3.0983866769659336 * o_residual_tmp) - 3.0983866769659336 *
        p_residual_tmp) - 2.2131333406899527 * q_residual_tmp) -
      2.2131333406899527 * r_residual_tmp) - 2.0655911179772892 *
     lambda_sample_12_3);
  residual[7] = (((((((4.47213595499958 * q[0] - 7.745966692414834 * q[4]) +
                      10.0 * q[7]) + 11.832159566199232 * q[9]) -
                    4.47213595499958 * q_past[0]) - 7.745966692414834 * q_past[4])
                  - 10.0 * q_past[7]) - 11.832159566199232 * q_past[9]) + nuv1 *
    ((((((((((((((((((((((((((((d5 + d4) - 2.0 * Uphi1_west_1 * Uphi1_west_3) -
    1.7566201313073597 * Uphi1_west_2 * q[9]) + x_residual_tmp) +
    0.87831006565367986 * Uphi2_west_2 * q[19]) + y_residual_tmp) +
    ab_residual_tmp) + bb_residual_tmp) + cb_residual_tmp) + db_residual_tmp) +
                      eb_residual_tmp) + fb_residual_tmp) - gb_residual_tmp) -
                   hb_residual_tmp) - ib_residual_tmp) - jb_residual_tmp) -
                kb_residual_tmp) - lb_residual_tmp) - mb_residual_tmp) +
             nb_residual_tmp) + ob_residual_tmp) + pb_residual_tmp) -
          0.89442719099991586 * k_residual_tmp) - 0.63887656499993994 *
         l_residual_tmp) - 0.59628479399994394 * lambda_sample_12_3) +
       qb_residual_tmp) + rb_residual_tmp) + qc_residual_tmp);
  residual[6] = (((6.0 * q[6] - 3.4641016151377544 * q[2]) + 3.4641016151377544 *
                  q_past[2]) + 6.0 * q_past[6]) + nuv1 *
    (((((((((((((((((((((((((((((tb_residual_tmp + 2.0 * Uphi1_east_2 *
    Uphi1_east_3) + 1.9639610121239315 * Uphi1_east_3 * q[9]) - 4.47213595499958
    * Uphi1_west_1 * Uphi1_west_2) - 4.0 * Uphi1_west_2 * Uphi1_west_3) -
    3.927922024247863 * Uphi1_west_3 * q[9]) + ub_residual_tmp) +
    vb_residual_tmp) + wb_residual_tmp) - 15.491933384829668 * q[0] * q[5]) -
    15.491933384829668 * q[1] * q[4]) - 13.856406460551018 * q[1] * q[6]) -
                      13.856406460551018 * q[2] * q[5]) - 13.60672102833218 * q
                     [3] * q[6]) - 13.856406460551018 * q[4] * q[8]) -
                   13.856406460551018 * q[5] * q[7]) - 12.393546707863734 * q[6]
                  * q[8]) - 13.60672102833218 * q[8] * q[9]) + xb_residual_tmp)
               + yb_residual_tmp) + ac_residual_tmp) + bc_residual_tmp) +
            cc_residual_tmp) + dc_residual_tmp) - ec_residual_tmp) -
         fc_residual_tmp) - gc_residual_tmp) - hc_residual_tmp) -
      ic_residual_tmp) - jc_residual_tmp);
  residual[5] = (((((6.0 * q[5] - 3.4641016151377544 * q[1]) + 7.745966692414834
                    * q[8]) + 3.4641016151377544 * q_past[1]) + 6.0 * q_past[5])
                 + 7.745966692414834 * q_past[8]) + nuv1 *
    ((((((((((((((((((((((((kc_residual_tmp + 3.4641016151377544 * Uphi1_west_1 *
    Uphi1_west_2) + 3.0983866769659336 * Uphi1_west_2 * Uphi1_west_3) +
    3.04255531702266 * Uphi1_west_3 * q[9]) + lc_residual_tmp) + mc_residual_tmp)
                       + nc_residual_tmp) - 6.9282032302755088 * q[0] * q[4]) -
                     6.9282032302755088 * q[1] * q[5]) - 6.9282032302755088 * q
                    [2] * q[6]) - 6.1967733539318672 * q[4] * q[7]) -
                  6.1967733539318672 * q[5] * q[8]) - 6.08511063404532 * q[7] *
                 q[9]) + oc_residual_tmp) + pc_residual_tmp) + rc_residual_tmp)
             + sc_residual_tmp) + tc_residual_tmp) + uc_residual_tmp) -
          vc_residual_tmp) - wc_residual_tmp) - xc_residual_tmp) -
       yc_residual_tmp) - ad_residual_tmp) - bd_residual_tmp);
  residual[4] = (((((((6.0 * q[4] - 3.4641016151377544 * q[0]) +
                      7.745966692414834 * q[7]) + 9.16515138991168 * q[9]) +
                    3.4641016151377544 * q_past[0]) + 6.0 * q_past[4]) +
                  7.745966692414834 * q_past[7]) + 9.16515138991168 * q_past[9])
    + nuv1 * ((((((((((((((((((((cd_residual_tmp + dd_residual_tmp) +
    ed_residual_tmp) - 2.0 * Uphi1_west_1 * Uphi1_west_2) - 1.7888543819998317 *
    Uphi1_west_2 * Uphi1_west_3) - 1.7566201313073597 * Uphi1_west_3 * q[9]) +
    fd_residual_tmp) + gd_residual_tmp) + hd_residual_tmp) + id_residual_tmp) +
                        jd_residual_tmp) + kd_residual_tmp) + ld_residual_tmp) +
                     md_residual_tmp) + nd_residual_tmp) - od_residual_tmp) -
                  pd_residual_tmp) - qd_residual_tmp) - rd_residual_tmp) -
               sd_residual_tmp) - td_residual_tmp);
  c_residual_tmp = Uphi1_west_1 * Uphi1_west_1;
  d_residual_tmp = q[0] * q[0];
  e_residual_tmp = q[1] * q[1];
  f_residual_tmp = q[2] * q[2];
  g_residual_tmp = q[3] * q[3];
  residual[3] = (2.0 * q[3] - 2.0 * q_past[3]) + nuv1 *
    ((((((((((((((((((((((((((((((((2.6457513110645907 * Uphi1_east_1 *
    lambda_12_1 - 20.784609690826528 * q[1] * q[3]) - 23.664319132398465 * q[4] *
    q[6]) - 23.664319132398465 * q[0] * q[2]) + be_residual_tmp) +
    ce_residual_tmp) + de_residual_tmp) - ee_residual_tmp) - fe_residual_tmp) -
    ge_residual_tmp) - he_residual_tmp) + ie_residual_tmp) + je_residual_tmp) +
    ke_residual_tmp) + le_residual_tmp) + 2.6457513110645907 * c_residual_tmp) +
                     2.6457513110645907 * k_residual_tmp) + 2.6457513110645907 *
                    l_residual_tmp) + 2.6457513110645907 * lambda_sample_12_3) +
                  me_residual_tmp) + ne_residual_tmp) + oe_residual_tmp) +
               qf_residual_tmp) - 5.2915026221291814 * d_residual_tmp) -
             15.874507866387544 * e_residual_tmp) - 12.850792082313726 *
            f_residual_tmp) - 12.346839451634755 * g_residual_tmp) -
          5.2915026221291814 * m_residual_tmp) - 15.874507866387544 *
         o_residual_tmp) - 12.850792082313726 * p_residual_tmp) -
       5.2915026221291814 * q_residual_tmp) - 15.874507866387544 *
      r_residual_tmp) - 5.2915026221291814 * lambda_sample_12_3);
  residual[2] = (((2.0 * q[2] + 3.4641016151377544 * q[6]) - 2.0 * q_past[2]) -
                 3.4641016151377544 * q_past[6]) + nuv1 *
    (((((((((((((((((((((((((pe_residual_tmp - 13.856406460551018 * q[1] * q[2])
    - 13.60672102833218 * q[2] * q[3]) - 15.491933384829668 * q[4] * q[5]) -
    13.856406460551018 * q[5] * q[6]) - 15.491933384829668 * q[7] * q[8]) -
    15.491933384829668 * q[0] * q[1]) + qe_residual_tmp) + re_residual_tmp) +
                     se_residual_tmp) - te_residual_tmp) - ue_residual_tmp) -
                  ve_residual_tmp) - we_residual_tmp) + ff_residual_tmp) +
               gf_residual_tmp) + hf_residual_tmp) + if_residual_tmp) -
            2.23606797749979 * c_residual_tmp) - 2.23606797749979 *
           k_residual_tmp) - 2.23606797749979 * l_residual_tmp) -
         2.23606797749979 * lambda_sample_12_3) + jf_residual_tmp) +
       kf_residual_tmp) + lf_residual_tmp) + hg_residual_tmp);
  residual[1] = (((((2.0 * q[1] + 3.4641016151377544 * q[5]) + 4.47213595499958 *
                    q[8]) - 2.0 * q_past[1]) - 3.4641016151377544 * q_past[5]) -
                 4.47213595499958 * q_past[8]) + nuv1 *
    ((((((((((((((((((((((mf_residual_tmp + nf_residual_tmp) + of_residual_tmp)
    + pf_residual_tmp) + rf_residual_tmp) + 1.7320508075688772 * c_residual_tmp)
                     + 1.7320508075688772 * k_residual_tmp) + 1.7320508075688772
                    * l_residual_tmp) + 1.7320508075688772 * lambda_sample_12_3)
                  + sf_residual_tmp) + tf_residual_tmp) + uf_residual_tmp) +
               xg_residual_tmp) - 3.4641016151377544 * d_residual_tmp) -
             3.4641016151377544 * e_residual_tmp) - 3.4641016151377544 *
            f_residual_tmp) - 3.4641016151377544 * g_residual_tmp) -
          3.4641016151377544 * m_residual_tmp) - 3.4641016151377544 *
         o_residual_tmp) - 3.4641016151377544 * p_residual_tmp) -
       3.4641016151377544 * q_residual_tmp) - 3.4641016151377544 *
      r_residual_tmp) - 3.4641016151377544 * lambda_sample_12_3);
  residual[0] = (((((((2.0 * q[0] + 3.4641016151377544 * q[4]) +
                      4.47213595499958 * q[7]) + 5.2915026221291814 * q[9]) -
                    2.0 * q_past[0]) - 3.4641016151377544 * q_past[4]) -
                  4.47213595499958 * q_past[7]) - 5.2915026221291814 * q_past[9])
    + nuv1 * (((((((((((((((((((lambda_sample_23_2 + lambda_sample_23_3) +
    lambda_sample_23_4) + residual_tmp) - b_residual_tmp) - n_residual_tmp) -
    wd_residual_tmp) - xd_residual_tmp) + ud_residual_tmp) + lambda_sample_12_1)
                       + lambda_sample_12_2) + yd_residual_tmp) - c_residual_tmp)
                    - k_residual_tmp) - l_residual_tmp) - lambda_sample_12_3) +
                 vd_residual_tmp) + lambda_sample_12_4) + lambda_sample_23_1) +
              ef_residual_tmp);
  memset(&Jacobian[0], 0, 900U * sizeof(double));
  Jacobian[899] = 14.0 + nuv1 * (((((2.0 * Uphi3_east_1 + 1.1925695879998879 *
    Uphi3_east_3) - Uphi3_west_1) - 0.59628479399994394 * Uphi3_west_3) +
    lambda_23_1) + 0.59628479399994394 * lambda_23_3);
  l_residual_tmp = ((3.04255531702266 * Uphi3_east_2 + 2.0655911179772892 * q[29])
                    + 1.52127765851133 * Uphi3_west_2) + 1.0327955589886444 * q
    [29];
  Jacobian[869] = nuv1 * ((l_residual_tmp - 1.52127765851133 * lambda_23_2) -
    1.0327955589886444 * lambda_23_4);
  m_residual_tmp = nuv1 * (((((1.7566201313073597 * Uphi3_east_2 +
    1.1925695879998879 * q[29]) - 0.87831006565367986 * Uphi3_west_2) -
    0.59628479399994394 * q[29]) + 0.87831006565367986 * lambda_23_2) +
    0.59628479399994394 * lambda_23_4);
  Jacobian[839] = -11.832159566199232 + m_residual_tmp;
  n_residual_tmp = 3.927922024247863 * Uphi3_east_3 - 1.9639610121239315 *
    Uphi3_west_3;
  Jacobian[809] = nuv1 * (n_residual_tmp + 1.9639610121239315 * lambda_23_3);
  o_residual_tmp = 3.04255531702266 * Uphi3_east_3 + 1.52127765851133 *
    Uphi3_west_3;
  Jacobian[779] = nuv1 * (o_residual_tmp - 1.52127765851133 * lambda_23_3);
  p_residual_tmp = 9.16515138991168 + nuv1 * ((1.7566201313073597 * Uphi3_east_3
    - 0.87831006565367986 * Uphi3_west_3) + 0.87831006565367986 * lambda_23_3);
  Jacobian[749] = p_residual_tmp;
  q_residual_tmp = 5.2915026221291814 * q[29] + 2.6457513110645907 * q[29];
  Jacobian[719] = nuv1 * (q_residual_tmp - 2.6457513110645907 * lambda_23_4);
  r_residual_tmp = nuv1 * ((4.47213595499958 * q[29] - 2.23606797749979 * q[29])
    + 2.23606797749979 * lambda_23_4);
  Jacobian[689] = r_residual_tmp;
  s_residual_tmp = 3.4641016151377544 * q[29] + 1.7320508075688772 * q[29];
  Jacobian[659] = nuv1 * (s_residual_tmp - 1.7320508075688772 * lambda_23_4);
  t_residual_tmp = nuv1 * ((2.0 * q[29] - q[29]) + lambda_23_4);
  Jacobian[629] = -5.2915026221291814 + t_residual_tmp;
  Jacobian[599] = nuv1 * (((-Uphi2_east_1 - 0.59628479399994394 * Uphi2_east_3)
    - lambda_23_1) - 0.59628479399994394 * lambda_23_3);
  Jacobian[569] = nuv1 * (((-1.52127765851133 * Uphi2_east_2 -
    1.0327955589886444 * q[19]) - 1.52127765851133 * lambda_23_2) -
    1.0327955589886444 * lambda_23_4);
  u_residual_tmp = nuv1 * (((-0.87831006565367986 * Uphi2_east_2 -
    0.59628479399994394 * q[19]) - 0.87831006565367986 * lambda_23_2) -
    0.59628479399994394 * lambda_23_4);
  Jacobian[539] = u_residual_tmp;
  v_residual_tmp = nuv1 * (-1.9639610121239315 * Uphi2_east_3 -
    1.9639610121239315 * lambda_23_3);
  Jacobian[509] = v_residual_tmp;
  Jacobian[479] = nuv1 * (-1.52127765851133 * Uphi2_east_3 - 1.52127765851133 *
    lambda_23_3);
  w_residual_tmp = nuv1 * (-0.87831006565367986 * Uphi2_east_3 -
    0.87831006565367986 * lambda_23_3);
  Jacobian[449] = w_residual_tmp;
  Jacobian[419] = nuv1 * (-2.6457513110645907 * q[19] - 2.6457513110645907 *
    lambda_23_4);
  x_residual_tmp = nuv1 * (-2.23606797749979 * q[19] - 2.23606797749979 *
    lambda_23_4);
  Jacobian[389] = x_residual_tmp;
  Jacobian[359] = nuv1 * (-1.7320508075688772 * q[19] - 1.7320508075688772 *
    lambda_23_4);
  y_residual_tmp = nuv1 * (-q[19] - lambda_23_4);
  Jacobian[329] = y_residual_tmp;
  Jacobian[898] = nuv1 * ((((l_residual_tmp - 6.08511063404532 * q[24]) -
    4.1311822359545785 * q[29]) - 1.52127765851133 * lambda_23_2) -
    1.0327955589886444 * lambda_23_4);
  Jacobian[868] = 10.0 + nuv1 * (((((((6.0 * Uphi3_east_1 + 3.8332593899996397 *
    Uphi3_east_3) - 3.0 * Uphi3_west_1) - 1.9166296949998198 * Uphi3_west_3) -
    6.9282032302755088 * q[21]) - 4.4262666813799054 * q[28]) + 3.0 *
    lambda_23_1) + 1.9166296949998198 * lambda_23_3);
  l_residual_tmp = ((3.4641016151377544 * Uphi3_east_1 + 2.2131333406899527 *
                     Uphi3_east_3) + 1.7320508075688772 * Uphi3_west_1) +
    1.1065666703449764 * Uphi3_west_3;
  Jacobian[838] = nuv1 * ((((l_residual_tmp - 6.9282032302755088 * q[20]) -
    4.4262666813799054 * q[27]) - 1.7320508075688772 * lambda_23_1) -
    1.1065666703449764 * lambda_23_3);
  ab_residual_tmp = ((6.9282032302755088 * Uphi3_east_2 + 6.80336051416609 * q
                      [29]) + 3.4641016151377544 * Uphi3_west_2) +
    3.4016802570830449 * q[29];
  Jacobian[808] = nuv1 * (((ab_residual_tmp - 6.1967733539318672 * q[26]) -
    3.4641016151377544 * lambda_23_2) - 3.4016802570830449 * lambda_23_4);
  bb_residual_tmp = nuv1 * ((((((5.3665631459994954 * Uphi3_east_2 +
    5.2698603939220794 * q[29]) - 2.6832815729997477 * Uphi3_west_2) -
    2.6349301969610397 * q[29]) - 6.1967733539318672 * q[25]) +
    2.6832815729997477 * lambda_23_2) + 2.6349301969610397 * lambda_23_4);
  Jacobian[778] = -7.745966692414834 + bb_residual_tmp;
  cb_residual_tmp = ((3.0983866769659336 * Uphi3_east_2 + 3.04255531702266 * q
                      [29]) + 1.5491933384829668 * Uphi3_west_2) +
    1.52127765851133 * q[29];
  db_residual_tmp = nuv1 * ((((cb_residual_tmp - 6.1967733539318672 * q[24]) -
    6.08511063404532 * q[29]) - 1.5491933384829668 * lambda_23_2) -
    1.52127765851133 * lambda_23_4);
  Jacobian[748] = db_residual_tmp;
  eb_residual_tmp = 9.16515138991168 * Uphi3_east_3 - 4.58257569495584 *
    Uphi3_west_3;
  Jacobian[718] = nuv1 * (eb_residual_tmp + 4.58257569495584 * lambda_23_3);
  fb_residual_tmp = 7.745966692414834 * Uphi3_east_3 + 3.872983346207417 *
    Uphi3_west_3;
  Jacobian[688] = nuv1 * (fb_residual_tmp - 3.872983346207417 * lambda_23_3);
  gb_residual_tmp = 4.47213595499958 + nuv1 * (((6.0 * Uphi3_east_3 - 3.0 *
    Uphi3_west_3) - 6.9282032302755088 * q[28]) + 3.0 * lambda_23_3);
  Jacobian[658] = gb_residual_tmp;
  hb_residual_tmp = 3.4641016151377544 * Uphi3_east_3 + 1.7320508075688772 *
    Uphi3_west_3;
  ib_residual_tmp = nuv1 * ((hb_residual_tmp - 6.9282032302755088 * q[27]) -
    1.7320508075688772 * lambda_23_3);
  Jacobian[628] = ib_residual_tmp;
  jb_residual_tmp = 1.52127765851133 * Uphi2_east_2 + 1.0327955589886444 * q[19];
  Jacobian[598] = nuv1 * ((jb_residual_tmp + 1.52127765851133 * lambda_23_2) +
    1.0327955589886444 * lambda_23_4);
  kb_residual_tmp = 3.0 * Uphi2_east_1 + 1.9166296949998198 * Uphi2_east_3;
  Jacobian[568] = nuv1 * ((kb_residual_tmp + 3.0 * lambda_23_1) +
    1.9166296949998198 * lambda_23_3);
  lb_residual_tmp = 1.7320508075688772 * Uphi2_east_1 + 1.1065666703449764 *
    Uphi2_east_3;
  Jacobian[538] = nuv1 * ((lb_residual_tmp + 1.7320508075688772 * lambda_23_1) +
    1.1065666703449764 * lambda_23_3);
  mb_residual_tmp = 3.4641016151377544 * Uphi2_east_2 + 3.4016802570830449 * q
    [19];
  Jacobian[508] = nuv1 * ((mb_residual_tmp + 3.4641016151377544 * lambda_23_2) +
    3.4016802570830449 * lambda_23_4);
  nb_residual_tmp = 2.6832815729997477 * Uphi2_east_2 + 2.6349301969610397 * q
    [19];
  ob_residual_tmp = nuv1 * ((nb_residual_tmp + 2.6832815729997477 * lambda_23_2)
    + 2.6349301969610397 * lambda_23_4);
  Jacobian[478] = ob_residual_tmp;
  pb_residual_tmp = 1.5491933384829668 * Uphi2_east_2 + 1.52127765851133 * q[19];
  qb_residual_tmp = nuv1 * ((pb_residual_tmp + 1.5491933384829668 * lambda_23_2)
    + 1.52127765851133 * lambda_23_4);
  Jacobian[448] = qb_residual_tmp;
  rb_residual_tmp = nuv1 * (4.58257569495584 * Uphi2_east_3 + 4.58257569495584 *
    lambda_23_3);
  Jacobian[418] = rb_residual_tmp;
  Jacobian[388] = nuv1 * (3.872983346207417 * Uphi2_east_3 + 3.872983346207417 *
    lambda_23_3);
  sb_residual_tmp = nuv1 * (3.0 * Uphi2_east_3 + 3.0 * lambda_23_3);
  Jacobian[358] = sb_residual_tmp;
  tb_residual_tmp = nuv1 * (1.7320508075688772 * Uphi2_east_3 +
    1.7320508075688772 * lambda_23_3);
  Jacobian[328] = tb_residual_tmp;
  Jacobian[897] = 11.832159566199232 + m_residual_tmp;
  Jacobian[867] = nuv1 * ((l_residual_tmp - 1.7320508075688772 * lambda_23_1) -
    1.1065666703449764 * lambda_23_3);
  Jacobian[837] = 10.0 + nuv1 * (((((2.0 * Uphi3_east_1 + 1.2777531299998799 *
    Uphi3_east_3) - Uphi3_west_1) - 0.63887656499993994 * Uphi3_west_3) +
    lambda_23_1) + 0.63887656499993994 * lambda_23_3);
  l_residual_tmp = ((4.0 * Uphi3_east_2 + 3.927922024247863 * q[29]) - 2.0 *
                    Uphi3_west_2) - 1.9639610121239315 * q[29];
  Jacobian[807] = nuv1 * ((l_residual_tmp + 2.0 * lambda_23_2) +
    1.9639610121239315 * lambda_23_4);
  m_residual_tmp = nuv1 * ((cb_residual_tmp - 1.5491933384829668 * lambda_23_2)
    - 1.52127765851133 * lambda_23_4);
  Jacobian[777] = m_residual_tmp;
  cb_residual_tmp = nuv1 * (((((1.7888543819998317 * Uphi3_east_2 +
    1.7566201313073597 * q[29]) - 0.89442719099991586 * Uphi3_west_2) -
    0.87831006565367986 * q[29]) + 0.89442719099991586 * lambda_23_2) +
    0.87831006565367986 * lambda_23_4);
  Jacobian[747] = -7.745966692414834 + cb_residual_tmp;
  ub_residual_tmp = 5.2915026221291814 * Uphi3_east_3 + 2.6457513110645907 *
    Uphi3_west_3;
  Jacobian[717] = nuv1 * (ub_residual_tmp - 2.6457513110645907 * lambda_23_3);
  vb_residual_tmp = 4.47213595499958 * Uphi3_east_3 - 2.23606797749979 *
    Uphi3_west_3;
  Jacobian[687] = nuv1 * (vb_residual_tmp + 2.23606797749979 * lambda_23_3);
  hb_residual_tmp = nuv1 * (hb_residual_tmp - 1.7320508075688772 * lambda_23_3);
  Jacobian[657] = hb_residual_tmp;
  wb_residual_tmp = 4.47213595499958 + nuv1 * ((2.0 * Uphi3_east_3 -
    Uphi3_west_3) + lambda_23_3);
  Jacobian[627] = wb_residual_tmp;
  Jacobian[597] = u_residual_tmp;
  Jacobian[567] = nuv1 * (((-1.7320508075688772 * Uphi2_east_1 -
    1.1065666703449764 * Uphi2_east_3) - 1.7320508075688772 * lambda_23_1) -
    1.1065666703449764 * lambda_23_3);
  Jacobian[537] = nuv1 * (((-Uphi2_east_1 - 0.63887656499993994 * Uphi2_east_3)
    - lambda_23_1) - 0.63887656499993994 * lambda_23_3);
  u_residual_tmp = nuv1 * (((-2.0 * Uphi2_east_2 - 1.9639610121239315 * q[19]) -
    2.0 * lambda_23_2) - 1.9639610121239315 * lambda_23_4);
  Jacobian[507] = u_residual_tmp;
  xb_residual_tmp = nuv1 * (((-1.5491933384829668 * Uphi2_east_2 -
    1.52127765851133 * q[19]) - 1.5491933384829668 * lambda_23_2) -
    1.52127765851133 * lambda_23_4);
  Jacobian[477] = xb_residual_tmp;
  yb_residual_tmp = nuv1 * (((-0.89442719099991586 * Uphi2_east_2 -
    0.87831006565367986 * q[19]) - 0.89442719099991586 * lambda_23_2) -
    0.87831006565367986 * lambda_23_4);
  Jacobian[447] = yb_residual_tmp;
  Jacobian[417] = nuv1 * (-2.6457513110645907 * Uphi2_east_3 -
    2.6457513110645907 * lambda_23_3);
  ac_residual_tmp = nuv1 * (-2.23606797749979 * Uphi2_east_3 - 2.23606797749979 *
    lambda_23_3);
  Jacobian[387] = ac_residual_tmp;
  bc_residual_tmp = nuv1 * (-1.7320508075688772 * Uphi2_east_3 -
    1.7320508075688772 * lambda_23_3);
  Jacobian[357] = bc_residual_tmp;
  cc_residual_tmp = nuv1 * (-Uphi2_east_3 - lambda_23_3);
  Jacobian[327] = cc_residual_tmp;
  Jacobian[896] = nuv1 * ((n_residual_tmp - 13.60672102833218 * q[28]) +
    1.9639610121239315 * lambda_23_3);
  Jacobian[866] = nuv1 * (((((ab_residual_tmp - 13.856406460551018 * q[24]) -
    12.393546707863734 * q[26]) - 13.60672102833218 * q[29]) -
    3.4641016151377544 * lambda_23_2) - 3.4016802570830449 * lambda_23_4);
  Jacobian[836] = nuv1 * (((l_residual_tmp - 13.856406460551018 * q[25]) + 2.0 *
    lambda_23_2) + 1.9639610121239315 * lambda_23_4);
  Jacobian[806] = 6.0 + nuv1 * ((((((((10.0 * Uphi3_east_1 + 8.94427190999916 *
    Uphi3_east_3) - 5.0 * Uphi3_west_1) - 4.47213595499958 * Uphi3_west_3) -
    13.856406460551018 * q[21]) - 13.60672102833218 * q[23]) -
    12.393546707863734 * q[28]) + 5.0 * lambda_23_1) + 4.47213595499958 *
    lambda_23_3);
  l_residual_tmp = ((7.745966692414834 * Uphi3_east_1 + 6.9282032302755088 *
                     Uphi3_east_3) + 3.872983346207417 * Uphi3_west_1) +
    3.4641016151377544 * Uphi3_west_3;
  Jacobian[776] = nuv1 * (((((l_residual_tmp - 15.491933384829668 * q[20]) -
    13.856406460551018 * q[22]) - 13.856406460551018 * q[27]) -
    3.872983346207417 * lambda_23_1) - 3.4641016151377544 * lambda_23_3);
  n_residual_tmp = ((4.47213595499958 * Uphi3_east_1 + 4.0 * Uphi3_east_3) -
                    2.23606797749979 * Uphi3_west_1) - 2.0 * Uphi3_west_3;
  Jacobian[746] = nuv1 * ((((n_residual_tmp - 15.491933384829668 * q[21]) -
    13.856406460551018 * q[28]) + 2.23606797749979 * lambda_23_1) + 2.0 *
    lambda_23_3);
  ab_residual_tmp = 11.832159566199232 * Uphi3_east_2 + 5.9160797830996161 *
    Uphi3_west_2;
  Jacobian[716] = nuv1 * ((ab_residual_tmp - 13.60672102833218 * q[26]) -
    5.9160797830996161 * lambda_23_2);
  dc_residual_tmp = nuv1 * (((10.0 * Uphi3_east_2 - 5.0 * Uphi3_west_2) -
    13.856406460551018 * q[25]) + 5.0 * lambda_23_2);
  Jacobian[686] = -3.4641016151377544 + dc_residual_tmp;
  ec_residual_tmp = 7.745966692414834 * Uphi3_east_2 + 3.872983346207417 *
    Uphi3_west_2;
  fc_residual_tmp = nuv1 * (((ec_residual_tmp - 15.491933384829668 * q[24]) -
    13.856406460551018 * q[26]) - 3.872983346207417 * lambda_23_2);
  Jacobian[656] = fc_residual_tmp;
  gc_residual_tmp = 4.47213595499958 * Uphi3_east_2 - 2.23606797749979 *
    Uphi3_west_2;
  hc_residual_tmp = nuv1 * ((gc_residual_tmp - 15.491933384829668 * q[25]) +
    2.23606797749979 * lambda_23_2);
  Jacobian[626] = hc_residual_tmp;
  Jacobian[596] = v_residual_tmp;
  Jacobian[566] = nuv1 * (((-3.4641016151377544 * Uphi2_east_2 -
    3.4016802570830449 * q[19]) - 3.4641016151377544 * lambda_23_2) -
    3.4016802570830449 * lambda_23_4);
  Jacobian[536] = u_residual_tmp;
  Jacobian[506] = nuv1 * (((-5.0 * Uphi2_east_1 - 4.47213595499958 *
    Uphi2_east_3) - 5.0 * lambda_23_1) - 4.47213595499958 * lambda_23_3);
  Jacobian[476] = nuv1 * (((-3.872983346207417 * Uphi2_east_1 -
    3.4641016151377544 * Uphi2_east_3) - 3.872983346207417 * lambda_23_1) -
    3.4641016151377544 * lambda_23_3);
  u_residual_tmp = nuv1 * (((-2.23606797749979 * Uphi2_east_1 - 2.0 *
    Uphi2_east_3) - 2.23606797749979 * lambda_23_1) - 2.0 * lambda_23_3);
  Jacobian[446] = u_residual_tmp;
  Jacobian[416] = nuv1 * (-5.9160797830996161 * Uphi2_east_2 -
    5.9160797830996161 * lambda_23_2);
  v_residual_tmp = nuv1 * (-5.0 * Uphi2_east_2 - 5.0 * lambda_23_2);
  Jacobian[386] = v_residual_tmp;
  ic_residual_tmp = nuv1 * (-3.872983346207417 * Uphi2_east_2 -
    3.872983346207417 * lambda_23_2);
  Jacobian[356] = ic_residual_tmp;
  jc_residual_tmp = nuv1 * (-2.23606797749979 * Uphi2_east_2 - 2.23606797749979 *
    lambda_23_2);
  Jacobian[326] = jc_residual_tmp;
  Jacobian[895] = nuv1 * ((o_residual_tmp - 6.08511063404532 * q[27]) -
    1.52127765851133 * lambda_23_3);
  Jacobian[865] = 7.745966692414834 + bb_residual_tmp;
  Jacobian[835] = db_residual_tmp;
  Jacobian[805] = nuv1 * (((l_residual_tmp - 6.9282032302755088 * q[22]) -
    3.872983346207417 * lambda_23_1) - 3.4641016151377544 * lambda_23_3);
  Jacobian[775] = 6.0 + nuv1 * (((((((6.0 * Uphi3_east_1 + 5.3665631459994954 *
    Uphi3_east_3) - 3.0 * Uphi3_west_1) - 2.6832815729997477 * Uphi3_west_3) -
    6.9282032302755088 * q[21]) - 6.1967733539318672 * q[28]) + 3.0 *
    lambda_23_1) + 2.6832815729997477 * lambda_23_3);
  l_residual_tmp = ((3.4641016151377544 * Uphi3_east_1 + 3.0983866769659336 *
                     Uphi3_east_3) + 1.7320508075688772 * Uphi3_west_1) +
    1.5491933384829668 * Uphi3_west_3;
  Jacobian[745] = nuv1 * ((((l_residual_tmp - 6.9282032302755088 * q[20]) -
    6.1967733539318672 * q[27]) - 1.7320508075688772 * lambda_23_1) -
    1.5491933384829668 * lambda_23_3);
  o_residual_tmp = 9.16515138991168 * Uphi3_east_2 - 4.58257569495584 *
    Uphi3_west_2;
  Jacobian[715] = nuv1 * (o_residual_tmp + 4.58257569495584 * lambda_23_2);
  bb_residual_tmp = nuv1 * ((ec_residual_tmp - 6.9282032302755088 * q[26]) -
    3.872983346207417 * lambda_23_2);
  Jacobian[685] = bb_residual_tmp;
  db_residual_tmp = nuv1 * (((6.0 * Uphi3_east_2 - 3.0 * Uphi3_west_2) -
    6.9282032302755088 * q[25]) + 3.0 * lambda_23_2);
  Jacobian[655] = -3.4641016151377544 + db_residual_tmp;
  ec_residual_tmp = 3.4641016151377544 * Uphi3_east_2 + 1.7320508075688772 *
    Uphi3_west_2;
  kc_residual_tmp = nuv1 * ((ec_residual_tmp - 6.9282032302755088 * q[24]) -
    1.7320508075688772 * lambda_23_2);
  Jacobian[625] = kc_residual_tmp;
  Jacobian[595] = nuv1 * (1.52127765851133 * Uphi2_east_3 + 1.52127765851133 *
    lambda_23_3);
  Jacobian[565] = ob_residual_tmp;
  Jacobian[535] = qb_residual_tmp;
  ob_residual_tmp = 3.872983346207417 * Uphi2_east_1 + 3.4641016151377544 *
    Uphi2_east_3;
  Jacobian[505] = nuv1 * ((ob_residual_tmp + 3.872983346207417 * lambda_23_1) +
    3.4641016151377544 * lambda_23_3);
  qb_residual_tmp = 3.0 * Uphi2_east_1 + 2.6832815729997477 * Uphi2_east_3;
  Jacobian[475] = nuv1 * ((qb_residual_tmp + 3.0 * lambda_23_1) +
    2.6832815729997477 * lambda_23_3);
  lc_residual_tmp = 1.7320508075688772 * Uphi2_east_1 + 1.5491933384829668 *
    Uphi2_east_3;
  Jacobian[445] = nuv1 * ((lc_residual_tmp + 1.7320508075688772 * lambda_23_1) +
    1.5491933384829668 * lambda_23_3);
  mc_residual_tmp = nuv1 * (4.58257569495584 * Uphi2_east_2 + 4.58257569495584 *
    lambda_23_2);
  Jacobian[415] = mc_residual_tmp;
  nc_residual_tmp = nuv1 * (3.872983346207417 * Uphi2_east_2 + 3.872983346207417
    * lambda_23_2);
  Jacobian[385] = nc_residual_tmp;
  oc_residual_tmp = nuv1 * (3.0 * Uphi2_east_2 + 3.0 * lambda_23_2);
  Jacobian[355] = oc_residual_tmp;
  pc_residual_tmp = nuv1 * (1.7320508075688772 * Uphi2_east_2 +
    1.7320508075688772 * lambda_23_2);
  Jacobian[325] = pc_residual_tmp;
  Jacobian[894] = p_residual_tmp;
  Jacobian[864] = m_residual_tmp;
  Jacobian[834] = 7.745966692414834 + cb_residual_tmp;
  Jacobian[804] = nuv1 * ((n_residual_tmp + 2.23606797749979 * lambda_23_1) +
    2.0 * lambda_23_3);
  Jacobian[774] = nuv1 * ((l_residual_tmp - 1.7320508075688772 * lambda_23_1) -
    1.5491933384829668 * lambda_23_3);
  Jacobian[744] = 6.0 + nuv1 * (((((2.0 * Uphi3_east_1 + 1.7888543819998317 *
    Uphi3_east_3) - Uphi3_west_1) - 0.89442719099991586 * Uphi3_west_3) +
    lambda_23_1) + 0.89442719099991586 * lambda_23_3);
  l_residual_tmp = 5.2915026221291814 * Uphi3_east_2 + 2.6457513110645907 *
    Uphi3_west_2;
  Jacobian[714] = nuv1 * (l_residual_tmp - 2.6457513110645907 * lambda_23_2);
  m_residual_tmp = nuv1 * (gc_residual_tmp + 2.23606797749979 * lambda_23_2);
  Jacobian[684] = m_residual_tmp;
  n_residual_tmp = nuv1 * (ec_residual_tmp - 1.7320508075688772 * lambda_23_2);
  Jacobian[654] = n_residual_tmp;
  p_residual_tmp = nuv1 * ((2.0 * Uphi3_east_2 - Uphi3_west_2) + lambda_23_2);
  Jacobian[624] = -3.4641016151377544 + p_residual_tmp;
  Jacobian[594] = w_residual_tmp;
  Jacobian[564] = xb_residual_tmp;
  Jacobian[534] = yb_residual_tmp;
  Jacobian[504] = u_residual_tmp;
  Jacobian[474] = nuv1 * (((-1.7320508075688772 * Uphi2_east_1 -
    1.5491933384829668 * Uphi2_east_3) - 1.7320508075688772 * lambda_23_1) -
    1.5491933384829668 * lambda_23_3);
  Jacobian[444] = nuv1 * (((-Uphi2_east_1 - 0.89442719099991586 * Uphi2_east_3)
    - lambda_23_1) - 0.89442719099991586 * lambda_23_3);
  Jacobian[414] = nuv1 * (-2.6457513110645907 * Uphi2_east_2 -
    2.6457513110645907 * lambda_23_2);
  Jacobian[384] = jc_residual_tmp;
  u_residual_tmp = nuv1 * (-1.7320508075688772 * Uphi2_east_2 -
    1.7320508075688772 * lambda_23_2);
  Jacobian[354] = u_residual_tmp;
  w_residual_tmp = nuv1 * (-Uphi2_east_2 - lambda_23_2);
  Jacobian[324] = w_residual_tmp;
  Jacobian[893] = nuv1 * ((q_residual_tmp - 10.583005244258363 * q[29]) -
    2.6457513110645907 * lambda_23_4);
  Jacobian[863] = nuv1 * ((eb_residual_tmp - 31.749015732775089 * q[28]) +
    4.58257569495584 * lambda_23_3);
  Jacobian[833] = nuv1 * ((ub_residual_tmp - 10.583005244258363 * q[27]) -
    2.6457513110645907 * lambda_23_3);
  Jacobian[803] = nuv1 * (((ab_residual_tmp - 23.664319132398465 * q[24]) -
    25.701584164627452 * q[26]) - 5.9160797830996161 * lambda_23_2);
  Jacobian[773] = nuv1 * ((o_residual_tmp - 31.749015732775089 * q[25]) +
    4.58257569495584 * lambda_23_2);
  Jacobian[743] = nuv1 * (((l_residual_tmp - 10.583005244258363 * q[24]) -
    23.664319132398465 * q[26]) - 2.6457513110645907 * lambda_23_2);
  Jacobian[713] = 2.0 + nuv1 * ((((14.0 * Uphi3_east_1 - 7.0 * Uphi3_west_1) -
    20.784609690826528 * q[21]) - 24.693678903269511 * q[23]) + 7.0 *
    lambda_23_1);
  l_residual_tmp = 11.832159566199232 * Uphi3_east_1 + 5.9160797830996161 *
    Uphi3_west_1;
  Jacobian[683] = nuv1 * (((l_residual_tmp - 23.664319132398465 * q[20]) -
    25.701584164627452 * q[22]) - 5.9160797830996161 * lambda_23_1);
  o_residual_tmp = 9.16515138991168 * Uphi3_east_1 - 4.58257569495584 *
    Uphi3_west_1;
  Jacobian[653] = nuv1 * (((o_residual_tmp - 31.749015732775089 * q[21]) -
    20.784609690826528 * q[23]) + 4.58257569495584 * lambda_23_1);
  q_residual_tmp = 5.2915026221291814 * Uphi3_east_1 + 2.6457513110645907 *
    Uphi3_west_1;
  Jacobian[623] = nuv1 * (((q_residual_tmp - 10.583005244258363 * q[20]) -
    23.664319132398465 * q[22]) - 2.6457513110645907 * lambda_23_1);
  Jacobian[593] = nuv1 * (2.6457513110645907 * q[19] + 2.6457513110645907 *
    lambda_23_4);
  Jacobian[563] = rb_residual_tmp;
  Jacobian[533] = nuv1 * (2.6457513110645907 * Uphi2_east_3 + 2.6457513110645907
    * lambda_23_3);
  Jacobian[503] = nuv1 * (5.9160797830996161 * Uphi2_east_2 + 5.9160797830996161
    * lambda_23_2);
  Jacobian[473] = mc_residual_tmp;
  Jacobian[443] = nuv1 * (2.6457513110645907 * Uphi2_east_2 + 2.6457513110645907
    * lambda_23_2);
  Jacobian[413] = nuv1 * (7.0 * Uphi2_east_1 + 7.0 * lambda_23_1);
  Jacobian[383] = nuv1 * (5.9160797830996161 * Uphi2_east_1 + 5.9160797830996161
    * lambda_23_1);
  ab_residual_tmp = nuv1 * (4.58257569495584 * Uphi2_east_1 + 4.58257569495584 *
    lambda_23_1);
  Jacobian[353] = ab_residual_tmp;
  Jacobian[323] = nuv1 * (2.6457513110645907 * Uphi2_east_1 + 2.6457513110645907
    * lambda_23_1);
  Jacobian[892] = r_residual_tmp;
  Jacobian[862] = nuv1 * ((fb_residual_tmp - 15.491933384829668 * q[27]) -
    3.872983346207417 * lambda_23_3);
  Jacobian[832] = nuv1 * ((vb_residual_tmp - 15.491933384829668 * q[28]) +
    2.23606797749979 * lambda_23_3);
  Jacobian[802] = 3.4641016151377544 + dc_residual_tmp;
  Jacobian[772] = fc_residual_tmp;
  Jacobian[742] = hc_residual_tmp;
  Jacobian[712] = nuv1 * ((l_residual_tmp - 13.60672102833218 * q[22]) -
    5.9160797830996161 * lambda_23_1);
  Jacobian[682] = 2.0 + nuv1 * ((((10.0 * Uphi3_east_1 - 5.0 * Uphi3_west_1) -
    13.856406460551018 * q[21]) - 13.60672102833218 * q[23]) + 5.0 * lambda_23_1);
  l_residual_tmp = 7.745966692414834 * Uphi3_east_1 + 3.872983346207417 *
    Uphi3_west_1;
  Jacobian[652] = nuv1 * (((l_residual_tmp - 15.491933384829668 * q[20]) -
    13.856406460551018 * q[22]) - 3.872983346207417 * lambda_23_1);
  r_residual_tmp = 4.47213595499958 * Uphi3_east_1 - 2.23606797749979 *
    Uphi3_west_1;
  Jacobian[622] = nuv1 * ((r_residual_tmp - 15.491933384829668 * q[21]) +
    2.23606797749979 * lambda_23_1);
  Jacobian[592] = x_residual_tmp;
  Jacobian[562] = nuv1 * (-3.872983346207417 * Uphi2_east_3 - 3.872983346207417 *
    lambda_23_3);
  Jacobian[532] = ac_residual_tmp;
  Jacobian[502] = v_residual_tmp;
  Jacobian[472] = ic_residual_tmp;
  Jacobian[442] = jc_residual_tmp;
  Jacobian[412] = nuv1 * (-5.9160797830996161 * Uphi2_east_1 -
    5.9160797830996161 * lambda_23_1);
  Jacobian[382] = nuv1 * (-5.0 * Uphi2_east_1 - 5.0 * lambda_23_1);
  Jacobian[352] = nuv1 * (-3.872983346207417 * Uphi2_east_1 - 3.872983346207417 *
    lambda_23_1);
  v_residual_tmp = nuv1 * (-2.23606797749979 * Uphi2_east_1 - 2.23606797749979 *
    lambda_23_1);
  Jacobian[322] = v_residual_tmp;
  Jacobian[891] = nuv1 * ((s_residual_tmp - 6.9282032302755088 * q[29]) -
    1.7320508075688772 * lambda_23_4);
  Jacobian[861] = gb_residual_tmp;
  Jacobian[831] = ib_residual_tmp;
  Jacobian[801] = bb_residual_tmp;
  Jacobian[771] = 3.4641016151377544 + db_residual_tmp;
  Jacobian[741] = kc_residual_tmp;
  Jacobian[711] = nuv1 * ((o_residual_tmp - 6.9282032302755088 * q[23]) +
    4.58257569495584 * lambda_23_1);
  Jacobian[681] = nuv1 * ((l_residual_tmp - 6.9282032302755088 * q[22]) -
    3.872983346207417 * lambda_23_1);
  Jacobian[651] = 2.0 + nuv1 * (((6.0 * Uphi3_east_1 - 3.0 * Uphi3_west_1) -
    6.9282032302755088 * q[21]) + 3.0 * lambda_23_1);
  l_residual_tmp = 3.4641016151377544 * Uphi3_east_1 + 1.7320508075688772 *
    Uphi3_west_1;
  Jacobian[621] = nuv1 * ((l_residual_tmp - 6.9282032302755088 * q[20]) -
    1.7320508075688772 * lambda_23_1);
  Jacobian[591] = nuv1 * (1.7320508075688772 * q[19] + 1.7320508075688772 *
    lambda_23_4);
  Jacobian[561] = sb_residual_tmp;
  Jacobian[531] = tb_residual_tmp;
  Jacobian[501] = nc_residual_tmp;
  Jacobian[471] = oc_residual_tmp;
  Jacobian[441] = pc_residual_tmp;
  Jacobian[411] = ab_residual_tmp;
  Jacobian[381] = nuv1 * (3.872983346207417 * Uphi2_east_1 + 3.872983346207417 *
    lambda_23_1);
  Jacobian[351] = nuv1 * (3.0 * Uphi2_east_1 + 3.0 * lambda_23_1);
  Jacobian[321] = nuv1 * (1.7320508075688772 * Uphi2_east_1 + 1.7320508075688772
    * lambda_23_1);
  Jacobian[890] = 5.2915026221291814 + t_residual_tmp;
  Jacobian[860] = hb_residual_tmp;
  Jacobian[830] = wb_residual_tmp;
  Jacobian[800] = m_residual_tmp;
  Jacobian[770] = n_residual_tmp;
  Jacobian[740] = 3.4641016151377544 + p_residual_tmp;
  Jacobian[710] = nuv1 * (q_residual_tmp - 2.6457513110645907 * lambda_23_1);
  Jacobian[680] = nuv1 * (r_residual_tmp + 2.23606797749979 * lambda_23_1);
  Jacobian[650] = nuv1 * (l_residual_tmp - 1.7320508075688772 * lambda_23_1);
  Jacobian[620] = 2.0 + nuv1 * ((2.0 * Uphi3_east_1 - Uphi3_west_1) +
    lambda_23_1);
  Jacobian[590] = y_residual_tmp;
  Jacobian[560] = bc_residual_tmp;
  Jacobian[530] = cc_residual_tmp;
  Jacobian[500] = jc_residual_tmp;
  Jacobian[470] = u_residual_tmp;
  Jacobian[440] = w_residual_tmp;
  Jacobian[410] = nuv1 * (-2.6457513110645907 * Uphi2_east_1 -
    2.6457513110645907 * lambda_23_1);
  Jacobian[380] = v_residual_tmp;
  Jacobian[350] = nuv1 * (-1.7320508075688772 * Uphi2_east_1 -
    1.7320508075688772 * lambda_23_1);
  Jacobian[320] = nuv1 * (-Uphi2_east_1 - lambda_23_1);
  Jacobian[889] = nuv1 * (((Uphi3_west_1 + 0.59628479399994394 * Uphi3_west_3) -
    lambda_23_1) - 0.59628479399994394 * lambda_23_3);
  Jacobian[859] = nuv1 * (((1.52127765851133 * lambda_23_2 - 1.0327955589886444 *
    q[29]) - 1.52127765851133 * Uphi3_west_2) + 1.0327955589886444 * lambda_23_4);
  l_residual_tmp = nuv1 * (((0.87831006565367986 * Uphi3_west_2 +
    0.59628479399994394 * q[29]) - 0.87831006565367986 * lambda_23_2) -
    0.59628479399994394 * lambda_23_4);
  Jacobian[829] = l_residual_tmp;
  m_residual_tmp = nuv1 * (1.9639610121239315 * Uphi3_west_3 -
    1.9639610121239315 * lambda_23_3);
  Jacobian[799] = m_residual_tmp;
  Jacobian[769] = nuv1 * (1.52127765851133 * lambda_23_3 - 1.52127765851133 *
    Uphi3_west_3);
  n_residual_tmp = nuv1 * (0.87831006565367986 * Uphi3_west_3 -
    0.87831006565367986 * lambda_23_3);
  Jacobian[739] = n_residual_tmp;
  Jacobian[709] = nuv1 * (2.6457513110645907 * lambda_23_4 - 2.6457513110645907 *
    q[29]);
  o_residual_tmp = nuv1 * (2.23606797749979 * q[29] - 2.23606797749979 *
    lambda_23_4);
  Jacobian[679] = o_residual_tmp;
  Jacobian[649] = nuv1 * (1.7320508075688772 * lambda_23_4 - 1.7320508075688772 *
    q[29]);
  p_residual_tmp = nuv1 * (q[29] - lambda_23_4);
  Jacobian[619] = p_residual_tmp;
  Jacobian[589] = 14.0 + nuv1 * (((((((Uphi2_east_1 + 0.59628479399994394 *
    Uphi2_east_3) - Uphi2_west_1) - 0.59628479399994394 * Uphi2_west_3) +
    lambda_12_1) + 0.59628479399994394 * lambda_12_3) + lambda_23_1) +
    0.59628479399994394 * lambda_23_3);
  q_residual_tmp = (jb_residual_tmp + 1.52127765851133 * Uphi2_west_2) +
    1.0327955589886444 * q[19];
  Jacobian[559] = nuv1 * ((((q_residual_tmp - 1.52127765851133 * lambda_12_2) -
    1.0327955589886444 * lambda_12_4) + 1.52127765851133 * lambda_23_2) +
    1.0327955589886444 * lambda_23_4);
  r_residual_tmp = nuv1 * (((((((0.87831006565367986 * Uphi2_east_2 +
    0.59628479399994394 * q[19]) - 0.87831006565367986 * Uphi2_west_2) -
    0.59628479399994394 * q[19]) + 0.87831006565367986 * lambda_12_2) +
    0.59628479399994394 * lambda_12_4) + 0.87831006565367986 * lambda_23_2) +
    0.59628479399994394 * lambda_23_4);
  Jacobian[529] = -11.832159566199232 + r_residual_tmp;
  s_residual_tmp = 1.9639610121239315 * Uphi2_east_3 - 1.9639610121239315 *
    Uphi2_west_3;
  Jacobian[499] = nuv1 * ((s_residual_tmp + 1.9639610121239315 * lambda_12_3) +
    1.9639610121239315 * lambda_23_3);
  t_residual_tmp = 1.52127765851133 * Uphi2_east_3 + 1.52127765851133 *
    Uphi2_west_3;
  Jacobian[469] = nuv1 * ((t_residual_tmp - 1.52127765851133 * lambda_12_3) +
    1.52127765851133 * lambda_23_3);
  u_residual_tmp = 9.16515138991168 + nuv1 * (((0.87831006565367986 *
    Uphi2_east_3 - 0.87831006565367986 * Uphi2_west_3) + 0.87831006565367986 *
    lambda_12_3) + 0.87831006565367986 * lambda_23_3);
  Jacobian[439] = u_residual_tmp;
  v_residual_tmp = 2.6457513110645907 * q[19] + 2.6457513110645907 * q[19];
  Jacobian[409] = nuv1 * ((v_residual_tmp - 2.6457513110645907 * lambda_12_4) +
    2.6457513110645907 * lambda_23_4);
  w_residual_tmp = nuv1 * (((2.23606797749979 * q[19] - 2.23606797749979 * q[19])
    + 2.23606797749979 * lambda_12_4) + 2.23606797749979 * lambda_23_4);
  Jacobian[379] = w_residual_tmp;
  x_residual_tmp = 1.7320508075688772 * q[19] + 1.7320508075688772 * q[19];
  Jacobian[349] = nuv1 * ((x_residual_tmp - 1.7320508075688772 * lambda_12_4) +
    1.7320508075688772 * lambda_23_4);
  y_residual_tmp = nuv1 * (((q[19] - q[19]) + lambda_12_4) + lambda_23_4);
  Jacobian[319] = -5.2915026221291814 + y_residual_tmp;
  Jacobian[289] = nuv1 * (((-Uphi1_east_1 - 0.59628479399994394 * Uphi1_east_3)
    - lambda_12_1) - 0.59628479399994394 * lambda_12_3);
  Jacobian[259] = nuv1 * (((-1.52127765851133 * Uphi1_east_2 -
    1.0327955589886444 * q[9]) - 1.52127765851133 * lambda_12_2) -
    1.0327955589886444 * lambda_12_4);
  ab_residual_tmp = nuv1 * (((-0.87831006565367986 * Uphi1_east_2 -
    0.59628479399994394 * q[9]) - 0.87831006565367986 * lambda_12_2) -
    0.59628479399994394 * lambda_12_4);
  Jacobian[229] = ab_residual_tmp;
  bb_residual_tmp = nuv1 * (-1.9639610121239315 * Uphi1_east_3 -
    1.9639610121239315 * lambda_12_3);
  Jacobian[199] = bb_residual_tmp;
  Jacobian[169] = nuv1 * (-1.52127765851133 * Uphi1_east_3 - 1.52127765851133 *
    lambda_12_3);
  cb_residual_tmp = nuv1 * (-0.87831006565367986 * Uphi1_east_3 -
    0.87831006565367986 * lambda_12_3);
  Jacobian[139] = cb_residual_tmp;
  Jacobian[109] = nuv1 * (-2.6457513110645907 * q[9] - 2.6457513110645907 *
    lambda_12_4);
  db_residual_tmp = nuv1 * (-2.23606797749979 * q[9] - 2.23606797749979 *
    lambda_12_4);
  Jacobian[79] = db_residual_tmp;
  Jacobian[49] = nuv1 * (-1.7320508075688772 * q[9] - 1.7320508075688772 *
    lambda_12_4);
  eb_residual_tmp = nuv1 * (-q[9] - lambda_12_4);
  Jacobian[19] = eb_residual_tmp;
  Jacobian[888] = nuv1 * (((1.52127765851133 * Uphi3_west_2 + 1.0327955589886444
    * q[29]) - 1.52127765851133 * lambda_23_2) - 1.0327955589886444 *
    lambda_23_4);
  Jacobian[858] = nuv1 * (((3.0 * lambda_23_1 - 1.9166296949998198 *
    Uphi3_west_3) - 3.0 * Uphi3_west_1) + 1.9166296949998198 * lambda_23_3);
  Jacobian[828] = nuv1 * (((1.7320508075688772 * Uphi3_west_1 +
    1.1065666703449764 * Uphi3_west_3) - 1.7320508075688772 * lambda_23_1) -
    1.1065666703449764 * lambda_23_3);
  Jacobian[798] = nuv1 * (((3.4641016151377544 * Uphi3_west_2 +
    3.4016802570830449 * q[29]) - 3.4641016151377544 * lambda_23_2) -
    3.4016802570830449 * lambda_23_4);
  fb_residual_tmp = nuv1 * (((2.6832815729997477 * lambda_23_2 -
    2.6349301969610397 * q[29]) - 2.6832815729997477 * Uphi3_west_2) +
    2.6349301969610397 * lambda_23_4);
  Jacobian[768] = fb_residual_tmp;
  gb_residual_tmp = nuv1 * (((1.5491933384829668 * Uphi3_west_2 +
    1.52127765851133 * q[29]) - 1.5491933384829668 * lambda_23_2) -
    1.52127765851133 * lambda_23_4);
  Jacobian[738] = gb_residual_tmp;
  hb_residual_tmp = nuv1 * (4.58257569495584 * lambda_23_3 - 4.58257569495584 *
    Uphi3_west_3);
  Jacobian[708] = hb_residual_tmp;
  Jacobian[678] = nuv1 * (3.872983346207417 * Uphi3_west_3 - 3.872983346207417 *
    lambda_23_3);
  ib_residual_tmp = nuv1 * (3.0 * lambda_23_3 - 3.0 * Uphi3_west_3);
  Jacobian[648] = ib_residual_tmp;
  jb_residual_tmp = nuv1 * (1.7320508075688772 * Uphi3_west_3 -
    1.7320508075688772 * lambda_23_3);
  Jacobian[618] = jb_residual_tmp;
  Jacobian[588] = nuv1 * ((((((q_residual_tmp - 6.08511063404532 * q[14]) -
    4.1311822359545785 * q[19]) - 1.52127765851133 * lambda_12_2) -
    1.0327955589886444 * lambda_12_4) + 1.52127765851133 * lambda_23_2) +
    1.0327955589886444 * lambda_23_4);
  Jacobian[558] = 10.0 + nuv1 * ((((((((kb_residual_tmp - 3.0 * Uphi2_west_1) -
    1.9166296949998198 * Uphi2_west_3) - 6.9282032302755088 * q[11]) -
    4.4262666813799054 * q[18]) + 3.0 * lambda_12_1) + 1.9166296949998198 *
    lambda_12_3) + 3.0 * lambda_23_1) + 1.9166296949998198 * lambda_23_3);
  q_residual_tmp = (lb_residual_tmp + 1.7320508075688772 * Uphi2_west_1) +
    1.1065666703449764 * Uphi2_west_3;
  Jacobian[528] = nuv1 * ((((((q_residual_tmp - 6.9282032302755088 * q[10]) -
    4.4262666813799054 * q[17]) - 1.7320508075688772 * lambda_12_1) -
    1.1065666703449764 * lambda_12_3) + 1.7320508075688772 * lambda_23_1) +
    1.1065666703449764 * lambda_23_3);
  kb_residual_tmp = (mb_residual_tmp + 3.4641016151377544 * Uphi2_west_2) +
    3.4016802570830449 * q[19];
  Jacobian[498] = nuv1 * (((((kb_residual_tmp - 6.1967733539318672 * q[16]) -
    3.4641016151377544 * lambda_12_2) - 3.4016802570830449 * lambda_12_4) +
    3.4641016151377544 * lambda_23_2) + 3.4016802570830449 * lambda_23_4);
  lb_residual_tmp = nuv1 * (((((((nb_residual_tmp - 2.6832815729997477 *
    Uphi2_west_2) - 2.6349301969610397 * q[19]) - 6.1967733539318672 * q[15]) +
    2.6832815729997477 * lambda_12_2) + 2.6349301969610397 * lambda_12_4) +
    2.6832815729997477 * lambda_23_2) + 2.6349301969610397 * lambda_23_4);
  Jacobian[468] = -7.745966692414834 + lb_residual_tmp;
  mb_residual_tmp = (pb_residual_tmp + 1.5491933384829668 * Uphi2_west_2) +
    1.52127765851133 * q[19];
  rb_residual_tmp = nuv1 * ((((((mb_residual_tmp - 6.1967733539318672 * q[14]) -
    6.08511063404532 * q[19]) - 1.5491933384829668 * lambda_12_2) -
    1.52127765851133 * lambda_12_4) + 1.5491933384829668 * lambda_23_2) +
    1.52127765851133 * lambda_23_4);
  Jacobian[438] = rb_residual_tmp;
  sb_residual_tmp = 4.58257569495584 * Uphi2_east_3 - 4.58257569495584 *
    Uphi2_west_3;
  Jacobian[408] = nuv1 * ((sb_residual_tmp + 4.58257569495584 * lambda_12_3) +
    4.58257569495584 * lambda_23_3);
  tb_residual_tmp = 3.872983346207417 * Uphi2_east_3 + 3.872983346207417 *
    Uphi2_west_3;
  Jacobian[378] = nuv1 * ((tb_residual_tmp - 3.872983346207417 * lambda_12_3) +
    3.872983346207417 * lambda_23_3);
  ub_residual_tmp = 4.47213595499958 + nuv1 * ((((3.0 * Uphi2_east_3 - 3.0 *
    Uphi2_west_3) - 6.9282032302755088 * q[18]) + 3.0 * lambda_12_3) + 3.0 *
    lambda_23_3);
  Jacobian[348] = ub_residual_tmp;
  vb_residual_tmp = 1.7320508075688772 * Uphi2_east_3 + 1.7320508075688772 *
    Uphi2_west_3;
  wb_residual_tmp = nuv1 * (((vb_residual_tmp - 6.9282032302755088 * q[17]) -
    1.7320508075688772 * lambda_12_3) + 1.7320508075688772 * lambda_23_3);
  Jacobian[318] = wb_residual_tmp;
  xb_residual_tmp = 1.52127765851133 * Uphi1_east_2 + 1.0327955589886444 * q[9];
  Jacobian[288] = nuv1 * ((xb_residual_tmp + 1.52127765851133 * lambda_12_2) +
    1.0327955589886444 * lambda_12_4);
  yb_residual_tmp = 3.0 * Uphi1_east_1 + 1.9166296949998198 * Uphi1_east_3;
  Jacobian[258] = nuv1 * ((yb_residual_tmp + 3.0 * lambda_12_1) +
    1.9166296949998198 * lambda_12_3);
  ac_residual_tmp = 1.7320508075688772 * Uphi1_east_1 + 1.1065666703449764 *
    Uphi1_east_3;
  Jacobian[228] = nuv1 * ((ac_residual_tmp + 1.7320508075688772 * lambda_12_1) +
    1.1065666703449764 * lambda_12_3);
  bc_residual_tmp = 3.4641016151377544 * Uphi1_east_2 + 3.4016802570830449 * q[9];
  Jacobian[198] = nuv1 * ((bc_residual_tmp + 3.4641016151377544 * lambda_12_2) +
    3.4016802570830449 * lambda_12_4);
  nb_residual_tmp = 2.6832815729997477 * Uphi1_east_2 + 2.6349301969610397 * q[9];
  cc_residual_tmp = nuv1 * ((nb_residual_tmp + 2.6832815729997477 * lambda_12_2)
    + 2.6349301969610397 * lambda_12_4);
  Jacobian[168] = cc_residual_tmp;
  pb_residual_tmp = 1.5491933384829668 * Uphi1_east_2 + 1.52127765851133 * q[9];
  dc_residual_tmp = nuv1 * ((pb_residual_tmp + 1.5491933384829668 * lambda_12_2)
    + 1.52127765851133 * lambda_12_4);
  Jacobian[138] = dc_residual_tmp;
  ec_residual_tmp = nuv1 * (4.58257569495584 * Uphi1_east_3 + 4.58257569495584 *
    lambda_12_3);
  Jacobian[108] = ec_residual_tmp;
  Jacobian[78] = nuv1 * (3.872983346207417 * Uphi1_east_3 + 3.872983346207417 *
    lambda_12_3);
  fc_residual_tmp = nuv1 * (3.0 * Uphi1_east_3 + 3.0 * lambda_12_3);
  Jacobian[48] = fc_residual_tmp;
  gc_residual_tmp = nuv1 * (1.7320508075688772 * Uphi1_east_3 +
    1.7320508075688772 * lambda_12_3);
  Jacobian[18] = gc_residual_tmp;
  Jacobian[887] = l_residual_tmp;
  Jacobian[857] = nuv1 * (((1.7320508075688772 * lambda_23_1 -
    1.1065666703449764 * Uphi3_west_3) - 1.7320508075688772 * Uphi3_west_1) +
    1.1065666703449764 * lambda_23_3);
  Jacobian[827] = nuv1 * (((Uphi3_west_1 + 0.63887656499993994 * Uphi3_west_3) -
    lambda_23_1) - 0.63887656499993994 * lambda_23_3);
  l_residual_tmp = nuv1 * (((2.0 * Uphi3_west_2 + 1.9639610121239315 * q[29]) -
    2.0 * lambda_23_2) - 1.9639610121239315 * lambda_23_4);
  Jacobian[797] = l_residual_tmp;
  hc_residual_tmp = nuv1 * (((1.5491933384829668 * lambda_23_2 -
    1.52127765851133 * q[29]) - 1.5491933384829668 * Uphi3_west_2) +
    1.52127765851133 * lambda_23_4);
  Jacobian[767] = hc_residual_tmp;
  ic_residual_tmp = nuv1 * (((0.89442719099991586 * Uphi3_west_2 +
    0.87831006565367986 * q[29]) - 0.89442719099991586 * lambda_23_2) -
    0.87831006565367986 * lambda_23_4);
  Jacobian[737] = ic_residual_tmp;
  Jacobian[707] = nuv1 * (2.6457513110645907 * lambda_23_3 - 2.6457513110645907 *
    Uphi3_west_3);
  jc_residual_tmp = nuv1 * (2.23606797749979 * Uphi3_west_3 - 2.23606797749979 *
    lambda_23_3);
  Jacobian[677] = jc_residual_tmp;
  kc_residual_tmp = nuv1 * (1.7320508075688772 * lambda_23_3 -
    1.7320508075688772 * Uphi3_west_3);
  Jacobian[647] = kc_residual_tmp;
  mc_residual_tmp = nuv1 * (Uphi3_west_3 - lambda_23_3);
  Jacobian[617] = mc_residual_tmp;
  Jacobian[587] = 11.832159566199232 + r_residual_tmp;
  Jacobian[557] = nuv1 * ((((q_residual_tmp - 1.7320508075688772 * lambda_12_1)
    - 1.1065666703449764 * lambda_12_3) + 1.7320508075688772 * lambda_23_1) +
    1.1065666703449764 * lambda_23_3);
  Jacobian[527] = 10.0 + nuv1 * (((((((Uphi2_east_1 + 0.63887656499993994 *
    Uphi2_east_3) - Uphi2_west_1) - 0.63887656499993994 * Uphi2_west_3) +
    lambda_12_1) + 0.63887656499993994 * lambda_12_3) + lambda_23_1) +
    0.63887656499993994 * lambda_23_3);
  q_residual_tmp = ((2.0 * Uphi2_east_2 + 1.9639610121239315 * q[19]) - 2.0 *
                    Uphi2_west_2) - 1.9639610121239315 * q[19];
  Jacobian[497] = nuv1 * ((((q_residual_tmp + 2.0 * lambda_12_2) +
    1.9639610121239315 * lambda_12_4) + 2.0 * lambda_23_2) + 1.9639610121239315 *
    lambda_23_4);
  r_residual_tmp = nuv1 * ((((mb_residual_tmp - 1.5491933384829668 * lambda_12_2)
    - 1.52127765851133 * lambda_12_4) + 1.5491933384829668 * lambda_23_2) +
    1.52127765851133 * lambda_23_4);
  Jacobian[467] = r_residual_tmp;
  mb_residual_tmp = nuv1 * (((((((0.89442719099991586 * Uphi2_east_2 +
    0.87831006565367986 * q[19]) - 0.89442719099991586 * Uphi2_west_2) -
    0.87831006565367986 * q[19]) + 0.89442719099991586 * lambda_12_2) +
    0.87831006565367986 * lambda_12_4) + 0.89442719099991586 * lambda_23_2) +
    0.87831006565367986 * lambda_23_4);
  Jacobian[437] = -7.745966692414834 + mb_residual_tmp;
  nc_residual_tmp = 2.6457513110645907 * Uphi2_east_3 + 2.6457513110645907 *
    Uphi2_west_3;
  Jacobian[407] = nuv1 * ((nc_residual_tmp - 2.6457513110645907 * lambda_12_3) +
    2.6457513110645907 * lambda_23_3);
  oc_residual_tmp = 2.23606797749979 * Uphi2_east_3 - 2.23606797749979 *
    Uphi2_west_3;
  Jacobian[377] = nuv1 * ((oc_residual_tmp + 2.23606797749979 * lambda_12_3) +
    2.23606797749979 * lambda_23_3);
  vb_residual_tmp = nuv1 * ((vb_residual_tmp - 1.7320508075688772 * lambda_12_3)
    + 1.7320508075688772 * lambda_23_3);
  Jacobian[347] = vb_residual_tmp;
  pc_residual_tmp = 4.47213595499958 + nuv1 * (((Uphi2_east_3 - Uphi2_west_3) +
    lambda_12_3) + lambda_23_3);
  Jacobian[317] = pc_residual_tmp;
  Jacobian[287] = ab_residual_tmp;
  Jacobian[257] = nuv1 * (((-1.7320508075688772 * Uphi1_east_1 -
    1.1065666703449764 * Uphi1_east_3) - 1.7320508075688772 * lambda_12_1) -
    1.1065666703449764 * lambda_12_3);
  Jacobian[227] = nuv1 * (((-Uphi1_east_1 - 0.63887656499993994 * Uphi1_east_3)
    - lambda_12_1) - 0.63887656499993994 * lambda_12_3);
  ab_residual_tmp = nuv1 * (((-2.0 * Uphi1_east_2 - 1.9639610121239315 * q[9]) -
    2.0 * lambda_12_2) - 1.9639610121239315 * lambda_12_4);
  Jacobian[197] = ab_residual_tmp;
  lambda_sample_12_1 = nuv1 * (((-1.5491933384829668 * Uphi1_east_2 -
    1.52127765851133 * q[9]) - 1.5491933384829668 * lambda_12_2) -
    1.52127765851133 * lambda_12_4);
  Jacobian[167] = lambda_sample_12_1;
  lambda_sample_12_2 = nuv1 * (((-0.89442719099991586 * Uphi1_east_2 -
    0.87831006565367986 * q[9]) - 0.89442719099991586 * lambda_12_2) -
    0.87831006565367986 * lambda_12_4);
  Jacobian[137] = lambda_sample_12_2;
  Jacobian[107] = nuv1 * (-2.6457513110645907 * Uphi1_east_3 -
    2.6457513110645907 * lambda_12_3);
  lambda_sample_12_3 = nuv1 * (-2.23606797749979 * Uphi1_east_3 -
    2.23606797749979 * lambda_12_3);
  Jacobian[77] = lambda_sample_12_3;
  lambda_sample_12_4 = nuv1 * (-1.7320508075688772 * Uphi1_east_3 -
    1.7320508075688772 * lambda_12_3);
  Jacobian[47] = lambda_sample_12_4;
  lambda_sample_23_1 = nuv1 * (-Uphi1_east_3 - lambda_12_3);
  Jacobian[17] = lambda_sample_23_1;
  Jacobian[886] = m_residual_tmp;
  Jacobian[856] = nuv1 * (((3.4641016151377544 * lambda_23_2 -
    3.4016802570830449 * q[29]) - 3.4641016151377544 * Uphi3_west_2) +
    3.4016802570830449 * lambda_23_4);
  Jacobian[826] = l_residual_tmp;
  Jacobian[796] = nuv1 * (((5.0 * Uphi3_west_1 + 4.47213595499958 * Uphi3_west_3)
    - 5.0 * lambda_23_1) - 4.47213595499958 * lambda_23_3);
  Jacobian[766] = nuv1 * (((3.872983346207417 * lambda_23_1 - 3.4641016151377544
    * Uphi3_west_3) - 3.872983346207417 * Uphi3_west_1) + 3.4641016151377544 *
    lambda_23_3);
  l_residual_tmp = nuv1 * (((2.23606797749979 * Uphi3_west_1 + 2.0 *
    Uphi3_west_3) - 2.23606797749979 * lambda_23_1) - 2.0 * lambda_23_3);
  Jacobian[736] = l_residual_tmp;
  Jacobian[706] = nuv1 * (5.9160797830996161 * lambda_23_2 - 5.9160797830996161 *
    Uphi3_west_2);
  m_residual_tmp = nuv1 * (5.0 * Uphi3_west_2 - 5.0 * lambda_23_2);
  Jacobian[676] = m_residual_tmp;
  lambda_sample_23_2 = nuv1 * (3.872983346207417 * lambda_23_2 -
    3.872983346207417 * Uphi3_west_2);
  Jacobian[646] = lambda_sample_23_2;
  lambda_sample_23_3 = nuv1 * (2.23606797749979 * Uphi3_west_2 -
    2.23606797749979 * lambda_23_2);
  Jacobian[616] = lambda_sample_23_3;
  Jacobian[586] = nuv1 * (((s_residual_tmp - 13.60672102833218 * q[18]) +
    1.9639610121239315 * lambda_12_3) + 1.9639610121239315 * lambda_23_3);
  Jacobian[556] = nuv1 * (((((((kb_residual_tmp - 13.856406460551018 * q[14]) -
    12.393546707863734 * q[16]) - 13.60672102833218 * q[19]) -
    3.4641016151377544 * lambda_12_2) - 3.4016802570830449 * lambda_12_4) +
    3.4641016151377544 * lambda_23_2) + 3.4016802570830449 * lambda_23_4);
  Jacobian[526] = nuv1 * (((((q_residual_tmp - 13.856406460551018 * q[15]) + 2.0
    * lambda_12_2) + 1.9639610121239315 * lambda_12_4) + 2.0 * lambda_23_2) +
    1.9639610121239315 * lambda_23_4);
  Jacobian[496] = 6.0 + nuv1 * ((((((((((5.0 * Uphi2_east_1 + 4.47213595499958 *
    Uphi2_east_3) - 5.0 * Uphi2_west_1) - 4.47213595499958 * Uphi2_west_3) -
    13.856406460551018 * q[11]) - 13.60672102833218 * q[13]) -
    12.393546707863734 * q[18]) + 5.0 * lambda_12_1) + 4.47213595499958 *
    lambda_12_3) + 5.0 * lambda_23_1) + 4.47213595499958 * lambda_23_3);
  q_residual_tmp = (ob_residual_tmp + 3.872983346207417 * Uphi2_west_1) +
    3.4641016151377544 * Uphi2_west_3;
  Jacobian[466] = nuv1 * (((((((q_residual_tmp - 15.491933384829668 * q[10]) -
    13.856406460551018 * q[12]) - 13.856406460551018 * q[17]) -
    3.872983346207417 * lambda_12_1) - 3.4641016151377544 * lambda_12_3) +
    3.872983346207417 * lambda_23_1) + 3.4641016151377544 * lambda_23_3);
  s_residual_tmp = ((2.23606797749979 * Uphi2_east_1 + 2.0 * Uphi2_east_3) -
                    2.23606797749979 * Uphi2_west_1) - 2.0 * Uphi2_west_3;
  Jacobian[436] = nuv1 * ((((((s_residual_tmp - 15.491933384829668 * q[11]) -
    13.856406460551018 * q[18]) + 2.23606797749979 * lambda_12_1) + 2.0 *
    lambda_12_3) + 2.23606797749979 * lambda_23_1) + 2.0 * lambda_23_3);
  kb_residual_tmp = 5.9160797830996161 * Uphi2_east_2 + 5.9160797830996161 *
    Uphi2_west_2;
  Jacobian[406] = nuv1 * (((kb_residual_tmp - 13.60672102833218 * q[16]) -
    5.9160797830996161 * lambda_12_2) + 5.9160797830996161 * lambda_23_2);
  ob_residual_tmp = nuv1 * ((((5.0 * Uphi2_east_2 - 5.0 * Uphi2_west_2) -
    13.856406460551018 * q[15]) + 5.0 * lambda_12_2) + 5.0 * lambda_23_2);
  Jacobian[376] = -3.4641016151377544 + ob_residual_tmp;
  lambda_sample_23_4 = 3.872983346207417 * Uphi2_east_2 + 3.872983346207417 *
    Uphi2_west_2;
  residual_tmp = nuv1 * ((((lambda_sample_23_4 - 15.491933384829668 * q[14]) -
    13.856406460551018 * q[16]) - 3.872983346207417 * lambda_12_2) +
    3.872983346207417 * lambda_23_2);
  Jacobian[346] = residual_tmp;
  b_residual_tmp = 2.23606797749979 * Uphi2_east_2 - 2.23606797749979 *
    Uphi2_west_2;
  c_residual_tmp = nuv1 * (((b_residual_tmp - 15.491933384829668 * q[15]) +
    2.23606797749979 * lambda_12_2) + 2.23606797749979 * lambda_23_2);
  Jacobian[316] = c_residual_tmp;
  Jacobian[286] = bb_residual_tmp;
  Jacobian[256] = nuv1 * (((-3.4641016151377544 * Uphi1_east_2 -
    3.4016802570830449 * q[9]) - 3.4641016151377544 * lambda_12_2) -
    3.4016802570830449 * lambda_12_4);
  Jacobian[226] = ab_residual_tmp;
  Jacobian[196] = nuv1 * (((-5.0 * Uphi1_east_1 - 4.47213595499958 *
    Uphi1_east_3) - 5.0 * lambda_12_1) - 4.47213595499958 * lambda_12_3);
  Jacobian[166] = nuv1 * (((-3.872983346207417 * Uphi1_east_1 -
    3.4641016151377544 * Uphi1_east_3) - 3.872983346207417 * lambda_12_1) -
    3.4641016151377544 * lambda_12_3);
  ab_residual_tmp = nuv1 * (((-2.23606797749979 * Uphi1_east_1 - 2.0 *
    Uphi1_east_3) - 2.23606797749979 * lambda_12_1) - 2.0 * lambda_12_3);
  Jacobian[136] = ab_residual_tmp;
  Jacobian[106] = nuv1 * (-5.9160797830996161 * Uphi1_east_2 -
    5.9160797830996161 * lambda_12_2);
  bb_residual_tmp = nuv1 * (-5.0 * Uphi1_east_2 - 5.0 * lambda_12_2);
  Jacobian[76] = bb_residual_tmp;
  d_residual_tmp = nuv1 * (-3.872983346207417 * Uphi1_east_2 - 3.872983346207417
    * lambda_12_2);
  Jacobian[46] = d_residual_tmp;
  e_residual_tmp = nuv1 * (-2.23606797749979 * Uphi1_east_2 - 2.23606797749979 *
    lambda_12_2);
  Jacobian[16] = e_residual_tmp;
  Jacobian[885] = nuv1 * (1.52127765851133 * Uphi3_west_3 - 1.52127765851133 *
    lambda_23_3);
  Jacobian[855] = fb_residual_tmp;
  Jacobian[825] = gb_residual_tmp;
  Jacobian[795] = nuv1 * (((3.872983346207417 * Uphi3_west_1 +
    3.4641016151377544 * Uphi3_west_3) - 3.872983346207417 * lambda_23_1) -
    3.4641016151377544 * lambda_23_3);
  Jacobian[765] = nuv1 * (((3.0 * lambda_23_1 - 2.6832815729997477 *
    Uphi3_west_3) - 3.0 * Uphi3_west_1) + 2.6832815729997477 * lambda_23_3);
  Jacobian[735] = nuv1 * (((1.7320508075688772 * Uphi3_west_1 +
    1.5491933384829668 * Uphi3_west_3) - 1.7320508075688772 * lambda_23_1) -
    1.5491933384829668 * lambda_23_3);
  fb_residual_tmp = nuv1 * (4.58257569495584 * lambda_23_2 - 4.58257569495584 *
    Uphi3_west_2);
  Jacobian[705] = fb_residual_tmp;
  gb_residual_tmp = nuv1 * (3.872983346207417 * Uphi3_west_2 - 3.872983346207417
    * lambda_23_2);
  Jacobian[675] = gb_residual_tmp;
  f_residual_tmp = nuv1 * (3.0 * lambda_23_2 - 3.0 * Uphi3_west_2);
  Jacobian[645] = f_residual_tmp;
  g_residual_tmp = nuv1 * (1.7320508075688772 * Uphi3_west_2 -
    1.7320508075688772 * lambda_23_2);
  Jacobian[615] = g_residual_tmp;
  Jacobian[585] = nuv1 * (((t_residual_tmp - 6.08511063404532 * q[17]) -
    1.52127765851133 * lambda_12_3) + 1.52127765851133 * lambda_23_3);
  Jacobian[555] = 7.745966692414834 + lb_residual_tmp;
  Jacobian[525] = rb_residual_tmp;
  Jacobian[495] = nuv1 * (((((q_residual_tmp - 6.9282032302755088 * q[12]) -
    3.872983346207417 * lambda_12_1) - 3.4641016151377544 * lambda_12_3) +
    3.872983346207417 * lambda_23_1) + 3.4641016151377544 * lambda_23_3);
  Jacobian[465] = 6.0 + nuv1 * ((((((((qb_residual_tmp - 3.0 * Uphi2_west_1) -
    2.6832815729997477 * Uphi2_west_3) - 6.9282032302755088 * q[11]) -
    6.1967733539318672 * q[18]) + 3.0 * lambda_12_1) + 2.6832815729997477 *
    lambda_12_3) + 3.0 * lambda_23_1) + 2.6832815729997477 * lambda_23_3);
  q_residual_tmp = (lc_residual_tmp + 1.7320508075688772 * Uphi2_west_1) +
    1.5491933384829668 * Uphi2_west_3;
  Jacobian[435] = nuv1 * ((((((q_residual_tmp - 6.9282032302755088 * q[10]) -
    6.1967733539318672 * q[17]) - 1.7320508075688772 * lambda_12_1) -
    1.5491933384829668 * lambda_12_3) + 1.7320508075688772 * lambda_23_1) +
    1.5491933384829668 * lambda_23_3);
  t_residual_tmp = 4.58257569495584 * Uphi2_east_2 - 4.58257569495584 *
    Uphi2_west_2;
  Jacobian[405] = nuv1 * ((t_residual_tmp + 4.58257569495584 * lambda_12_2) +
    4.58257569495584 * lambda_23_2);
  lb_residual_tmp = nuv1 * (((lambda_sample_23_4 - 6.9282032302755088 * q[16]) -
    3.872983346207417 * lambda_12_2) + 3.872983346207417 * lambda_23_2);
  Jacobian[375] = lb_residual_tmp;
  qb_residual_tmp = nuv1 * ((((3.0 * Uphi2_east_2 - 3.0 * Uphi2_west_2) -
    6.9282032302755088 * q[15]) + 3.0 * lambda_12_2) + 3.0 * lambda_23_2);
  Jacobian[345] = -3.4641016151377544 + qb_residual_tmp;
  rb_residual_tmp = 1.7320508075688772 * Uphi2_east_2 + 1.7320508075688772 *
    Uphi2_west_2;
  lc_residual_tmp = nuv1 * (((rb_residual_tmp - 6.9282032302755088 * q[14]) -
    1.7320508075688772 * lambda_12_2) + 1.7320508075688772 * lambda_23_2);
  Jacobian[315] = lc_residual_tmp;
  Jacobian[285] = nuv1 * (1.52127765851133 * Uphi1_east_3 + 1.52127765851133 *
    lambda_12_3);
  Jacobian[255] = cc_residual_tmp;
  Jacobian[225] = dc_residual_tmp;
  cc_residual_tmp = 3.872983346207417 * Uphi1_east_1 + 3.4641016151377544 *
    Uphi1_east_3;
  Jacobian[195] = nuv1 * ((cc_residual_tmp + 3.872983346207417 * lambda_12_1) +
    3.4641016151377544 * lambda_12_3);
  dc_residual_tmp = 3.0 * Uphi1_east_1 + 2.6832815729997477 * Uphi1_east_3;
  Jacobian[165] = nuv1 * ((dc_residual_tmp + 3.0 * lambda_12_1) +
    2.6832815729997477 * lambda_12_3);
  lambda_sample_23_4 = 1.7320508075688772 * Uphi1_east_1 + 1.5491933384829668 *
    Uphi1_east_3;
  Jacobian[135] = nuv1 * ((lambda_sample_23_4 + 1.7320508075688772 * lambda_12_1)
    + 1.5491933384829668 * lambda_12_3);
  h_residual_tmp = nuv1 * (4.58257569495584 * Uphi1_east_2 + 4.58257569495584 *
    lambda_12_2);
  Jacobian[105] = h_residual_tmp;
  i_residual_tmp = nuv1 * (3.872983346207417 * Uphi1_east_2 + 3.872983346207417 *
    lambda_12_2);
  Jacobian[75] = i_residual_tmp;
  j_residual_tmp = nuv1 * (3.0 * Uphi1_east_2 + 3.0 * lambda_12_2);
  Jacobian[45] = j_residual_tmp;
  k_residual_tmp = nuv1 * (1.7320508075688772 * Uphi1_east_2 +
    1.7320508075688772 * lambda_12_2);
  Jacobian[15] = k_residual_tmp;
  Jacobian[884] = n_residual_tmp;
  Jacobian[854] = hc_residual_tmp;
  Jacobian[824] = ic_residual_tmp;
  Jacobian[794] = l_residual_tmp;
  Jacobian[764] = nuv1 * (((1.7320508075688772 * lambda_23_1 -
    1.5491933384829668 * Uphi3_west_3) - 1.7320508075688772 * Uphi3_west_1) +
    1.5491933384829668 * lambda_23_3);
  Jacobian[734] = nuv1 * (((Uphi3_west_1 + 0.89442719099991586 * Uphi3_west_3) -
    lambda_23_1) - 0.89442719099991586 * lambda_23_3);
  Jacobian[704] = nuv1 * (2.6457513110645907 * lambda_23_2 - 2.6457513110645907 *
    Uphi3_west_2);
  Jacobian[674] = lambda_sample_23_3;
  l_residual_tmp = nuv1 * (1.7320508075688772 * lambda_23_2 - 1.7320508075688772
    * Uphi3_west_2);
  Jacobian[644] = l_residual_tmp;
  n_residual_tmp = nuv1 * (Uphi3_west_2 - lambda_23_2);
  Jacobian[614] = n_residual_tmp;
  Jacobian[584] = u_residual_tmp;
  Jacobian[554] = r_residual_tmp;
  Jacobian[524] = 7.745966692414834 + mb_residual_tmp;
  Jacobian[494] = nuv1 * ((((s_residual_tmp + 2.23606797749979 * lambda_12_1) +
    2.0 * lambda_12_3) + 2.23606797749979 * lambda_23_1) + 2.0 * lambda_23_3);
  Jacobian[464] = nuv1 * ((((q_residual_tmp - 1.7320508075688772 * lambda_12_1)
    - 1.5491933384829668 * lambda_12_3) + 1.7320508075688772 * lambda_23_1) +
    1.5491933384829668 * lambda_23_3);
  Jacobian[434] = 6.0 + nuv1 * (((((((Uphi2_east_1 + 0.89442719099991586 *
    Uphi2_east_3) - Uphi2_west_1) - 0.89442719099991586 * Uphi2_west_3) +
    lambda_12_1) + 0.89442719099991586 * lambda_12_3) + lambda_23_1) +
    0.89442719099991586 * lambda_23_3);
  q_residual_tmp = 2.6457513110645907 * Uphi2_east_2 + 2.6457513110645907 *
    Uphi2_west_2;
  Jacobian[404] = nuv1 * ((q_residual_tmp - 2.6457513110645907 * lambda_12_2) +
    2.6457513110645907 * lambda_23_2);
  r_residual_tmp = nuv1 * ((b_residual_tmp + 2.23606797749979 * lambda_12_2) +
    2.23606797749979 * lambda_23_2);
  Jacobian[374] = r_residual_tmp;
  s_residual_tmp = nuv1 * ((rb_residual_tmp - 1.7320508075688772 * lambda_12_2)
    + 1.7320508075688772 * lambda_23_2);
  Jacobian[344] = s_residual_tmp;
  u_residual_tmp = nuv1 * (((Uphi2_east_2 - Uphi2_west_2) + lambda_12_2) +
    lambda_23_2);
  Jacobian[314] = -3.4641016151377544 + u_residual_tmp;
  Jacobian[284] = cb_residual_tmp;
  Jacobian[254] = lambda_sample_12_1;
  Jacobian[224] = lambda_sample_12_2;
  Jacobian[194] = ab_residual_tmp;
  Jacobian[164] = nuv1 * (((-1.7320508075688772 * Uphi1_east_1 -
    1.5491933384829668 * Uphi1_east_3) - 1.7320508075688772 * lambda_12_1) -
    1.5491933384829668 * lambda_12_3);
  Jacobian[134] = nuv1 * (((-Uphi1_east_1 - 0.89442719099991586 * Uphi1_east_3)
    - lambda_12_1) - 0.89442719099991586 * lambda_12_3);
  Jacobian[104] = nuv1 * (-2.6457513110645907 * Uphi1_east_2 -
    2.6457513110645907 * lambda_12_2);
  Jacobian[74] = e_residual_tmp;
  ab_residual_tmp = nuv1 * (-1.7320508075688772 * Uphi1_east_2 -
    1.7320508075688772 * lambda_12_2);
  Jacobian[44] = ab_residual_tmp;
  cb_residual_tmp = nuv1 * (-Uphi1_east_2 - lambda_12_2);
  Jacobian[14] = cb_residual_tmp;
  Jacobian[883] = nuv1 * (2.6457513110645907 * q[29] - 2.6457513110645907 *
    lambda_23_4);
  Jacobian[853] = hb_residual_tmp;
  Jacobian[823] = nuv1 * (2.6457513110645907 * Uphi3_west_3 - 2.6457513110645907
    * lambda_23_3);
  Jacobian[793] = nuv1 * (5.9160797830996161 * Uphi3_west_2 - 5.9160797830996161
    * lambda_23_2);
  Jacobian[763] = fb_residual_tmp;
  Jacobian[733] = nuv1 * (2.6457513110645907 * Uphi3_west_2 - 2.6457513110645907
    * lambda_23_2);
  Jacobian[703] = nuv1 * (7.0 * lambda_23_1 - 7.0 * Uphi3_west_1);
  Jacobian[673] = nuv1 * (5.9160797830996161 * Uphi3_west_1 - 5.9160797830996161
    * lambda_23_1);
  fb_residual_tmp = nuv1 * (4.58257569495584 * lambda_23_1 - 4.58257569495584 *
    Uphi3_west_1);
  Jacobian[643] = fb_residual_tmp;
  Jacobian[613] = nuv1 * (2.6457513110645907 * Uphi3_west_1 - 2.6457513110645907
    * lambda_23_1);
  Jacobian[583] = nuv1 * (((v_residual_tmp - 10.583005244258363 * q[19]) -
    2.6457513110645907 * lambda_12_4) + 2.6457513110645907 * lambda_23_4);
  Jacobian[553] = nuv1 * (((sb_residual_tmp - 31.749015732775089 * q[18]) +
    4.58257569495584 * lambda_12_3) + 4.58257569495584 * lambda_23_3);
  Jacobian[523] = nuv1 * (((nc_residual_tmp - 10.583005244258363 * q[17]) -
    2.6457513110645907 * lambda_12_3) + 2.6457513110645907 * lambda_23_3);
  Jacobian[493] = nuv1 * ((((kb_residual_tmp - 23.664319132398465 * q[14]) -
    25.701584164627452 * q[16]) - 5.9160797830996161 * lambda_12_2) +
    5.9160797830996161 * lambda_23_2);
  Jacobian[463] = nuv1 * (((t_residual_tmp - 31.749015732775089 * q[15]) +
    4.58257569495584 * lambda_12_2) + 4.58257569495584 * lambda_23_2);
  Jacobian[433] = nuv1 * ((((q_residual_tmp - 10.583005244258363 * q[14]) -
    23.664319132398465 * q[16]) - 2.6457513110645907 * lambda_12_2) +
    2.6457513110645907 * lambda_23_2);
  Jacobian[403] = 2.0 + nuv1 * (((((7.0 * Uphi2_east_1 - 7.0 * Uphi2_west_1) -
    20.784609690826528 * q[11]) - 24.693678903269511 * q[13]) + 7.0 *
    lambda_12_1) + 7.0 * lambda_23_1);
  q_residual_tmp = 5.9160797830996161 * Uphi2_east_1 + 5.9160797830996161 *
    Uphi2_west_1;
  Jacobian[373] = nuv1 * ((((q_residual_tmp - 23.664319132398465 * q[10]) -
    25.701584164627452 * q[12]) - 5.9160797830996161 * lambda_12_1) +
    5.9160797830996161 * lambda_23_1);
  t_residual_tmp = 4.58257569495584 * Uphi2_east_1 - 4.58257569495584 *
    Uphi2_west_1;
  Jacobian[343] = nuv1 * ((((t_residual_tmp - 31.749015732775089 * q[11]) -
    20.784609690826528 * q[13]) + 4.58257569495584 * lambda_12_1) +
    4.58257569495584 * lambda_23_1);
  v_residual_tmp = 2.6457513110645907 * Uphi2_east_1 + 2.6457513110645907 *
    Uphi2_west_1;
  Jacobian[313] = nuv1 * ((((v_residual_tmp - 10.583005244258363 * q[10]) -
    23.664319132398465 * q[12]) - 2.6457513110645907 * lambda_12_1) +
    2.6457513110645907 * lambda_23_1);
  Jacobian[283] = nuv1 * (2.6457513110645907 * q[9] + 2.6457513110645907 *
    lambda_12_4);
  Jacobian[253] = ec_residual_tmp;
  Jacobian[223] = nuv1 * (2.6457513110645907 * Uphi1_east_3 + 2.6457513110645907
    * lambda_12_3);
  Jacobian[193] = nuv1 * (5.9160797830996161 * Uphi1_east_2 + 5.9160797830996161
    * lambda_12_2);
  Jacobian[163] = h_residual_tmp;
  Jacobian[133] = nuv1 * (2.6457513110645907 * Uphi1_east_2 + 2.6457513110645907
    * lambda_12_2);
  Jacobian[103] = nuv1 * (7.0 * Uphi1_east_1 + 7.0 * lambda_12_1);
  Jacobian[73] = nuv1 * (5.9160797830996161 * Uphi1_east_1 + 5.9160797830996161 *
    lambda_12_1);
  hb_residual_tmp = nuv1 * (4.58257569495584 * Uphi1_east_1 + 4.58257569495584 *
    lambda_12_1);
  Jacobian[43] = hb_residual_tmp;
  Jacobian[13] = nuv1 * (2.6457513110645907 * Uphi1_east_1 + 2.6457513110645907 *
    lambda_12_1);
  Jacobian[882] = o_residual_tmp;
  Jacobian[852] = nuv1 * (3.872983346207417 * lambda_23_3 - 3.872983346207417 *
    Uphi3_west_3);
  Jacobian[822] = jc_residual_tmp;
  Jacobian[792] = m_residual_tmp;
  Jacobian[762] = lambda_sample_23_2;
  Jacobian[732] = lambda_sample_23_3;
  Jacobian[702] = nuv1 * (5.9160797830996161 * lambda_23_1 - 5.9160797830996161 *
    Uphi3_west_1);
  Jacobian[672] = nuv1 * (5.0 * Uphi3_west_1 - 5.0 * lambda_23_1);
  Jacobian[642] = nuv1 * (3.872983346207417 * lambda_23_1 - 3.872983346207417 *
    Uphi3_west_1);
  m_residual_tmp = nuv1 * (2.23606797749979 * Uphi3_west_1 - 2.23606797749979 *
    lambda_23_1);
  Jacobian[612] = m_residual_tmp;
  Jacobian[582] = w_residual_tmp;
  Jacobian[552] = nuv1 * (((tb_residual_tmp - 15.491933384829668 * q[17]) -
    3.872983346207417 * lambda_12_3) + 3.872983346207417 * lambda_23_3);
  Jacobian[522] = nuv1 * (((oc_residual_tmp - 15.491933384829668 * q[18]) +
    2.23606797749979 * lambda_12_3) + 2.23606797749979 * lambda_23_3);
  Jacobian[492] = 3.4641016151377544 + ob_residual_tmp;
  Jacobian[462] = residual_tmp;
  Jacobian[432] = c_residual_tmp;
  Jacobian[402] = nuv1 * (((q_residual_tmp - 13.60672102833218 * q[12]) -
    5.9160797830996161 * lambda_12_1) + 5.9160797830996161 * lambda_23_1);
  Jacobian[372] = 2.0 + nuv1 * (((((5.0 * Uphi2_east_1 - 5.0 * Uphi2_west_1) -
    13.856406460551018 * q[11]) - 13.60672102833218 * q[13]) + 5.0 * lambda_12_1)
    + 5.0 * lambda_23_1);
  o_residual_tmp = 3.872983346207417 * Uphi2_east_1 + 3.872983346207417 *
    Uphi2_west_1;
  Jacobian[342] = nuv1 * ((((o_residual_tmp - 15.491933384829668 * q[10]) -
    13.856406460551018 * q[12]) - 3.872983346207417 * lambda_12_1) +
    3.872983346207417 * lambda_23_1);
  q_residual_tmp = 2.23606797749979 * Uphi2_east_1 - 2.23606797749979 *
    Uphi2_west_1;
  Jacobian[312] = nuv1 * (((q_residual_tmp - 15.491933384829668 * q[11]) +
    2.23606797749979 * lambda_12_1) + 2.23606797749979 * lambda_23_1);
  Jacobian[282] = db_residual_tmp;
  Jacobian[252] = nuv1 * (-3.872983346207417 * Uphi1_east_3 - 3.872983346207417 *
    lambda_12_3);
  Jacobian[222] = lambda_sample_12_3;
  Jacobian[192] = bb_residual_tmp;
  Jacobian[162] = d_residual_tmp;
  Jacobian[132] = e_residual_tmp;
  Jacobian[102] = nuv1 * (-5.9160797830996161 * Uphi1_east_1 -
    5.9160797830996161 * lambda_12_1);
  Jacobian[72] = nuv1 * (-5.0 * Uphi1_east_1 - 5.0 * lambda_12_1);
  Jacobian[42] = nuv1 * (-3.872983346207417 * Uphi1_east_1 - 3.872983346207417 *
    lambda_12_1);
  w_residual_tmp = nuv1 * (-2.23606797749979 * Uphi1_east_1 - 2.23606797749979 *
    lambda_12_1);
  Jacobian[12] = w_residual_tmp;
  Jacobian[881] = nuv1 * (1.7320508075688772 * q[29] - 1.7320508075688772 *
    lambda_23_4);
  Jacobian[851] = ib_residual_tmp;
  Jacobian[821] = jb_residual_tmp;
  Jacobian[791] = gb_residual_tmp;
  Jacobian[761] = f_residual_tmp;
  Jacobian[731] = g_residual_tmp;
  Jacobian[701] = fb_residual_tmp;
  Jacobian[671] = nuv1 * (3.872983346207417 * Uphi3_west_1 - 3.872983346207417 *
    lambda_23_1);
  Jacobian[641] = nuv1 * (3.0 * lambda_23_1 - 3.0 * Uphi3_west_1);
  Jacobian[611] = nuv1 * (1.7320508075688772 * Uphi3_west_1 - 1.7320508075688772
    * lambda_23_1);
  Jacobian[581] = nuv1 * (((x_residual_tmp - 6.9282032302755088 * q[19]) -
    1.7320508075688772 * lambda_12_4) + 1.7320508075688772 * lambda_23_4);
  Jacobian[551] = ub_residual_tmp;
  Jacobian[521] = wb_residual_tmp;
  Jacobian[491] = lb_residual_tmp;
  Jacobian[461] = 3.4641016151377544 + qb_residual_tmp;
  Jacobian[431] = lc_residual_tmp;
  Jacobian[401] = nuv1 * (((t_residual_tmp - 6.9282032302755088 * q[13]) +
    4.58257569495584 * lambda_12_1) + 4.58257569495584 * lambda_23_1);
  Jacobian[371] = nuv1 * (((o_residual_tmp - 6.9282032302755088 * q[12]) -
    3.872983346207417 * lambda_12_1) + 3.872983346207417 * lambda_23_1);
  Jacobian[341] = 2.0 + nuv1 * ((((3.0 * Uphi2_east_1 - 3.0 * Uphi2_west_1) -
    6.9282032302755088 * q[11]) + 3.0 * lambda_12_1) + 3.0 * lambda_23_1);
  o_residual_tmp = 1.7320508075688772 * Uphi2_east_1 + 1.7320508075688772 *
    Uphi2_west_1;
  Jacobian[311] = nuv1 * (((o_residual_tmp - 6.9282032302755088 * q[10]) -
    1.7320508075688772 * lambda_12_1) + 1.7320508075688772 * lambda_23_1);
  Jacobian[281] = nuv1 * (1.7320508075688772 * q[9] + 1.7320508075688772 *
    lambda_12_4);
  Jacobian[251] = fc_residual_tmp;
  Jacobian[221] = gc_residual_tmp;
  Jacobian[191] = i_residual_tmp;
  Jacobian[161] = j_residual_tmp;
  Jacobian[131] = k_residual_tmp;
  Jacobian[101] = hb_residual_tmp;
  Jacobian[71] = nuv1 * (3.872983346207417 * Uphi1_east_1 + 3.872983346207417 *
    lambda_12_1);
  Jacobian[41] = nuv1 * (3.0 * Uphi1_east_1 + 3.0 * lambda_12_1);
  Jacobian[11] = nuv1 * (1.7320508075688772 * Uphi1_east_1 + 1.7320508075688772 *
    lambda_12_1);
  Jacobian[880] = p_residual_tmp;
  Jacobian[850] = kc_residual_tmp;
  Jacobian[820] = mc_residual_tmp;
  Jacobian[790] = lambda_sample_23_3;
  Jacobian[760] = l_residual_tmp;
  Jacobian[730] = n_residual_tmp;
  Jacobian[700] = nuv1 * (2.6457513110645907 * lambda_23_1 - 2.6457513110645907 *
    Uphi3_west_1);
  Jacobian[670] = m_residual_tmp;
  Jacobian[640] = nuv1 * (1.7320508075688772 * lambda_23_1 - 1.7320508075688772 *
    Uphi3_west_1);
  Jacobian[610] = nuv1 * (Uphi3_west_1 - lambda_23_1);
  Jacobian[580] = 5.2915026221291814 + y_residual_tmp;
  Jacobian[550] = vb_residual_tmp;
  Jacobian[520] = pc_residual_tmp;
  Jacobian[490] = r_residual_tmp;
  Jacobian[460] = s_residual_tmp;
  Jacobian[430] = 3.4641016151377544 + u_residual_tmp;
  Jacobian[400] = nuv1 * ((v_residual_tmp - 2.6457513110645907 * lambda_12_1) +
    2.6457513110645907 * lambda_23_1);
  Jacobian[370] = nuv1 * ((q_residual_tmp + 2.23606797749979 * lambda_12_1) +
    2.23606797749979 * lambda_23_1);
  Jacobian[340] = nuv1 * ((o_residual_tmp - 1.7320508075688772 * lambda_12_1) +
    1.7320508075688772 * lambda_23_1);
  Jacobian[310] = 2.0 + nuv1 * (((Uphi2_east_1 - Uphi2_west_1) + lambda_12_1) +
    lambda_23_1);
  Jacobian[280] = eb_residual_tmp;
  Jacobian[250] = lambda_sample_12_4;
  Jacobian[220] = lambda_sample_23_1;
  Jacobian[190] = e_residual_tmp;
  Jacobian[160] = ab_residual_tmp;
  Jacobian[130] = cb_residual_tmp;
  Jacobian[100] = nuv1 * (-2.6457513110645907 * Uphi1_east_1 -
    2.6457513110645907 * lambda_12_1);
  Jacobian[70] = w_residual_tmp;
  Jacobian[40] = nuv1 * (-1.7320508075688772 * Uphi1_east_1 - 1.7320508075688772
    * lambda_12_1);
  Jacobian[10] = nuv1 * (-Uphi1_east_1 - lambda_12_1);
  Jacobian[579] = nuv1 * (((Uphi2_west_1 + 0.59628479399994394 * Uphi2_west_3) -
    lambda_12_1) - 0.59628479399994394 * lambda_12_3);
  Jacobian[549] = nuv1 * (((1.52127765851133 * lambda_12_2 - 1.0327955589886444 *
    q[19]) - 1.52127765851133 * Uphi2_west_2) + 1.0327955589886444 * lambda_12_4);
  l_residual_tmp = nuv1 * (((0.87831006565367986 * Uphi2_west_2 +
    0.59628479399994394 * q[19]) - 0.87831006565367986 * lambda_12_2) -
    0.59628479399994394 * lambda_12_4);
  Jacobian[519] = l_residual_tmp;
  m_residual_tmp = nuv1 * (1.9639610121239315 * Uphi2_west_3 -
    1.9639610121239315 * lambda_12_3);
  Jacobian[489] = m_residual_tmp;
  Jacobian[459] = nuv1 * (1.52127765851133 * lambda_12_3 - 1.52127765851133 *
    Uphi2_west_3);
  n_residual_tmp = nuv1 * (0.87831006565367986 * Uphi2_west_3 -
    0.87831006565367986 * lambda_12_3);
  Jacobian[429] = n_residual_tmp;
  Jacobian[399] = nuv1 * (2.6457513110645907 * lambda_12_4 - 2.6457513110645907 *
    q[19]);
  o_residual_tmp = nuv1 * (2.23606797749979 * q[19] - 2.23606797749979 *
    lambda_12_4);
  Jacobian[369] = o_residual_tmp;
  Jacobian[339] = nuv1 * (1.7320508075688772 * lambda_12_4 - 1.7320508075688772 *
    q[19]);
  p_residual_tmp = nuv1 * (q[19] - lambda_12_4);
  Jacobian[309] = p_residual_tmp;
  Jacobian[279] = 14.0 + nuv1 * (((((Uphi1_east_1 + 0.59628479399994394 *
    Uphi1_east_3) - 2.0 * Uphi1_west_1) - 1.1925695879998879 * Uphi1_west_3) +
    lambda_12_1) + 0.59628479399994394 * lambda_12_3);
  q_residual_tmp = (xb_residual_tmp + 3.04255531702266 * Uphi1_west_2) +
    2.0655911179772892 * q[9];
  Jacobian[249] = nuv1 * ((q_residual_tmp + 1.52127765851133 * lambda_12_2) +
    1.0327955589886444 * lambda_12_4);
  r_residual_tmp = nuv1 * (((((0.87831006565367986 * Uphi1_east_2 +
    0.59628479399994394 * q[9]) - 1.7566201313073597 * Uphi1_west_2) -
    1.1925695879998879 * q[9]) + 0.87831006565367986 * lambda_12_2) +
    0.59628479399994394 * lambda_12_4);
  Jacobian[219] = -11.832159566199232 + r_residual_tmp;
  s_residual_tmp = 1.9639610121239315 * Uphi1_east_3 - 3.927922024247863 *
    Uphi1_west_3;
  Jacobian[189] = nuv1 * (s_residual_tmp + 1.9639610121239315 * lambda_12_3);
  t_residual_tmp = 1.52127765851133 * Uphi1_east_3 + 3.04255531702266 *
    Uphi1_west_3;
  Jacobian[159] = nuv1 * (t_residual_tmp + 1.52127765851133 * lambda_12_3);
  u_residual_tmp = 9.16515138991168 + nuv1 * ((0.87831006565367986 *
    Uphi1_east_3 - 1.7566201313073597 * Uphi1_west_3) + 0.87831006565367986 *
    lambda_12_3);
  Jacobian[129] = u_residual_tmp;
  v_residual_tmp = 2.6457513110645907 * q[9] + 5.2915026221291814 * q[9];
  Jacobian[99] = nuv1 * (v_residual_tmp + 2.6457513110645907 * lambda_12_4);
  w_residual_tmp = nuv1 * ((2.23606797749979 * q[9] - 4.47213595499958 * q[9]) +
    2.23606797749979 * lambda_12_4);
  Jacobian[69] = w_residual_tmp;
  x_residual_tmp = 1.7320508075688772 * q[9] + 3.4641016151377544 * q[9];
  Jacobian[39] = nuv1 * (x_residual_tmp + 1.7320508075688772 * lambda_12_4);
  y_residual_tmp = nuv1 * ((q[9] - 2.0 * q[9]) + lambda_12_4);
  Jacobian[9] = -5.2915026221291814 + y_residual_tmp;
  Jacobian[578] = nuv1 * (((1.52127765851133 * Uphi2_west_2 + 1.0327955589886444
    * q[19]) - 1.52127765851133 * lambda_12_2) - 1.0327955589886444 *
    lambda_12_4);
  Jacobian[548] = nuv1 * (((3.0 * lambda_12_1 - 1.9166296949998198 *
    Uphi2_west_3) - 3.0 * Uphi2_west_1) + 1.9166296949998198 * lambda_12_3);
  Jacobian[518] = nuv1 * (((1.7320508075688772 * Uphi2_west_1 +
    1.1065666703449764 * Uphi2_west_3) - 1.7320508075688772 * lambda_12_1) -
    1.1065666703449764 * lambda_12_3);
  Jacobian[488] = nuv1 * (((3.4641016151377544 * Uphi2_west_2 +
    3.4016802570830449 * q[19]) - 3.4641016151377544 * lambda_12_2) -
    3.4016802570830449 * lambda_12_4);
  ab_residual_tmp = nuv1 * (((2.6832815729997477 * lambda_12_2 -
    2.6349301969610397 * q[19]) - 2.6832815729997477 * Uphi2_west_2) +
    2.6349301969610397 * lambda_12_4);
  Jacobian[458] = ab_residual_tmp;
  bb_residual_tmp = nuv1 * (((1.5491933384829668 * Uphi2_west_2 +
    1.52127765851133 * q[19]) - 1.5491933384829668 * lambda_12_2) -
    1.52127765851133 * lambda_12_4);
  Jacobian[428] = bb_residual_tmp;
  cb_residual_tmp = nuv1 * (4.58257569495584 * lambda_12_3 - 4.58257569495584 *
    Uphi2_west_3);
  Jacobian[398] = cb_residual_tmp;
  Jacobian[368] = nuv1 * (3.872983346207417 * Uphi2_west_3 - 3.872983346207417 *
    lambda_12_3);
  db_residual_tmp = nuv1 * (3.0 * lambda_12_3 - 3.0 * Uphi2_west_3);
  Jacobian[338] = db_residual_tmp;
  eb_residual_tmp = nuv1 * (1.7320508075688772 * Uphi2_west_3 -
    1.7320508075688772 * lambda_12_3);
  Jacobian[308] = eb_residual_tmp;
  Jacobian[278] = nuv1 * ((((q_residual_tmp - 6.08511063404532 * q[4]) -
    4.1311822359545785 * q[9]) + 1.52127765851133 * lambda_12_2) +
    1.0327955589886444 * lambda_12_4);
  Jacobian[248] = 10.0 + nuv1 * ((((((yb_residual_tmp - 6.0 * Uphi1_west_1) -
    3.8332593899996397 * Uphi1_west_3) - 6.9282032302755088 * q[1]) -
    4.4262666813799054 * q[8]) + 3.0 * lambda_12_1) + 1.9166296949998198 *
    lambda_12_3);
  q_residual_tmp = (ac_residual_tmp + 3.4641016151377544 * Uphi1_west_1) +
    2.2131333406899527 * Uphi1_west_3;
  Jacobian[218] = nuv1 * ((((q_residual_tmp - 6.9282032302755088 * q[0]) -
    4.4262666813799054 * q[7]) + 1.7320508075688772 * lambda_12_1) +
    1.1065666703449764 * lambda_12_3);
  fb_residual_tmp = (bc_residual_tmp + 6.9282032302755088 * Uphi1_west_2) +
    6.80336051416609 * q[9];
  Jacobian[188] = nuv1 * (((fb_residual_tmp - 6.1967733539318672 * q[6]) +
    3.4641016151377544 * lambda_12_2) + 3.4016802570830449 * lambda_12_4);
  gb_residual_tmp = nuv1 * (((((nb_residual_tmp - 5.3665631459994954 *
    Uphi1_west_2) - 5.2698603939220794 * q[9]) - 6.1967733539318672 * q[5]) +
    2.6832815729997477 * lambda_12_2) + 2.6349301969610397 * lambda_12_4);
  Jacobian[158] = -7.745966692414834 + gb_residual_tmp;
  hb_residual_tmp = (pb_residual_tmp + 3.0983866769659336 * Uphi1_west_2) +
    3.04255531702266 * q[9];
  ib_residual_tmp = nuv1 * ((((hb_residual_tmp - 6.1967733539318672 * q[4]) -
    6.08511063404532 * q[9]) + 1.5491933384829668 * lambda_12_2) +
    1.52127765851133 * lambda_12_4);
  Jacobian[128] = ib_residual_tmp;
  jb_residual_tmp = 4.58257569495584 * Uphi1_east_3 - 9.16515138991168 *
    Uphi1_west_3;
  Jacobian[98] = nuv1 * (jb_residual_tmp + 4.58257569495584 * lambda_12_3);
  kb_residual_tmp = 3.872983346207417 * Uphi1_east_3 + 7.745966692414834 *
    Uphi1_west_3;
  Jacobian[68] = nuv1 * (kb_residual_tmp + 3.872983346207417 * lambda_12_3);
  lb_residual_tmp = 4.47213595499958 + nuv1 * (((3.0 * Uphi1_east_3 - 6.0 *
    Uphi1_west_3) - 6.9282032302755088 * q[8]) + 3.0 * lambda_12_3);
  Jacobian[38] = lb_residual_tmp;
  mb_residual_tmp = 1.7320508075688772 * Uphi1_east_3 + 3.4641016151377544 *
    Uphi1_west_3;
  ob_residual_tmp = nuv1 * ((mb_residual_tmp - 6.9282032302755088 * q[7]) +
    1.7320508075688772 * lambda_12_3);
  Jacobian[8] = ob_residual_tmp;
  Jacobian[577] = l_residual_tmp;
  Jacobian[547] = nuv1 * (((1.7320508075688772 * lambda_12_1 -
    1.1065666703449764 * Uphi2_west_3) - 1.7320508075688772 * Uphi2_west_1) +
    1.1065666703449764 * lambda_12_3);
  Jacobian[517] = nuv1 * (((Uphi2_west_1 + 0.63887656499993994 * Uphi2_west_3) -
    lambda_12_1) - 0.63887656499993994 * lambda_12_3);
  l_residual_tmp = nuv1 * (((2.0 * Uphi2_west_2 + 1.9639610121239315 * q[19]) -
    2.0 * lambda_12_2) - 1.9639610121239315 * lambda_12_4);
  Jacobian[487] = l_residual_tmp;
  qb_residual_tmp = nuv1 * (((1.5491933384829668 * lambda_12_2 -
    1.52127765851133 * q[19]) - 1.5491933384829668 * Uphi2_west_2) +
    1.52127765851133 * lambda_12_4);
  Jacobian[457] = qb_residual_tmp;
  rb_residual_tmp = nuv1 * (((0.89442719099991586 * Uphi2_west_2 +
    0.87831006565367986 * q[19]) - 0.89442719099991586 * lambda_12_2) -
    0.87831006565367986 * lambda_12_4);
  Jacobian[427] = rb_residual_tmp;
  Jacobian[397] = nuv1 * (2.6457513110645907 * lambda_12_3 - 2.6457513110645907 *
    Uphi2_west_3);
  sb_residual_tmp = nuv1 * (2.23606797749979 * Uphi2_west_3 - 2.23606797749979 *
    lambda_12_3);
  Jacobian[367] = sb_residual_tmp;
  tb_residual_tmp = nuv1 * (1.7320508075688772 * lambda_12_3 -
    1.7320508075688772 * Uphi2_west_3);
  Jacobian[337] = tb_residual_tmp;
  ub_residual_tmp = nuv1 * (Uphi2_west_3 - lambda_12_3);
  Jacobian[307] = ub_residual_tmp;
  Jacobian[277] = 11.832159566199232 + r_residual_tmp;
  Jacobian[247] = nuv1 * ((q_residual_tmp + 1.7320508075688772 * lambda_12_1) +
    1.1065666703449764 * lambda_12_3);
  Jacobian[217] = 10.0 + nuv1 * (((((Uphi1_east_1 + 0.63887656499993994 *
    Uphi1_east_3) - 2.0 * Uphi1_west_1) - 1.2777531299998799 * Uphi1_west_3) +
    lambda_12_1) + 0.63887656499993994 * lambda_12_3);
  q_residual_tmp = ((2.0 * Uphi1_east_2 + 1.9639610121239315 * q[9]) - 4.0 *
                    Uphi1_west_2) - 3.927922024247863 * q[9];
  Jacobian[187] = nuv1 * ((q_residual_tmp + 2.0 * lambda_12_2) +
    1.9639610121239315 * lambda_12_4);
  r_residual_tmp = nuv1 * ((hb_residual_tmp + 1.5491933384829668 * lambda_12_2)
    + 1.52127765851133 * lambda_12_4);
  Jacobian[157] = r_residual_tmp;
  hb_residual_tmp = nuv1 * (((((0.89442719099991586 * Uphi1_east_2 +
    0.87831006565367986 * q[9]) - 1.7888543819998317 * Uphi1_west_2) -
    1.7566201313073597 * q[9]) + 0.89442719099991586 * lambda_12_2) +
    0.87831006565367986 * lambda_12_4);
  Jacobian[127] = -7.745966692414834 + hb_residual_tmp;
  vb_residual_tmp = 2.6457513110645907 * Uphi1_east_3 + 5.2915026221291814 *
    Uphi1_west_3;
  Jacobian[97] = nuv1 * (vb_residual_tmp + 2.6457513110645907 * lambda_12_3);
  wb_residual_tmp = 2.23606797749979 * Uphi1_east_3 - 4.47213595499958 *
    Uphi1_west_3;
  Jacobian[67] = nuv1 * (wb_residual_tmp + 2.23606797749979 * lambda_12_3);
  mb_residual_tmp = nuv1 * (mb_residual_tmp + 1.7320508075688772 * lambda_12_3);
  Jacobian[37] = mb_residual_tmp;
  xb_residual_tmp = 4.47213595499958 + nuv1 * ((Uphi1_east_3 - 2.0 *
    Uphi1_west_3) + lambda_12_3);
  Jacobian[7] = xb_residual_tmp;
  Jacobian[576] = m_residual_tmp;
  Jacobian[546] = nuv1 * (((3.4641016151377544 * lambda_12_2 -
    3.4016802570830449 * q[19]) - 3.4641016151377544 * Uphi2_west_2) +
    3.4016802570830449 * lambda_12_4);
  Jacobian[516] = l_residual_tmp;
  Jacobian[486] = nuv1 * (((5.0 * Uphi2_west_1 + 4.47213595499958 * Uphi2_west_3)
    - 5.0 * lambda_12_1) - 4.47213595499958 * lambda_12_3);
  Jacobian[456] = nuv1 * (((3.872983346207417 * lambda_12_1 - 3.4641016151377544
    * Uphi2_west_3) - 3.872983346207417 * Uphi2_west_1) + 3.4641016151377544 *
    lambda_12_3);
  l_residual_tmp = nuv1 * (((2.23606797749979 * Uphi2_west_1 + 2.0 *
    Uphi2_west_3) - 2.23606797749979 * lambda_12_1) - 2.0 * lambda_12_3);
  Jacobian[426] = l_residual_tmp;
  Jacobian[396] = nuv1 * (5.9160797830996161 * lambda_12_2 - 5.9160797830996161 *
    Uphi2_west_2);
  m_residual_tmp = nuv1 * (5.0 * Uphi2_west_2 - 5.0 * lambda_12_2);
  Jacobian[366] = m_residual_tmp;
  yb_residual_tmp = nuv1 * (3.872983346207417 * lambda_12_2 - 3.872983346207417 *
    Uphi2_west_2);
  Jacobian[336] = yb_residual_tmp;
  ac_residual_tmp = nuv1 * (2.23606797749979 * Uphi2_west_2 - 2.23606797749979 *
    lambda_12_2);
  Jacobian[306] = ac_residual_tmp;
  Jacobian[276] = nuv1 * ((s_residual_tmp - 13.60672102833218 * q[8]) +
    1.9639610121239315 * lambda_12_3);
  Jacobian[246] = nuv1 * (((((fb_residual_tmp - 13.856406460551018 * q[4]) -
    12.393546707863734 * q[6]) - 13.60672102833218 * q[9]) + 3.4641016151377544 *
    lambda_12_2) + 3.4016802570830449 * lambda_12_4);
  Jacobian[216] = nuv1 * (((q_residual_tmp - 13.856406460551018 * q[5]) + 2.0 *
    lambda_12_2) + 1.9639610121239315 * lambda_12_4);
  Jacobian[186] = 6.0 + nuv1 * ((((((((5.0 * Uphi1_east_1 + 4.47213595499958 *
    Uphi1_east_3) - 10.0 * Uphi1_west_1) - 8.94427190999916 * Uphi1_west_3) -
    13.856406460551018 * q[1]) - 13.60672102833218 * q[3]) - 12.393546707863734 *
    q[8]) + 5.0 * lambda_12_1) + 4.47213595499958 * lambda_12_3);
  q_residual_tmp = (cc_residual_tmp + 7.745966692414834 * Uphi1_west_1) +
    6.9282032302755088 * Uphi1_west_3;
  Jacobian[156] = nuv1 * (((((q_residual_tmp - 15.491933384829668 * q[0]) -
    13.856406460551018 * q[2]) - 13.856406460551018 * q[7]) + 3.872983346207417 *
    lambda_12_1) + 3.4641016151377544 * lambda_12_3);
  s_residual_tmp = ((2.23606797749979 * Uphi1_east_1 + 2.0 * Uphi1_east_3) -
                    4.47213595499958 * Uphi1_west_1) - 4.0 * Uphi1_west_3;
  Jacobian[126] = nuv1 * ((((s_residual_tmp - 15.491933384829668 * q[1]) -
    13.856406460551018 * q[8]) + 2.23606797749979 * lambda_12_1) + 2.0 *
    lambda_12_3);
  fb_residual_tmp = 5.9160797830996161 * Uphi1_east_2 + 11.832159566199232 *
    Uphi1_west_2;
  Jacobian[96] = nuv1 * ((fb_residual_tmp - 13.60672102833218 * q[6]) +
    5.9160797830996161 * lambda_12_2);
  bc_residual_tmp = nuv1 * (((5.0 * Uphi1_east_2 - 10.0 * Uphi1_west_2) -
    13.856406460551018 * q[5]) + 5.0 * lambda_12_2);
  Jacobian[66] = -3.4641016151377544 + bc_residual_tmp;
  cc_residual_tmp = 3.872983346207417 * Uphi1_east_2 + 7.745966692414834 *
    Uphi1_west_2;
  ec_residual_tmp = nuv1 * (((cc_residual_tmp - 15.491933384829668 * q[4]) -
    13.856406460551018 * q[6]) + 3.872983346207417 * lambda_12_2);
  Jacobian[36] = ec_residual_tmp;
  fc_residual_tmp = 2.23606797749979 * Uphi1_east_2 - 4.47213595499958 *
    Uphi1_west_2;
  gc_residual_tmp = nuv1 * ((fc_residual_tmp - 15.491933384829668 * q[5]) +
    2.23606797749979 * lambda_12_2);
  Jacobian[6] = gc_residual_tmp;
  Jacobian[575] = nuv1 * (1.52127765851133 * Uphi2_west_3 - 1.52127765851133 *
    lambda_12_3);
  Jacobian[545] = ab_residual_tmp;
  Jacobian[515] = bb_residual_tmp;
  Jacobian[485] = nuv1 * (((3.872983346207417 * Uphi2_west_1 +
    3.4641016151377544 * Uphi2_west_3) - 3.872983346207417 * lambda_12_1) -
    3.4641016151377544 * lambda_12_3);
  Jacobian[455] = nuv1 * (((3.0 * lambda_12_1 - 2.6832815729997477 *
    Uphi2_west_3) - 3.0 * Uphi2_west_1) + 2.6832815729997477 * lambda_12_3);
  Jacobian[425] = nuv1 * (((1.7320508075688772 * Uphi2_west_1 +
    1.5491933384829668 * Uphi2_west_3) - 1.7320508075688772 * lambda_12_1) -
    1.5491933384829668 * lambda_12_3);
  ab_residual_tmp = nuv1 * (4.58257569495584 * lambda_12_2 - 4.58257569495584 *
    Uphi2_west_2);
  Jacobian[395] = ab_residual_tmp;
  bb_residual_tmp = nuv1 * (3.872983346207417 * Uphi2_west_2 - 3.872983346207417
    * lambda_12_2);
  Jacobian[365] = bb_residual_tmp;
  hc_residual_tmp = nuv1 * (3.0 * lambda_12_2 - 3.0 * Uphi2_west_2);
  Jacobian[335] = hc_residual_tmp;
  ic_residual_tmp = nuv1 * (1.7320508075688772 * Uphi2_west_2 -
    1.7320508075688772 * lambda_12_2);
  Jacobian[305] = ic_residual_tmp;
  Jacobian[275] = nuv1 * ((t_residual_tmp - 6.08511063404532 * q[7]) +
    1.52127765851133 * lambda_12_3);
  Jacobian[245] = 7.745966692414834 + gb_residual_tmp;
  Jacobian[215] = ib_residual_tmp;
  Jacobian[185] = nuv1 * (((q_residual_tmp - 6.9282032302755088 * q[2]) +
    3.872983346207417 * lambda_12_1) + 3.4641016151377544 * lambda_12_3);
  Jacobian[155] = 6.0 + nuv1 * ((((((dc_residual_tmp - 6.0 * Uphi1_west_1) -
    5.3665631459994954 * Uphi1_west_3) - 6.9282032302755088 * q[1]) -
    6.1967733539318672 * q[8]) + 3.0 * lambda_12_1) + 2.6832815729997477 *
    lambda_12_3);
  q_residual_tmp = (lambda_sample_23_4 + 3.4641016151377544 * Uphi1_west_1) +
    3.0983866769659336 * Uphi1_west_3;
  Jacobian[125] = nuv1 * ((((q_residual_tmp - 6.9282032302755088 * q[0]) -
    6.1967733539318672 * q[7]) + 1.7320508075688772 * lambda_12_1) +
    1.5491933384829668 * lambda_12_3);
  t_residual_tmp = 4.58257569495584 * Uphi1_east_2 - 9.16515138991168 *
    Uphi1_west_2;
  Jacobian[95] = nuv1 * (t_residual_tmp + 4.58257569495584 * lambda_12_2);
  gb_residual_tmp = nuv1 * ((cc_residual_tmp - 6.9282032302755088 * q[6]) +
    3.872983346207417 * lambda_12_2);
  Jacobian[65] = gb_residual_tmp;
  ib_residual_tmp = nuv1 * (((3.0 * Uphi1_east_2 - 6.0 * Uphi1_west_2) -
    6.9282032302755088 * q[5]) + 3.0 * lambda_12_2);
  Jacobian[35] = -3.4641016151377544 + ib_residual_tmp;
  cc_residual_tmp = 1.7320508075688772 * Uphi1_east_2 + 3.4641016151377544 *
    Uphi1_west_2;
  dc_residual_tmp = nuv1 * ((cc_residual_tmp - 6.9282032302755088 * q[4]) +
    1.7320508075688772 * lambda_12_2);
  Jacobian[5] = dc_residual_tmp;
  Jacobian[574] = n_residual_tmp;
  Jacobian[544] = qb_residual_tmp;
  Jacobian[514] = rb_residual_tmp;
  Jacobian[484] = l_residual_tmp;
  Jacobian[454] = nuv1 * (((1.7320508075688772 * lambda_12_1 -
    1.5491933384829668 * Uphi2_west_3) - 1.7320508075688772 * Uphi2_west_1) +
    1.5491933384829668 * lambda_12_3);
  Jacobian[424] = nuv1 * (((Uphi2_west_1 + 0.89442719099991586 * Uphi2_west_3) -
    lambda_12_1) - 0.89442719099991586 * lambda_12_3);
  Jacobian[394] = nuv1 * (2.6457513110645907 * lambda_12_2 - 2.6457513110645907 *
    Uphi2_west_2);
  Jacobian[364] = ac_residual_tmp;
  l_residual_tmp = nuv1 * (1.7320508075688772 * lambda_12_2 - 1.7320508075688772
    * Uphi2_west_2);
  Jacobian[334] = l_residual_tmp;
  n_residual_tmp = nuv1 * (Uphi2_west_2 - lambda_12_2);
  Jacobian[304] = n_residual_tmp;
  Jacobian[274] = u_residual_tmp;
  Jacobian[244] = r_residual_tmp;
  Jacobian[214] = 7.745966692414834 + hb_residual_tmp;
  Jacobian[184] = nuv1 * ((s_residual_tmp + 2.23606797749979 * lambda_12_1) +
    2.0 * lambda_12_3);
  Jacobian[154] = nuv1 * ((q_residual_tmp + 1.7320508075688772 * lambda_12_1) +
    1.5491933384829668 * lambda_12_3);
  Jacobian[124] = 6.0 + nuv1 * (((((Uphi1_east_1 + 0.89442719099991586 *
    Uphi1_east_3) - 2.0 * Uphi1_west_1) - 1.7888543819998317 * Uphi1_west_3) +
    lambda_12_1) + 0.89442719099991586 * lambda_12_3);
  q_residual_tmp = 2.6457513110645907 * Uphi1_east_2 + 5.2915026221291814 *
    Uphi1_west_2;
  Jacobian[94] = nuv1 * (q_residual_tmp + 2.6457513110645907 * lambda_12_2);
  r_residual_tmp = nuv1 * (fc_residual_tmp + 2.23606797749979 * lambda_12_2);
  Jacobian[64] = r_residual_tmp;
  s_residual_tmp = nuv1 * (cc_residual_tmp + 1.7320508075688772 * lambda_12_2);
  Jacobian[34] = s_residual_tmp;
  u_residual_tmp = nuv1 * ((Uphi1_east_2 - 2.0 * Uphi1_west_2) + lambda_12_2);
  Jacobian[4] = -3.4641016151377544 + u_residual_tmp;
  Jacobian[573] = nuv1 * (2.6457513110645907 * q[19] - 2.6457513110645907 *
    lambda_12_4);
  Jacobian[543] = cb_residual_tmp;
  Jacobian[513] = nuv1 * (2.6457513110645907 * Uphi2_west_3 - 2.6457513110645907
    * lambda_12_3);
  Jacobian[483] = nuv1 * (5.9160797830996161 * Uphi2_west_2 - 5.9160797830996161
    * lambda_12_2);
  Jacobian[453] = ab_residual_tmp;
  Jacobian[423] = nuv1 * (2.6457513110645907 * Uphi2_west_2 - 2.6457513110645907
    * lambda_12_2);
  Jacobian[393] = nuv1 * (7.0 * lambda_12_1 - 7.0 * Uphi2_west_1);
  Jacobian[363] = nuv1 * (5.9160797830996161 * Uphi2_west_1 - 5.9160797830996161
    * lambda_12_1);
  ab_residual_tmp = nuv1 * (4.58257569495584 * lambda_12_1 - 4.58257569495584 *
    Uphi2_west_1);
  Jacobian[333] = ab_residual_tmp;
  Jacobian[303] = nuv1 * (2.6457513110645907 * Uphi2_west_1 - 2.6457513110645907
    * lambda_12_1);
  Jacobian[273] = nuv1 * ((v_residual_tmp - 10.583005244258363 * q[9]) +
    2.6457513110645907 * lambda_12_4);
  Jacobian[243] = nuv1 * ((jb_residual_tmp - 31.749015732775089 * q[8]) +
    4.58257569495584 * lambda_12_3);
  Jacobian[213] = nuv1 * ((vb_residual_tmp - 10.583005244258363 * q[7]) +
    2.6457513110645907 * lambda_12_3);
  Jacobian[183] = nuv1 * (((fb_residual_tmp - 23.664319132398465 * q[4]) -
    25.701584164627452 * q[6]) + 5.9160797830996161 * lambda_12_2);
  Jacobian[153] = nuv1 * ((t_residual_tmp - 31.749015732775089 * q[5]) +
    4.58257569495584 * lambda_12_2);
  Jacobian[123] = nuv1 * (((q_residual_tmp - 10.583005244258363 * q[4]) -
    23.664319132398465 * q[6]) + 2.6457513110645907 * lambda_12_2);
  Jacobian[93] = 2.0 + nuv1 * ((((7.0 * Uphi1_east_1 - 14.0 * Uphi1_west_1) -
    20.784609690826528 * q[1]) - 24.693678903269511 * q[3]) + 7.0 * lambda_12_1);
  q_residual_tmp = 5.9160797830996161 * Uphi1_east_1 + 11.832159566199232 *
    Uphi1_west_1;
  Jacobian[63] = nuv1 * (((q_residual_tmp - 23.664319132398465 * q[0]) -
    25.701584164627452 * q[2]) + 5.9160797830996161 * lambda_12_1);
  t_residual_tmp = 4.58257569495584 * Uphi1_east_1 - 9.16515138991168 *
    Uphi1_west_1;
  Jacobian[33] = nuv1 * (((t_residual_tmp - 31.749015732775089 * q[1]) -
    20.784609690826528 * q[3]) + 4.58257569495584 * lambda_12_1);
  v_residual_tmp = 2.6457513110645907 * Uphi1_east_1 + 5.2915026221291814 *
    Uphi1_west_1;
  Jacobian[3] = nuv1 * (((v_residual_tmp - 10.583005244258363 * q[0]) -
    23.664319132398465 * q[2]) + 2.6457513110645907 * lambda_12_1);
  Jacobian[572] = o_residual_tmp;
  Jacobian[542] = nuv1 * (3.872983346207417 * lambda_12_3 - 3.872983346207417 *
    Uphi2_west_3);
  Jacobian[512] = sb_residual_tmp;
  Jacobian[482] = m_residual_tmp;
  Jacobian[452] = yb_residual_tmp;
  Jacobian[422] = ac_residual_tmp;
  Jacobian[392] = nuv1 * (5.9160797830996161 * lambda_12_1 - 5.9160797830996161 *
    Uphi2_west_1);
  Jacobian[362] = nuv1 * (5.0 * Uphi2_west_1 - 5.0 * lambda_12_1);
  Jacobian[332] = nuv1 * (3.872983346207417 * lambda_12_1 - 3.872983346207417 *
    Uphi2_west_1);
  m_residual_tmp = nuv1 * (2.23606797749979 * Uphi2_west_1 - 2.23606797749979 *
    lambda_12_1);
  Jacobian[302] = m_residual_tmp;
  Jacobian[272] = w_residual_tmp;
  Jacobian[242] = nuv1 * ((kb_residual_tmp - 15.491933384829668 * q[7]) +
    3.872983346207417 * lambda_12_3);
  Jacobian[212] = nuv1 * ((wb_residual_tmp - 15.491933384829668 * q[8]) +
    2.23606797749979 * lambda_12_3);
  Jacobian[182] = 3.4641016151377544 + bc_residual_tmp;
  Jacobian[152] = ec_residual_tmp;
  Jacobian[122] = gc_residual_tmp;
  Jacobian[92] = nuv1 * ((q_residual_tmp - 13.60672102833218 * q[2]) +
    5.9160797830996161 * lambda_12_1);
  Jacobian[62] = 2.0 + nuv1 * ((((5.0 * Uphi1_east_1 - 10.0 * Uphi1_west_1) -
    13.856406460551018 * q[1]) - 13.60672102833218 * q[3]) + 5.0 * lambda_12_1);
  o_residual_tmp = 3.872983346207417 * Uphi1_east_1 + 7.745966692414834 *
    Uphi1_west_1;
  Jacobian[32] = nuv1 * (((o_residual_tmp - 15.491933384829668 * q[0]) -
    13.856406460551018 * q[2]) + 3.872983346207417 * lambda_12_1);
  q_residual_tmp = 2.23606797749979 * Uphi1_east_1 - 4.47213595499958 *
    Uphi1_west_1;
  Jacobian[2] = nuv1 * ((q_residual_tmp - 15.491933384829668 * q[1]) +
                        2.23606797749979 * lambda_12_1);
  Jacobian[571] = nuv1 * (1.7320508075688772 * q[19] - 1.7320508075688772 *
    lambda_12_4);
  Jacobian[541] = db_residual_tmp;
  Jacobian[511] = eb_residual_tmp;
  Jacobian[481] = bb_residual_tmp;
  Jacobian[451] = hc_residual_tmp;
  Jacobian[421] = ic_residual_tmp;
  Jacobian[391] = ab_residual_tmp;
  Jacobian[361] = nuv1 * (3.872983346207417 * Uphi2_west_1 - 3.872983346207417 *
    lambda_12_1);
  Jacobian[331] = nuv1 * (3.0 * lambda_12_1 - 3.0 * Uphi2_west_1);
  Jacobian[301] = nuv1 * (1.7320508075688772 * Uphi2_west_1 - 1.7320508075688772
    * lambda_12_1);
  Jacobian[271] = nuv1 * ((x_residual_tmp - 6.9282032302755088 * q[9]) +
    1.7320508075688772 * lambda_12_4);
  Jacobian[241] = lb_residual_tmp;
  Jacobian[211] = ob_residual_tmp;
  Jacobian[181] = gb_residual_tmp;
  Jacobian[151] = 3.4641016151377544 + ib_residual_tmp;
  Jacobian[121] = dc_residual_tmp;
  Jacobian[91] = nuv1 * ((t_residual_tmp - 6.9282032302755088 * q[3]) +
    4.58257569495584 * lambda_12_1);
  Jacobian[61] = nuv1 * ((o_residual_tmp - 6.9282032302755088 * q[2]) +
    3.872983346207417 * lambda_12_1);
  Jacobian[31] = 2.0 + nuv1 * (((3.0 * Uphi1_east_1 - 6.0 * Uphi1_west_1) -
    6.9282032302755088 * q[1]) + 3.0 * lambda_12_1);
  o_residual_tmp = 1.7320508075688772 * Uphi1_east_1 + 3.4641016151377544 *
    Uphi1_west_1;
  Jacobian[1] = nuv1 * ((o_residual_tmp - 6.9282032302755088 * q[0]) +
                        1.7320508075688772 * lambda_12_1);
  Jacobian[570] = p_residual_tmp;
  Jacobian[540] = tb_residual_tmp;
  Jacobian[510] = ub_residual_tmp;
  Jacobian[480] = ac_residual_tmp;
  Jacobian[450] = l_residual_tmp;
  Jacobian[420] = n_residual_tmp;
  Jacobian[390] = nuv1 * (2.6457513110645907 * lambda_12_1 - 2.6457513110645907 *
    Uphi2_west_1);
  Jacobian[360] = m_residual_tmp;
  Jacobian[330] = nuv1 * (1.7320508075688772 * lambda_12_1 - 1.7320508075688772 *
    Uphi2_west_1);
  Jacobian[300] = nuv1 * (Uphi2_west_1 - lambda_12_1);
  Jacobian[270] = 5.2915026221291814 + y_residual_tmp;
  Jacobian[240] = mb_residual_tmp;
  Jacobian[210] = xb_residual_tmp;
  Jacobian[180] = r_residual_tmp;
  Jacobian[150] = s_residual_tmp;
  Jacobian[120] = 3.4641016151377544 + u_residual_tmp;
  Jacobian[90] = nuv1 * (v_residual_tmp + 2.6457513110645907 * lambda_12_1);
  Jacobian[60] = nuv1 * (q_residual_tmp + 2.23606797749979 * lambda_12_1);
  Jacobian[30] = nuv1 * (o_residual_tmp + 1.7320508075688772 * lambda_12_1);
  Jacobian[0] = 2.0 + nuv1 * ((Uphi1_east_1 - 2.0 * Uphi1_west_1) + lambda_12_1);
}

/* End of code generation (problem_region_Jacobian_residual_M4_P.c) */
