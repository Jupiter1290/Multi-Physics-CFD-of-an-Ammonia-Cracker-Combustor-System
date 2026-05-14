/*==========================================================================
 * UDF: Zone-Dependent Species Diffusivity — v4
 * System: Coupled Ammonia Decomposer – Hydrogen Combustor
 *
 * Confirmed working:
 *   - Zone ID 6 = fff-porous_decomposer (verified via minimal test)
 *   - Zone ID logic is correct
 *   - Flat constant 2e-4 works → Chapman-Enskog needed better T/P guards
 *
 * Fix in v4:
 *   - Aggressive temperature and pressure guards during initialization
 *   - Fallback to 2e-4 if T or P are unphysical (e.g. during first iter)
 *   - Chapman-Enskog only computed when T and P are trustworthy
 *
 * Species index order (confirmed from Selected Species panel):
 *   0: nh3,  1: n2,  2: h2,  3: o2,  4: h2o (bulk — not called)
 *
 * Operating pressure: 1e5 Pa
 *==========================================================================*/

#include "udf.h"

#define DECOMPOSER_ZONE_ID   6
#define D_SUPPRESSED         1e-20   /* m²/s — PFR behavior in decomposer */
#define D_FALLBACK           2e-4    /* m²/s — safe fallback during init  */
#define T_MIN                300.0   /* K — below this, use fallback       */
#define T_MAX                5000.0  /* K — above this, use fallback       */
#define P_MIN                1e4     /* Pa — below this, use fallback      */

/*--------------------------------------------------------------------------
 * Neufeld (1972) collision integral
 *--------------------------------------------------------------------------*/
static real omega_D(real T_star)
{
    real A=1.06036, B=0.15610, C=0.19300, D_=0.47635;
    real E=1.03587, F=1.52996, G=1.76474, H=3.89411;
    return A/pow(T_star,B) + C/exp(D_*T_star)
         + E/exp(F*T_star) + G/exp(H*T_star);
}

/*--------------------------------------------------------------------------
 * Binary Chapman-Enskog diffusivity [m²/s]
 *--------------------------------------------------------------------------*/
static real D_CE(real T, real P,
                 real sigAB, real epsAB,
                 real MA,    real MB)
{
    real T_star = T / epsAB;
    real OmD    = omega_D(T_star);
    real P_bar  = P / 1e5;
    real D_cm2s = 1.858e-3 * pow(T, 1.5)
                  * sqrt(1.0/MA + 1.0/MB)
                  / (P_bar * sigAB*sigAB * OmD);
    return D_cm2s * 1.0e-4;
}

/*--------------------------------------------------------------------------
 * Species-specific Chapman-Enskog, each paired against N2 bath gas
 *--------------------------------------------------------------------------*/
static real kinetic_theory_diffusivity(cell_t c, Thread *t, int i)
{
    real T = C_T(c, t);
    real P = C_P(c, t) + RP_Get_Real("operating-pressure");

    /* Aggressive guard — if T or P are unphysical (uninitialized cells,
     * first iteration, or boundary cells), return safe fallback.
     * This was the root cause of v3 instability.                          */
    if (T < T_MIN || T > T_MAX || P < P_MIN)
        return D_FALLBACK;

    /* N2 bath gas parameters */
    real sig_N2 = 3.798, eps_N2 = 71.4, M_N2 = 28.014;
    real sigAB, epsAB, MA;

    switch(i)
    {
        case 0:  /* nh3 */
            MA=17.031; sigAB=0.5*(3.758+sig_N2); epsAB=sqrt(481.0*eps_N2);
            break;
        case 1:  /* n2 — self-diffusion */
            MA=28.014; sigAB=sig_N2; epsAB=eps_N2;
            break;
        case 2:  /* h2 */
            MA=2.016;  sigAB=0.5*(2.827+sig_N2); epsAB=sqrt(59.7*eps_N2);
            break;
        case 3:  /* o2 */
            MA=31.999; sigAB=0.5*(3.467+sig_N2); epsAB=sqrt(106.7*eps_N2);
            break;
        default:
            return D_FALLBACK;
    }

    return D_CE(T, P, sigAB, epsAB, MA, M_N2);
}

/*==========================================================================
 * MAIN HOOK
 *==========================================================================*/
DEFINE_DIFFUSIVITY(zone_dependent_diffusivity, c, t, i)
{
    if (THREAD_ID(t) == DECOMPOSER_ZONE_ID)
        return D_SUPPRESSED;
    else
        return kinetic_theory_diffusivity(c, t, i);
}

/*==========================================================================
 * DIAGNOSTIC — uncomment for first 1-2 iterations only
 *
DEFINE_ADJUST(debug_diffusivity, domain)
{
    Thread *t;
    cell_t  c;
    static int done = 0;
    if (done) return;
    done = 1;

    Message("\n--- DIFFUSIVITY DIAGNOSTIC (v4) ---\n");
    Message("  ZoneID | T [K]  |  P [Pa]  | D_h2 [m2/s]      \n");
    Message("--------------------------------------------------\n");

    thread_loop_c(t, domain)
    {
        if (FLUID_THREAD_P(t))
        {
            begin_c_loop(c, t)
            {
                int  zid = THREAD_ID(t);
                real T   = C_T(c,t);
                real P   = C_P(c,t) + RP_Get_Real("operating-pressure");
                real D   = (zid == DECOMPOSER_ZONE_ID)
                           ? D_SUPPRESSED
                           : kinetic_theory_diffusivity(c, t, 2);
                Message("  %6d | %6.1f | %8.0f | %.4e  %s\n",
                        zid, T, P, D,
                        (zid==DECOMPOSER_ZONE_ID) ? "<-- SUPPRESSED" :
                        (T < T_MIN || P < P_MIN)  ? "<-- FALLBACK"   :
                                                    "<-- Chapman-Enskog");
                break;
            }
            end_c_loop(c, t)
        }
    }
    Message("--------------------------------------------------\n\n");
}
 *==========================================================================*/