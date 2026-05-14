/*==========================================================================
 * UDF: Flue Gas Recirculation Coupling - v2
 * System: Coupled Ammonia Decomposer - Hydrogen Combustor
 *
 * Fix from v1: INLET_T(), MASS_FLOW(), INLET_YI() are read-only macros
 * in Fluent's API and cannot be assigned to directly.
 * Correct approach: use F_PROFILE via DEFINE_PROFILE, or use the
 * Cortex RP variable system to push values to the BC thread.
 *
 * Strategy used here:
 *   DEFINE_ADJUST  - computes averaged outlet conditions each iteration
 *   DEFINE_PROFILE - called by Fluent when it evaluates the inlet BC,
 *                    reads the stored values and applies them face-by-face
 *
 * This is the standard Fluent pattern for dynamic inlet BCs.
 *
 * Recirculation ratio: 0.4327
 * Combustor outlet: combustor_outlet (pressure-outlet, id=22)
 * Flue inlet:       flue_inlet       (mass-flow-inlet, id=14)
 *
 * Species index order (confirmed):
 *   0: nh3,  1: n2,  2: h2,  3: o2,  4: h2o (bulk)
 *
 * Fluent hookup:
 *   DEFINE_ADJUST  - Define > UDF > Function Hooks > Adjust
 *   DEFINE_PROFILE - flue_inlet BC panel:
 *                    Temperature    - udf flue_inlet_temperature
 *                    Mass Flux      - udf flue_inlet_massflux
 *                    Species Yi     - udf flue_inlet_yi_0 .. yi_3
 *==========================================================================*/

#include "udf.h"

/*--------------------------------------------------------------------------
 * CONFIGURATION
 *--------------------------------------------------------------------------*/
#define RECIRCULATION_RATIO   0.2635   /* stoichiometric case: MATLAB phi_sim_titanium */
#define COMBUSTOR_OUTLET_ID   22
#define FLUE_INLET_ID         14
#define N_SPECIES             4       /* nh3=0, n2=1, h2=2, o2=3          */
#define URF                   0.5     /* under-relaxation, reduce if noisy */

/*--------------------------------------------------------------------------
 * Global storage - computed in DEFINE_ADJUST, read in DEFINE_PROFILE
 * Initialised to sensible values (close to expected converged state)
 *--------------------------------------------------------------------------*/
static real g_T_flue    = 2490.8;    /* K   - from converged 26.35% RR run  */
static real g_mdot_flue = 1.411e-6;  /* kg/s - 5.3534e-6 x 0.2635         */
static real g_Yi_flue[N_SPECIES] = {0.0015, 0.6888, 0.0, 0.0021};
                                      /* nh3, n2, h2, o2 - from converged 26.35% run */

/*==========================================================================
 * DEFINE_ADJUST: runs every iteration, computes outlet averages,
 * stores into global variables for DEFINE_PROFILE to consume.
 *==========================================================================*/
DEFINE_ADJUST(flue_recirculation_coupling, domain)
{
    Thread *t_outlet = Lookup_Thread(domain, COMBUSTOR_OUTLET_ID);

    if (!t_outlet)
    {
        Message("WARNING: flue_recirculation_coupling - "
                "combustor_outlet (id=%d) not found.\n", COMBUSTOR_OUTLET_ID);
        return;
    }

    /* --- Accumulators --- */
    real sum_area = 0.0;
    real sum_mdot = 0.0;
    real sum_T    = 0.0;
    real sum_Yi[N_SPECIES];
    face_t f;
    int i;

    for (i = 0; i < N_SPECIES; i++) sum_Yi[i] = 0.0;

    /* --- Loop over outlet faces ---
     * F_FLUX(f,t) returns the mass flux [kg/s] through each face,
     * already accounting for 2D axisymmetric geometry internally.
     * Sign convention: negative = leaving domain (outlet).                */
    begin_f_loop(f, t_outlet)
    {
        real area_vec[ND_ND];
        F_AREA(area_vec, f, t_outlet);
        real A = NV_MAG(area_vec);

        real T_f    = F_T(f, t_outlet);
        real flux_f = F_FLUX(f, t_outlet);          /* kg/s, negative at outlet */
        if (flux_f < 0.0) flux_f = -flux_f;         /* take absolute value      */

        sum_area += A;
        sum_T    += T_f * A;
        sum_mdot += flux_f * 2.0 * M_PI;           /* 2D axi: F_FLUX is per-radian, multiply by 2pi */

        for (i = 0; i < N_SPECIES; i++)
            sum_Yi[i] += F_YI(f, t_outlet, i) * A;
    }
    end_f_loop(f, t_outlet)

    if (sum_area < 1.0e-20) return;

    real mdot_out = sum_mdot;

    /* --- Averaged outlet values --- */
    real T_out = sum_T / sum_area;
    real Yi_out[N_SPECIES];
    for (i = 0; i < N_SPECIES; i++)
        Yi_out[i] = sum_Yi[i] / sum_area;

    /* --- Sanity guard --- */
    if (T_out < 300.0 || T_out > 5000.0) return;  /* raised for stoichiometric flame */

    /* --- Under-relaxed update of global storage ---
     * Temperature and species: intensive - no RR scaling
     * Mass flow: extensive - scale by RR                                  */
    g_T_flue    = URF * T_out
                + (1.0 - URF) * g_T_flue;

    g_mdot_flue = URF * (mdot_out * RECIRCULATION_RATIO)
                + (1.0 - URF) * g_mdot_flue;

    for (i = 0; i < N_SPECIES; i++)
        g_Yi_flue[i] = URF * Yi_out[i]
                     + (1.0 - URF) * g_Yi_flue[i];

    /* --- Monitor every 10 iterations --- */
    if ((N_ITER % 10) == 0)
    {
        Message("\n--- Flue Recirculation (iter %d) ---\n", N_ITER);
        Message("  Outlet:    mdot=%.4e kg/s  T=%.1f K\n", mdot_out, T_out);
        Message("  Flue inlet: mdot=%.4e kg/s  T=%.1f K\n",
                g_mdot_flue, g_T_flue);
        for (i = 0; i < N_SPECIES; i++)
            Message("  Yi[%d]=%.4f\n", i, g_Yi_flue[i]);
        Message("------------------------------------\n\n");
    }
}

/*==========================================================================
 * DEFINE_PROFILE: called by Fluent each iteration when evaluating the
 * flue_inlet boundary. Reads from global storage set by DEFINE_ADJUST.
 *
 * Hook these in the flue_inlet BC panel under:
 *   Temperature - udf flue_inlet_temperature
 *   Mass Flux   - udf flue_inlet_massflux
 *   Yi for each species - udf flue_inlet_yi_0 .. flue_inlet_yi_3
 *==========================================================================*/
DEFINE_PROFILE(flue_inlet_temperature, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = g_T_flue;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(flue_inlet_massflux, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = g_mdot_flue;
    }
    end_f_loop(f, t)
}

/* Species profiles - one per transported species */
DEFINE_PROFILE(flue_inlet_yi_0, t, i)   /* nh3 */
{
    face_t f;
    begin_f_loop(f, t)
    { F_PROFILE(f, t, i) = g_Yi_flue[0]; }
    end_f_loop(f, t)
}

DEFINE_PROFILE(flue_inlet_yi_1, t, i)   /* n2 */
{
    face_t f;
    begin_f_loop(f, t)
    { F_PROFILE(f, t, i) = g_Yi_flue[1]; }
    end_f_loop(f, t)
}

DEFINE_PROFILE(flue_inlet_yi_2, t, i)   /* h2 */
{
    face_t f;
    begin_f_loop(f, t)
    { F_PROFILE(f, t, i) = g_Yi_flue[2]; }
    end_f_loop(f, t)
}

DEFINE_PROFILE(flue_inlet_yi_3, t, i)   /* o2 */
{
    face_t f;
    begin_f_loop(f, t)
    { F_PROFILE(f, t, i) = g_Yi_flue[3]; }
    end_f_loop(f, t)
}
