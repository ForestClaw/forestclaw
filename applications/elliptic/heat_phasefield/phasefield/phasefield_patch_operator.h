#include <ThunderEgg.h>
#include "phasefield_options.h"
#include "fc2d_thunderegg_options.h"

class phasefield : public ThunderEgg::PatchOperator<2>
{
    private:
    static double lambda;

    public:
    ThunderEgg::Vector<2> phi_n;

    const phasefield_options_t *phase_opt;

    phasefield(fclaw_global_t *glob,
               const ThunderEgg::Vector<2>& phi_n_in,
               const ThunderEgg::Domain<2>& domain,
               const ThunderEgg::GhostFiller<2>& ghost_filler);

    phasefield(const fc2d_thunderegg_options *mg_opt,const phasefield_options* phase_opt,
               const ThunderEgg::Vector<2>& phi_n_in,
               const ThunderEgg::Domain<2>& domain,
               const ThunderEgg::GhostFiller<2>& ghost_filler);

    phasefield* clone() const override;

    void applySinglePatch(const ThunderEgg::PatchInfo<2>& pinfo,
                          const ThunderEgg::PatchView<const double, 2>& us,
                          const ThunderEgg::PatchView<double, 2>& fs) const override;

    void applySinglePatchWithInternalBoundaryConditions(const ThunderEgg::PatchInfo<2>& pinfo,
                                                        const ThunderEgg::PatchView<const double, 2>& us,
                                                        const ThunderEgg::PatchView<double, 2>& fs) const override;

    void modifyRHSForInternalBoundaryConditions(const ThunderEgg::PatchInfo<2> &pinfo,
	                                            const ThunderEgg::PatchView<const double, 2> &u_view,
	                                            const ThunderEgg::PatchView<double, 2> &f_view) const override;


    int s[4]; /* Determines sign when applying BCs */

    static void setLambda(double lambda);

    static double getLambda();

    const int * getS() const;
};