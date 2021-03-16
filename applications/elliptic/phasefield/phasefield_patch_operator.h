#include <ThunderEgg.h>
#include "phasefield_options.h"
#include "fc2d_thunderegg_options.h"

class phasefield : public ThunderEgg::PatchOperator<2>
{
    private:
    static double lambda;

    public:
    std::shared_ptr<const ThunderEgg::Vector<2>> phi_n;

    const phasefield_options_t *phase_opt;

    phasefield(fclaw2d_global_t *glob,
               std::shared_ptr<const ThunderEgg::Vector<2>> phi_n_in,
               std::shared_ptr<const ThunderEgg::Domain<2>> domain,
               std::shared_ptr<const ThunderEgg::GhostFiller<2>> ghost_filler);

    phasefield(const fc2d_thunderegg_options *mg_opt,const phasefield_options* phase_opt,
               std::shared_ptr<const ThunderEgg::Vector<2>> phi_n_in,
               std::shared_ptr<const ThunderEgg::Domain<2>> domain,
               std::shared_ptr<const ThunderEgg::GhostFiller<2>> ghost_filler);

    void applySinglePatch(std::shared_ptr<const ThunderEgg::PatchInfo<2>> pinfo,
                          const std::vector<ThunderEgg::LocalData<2>> &us,
                          std::vector<ThunderEgg::LocalData<2>> &fs,
                          bool interior_dirichlet) const override;
    void addGhostToRHS(std::shared_ptr<const ThunderEgg::PatchInfo<2>> pinfo,
                       const std::vector<ThunderEgg::LocalData<2>> &us,
                       std::vector<ThunderEgg::LocalData<2>> &fs) const override;

    int s[4]; /* Determines sign when applying BCs */

    static void setLambda(double lambda);

    static double getLambda();

    const int * getS() const;
};