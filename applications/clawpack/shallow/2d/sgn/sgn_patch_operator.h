#include "sgn_options.h"

#include <ThunderEgg/PatchOperator.h>
#include "fc2d_thunderegg_options.h"

class sgn : public ThunderEgg::PatchOperator<2>
{
    private:

    public:
    std::shared_ptr<const ThunderEgg::Vector<2>> q_n;

    const sgn_options_t *sgn_opt;

    sgn(fclaw2d_global_t *glob,
               std::shared_ptr<const ThunderEgg::Vector<2>> q_n_in,
               std::shared_ptr<const ThunderEgg::Domain<2>> domain,
               std::shared_ptr<const ThunderEgg::GhostFiller<2>> ghost_filler);

    sgn(const fc2d_thunderegg_options *mg_opt,const sgn_options_t* sgn_opt,
               std::shared_ptr<const ThunderEgg::Vector<2>> q_n_in,
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

    const int * getS() const;
};