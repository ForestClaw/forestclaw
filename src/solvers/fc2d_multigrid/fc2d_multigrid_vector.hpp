#include <ThunderEgg/Vector.h>
#include <fclaw2d_global.h>
/**
 * @brief Wrapper class for forestclaw data acess
 */
class fc2d_multigrid_vector : public ThunderEgg::Vector<2> {
   private:
    /**
     * @brief the number of cells in each direction on each patch
     */
    std::array<int, 2> ns;
    /**
     * @brief number of ghost cells
     */
    int mbc;
    /**
     * @brief strides in each direction
     */
    std::array<int, 2> strides;
    /**
     * @brief the number of equations
     */
    int meqn;
    /**
     * @brief stride to next eqn in patch
     */
    int eqn_stride;
    /**
     * @brief pointers to patch data, index corresponds to local_num in
     * forestclaw
     */
    std::vector<double *> patch_data;
    /**
     * @brief fclaw iterator to initialize path_data vector
     */
    static void enumeratePatchData(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t *patch, int blockno,
                                   int patchno, void *user);
    ThunderEgg::LocalData<2> getLocalDataPriv(int component_index, int local_patch_id) const;

   public:
    fc2d_multigrid_vector(fclaw2d_global_t *glob);

    ThunderEgg::LocalData<2> getLocalData(int component_index, int local_patch_id) override;
    const ThunderEgg::LocalData<2> getLocalData(int component_index, int local_patch_id) const override;
};
