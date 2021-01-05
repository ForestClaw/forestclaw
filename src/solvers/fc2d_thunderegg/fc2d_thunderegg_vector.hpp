/*
  Copyright (c) 2019-2020 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <ThunderEgg/Vector.h>

// #include <fclaw2d_global.h>

/* Avoid circular dependencies */
struct fclaw2d_patch;
struct fclaw2d_domain;
struct fclaw2d_global;

typedef enum fc2d_thunderegg_data_choice
{
    RHS=0,
    SOLN,
    STORE_STATE,
}  fc2d_thunderegg_data_choice_t;

/**
 * @brief Wrapper class for forestclaw data acess
 */
class fc2d_thunderegg_vector : public ThunderEgg::Vector<2> {
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
    int mfields;    

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
    static void enumeratePatchData(struct fclaw2d_domain *domain,
                                   struct fclaw2d_patch *patch, int blockno,
                                   int patchno, void *user);

    static void set_start_value(struct fclaw2d_domain *domain,
                                struct fclaw2d_patch *patch, int blockno,
                                int patchno, void *user);

#if 1
    static void store_state(struct fclaw2d_domain *domain,
                            struct fclaw2d_patch *patch, int blockno,
                            int patchno, void *user);
#endif                            

    ThunderEgg::LocalData<2> getLocalDataPriv(int component_index, int local_patch_id) const;

   public:
    fc2d_thunderegg_vector(struct fclaw2d_global *glob, int data_choice, int num_components);

    ThunderEgg::LocalData<2> getLocalData(int component_index, int local_patch_id) override;
    const ThunderEgg::LocalData<2> getLocalData(int component_index, int local_patch_id) const override;
};
