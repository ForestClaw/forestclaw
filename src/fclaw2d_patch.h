/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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


#ifndef FCLAW2D_PATCH_H
#define FCLAW2D_PATCH_H

#include <forestclaw2d.h>  /* Contains definition of patch-iterator callback */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fclaw2d_patch_vtable          fclaw2d_patch_vtable_t;
typedef struct fclaw2d_patch_data            fclaw2d_patch_data_t;
typedef struct fclaw2d_patch_transform_data  fclaw2d_patch_transform_data_t;

typedef enum
{
	FCLAW2D_BUILD_FOR_GHOST_AREA_COMPUTED = 0,
	FCLAW2D_BUILD_FOR_GHOST_AREA_PACKED,
	FCLAW2D_BUILD_FOR_UPDATE,
	FCLAW2D_BUILD_CUSTOM
} fclaw2d_build_mode_t;


/* The user patch (clawpatch, for example) is stored as a 'user_patch', below. */
struct fclaw2d_patch_data
{
    fclaw2d_patch_relation_t face_neighbors[4];
    fclaw2d_patch_relation_t corner_neighbors[4];
    int corners[4];
    int block_corner_count[4];
    int on_coarsefine_interface;
    int has_finegrid_neighbors;
    int neighbors_set;
    int blockno;

    int patch_idx;    /* for debugging!! */
    int block_idx;

    void *user_patch; /* User patch is virualized */
    void *user_data;  /* Data the user may want to attach to each patch */

};

struct fclaw2d_patch_transform_data
{
	struct fclaw2d_patch *this_patch;
	struct fclaw2d_patch *neighbor_patch;
	int transform[9];
	int icorner;
	int based;      /* 1 for cell-centered (1 .. mx); 0 for nodes (0 .. mx) */
	int is_block_corner;
	int block_iface;   /* -1 for interior faces or block corners */

	struct fclaw2d_global *glob;
	void* user;  /* Used by individual patches */
};


struct fclaw2d_global;
struct fclaw2d_domain;
struct fclaw2d_patch;


/* ------------------------------ Creating/deleting patches --------------------------- */

#if 0
void fclaw2d_patch_data_new(struct fclaw2d_global* glob,
							struct fclaw2d_patch* this_patch);
#endif                            

void fclaw2d_patch_reset_data(struct fclaw2d_global* glob,
							  struct fclaw2d_patch* old_patch,
							  struct fclaw2d_patch* new_patch,
							  int blockno,int old_patchno, int new_patchno);


void fclaw2d_patch_data_delete(struct fclaw2d_global *glob,
							   struct fclaw2d_patch *patch);

void fclaw2d_patch_build(struct fclaw2d_global *glob,
						 struct fclaw2d_patch *this_patch,
						 int blockno,
						 int patchno,
						 void *user);

void fclaw2d_patch_build_from_fine(struct fclaw2d_global *glob,
								   struct fclaw2d_patch *fine_patches,
								   struct fclaw2d_patch *coarse_patch,
								   int blockno,
								   int coarse_patchno,
								   int fine0_patchno,
								   fclaw2d_build_mode_t build_mode);


/* ---------------------------- Solver specific functions ----------------------------- */

void fclaw2d_patch_initialize(struct fclaw2d_global *glob,
							  struct fclaw2d_patch *this_patch,
							  int this_block_idx,
							  int this_patch_idx);

void fclaw2d_patch_physical_bc(struct fclaw2d_global *glob,
							   struct fclaw2d_patch *this_patch,
							   int this_block_idx,
							   int this_patch_idx,
							   double t,
							   double dt,
							   int *intersects_bc,
							   int time_interp);

double fclaw2d_patch_single_step_update(struct fclaw2d_global *glob,
                                        struct fclaw2d_patch *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        double t,
                                        double dt, void* user);

void fclaw2d_patch_set_rhs(struct fclaw2d_global *glob,
                           struct fclaw2d_patch *patch,
                           int blockno,
                           int patchno);


/* -------------------------------- time stepping ------------------------------------- */

void fclaw2d_patch_restore_step(struct fclaw2d_global* glob,
								struct fclaw2d_patch* this_patch);

void fclaw2d_patch_save_step(struct fclaw2d_global* glob,
							 struct fclaw2d_patch* this_patch);


void fclaw2d_patch_setup_timeinterp(struct fclaw2d_global *glob,
									struct fclaw2d_patch *this_patch,
									double alpha);


/* ---------------------------- Ghost filling - patch specific ------------------------ */

void fclaw2d_patch_copy_face(struct fclaw2d_global* glob,
							 struct fclaw2d_patch *this_patch,
							 struct fclaw2d_patch *neighbor_patch,
							 int iface,
							 int time_interp,
							 struct fclaw2d_patch_transform_data *transform_data);

void fclaw2d_patch_average_face(struct fclaw2d_global* glob,
								struct fclaw2d_patch *coarse_patch,
								struct fclaw2d_patch *fine_patch,
								int idir,
								int iface_coarse,
								int RefineFactor,
								int refratio,
								int time_interp,
								int igrid,
								struct fclaw2d_patch_transform_data* transform_data);

void fclaw2d_patch_interpolate_face(struct fclaw2d_global* glob,
									struct fclaw2d_patch *coarse_patch,
									struct fclaw2d_patch *fine_patch,
									int idir,
									int iside,
									int RefineFactor,
									int refratio,
									int time_interp,
									int igrid,
									struct fclaw2d_patch_transform_data* transform_data);


void fclaw2d_patch_copy_corner(struct fclaw2d_global* glob,
							   struct fclaw2d_patch *this_patch,
							   struct fclaw2d_patch *corner_patch,
							   int coarse_blockno,
							   int fine_blockno,
							   int is_block_corner,
							   int icorner,
							   int time_interp,
							   struct fclaw2d_patch_transform_data *transform_data);

void fclaw2d_patch_average_corner(struct fclaw2d_global* glob,
								  struct fclaw2d_patch *coarse_patch,
								  struct fclaw2d_patch *fine_patch,
								  int coarse_blockno,
								  int fine_blockno,
								  int is_block_corner,
								  int coarse_corner,
								  int time_interp,
								  struct fclaw2d_patch_transform_data* transform_data);



void fclaw2d_patch_interpolate_corner(struct fclaw2d_global* glob,
									  struct fclaw2d_patch* coarse_patch,
									  struct fclaw2d_patch* fine_patch,
									  int coarse_blockno,
									  int fine_blockno,
									  int is_block_corner,
									  int coarse_corner,
									  int time_interp,
									  struct fclaw2d_patch_transform_data* transform_data);

/* -------------------------- Transform functions  ------------------------------------ */
void fclaw2d_patch_create_user_data(struct fclaw2d_global* glob,
                                    struct fclaw2d_patch* patch);

void fclaw2d_patch_destroy_user_data(struct fclaw2d_global* glob,
                                     struct fclaw2d_patch* patch);

/* -------------------------- Transform functions (typedefs) -------------------------- */

void fclaw2d_patch_transform_init_data(struct fclaw2d_global* glob,
									   struct fclaw2d_patch* patch,
									   int blockno, int patchno,
									   struct fclaw2d_patch_transform_data *tdata);

/* Had to go with name "blockface", since "fclaw2d_patch_transform_face" is already 
used in forestclaw2d */
void fclaw2d_patch_transform_blockface(int faceno, int rfaceno,
									   int ftransform[]);

void fclaw2d_patch_transform_blockface_intra(int ftransform[]);
  

/* ------------------------------- Regridding functions ------------------------------- */

int fclaw2d_patch_tag4refinement(struct fclaw2d_global *glob,
								 struct fclaw2d_patch *this_patch,
								 int blockno, int patchno,
								 int initflag);

int fclaw2d_patch_tag4coarsening(struct fclaw2d_global *glob,
								 struct fclaw2d_patch *fine_patches,
								 int blockno,
								 int patchno,
                                 int initflag);

void fclaw2d_patch_interpolate2fine(struct fclaw2d_global *glob,
									struct fclaw2d_patch* coarse_patch,
									struct fclaw2d_patch* fine_patches,
									int this_blockno, int coarse_patchno,
									int fine0_patchno);

void fclaw2d_patch_average2coarse(struct fclaw2d_global *glob,
								  struct fclaw2d_patch *fine_patches,
								  struct fclaw2d_patch *coarse_patch,
								  int blockno, int fine0_patchno,
								  int coarse_patchno);

/* ----------------------------- Parallel ghost patches ------------------------------- */

size_t fclaw2d_patch_ghost_packsize(struct fclaw2d_global* glob);


void fclaw2d_patch_local_ghost_alloc(struct fclaw2d_global* glob,
									 void** q);

void fclaw2d_patch_local_ghost_free(struct fclaw2d_global* glob,
									void **q);

void fclaw2d_patch_local_ghost_pack(struct fclaw2d_global *glob,
									struct fclaw2d_patch *this_patch,
									void *patch_data,
									int time_interp);

void fclaw2d_patch_remote_ghost_build(struct fclaw2d_global *glob,
									  struct fclaw2d_patch *this_patch,
									  int blockno,
									  int patchno,
									  void *user);

void fclaw2d_patch_remote_ghost_unpack(struct fclaw2d_global* glob,
									   struct fclaw2d_patch* this_patch,
									   int this_block_idx, int this_patch_idx,
									   void *qdata, int time_interp);


void fclaw2d_patch_remote_ghost_delete(struct fclaw2d_global *glob,
									   struct fclaw2d_patch *ghost_patch);

/* -------------------------------- Parallel partitioning ----------------------------- */


void fclaw2d_patch_partition_pack(struct fclaw2d_global *glob,
								  struct fclaw2d_patch *this_patch,
								  int this_block_idx,
								  int this_patch_idx,
								  void *pack_data_here);

void fclaw2d_patch_partition_unpack(struct fclaw2d_global *glob,  /* contains old domain */
									struct fclaw2d_domain *new_domain,  
									struct fclaw2d_patch *this_patch,
									int this_block_idx,
									int this_patch_idx,
									void *packed_data);

size_t fclaw2d_patch_partition_packsize(struct fclaw2d_global* glob);


/* ------------------------------ Time syncing funtions ------------------------------- */
void fclaw2d_patch_time_sync_f2c(struct fclaw2d_global* glob,
								 struct fclaw2d_patch *coarse_patch,
								 struct fclaw2d_patch *fine_patch,
								 int coarse_blockno, int fine_blockno,
								 int coarse_patchno, 
								 int idir,
								 int igrid,
								 int iface_coarse,
								 int time_interp,
								 struct fclaw2d_patch_transform_data* transform_data);

void fclaw2d_patch_time_sync_samesize(struct fclaw2d_global* glob,
                                      struct fclaw2d_patch *this_patch,
                                      struct fclaw2d_patch *neighbor_patch,
                                      int iface, int idir,
                                      struct fclaw2d_patch_transform_data *transform_data);


void fclaw2d_patch_time_sync_reset(struct fclaw2d_global* glob,
                                   struct fclaw2d_patch* this_patch,
                                   int coarse_level,
                                   int reset_mode);

#if 0
void fclaw2d_patch_time_sync_reset_f2c(struct fclaw2d_global* glob,
                                       struct fclaw2d_patch *patch,
                                       int coarse_level);

void fclaw2d_patch_time_sync_reset_samesize(struct fclaw2d_global* glob, 
                                            struct fclaw2d_patch *patch);
#endif                                            

/* ------------------------------ Misc access functions ------------------------------- */

/* I don't completely trust this routine */
void fclaw2d_patch_get_info2(struct fclaw2d_domain * domain,
							struct fclaw2d_patch * this_patch,
							int *this_block_idx, int *this_patch_idx,
							int *global_num, int *level);

void*
fclaw2d_patch_get_user_patch(struct fclaw2d_patch* patch);

#if 0
struct fclaw2d_patch_data*
fclaw2d_patch_get_user_data(struct fclaw2d_patch* patch);
#endif


void* fclaw2d_patch_get_user_patch(struct fclaw2d_patch* patch);


/* ---------------------- Creating/deleting patches (typedefs) ------------------------ */

typedef void* (*fclaw2d_patch_new_t)();

typedef void (*fclaw2d_patch_delete_t)(void *user_patch);

typedef void (*fclaw2d_patch_build_t)(struct fclaw2d_global *glob,
									  struct fclaw2d_patch *this_patch,
									  int blockno,
									  int patchno,
									  void *user);

typedef void (*fclaw2d_patch_build_from_fine_t)(struct fclaw2d_global *glob,
												struct fclaw2d_patch *fine_patches,
												struct fclaw2d_patch *coarse_patch,
												int blockno,
												int coarse_patchno,
												int fine0_patchno,
												fclaw2d_build_mode_t build_mode);

typedef void (*fclaw2d_patch_setup_t)(struct fclaw2d_global *glob,
									  struct fclaw2d_patch *this_patch,
									  int this_block_idx,
									  int this_patch_idx);

/* --------------------- Solver specific functions (typedefs) ------------------------- */

typedef void (*fclaw2d_patch_initialize_t)(struct fclaw2d_global *glob,
										   struct fclaw2d_patch *this_patch,
										   int this_block_idx,
										   int this_patch_idx);

typedef void (*fclaw2d_patch_physical_bc_t)(struct fclaw2d_global *glob,
											struct fclaw2d_patch *this_patch,
											int this_block_idx,
											int this_patch_idx,
											double t,
											double dt,
											int *intersects_bc,
											int time_interp);

typedef double (*fclaw2d_patch_single_step_update_t)(struct fclaw2d_global *glob,
                                                     struct fclaw2d_patch *this_patch,
                                                     int this_block_idx,
                                                     int this_patch_idx,
                                                     double t,
                                                     double dt,
                                                     void* user);


typedef void (*fclaw2d_patch_rhs_t)(struct fclaw2d_global *glob,
                                    struct fclaw2d_patch *patch,
                                    int blockno,
                                    int patchno);


/* ----------------------------- Time stepping (typedefs) ----------------------------- */

typedef void (*fclaw2d_patch_setup_timeinterp_t)(struct fclaw2d_global *glob,
												 struct fclaw2d_patch *this_patch,
												 double alpha);

typedef void (*fclaw2d_patch_restore_step_t)(struct fclaw2d_global *glob,
											 struct fclaw2d_patch* this_patch);

typedef void (*fclaw2d_patch_save_step_t)(struct fclaw2d_global *glob,
										  struct fclaw2d_patch* this_patch);


/* --------------------- Ghost filling - patch specific (typedefs) -------------------- */

typedef void (*fclaw2d_patch_copy_face_t)(struct fclaw2d_global* glob,
										  struct fclaw2d_patch *this_patch,
										  struct fclaw2d_patch *neighbor_patch,
										  int iface,
										  int time_interp,
										  struct fclaw2d_patch_transform_data 
										  *transform_data);

typedef void (*fclaw2d_patch_average_face_t)(struct fclaw2d_global* glob,
											 struct fclaw2d_patch *coarse_patch,
											 struct fclaw2d_patch *fine_patch,
											 int idir,
											 int iface_coarse,
											 int RefineFactor,
											 int refratio,
											 int time_interp,
											 int igrid,
											 struct fclaw2d_patch_transform_data
											 *transform_data);

typedef void (*fclaw2d_patch_interpolate_face_t)(struct fclaw2d_global* glob,
												 struct fclaw2d_patch *coarse_patch,
												 struct fclaw2d_patch *fine_patch,
												 int idir,
												 int iside,
												 int RefineFactor,
												 int refratio,
												 int a_time_interp,
												 int igrid,
												 struct fclaw2d_patch_transform_data
												 *transform_data);

/* These three are identical;  do we need three separate functions? */

typedef void (*fclaw2d_patch_corner_t)(struct fclaw2d_global* glob,
									   struct fclaw2d_patch *this_patch,
									   struct fclaw2d_patch *corner_patch,
									   int coarse_blockno,
									   int fine_blockno,
									   int icorner,
									   int time_interp,
									   struct fclaw2d_patch_transform_data 
									   *transform_data);
	

/* -------------------------- Transform functions (typedefs) -------------------------- */

typedef void (*fclaw2d_patch_transform_init_data_t)(struct fclaw2d_global* glob,
													struct fclaw2d_patch* patch,
													int blockno, int patchno,
													struct fclaw2d_patch_transform_data *tdata);
  
typedef void (*fclaw2d_patch_transform_blockface_t)(int faceno, int rfaceno,
											   int ftransform[]);

typedef void (*fclaw2d_patch_transform_blockface_intra_t)(int ftransform[]);

/* ------------------------- Regridding functions (typedefs) -------------------------- */

typedef int (*fclaw2d_patch_tag4refinement_t)(struct fclaw2d_global *glob,
											  struct fclaw2d_patch *this_patch,
											  int this_block_idx, int this_patch_idx,
											  int initflag);

typedef int (*fclaw2d_patch_tag4coarsening_t)(struct fclaw2d_global *glob,
											   struct fclaw2d_patch *this_patch,
											   int this_blockno,
											   int this_patchno,
                                               int initflag);

typedef void (*fclaw2d_patch_interpolate2fine_t)(struct fclaw2d_global *glob,
												 struct fclaw2d_patch *coarse_patch,
												 struct fclaw2d_patch* fine_patches,
												 int this_blockno, int coarse_patchno,
												 int fine_patchno);

typedef void (*fclaw2d_patch_average2coarse_t)(struct fclaw2d_global *glob,
											   struct fclaw2d_patch *fine_siblings,
											   struct fclaw2d_patch *coarse_patch,
											   int blockno, int fine_patchno,
											   int coarse_patchno);

/* -------------------------- Parallel ghost patches (typedefs) ----------------------- */

typedef size_t (*fclaw2d_patch_ghost_packsize_t)(struct fclaw2d_global* glob);

typedef void (*fclaw2d_patch_local_ghost_pack_t)(struct fclaw2d_global *glob,
												 struct fclaw2d_patch *this_patch,
												 void *patch_data,
												 int time_interp);

typedef void (*fclaw2d_patch_local_ghost_alloc_t)(struct fclaw2d_global* glob,
												 void** q);

typedef void (*fclaw2d_patch_local_ghost_free_t)(struct fclaw2d_global* glob,
												 void **q);

typedef void (*fclaw2d_patch_remote_ghost_build_t)(struct fclaw2d_global *glob,
												   struct fclaw2d_patch *this_patch,
												   int blockno,
												   int patchno,
												   void *user);

typedef void (*fclaw2d_patch_remote_ghost_setup_t)(struct fclaw2d_global *glob,
												   struct fclaw2d_patch *this_patch,
												   int this_block_idx,
												   int this_patch_idx);


typedef void (*fclaw2d_patch_remote_ghost_unpack_t)(struct fclaw2d_global *glob,
													struct fclaw2d_patch* this_patch,
													int this_block_idx, int this_patch_idx,
													void *qdata, int time_interp);

typedef void (*fclaw2d_patch_remote_ghost_delete_t)(void *user_patch);



/* ----------------------------- Partitioning (typedefs) ------------------------------ */

/* Returns size, in bytes, i.e. psize = mx*my*sizeof(double) */ 
typedef size_t (*fclaw2d_patch_partition_packsize_t)(struct fclaw2d_global* glob);


typedef void (*fclaw2d_patch_partition_pack_t)(struct fclaw2d_global *glob,
											   struct fclaw2d_patch *this_patch,
											   int this_block_idx,
											   int this_patch_idx,
											   void *pack_data_here);

typedef void (*fclaw2d_patch_partition_unpack_t)(struct fclaw2d_global *glob,
												 struct fclaw2d_domain *new_domain,
												 struct fclaw2d_patch *this_patch,
												 int this_block_idx,
												 int this_patch_idx,
												 void *unpack_data_from_here);

/* ----------------------------- Conservative updates --------------------------------- */

typedef void (*fclaw2d_patch_time_sync_f2c_t)(struct fclaw2d_global* glob,
											  struct fclaw2d_patch *coarse_patch,
											  struct fclaw2d_patch *fine_patch,
											  int coarse_blockno, int fine_blockno,
											  int coarse_patchno, 
											  int idir,
											  int igrid,
											  int iface_coarse,
											  int time_interp,
											  struct fclaw2d_patch_transform_data
											  *transform_data);

typedef void (*fclaw2d_patch_time_sync_samesize_t)(struct fclaw2d_global* glob,
                                                   struct fclaw2d_patch* this_patch,
                                                   struct fclaw2d_patch* neighbor_patch,
                                                   int iface, int idir,
                                                   struct fclaw2d_patch_transform_data 
                                                   *transform_data);

typedef void (*fclaw2d_patch_time_sync_reset_t)(struct fclaw2d_global *glob, 
                                                struct fclaw2d_patch *this_patch,
                                                int coarse_level,
                                                int reset_mode);

typedef void (*fclaw2d_patch_time_sync_reset_f2c_t)(struct fclaw2d_global *glob, 
                                                    struct fclaw2d_patch *this_patch,
                                                    int coarse_level);

typedef void (*fclaw2d_patch_time_sync_reset_samesize_t)(struct fclaw2d_global *glob,
                                                         struct fclaw2d_patch *patch);

/* ------------------------------  User data functions -------------------------------- */

typedef void (*fclaw2d_patch_create_user_data_t)(struct fclaw2d_global *glob, 
                                              struct fclaw2d_patch *patch);

typedef void (*fclaw2d_patch_destroy_user_data_t)(struct fclaw2d_global* glob,
                                                  struct fclaw2d_patch* patch);

/* ---------------------------------  Access functions -------------------------------- */

typedef void* (*fclaw2d_patch_metric_patch_t)(struct fclaw2d_patch *patch);

/* ----------------------------------- Virtual table  --------------------------------- */
struct fclaw2d_patch_vtable
{
    /* Creating/deleting/building patches */
    fclaw2d_patch_new_t                   patch_new;
    fclaw2d_patch_delete_t                patch_delete;
    fclaw2d_patch_build_t                 build;
    fclaw2d_patch_build_from_fine_t       build_from_fine;
    fclaw2d_patch_setup_t                 setup;

    /* Return metric patch from patch struct. */
    fclaw2d_patch_metric_patch_t          metric_patch;

    /* Set user data */
    fclaw2d_patch_create_user_data_t      create_user_data;
    fclaw2d_patch_destroy_user_data_t     destroy_user_data;

    /* Solver functions */
    fclaw2d_patch_initialize_t            initialize;
    fclaw2d_patch_physical_bc_t           physical_bc;
    fclaw2d_patch_single_step_update_t    single_step_update;
    fclaw2d_patch_rhs_t                   rhs;

    /* Time stepping */
    fclaw2d_patch_restore_step_t          restore_step;
    fclaw2d_patch_save_step_t             save_step;
    fclaw2d_patch_setup_timeinterp_t      setup_timeinterp;


    /* regridding functions */
    fclaw2d_patch_tag4refinement_t        tag4refinement;
    fclaw2d_patch_tag4coarsening_t        tag4coarsening;
    fclaw2d_patch_average2coarse_t        average2coarse;
    fclaw2d_patch_interpolate2fine_t      interpolate2fine;

    /* Time syncing functions for conservation */
    fclaw2d_patch_time_sync_f2c_t         time_sync_f2c;      /* f2c = fine to coarse */
    fclaw2d_patch_time_sync_samesize_t    time_sync_samesize;
    fclaw2d_patch_time_sync_reset_t       time_sync_reset;    /* Virtualized for each patch */


    /* ghost filling functions */
    fclaw2d_patch_copy_face_t             copy_face;
    fclaw2d_patch_average_face_t          average_face;
    fclaw2d_patch_interpolate_face_t      interpolate_face;

    /* Block face and interior corners */
    fclaw2d_patch_corner_t                copy_corner;
    fclaw2d_patch_corner_t                average_corner;
    fclaw2d_patch_corner_t                interpolate_corner;

    /* Block corners */
    fclaw2d_patch_corner_t                copy_block_corner;
    fclaw2d_patch_corner_t                average_block_corner;
    fclaw2d_patch_corner_t                interpolate_block_corner;

    /* Transform functions */
    fclaw2d_patch_transform_init_data_t        transform_init_data;
    fclaw2d_patch_transform_blockface_t        transform_face;
    fclaw2d_patch_transform_blockface_intra_t  transform_face_intra;

    /* Ghost packing functions (for parallel use) */
    fclaw2d_patch_ghost_packsize_t        ghost_packsize;
    fclaw2d_patch_local_ghost_pack_t      local_ghost_pack;
    fclaw2d_patch_local_ghost_alloc_t     local_ghost_alloc;
    fclaw2d_patch_local_ghost_free_t      local_ghost_free;

    fclaw2d_patch_remote_ghost_build_t    remote_ghost_build;
    fclaw2d_patch_remote_ghost_setup_t    remote_ghost_setup;   /* Remote ghost patches */
    fclaw2d_patch_remote_ghost_unpack_t   remote_ghost_unpack;
    fclaw2d_patch_remote_ghost_delete_t   remote_ghost_delete;  /* Delete remote ghosts */

    /* Parallel load balancing (partitioning) */
    fclaw2d_patch_partition_pack_t         partition_pack;
    fclaw2d_patch_partition_unpack_t       partition_unpack;
    fclaw2d_patch_partition_packsize_t     partition_packsize;

    int is_set;
};

fclaw2d_patch_vtable_t* fclaw2d_patch_vt();

void fclaw2d_patch_vtable_initialize();


/* ------------------------------ Misc access functions ------------------------------- */
void fclaw2d_patch_get_info(struct fclaw2d_domain * domain,
                            struct fclaw2d_patch * patch,
                            int blockno, int patchno,
                            int *global_num, int* local_num, 
                            int *level);

void*
fclaw2d_patch_get_user_patch(struct fclaw2d_patch* patch);

struct fclaw2d_patch_data*
fclaw2d_patch_get_patch_data(struct fclaw2d_patch* patch);

void* fclaw2d_patch_get_user_data(struct fclaw2d_global* glob,
                                  struct fclaw2d_patch* this_patch);


void* fclaw2d_patch_get_user_patch(struct fclaw2d_patch* patch);

void* fclaw2d_patch_metric_patch(struct fclaw2d_patch *patch);

int fclaw2d_patch_get_blockno(struct fclaw2d_patch* this_patch);

/* Misc. user data */
void* fclaw2d_patch_user_data(struct fclaw2d_global* glob,
                              struct fclaw2d_patch* this_patch);

void fclaw2d_patch_set_user_data(struct fclaw2d_global* glob,
                                 struct fclaw2d_patch* this_patch, 
                                 void* user);
//int fclaw2d_patch_get_patchno(fclaw2d_patch_t* this_patch);


/* ------------------  Miscellaneous functions (mostly internal) ---------------------- */
typedef void (*fclaw2d_patch_iterator_t) (struct fclaw2d_global * glob, int level,
										  fclaw2d_patch_callback_t pcb, void *user);

int fclaw2d_patch_on_parallel_boundary (const struct fclaw2d_patch * patch);


void fclaw2d_patch_set_face_type(struct fclaw2d_patch *patch, int iface,
								 fclaw2d_patch_relation_t face_type);

void fclaw2d_patch_set_corner_type(struct fclaw2d_patch *patch, int icorner,
								   fclaw2d_patch_relation_t corner_type);

void fclaw2d_patch_set_missing_corner(struct fclaw2d_patch *patch, int icorner);

fclaw2d_patch_relation_t fclaw2d_patch_get_face_type(struct fclaw2d_patch* patch,
														int iface);
fclaw2d_patch_relation_t fclaw2d_patch_get_corner_type(struct fclaw2d_patch* patch,
														  int icorner);

int fclaw2d_patch_corner_is_missing(struct fclaw2d_patch* patch,
									int icorner);

void fclaw2d_patch_neighbors_set(struct fclaw2d_patch* patch);

void fclaw2d_patch_neighbors_reset(struct fclaw2d_patch* patch);

int fclaw2d_patch_has_finegrid_neighbors(struct fclaw2d_patch *patch);

int fclaw2d_patch_on_coarsefine_interface(struct fclaw2d_patch *patch);

int* fclaw2d_patch_block_corner_count(struct fclaw2d_global *glob,
									  struct fclaw2d_patch* this_patch);

void fclaw2d_patch_set_block_corner_count(struct fclaw2d_global *glob,
										  struct fclaw2d_patch* this_patch,
										  int icorner, int block_corner_count);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
