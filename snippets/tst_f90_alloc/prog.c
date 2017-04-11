#include <stdio.h>
#include <stdlib.h>

void use_array_(double *v1, double *v2);
void store_ptrs_(double ** ptr1, double** ptr2);
void copy_ptrs2mod_(double ** ptr1, double **ptr2);
void deallocate_arrays_(double** ptr1, double ** ptr2);

typedef struct fc_patch
{
    double * x_comp_sp;
    double * concen;
} fc_patch_t;


int main()
{
    double *fc_x_comp_sp_1, *fc_concen_1;  // patch 1
    double *fc_x_comp_sp_2, *fc_concen_2;  // patch 2
    fc_patch_t p1, p2;
    double value1, value2;

    value1 = 1.0;
    value2 = 5.0;
    use_array_(&value1,&value2);    
    store_ptrs_(&p1.x_comp_sp, &p1.concen);

    value1 = 2.0;
    value2 = 11.0;
    use_array_(&value1, &value2);
    store_ptrs_(&p2.x_comp_sp, &p2.concen);

    printf("Print values of x_comp_sp after allocation and usage\n");
    for(int i = 0; i < 5; i++)
    {
    	printf("%f %f\n",p1.x_comp_sp[i], p2.x_comp_sp[i]);
        p1.x_comp_sp[i] = 47;
    }

    printf("\n");
    printf("Copying ptrs back to module (after x_comp_sp_1 set to 47)\n");
    copy_ptrs2mod_(&p1.x_comp_sp, &p1.concen);
    copy_ptrs2mod_(&p2.x_comp_sp, &p2.concen);

    printf("\n");
    printf("Deallocating arrays\n\n");
    deallocate_arrays_(&p1.x_comp_sp, &p1.concen);
    deallocate_arrays_(&p2.x_comp_sp, &p2.concen);

    for(int i = 0; i < 5; i++)
    { 
        /* This should seg-fault */
        // fc_x_comp_sp_1[i] = 47;
    }

    return 0;
}
