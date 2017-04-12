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
    fc_patch_t p1, p2;
    double value1, value2;

    value1 = 1.0;
    value2 = 5.0;

    use_array_(&value1,&value2);    
    store_ptrs_(&p1.x_comp_sp, &p1.concen);

    use_array_(&value1, &value2);
    store_ptrs_(&p2.x_comp_sp, &p2.concen);

    for(int i = 0; i < 5; i++)
    {
    	printf("%5.0f %5.0f %5.0f %5.0f\n",p1.x_comp_sp[i], p1.concen[i],p2.x_comp_sp[i],p2.concen[i]);
        p1.x_comp_sp[i] = 47;
    }

    copy_ptrs2mod_(&p1.x_comp_sp, &p1.concen);
    copy_ptrs2mod_(&p2.x_comp_sp, &p2.concen);

    deallocate_arrays_(&p1.x_comp_sp, &p1.concen);
    deallocate_arrays_(&p2.x_comp_sp, &p2.concen);

    for(int i = 0; i < 5; i++)
    { 
        /* This should seg-fault */
        // fc_x_comp_sp_1[i] = 47;
    }

    return 0;
}
