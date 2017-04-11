#include <stdio.h>
#include <stdlib.h>

void use_array_(double *v);
void store_ptrs_(double ** ptr);
void copy_ptrs2mod_(double ** ptr);

int main()
{
    double *fc_x_comp_sp_1;  // patch 1
    double *fc_x_comp_sp_2;  // patch 2
    double value;

    value = 1.0;
    use_array_(&value);    
    store_ptrs_(&fc_x_comp_sp_1);

    
    value = 5.0;
    use_array_(&value);
    store_ptrs_(&fc_x_comp_sp_2);



    for(int i = 0; i< 5; i++)
    {
    	printf("%f %f\n",fc_x_comp_sp_1[i],fc_x_comp_sp_2[i]);
        fc_x_comp_sp_1[i] = 47;
    }

    copy_ptrs2mod_(&fc_x_comp_sp_1);
    copy_ptrs2mod_(&fc_x_comp_sp_2);

    return 0;
}
