#include <stdio.h>
#include <stdlib.h>

void use_array_(double *v);
void assign_ptrs_(double ** ptr);

int main()
{
    double *fc_x_comp_sp_1;  // patch 1
    double *fc_x_comp_sp_2;  // patch 2
    double value;

    assign_ptrs_(&fc_x_comp_sp_1);
    value = 1.0;
    use_array_(&value);

    /* 
    assign_ptrs_(&fc_x_comp_sp_2);
    value = 5.0;
    use_array_(&value);
    */

    for(int i = 0; i< 5; i++)
    {
    	printf("%f %f\n",fc_x_comp_sp_1[i],fc_x_comp_sp_1[i]);
    }


    return 0;
}
