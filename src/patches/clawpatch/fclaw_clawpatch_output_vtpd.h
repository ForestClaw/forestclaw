#ifndef MYVTKWRITER_H
#define MYVTKWRITER_H
void fclaw_clawpatch_output_vtpd_to_file (struct fclaw_global * glob, const char* basename);
void fclaw_clawpatch_output_vtpd (struct fclaw_global * glob, int iframe);

#endif
