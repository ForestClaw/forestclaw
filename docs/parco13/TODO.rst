Stuff To Address
================
 - Figure 1 ".. between five MPI processes (the three middle ones are shown from black to white)" - Some annotations in the graphic would be helpful, actually the illustration has four colors.

 - ".. as one uniform patch in $mathbb{R}^d$ of dimension $m^d$, where $m$ is ..", The word "dimension" is misleading in this case, because $d$ is already the spatial dimension: "size"?

 - "To this end.." TODO
 
 - TODO: Edges?

 * End of chapter 2, "The exchange of ghost.." - This part should be explained in more detail. How is AMR working in detail. For example: Figure 1 shows a 2:1 refinement, thus refinement of a single patch could affect a whole region of the computational domain. + The underlying equations are known to Forestclaw: How is is the information shared with p4est (equation dependent refinement criteria, refinement and coarsening strategies (interpolation, integration).

 * "Forestclaw uses iterators over all .. neighbor lookups" - That's not obvious for me: Explain in more detail