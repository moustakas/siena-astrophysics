
void force_treeallocate(int maxnodes, int maxpart);
void force_treefree(void);


int force_treeevaluate(double * pos, double softening, double * acc);

void force_update_node_recursive(int no, int sib);
int force_treebuild(void);

int force_treeevaluate_direct(double * pos, double softening, double * acc);
