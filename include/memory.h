#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED

void** calloc_2D(int dim1, int dim2, size_t si);
void free_2D(int dim1, ...);

void*** calloc_3D(int dim1, int dim2, int dim3, size_t si);
void free_3D(int dim1, int dim2, ...);


#endif // MEMORY_H_INCLUDED
