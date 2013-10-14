#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED

#include <stdint.h>
#include <inttypes.h>

void** calloc_2D(uint32_t dim1, uint32_t dim2, size_t si);
void free_2D(uint32_t dim1, ...);

void*** calloc_3D(uint32_t dim1, uint32_t dim2, uint32_t dim3, size_t si);
void free_3D(uint32_t dim1, uint32_t dim2, ...);


#endif // MEMORY_H_INCLUDED
