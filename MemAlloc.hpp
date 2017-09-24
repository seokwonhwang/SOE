#ifndef MEMALLOC_H
#define MEMALLOC_H

#include <iostream>
#include <new>

template <typename T>
T* memAlloc(size_t n1)
{
    T *v;
    try { 
        v = new T[n1];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nVector Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }
 
    return v;
}
template <typename T>
T* memAllocInit(size_t n1, T init)
{
    T *v;
    try { 
        v = new T[n1];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nVector Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }
    for(size_t i=0; i<n1; i++)
	v[i] = init;
 
    return v;
}
template <typename T>
T** memAlloc(size_t n1, size_t n2)
{   
    T **m;
    
    try {
        m = new T*[n1];
    }   
    catch(std::bad_alloc &ba) {
        std::cerr << "\nMatrix Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }   
    try {
        m[0] = new T[n1*n2];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nMatrix Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }
    for(size_t j=1; j<n1; j++) {
        m[j] = &m[0][n2*j];
    }

    return m;
}
template <typename T>
T** memAllocInit(size_t n1, size_t n2, T init)
{   
    T **m;
    
    try {
        m = new T*[n1];
    }   
    catch(std::bad_alloc &ba) {
        std::cerr << "\nMatrix Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }   
    try {
        m[0] = new T[n1*n2];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nMatrix Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }
    for(size_t j=1; j<n1; j++) {
        m[j] = &m[0][n2*j];
    }

    for(size_t i=0; i<n1; i++) {
	for(size_t j=0; j<n2; j++) {
	    m[i][j] = init;
	}
    }

    return m;
}
template <typename T>
T*** memAlloc(size_t n1, size_t n2, size_t n3)
{
    T ***t;
    
    try {
        t = new T**[n1];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nTensor Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }
    try {
        t[0] = new T*[n1*n2];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nTensor Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }
    try {
        t[0][0] = new T[n1*n2*n3];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nTensor Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }

    for(size_t j=1; j<n1; j++) {
        t[j] = &t[0][n2*j];
    }
    for(size_t j=0; j<n1; j++) {
        for(size_t k=0; k<n2; k++) {
            t[j][k] = &t[0][0][n2*n3*j+n3*k];
        }
    }
    return t;
}
template <typename T>
T*** memAllocInit(size_t n1, size_t n2, size_t n3, T init)
{
    T ***t;
    
    try {
        t = new T**[n1];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nTensor Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }
    try {
        t[0] = new T*[n1*n2];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nTensor Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }
    try {
        t[0][0] = new T[n1*n2*n3];
    }
    catch(std::bad_alloc &ba) {
        std::cerr << "\nTensor Memory allocation Error: "
                  << ba.what() << std::endl;
        return NULL;
    }

    for(size_t j=1; j<n1; j++) {
        t[j] = &t[0][n2*j];
    }
    for(size_t j=0; j<n1; j++) {
        for(size_t k=0; k<n2; k++) {
            t[j][k] = &t[0][0][n2*n3*j+n3*k];
        }
    }

    for(size_t j=0; j<n1; j++) {
	for(size_t k=0; k<n2; k++) {
	    for(size_t l=0; l<n3; l++) {
		t[j][k][l] = init;
	    }
	}
    }

    return t;
}
template <typename T>
void memFree(T* v)
{
    delete[] v;
}
template <typename T>
void memFree(T** m)
{
    delete[] (*m);
    delete[] m;
}
template <typename T>
void memFree(T*** t)
{
    delete[] (**t);
    delete[] (*t);
    delete[] t;
}

#endif

