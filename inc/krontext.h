// Copyright: Tom Nickson 2015

template <typename T>
class Krontext { // Context, but in my KRONMAT library yeah?
    T * hzero;
    T * hone;
    T * dzero;
    T * done;
    
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    cusolverDnHandle_t solverHandle;
    cublasHandle_t blasHandle;
#else 
    int solverHandle = 0;
    int blasHandle = 0;
#endif
    Krontext() {
           
    }
    
}


