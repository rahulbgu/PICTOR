#include <thrust/device_vector.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>

extern "C" {
//Scan for integer arrays
void scan_int_wrapper( int *data, int N)
{
// Wrap raw pointer with a device_ptr
thrust::device_ptr <int> dev_ptr(data);
// Use device_ptr in Thrust sort algorithm
thrust::inclusive_scan(dev_ptr, dev_ptr+N, dev_ptr);
}

}