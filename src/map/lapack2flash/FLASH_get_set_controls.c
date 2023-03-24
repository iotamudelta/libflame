/*

    Copyright (C) 2023, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLASH // Start lapack2flash

static dim_t        flash_blocksize = 1024;
static unsigned int flash_n_threads = 1;
static dim_t        flash_depth = 1;

/*---------------------Get/Set Helper---------------------*/

FLA_Error FLASH_set_preferred_blocksize( dim_t blocksize )
{
    flash_blocksize = blocksize;

    return FLA_SUCCESS;
}

dim_t FLASH_get_preferred_blocksize( void )
{
    return flash_blocksize;
}

FLA_Error FLASH_set_n_preferred_threads( unsigned int threads )
{
    flash_n_threads = threads;

    #ifdef FLA_ENABLE_HIP
        int device_count = 1;
        FLASH_Queue_available_devices_hip( &device_count );

        // Pass number of threads unless it's more than available devices
        FLASH_Queue_set_num_threads( min( threads, device_count) );
    #else
        FLASH_Queue_set_num_threads( threads );
    #endif

    return FLA_SUCCESS;
}

unsigned int FLASH_get_n_preferred_threads( void )
{
    // Returning the number of preferred threads, not actual number.
    return flash_n_threads;
}

FLA_Error FLASH_set_depth( dim_t depth )
{
    flash_depth = depth;

    return FLA_SUCCESS;
}

dim_t FLASH_get_depth( void )
{
    return flash_depth;
}

#endif // End lapack2flash