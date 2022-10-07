/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


FLA_Error FLASH_Svd_ext_u_unb_var1( FLA_Svd_type jobu, FLA_Svd_type jobv,
                                    dim_t n_iter_max,
                                    FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                                    dim_t k_accum,
                                    dim_t b_alg,
                                    dim_t depth,
                                    dim_t* blocksize  )
{
    FLA_Error    r_val = FLA_SUCCESS;
    FLA_Datatype dt;
    FLA_Datatype dt_real;
    FLA_Datatype dt_comp;
    FLA_Obj      scale, T, S, rL, rR, d, e, G, H, C; // C is dummy.
    dim_t        m_A, n_A, min_m_n;
    dim_t        n_GH;
    double       crossover_ratio = 17.0 / 9.0;
    FLA_Bool     u_is_formed = FALSE, 
                 v_is_formed = FALSE;
    int          apply_scale;

    n_GH    = k_accum;

    m_A     = FLA_Obj_length( A );
    n_A     = FLA_Obj_width( A );
    min_m_n = min( m_A, n_A );
    dt      = FLA_Obj_datatype( A );
    dt_real = FLA_Obj_datatype_proj_to_real( A );
    dt_comp = FLA_Obj_datatype_proj_to_complex( A );

    // Create matrices to hold block Householder transformations.
    FLASH_Bidiag_UT_create_T( A, &T, &S, depth, blocksize );

    // Create vectors to hold the realifying scalars.
    if ( FLA_Obj_is_complex( A ) )
    {
        FLASH_Obj_create( dt, min_m_n, 1, depth, blocksize, &rL );
        FLASH_Obj_create( dt, min_m_n, 1, depth, blocksize, &rR );
    }

    // Create vectors to hold the diagonal and sub-diagonal.
    FLASH_Obj_create( dt_real, min_m_n,   1, depth, blocksize, &d );
    FLASH_Obj_create( dt_real, min_m_n-1, 1, depth, blocksize, &e );

    // Create matrices to hold the left and right Givens scalars.
    FLASH_Obj_create( dt_comp, min_m_n-1, n_GH, depth, blocksize, &G );
    FLASH_Obj_create( dt_comp, min_m_n-1, n_GH, depth, blocksize, &H );

    // Create a real scaling factor.
    FLASH_Obj_create( dt_real, 1, 1, depth, blocksize, &scale );

    // Scale matrix A if necessary. 
    FLA_Max_abs_value( A, scale );
    apply_scale =
      ( FLA_Obj_gt( scale, FLA_OVERFLOW_SQUARE_THRES  ) == TRUE ) -     
      ( FLA_Obj_lt( scale, FLA_UNDERFLOW_SQUARE_THRES ) == TRUE ); 
    
    if ( apply_scale )
      FLASH_Scal( apply_scale > 0 ? FLA_SAFE_MIN : FLA_SAFE_INV_MIN, A );   

    if ( m_A < crossover_ratio * n_A )
    {
        // Reduce the matrix to bidiagonal form.
        // Apply scalars to rotate elements on the superdiagonal to the real domain.
        // Extract the diagonal and superdiagonal from A.
        FLA_Bidiag_UT( A, T, S );
        if ( FLA_Obj_is_complex( A ) )
            FLA_Bidiag_UT_realify( A, rL, rR );
        FLA_Bidiag_UT_extract_real_diagonals( A, d, e );

        // Form U and V.
        if ( u_is_formed == FALSE )
        {
            switch ( jobu )
            {
            case FLA_SVD_VECTORS_MIN_OVERWRITE:
                if ( jobv != FLA_SVD_VECTORS_NONE )
                    FLASH_Bidiag_UT_form_V_ext( FLA_UPPER_TRIANGULAR, A, S, FLA_NO_TRANSPOSE, V, depth, blocksize );
                v_is_formed = TRUE; // For this case, V should be formed here.
                U = A;
            case FLA_SVD_VECTORS_ALL:
            case FLA_SVD_VECTORS_MIN_COPY:
                FLASH_Bidiag_UT_form_U_ext( FLA_UPPER_TRIANGULAR, A, T, FLA_NO_TRANSPOSE, U, depth, blocksize );
                u_is_formed = TRUE;
                break;
            case FLA_SVD_VECTORS_NONE:
                // Do nothing
                break;
            }
        }
        if ( v_is_formed == FALSE )
        {
            if ( jobv == FLA_SVD_VECTORS_MIN_OVERWRITE )
            {
                FLASH_Bidiag_UT_form_V_ext( FLA_UPPER_TRIANGULAR, A, S, FLA_CONJ_TRANSPOSE, A, depth, blocksize );
                v_is_formed = TRUE; /* and */
                V = A; // This V is actually V^H.

                // V^H -> V
                FLA_Obj_flip_base( &V );
                FLA_Obj_flip_view( &V );
                if ( FLA_Obj_is_complex( A ) )
                    FLA_Conjugate( V );
            }
            else if ( jobv != FLA_SVD_VECTORS_NONE )
            {
                FLASH_Bidiag_UT_form_V_ext( FLA_UPPER_TRIANGULAR, A, S, FLA_NO_TRANSPOSE, V, depth, blocksize );
                v_is_formed = TRUE;
            }
        }

        // For complex matrices, apply realification transformation.
        if ( FLA_Obj_is_complex( A ) && jobu != FLA_SVD_VECTORS_NONE )
        {
            FLA_Obj UL, UR;
            FLA_Part_1x2( U,   &UL, &UR,   min_m_n, FLA_LEFT );
            FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, rL, UL );
        }
        if ( FLA_Obj_is_complex( A ) && jobv != FLA_SVD_VECTORS_NONE )
        {
            FLA_Obj VL, VR;
            FLA_Part_1x2( V,   &VL, &VR,   min_m_n, FLA_LEFT );
            FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, rR, VL );
        }

        // Perform a singular value decomposition on the upper bidiagonal matrix.
        r_val = FLA_Bsvd_ext_opt_var1( n_iter_max,
                                       d, e, G, H,
                                       jobu, U, jobv, V,
                                       FALSE, C, // C is not referenced
                                       b_alg );
    }
    else // if ( crossover_ratio * n_A <= m_A )
    {
        FLA_Obj TQ, R;
        FLA_Obj AT,
                AB;

        // Perform a QR factorization on A.
        FLASH_QR_UT_create_T( A, &TQ, depth, blocksize );
        FLASH_QR_UT( A, TQ );

        // Set the lower triangle of R to zero and then copy the upper
        // triangle of A to R.
        FLA_Part_2x1( A,   &AT,
                           &AB,   n_A, FLA_TOP );
        FLASH_Obj_create( dt, n_A, n_A, depth, blocksize, &R );
        FLA_Setr( FLA_LOWER_TRIANGULAR, FLA_ZERO, R );
        FLASH_Copyr( FLA_UPPER_TRIANGULAR, AT, R );

        // Form U; if necessary overwrite on A.
        if ( u_is_formed == FALSE )
        {
            switch ( jobu )
            {
            case FLA_SVD_VECTORS_MIN_OVERWRITE:
                U = A;
            case FLA_SVD_VECTORS_ALL:
            case FLA_SVD_VECTORS_MIN_COPY:
                FLASH_QR_UT_form_Q( A, TQ, U );
                u_is_formed = TRUE;
                break;
            case FLA_SVD_VECTORS_NONE:
                // Do nothing
                break;
            }
        }
        FLASH_Obj_free( &TQ );

        // Reduce the matrix to bidiagonal form.
        // Apply scalars to rotate elements on the superdiagonal to the real domain.
        // Extract the diagonal and superdiagonal from A.
        FLA_Bidiag_UT( R, T, S );
        if ( FLA_Obj_is_complex( R ) )
            FLA_Bidiag_UT_realify( R, rL, rR );
        FLA_Bidiag_UT_extract_real_diagonals( R, d, e );

        if ( v_is_formed == FALSE )
        {
            if ( jobv == FLA_SVD_VECTORS_MIN_OVERWRITE )
            {
                FLASH_Bidiag_UT_form_V_ext( FLA_UPPER_TRIANGULAR, R, S, FLA_CONJ_TRANSPOSE, AT, depth, blocksize );
                v_is_formed = TRUE; /* and */
                V = AT; // This V is actually V^H.

                // V^H -> V
                FLA_Obj_flip_base( &V );
                FLA_Obj_flip_view( &V );
                if ( FLA_Obj_is_complex( A ) )
                    FLA_Conjugate( V );
            }
            else if ( jobv != FLA_SVD_VECTORS_NONE )
            {
                FLASH_Bidiag_UT_form_V_ext( FLA_UPPER_TRIANGULAR, R, S, FLA_NO_TRANSPOSE, V, depth, blocksize );
                v_is_formed = TRUE;
            }
        }

        // Apply householder vectors U in R.
        FLASH_Bidiag_UT_form_U_ext( FLA_UPPER_TRIANGULAR, R, T, FLA_NO_TRANSPOSE, R, depth, blocksize );

        // Apply the realifying scalars in rL and rR to U and V, respectively.
        if ( FLA_Obj_is_complex( A ) && jobu != FLA_SVD_VECTORS_NONE )
        {
            FLA_Obj RL, RR;
            FLA_Part_1x2( R,   &RL, &RR,   min_m_n, FLA_LEFT );
            FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, rL, RL );
        }
        if ( FLA_Obj_is_complex( A ) && jobv != FLA_SVD_VECTORS_NONE )
        {
            FLA_Obj VL, VR;
            FLA_Part_1x2( V,   &VL, &VR,   min_m_n, FLA_LEFT );
            FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, rR, VL );
        }

        // Perform a singular value decomposition on the bidiagonal matrix.
        r_val = FLA_Bsvd_ext_opt_var1( n_iter_max,
                                       d, e, G, H,
                                       jobu, R, jobv, V,
                                       FALSE, C,
                                       b_alg );

        // Multiply R into U, storing the result in A and then copying back
        // to U.
        if ( jobu != FLA_SVD_VECTORS_NONE )
        {
            FLA_Obj UL, UR;
            FLA_Part_1x2( U,   &UL, &UR,   min_m_n, FLA_LEFT );

            if ( jobu == FLA_SVD_VECTORS_MIN_OVERWRITE || 
                 jobv == FLA_SVD_VECTORS_MIN_OVERWRITE )
            {
                FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, UL, &C );
                FLASH_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                          FLA_ONE, UL, R, FLA_ZERO, C );
                FLASH_Copy( C, UL );
                FLASH_Obj_free( &C );
            }
            else
            {
                FLASH_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                          FLA_ONE, UL, R, FLA_ZERO, A );
                FLASH_Copy( A, UL );
            }
        }
        FLASH_Obj_free( &R );
    }

    // Copy the converged eigenvalues to the output vector.
    FLASH_Copy( d, s );

    // No sort is required as it is applied on FLA_Bsvd.

    if ( apply_scale ) 
      FLASH_Scal( apply_scale < 0 ? FLA_SAFE_MIN : FLA_SAFE_INV_MIN, s ); 

    // When V is overwritten, flip it again.
    if ( jobv == FLA_SVD_VECTORS_MIN_OVERWRITE )
    {
        // Always apply conjugation first wrt dimensions used; then, flip base.
        if ( FLA_Obj_is_complex( V ) )
            FLA_Conjugate( V );
        FLA_Obj_flip_base( &V );
    }

    FLASH_Obj_free( &scale );
    FLASH_Obj_free( &T );
    FLASH_Obj_free( &S );

    if ( FLA_Obj_is_complex( A ) )
    {
        FLASH_Obj_free( &rL );
        FLASH_Obj_free( &rR );
    }

    FLASH_Obj_free( &d );
    FLASH_Obj_free( &e );
    FLASH_Obj_free( &G );
    FLASH_Obj_free( &H );

    return r_val;
}

