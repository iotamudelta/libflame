/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_QR_UT_form_Q( FLA_Obj A, FLA_Obj T, FLA_Obj Q )
{
    FLA_Error r_val = FLA_SUCCESS;
    FLA_Obj   QTL, QTR,
              QBL, QBR;
    FLA_Obj   W;
    dim_t     b;

    if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
        FLA_QR_UT_form_Q_check( A, T, Q );

    if ( FLA_Obj_is_overlapped( A, Q ) == FALSE )
    {
        // If A and Q are different objects, Q is explicitly formed with A.

        // Set Q identify
        FLA_Set_to_identity( Q );

        // Q = H_{0} H_{1} ... H_{k-1}
        FLASH_Apply_Q_UT_create_workspace_side( FLA_LEFT, T, Q, &W );
        r_val = FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE,
                                FLA_FORWARD, FLA_COLUMNWISE,
                                A, T, W, Q );
        //FLASH_Obj_free( &W );  Need to see if it's alright to sub these out, this might cause a mem leak
        FLASH_Obj_free_without_buffer( &W );

    }
    else
    {
        // If A and Q are the same objects, Q is formed in-place.
        // - even if A and Q has the same base, they may have different
        //   dimensions.
        // - width of T controls the loop in FLA_QR_UT_form_Q_blk_var1.

        // Zero out the upper triangle of Q.
        FLA_Setr( FLA_UPPER_TRIANGULAR, FLA_ZERO, Q );

        // Adjust T w.r.t A; W is a place holder.
        if ( FLA_Obj_width( T ) > FLA_Obj_width( A ) )
            FLA_Part_1x2( T,    &T,   &W,
                          FLA_Obj_width( A ),
                          FLA_LEFT );

        // Zero out the lower triangle of QBR
        if ( FLA_Obj_width( Q ) > FLA_Obj_width( T ) )
        {
            b = FLA_Obj_width( T );
            FLA_Part_2x2( Q, &QTL, &QTR,
                             &QBL, &QBR, b, b, FLA_TL );
            FLA_Setr( FLA_LOWER_TRIANGULAR, FLA_ZERO, QBR );
        }

        // Set the digaonal to one.
        FLA_Set_diag( FLA_ONE, Q );

        // Create workspace for applying the block Householder transforms.
        FLASH_Apply_Q_UT_create_workspace_side( FLA_LEFT, T, Q, &W );

        // Overwrite Q, which currently contains Householder vectors in the
        // strictly lower triangle and identity in the upper triangle, with
        // the unitary matrix associated with those Householder transforms.
        r_val = FLA_QR_UT_form_Q_blk_var1( Q, T, W );

        // Free the temporary workspace.
        FLASH_Obj_free( &W );
    }

    return r_val;
}
