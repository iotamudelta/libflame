/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_LQ_UT_form_Q( FLA_Obj A, FLA_Obj T, FLA_Obj Q )
{
  FLA_Error r_val = FLA_SUCCESS;

  // Flip a base once.
  FLA_Obj_flip_base( &A );
  if ( FLA_Obj_is( A, Q ) == FALSE )
    FLA_Obj_flip_base( &Q );

  // Dimensions of the both matrices should be flipped.
  FLA_Obj_flip_view( &A );
  FLA_Obj_flip_view( &Q );

  // Run the QR utility function.
  r_val = FLASH_QR_UT_form_Q( A, T, Q );

  // Apply conjugation on Q as we use QR_UT_form_Q.
  if ( FLA_Obj_is_complex( Q ) )
    FLA_Conjugate( Q );
  
  // Recover the base object.
  if ( FLA_Obj_is( A, Q ) == FALSE )
    FLA_Obj_flip_base( &Q );
  FLA_Obj_flip_base( &A );

  return r_val;
}
