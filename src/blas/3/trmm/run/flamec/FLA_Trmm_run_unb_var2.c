/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trmm_run_unb_var2( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Scal_external( alpha, B );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_BR );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

  while ( FLA_Obj_length( ABR ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  &a01,     /**/ &A02,
                                                &a10t, &alpha11, /**/ &a12t,
                        /* ************* */   /* ************************** */
                           ABL, /**/ ABR,       &A20,  &a21,     /**/ &A22,
                           1, 1, FLA_TL );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &b1, /**/ &B2,
                           1, FLA_LEFT );

    /*------------------------------------------------------------*/

    /* B2 = B2 + b1 * a12t; */
    FLA_Ger_external( FLA_ONE, b1, a12t, B2 );

    /* b1 = b1 * alpha11; */
    if ( diagA != FLA_UNIT_DIAG )
      FLA_Scal_external( alpha11, b1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  /**/ a01,     A02,
                            /* ************** */  /* ************************ */
                                                     a10t, /**/ alpha11, a12t,
                              &ABL, /**/ &ABR,       A20,  /**/ a21,     A22,
                              FLA_BR );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ b1, B2,
                              FLA_RIGHT );

  }

  return FLA_SUCCESS;
}

#endif
