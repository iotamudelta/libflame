/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Absolute_square( FLA_Obj alpha )
{
  FLA_Datatype datatype;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Absolute_square_check( alpha );

  datatype = FLA_Obj_datatype( alpha );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    *buff_alpha = (*buff_alpha) * (*buff_alpha);

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    *buff_alpha = (*buff_alpha) * (*buff_alpha);

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );

    buff_alpha->real = buff_alpha->real * buff_alpha->real + 
                       buff_alpha->imag * buff_alpha->imag; 
    buff_alpha->imag = 0.0F; 

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    buff_alpha->real = buff_alpha->real * buff_alpha->real + 
                       buff_alpha->imag * buff_alpha->imag; 
    buff_alpha->imag = 0.0; 

    break;
  }

  }

  return FLA_SUCCESS;
}

