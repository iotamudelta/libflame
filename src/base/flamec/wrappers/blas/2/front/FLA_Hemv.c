
#include "FLAME.h"

FLA_Error FLA_Hemv( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
	return FLA_Hemv_external( uplo, alpha, A, x, beta, y );
}
