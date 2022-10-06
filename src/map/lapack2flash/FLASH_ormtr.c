/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


//  This function calls ormql (upper) and ormqr (lower). As FLA_hetrd
//  is only implemented for the lower triangular, we do not need to
//  directly interface this.
#ifdef FLA_ENABLE_LAPACK2FLASH

#include "FLASH_lapack2flash_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  ORMTR overwrites the general real M-by-N matrix C with
                 SIDE = 'L'     SIDE = 'R'
  TRANS = 'N':      Q * C          C * Q
  TRANS = 'T':      Q**T * C       C * Q**T

  where Q is a real orthogonal matrix of order nq, with nq = m if
  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
  nq-1 elementary reflectors, as returned by SYTRD:

  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).

  TODO:: Complete FLA_Accum_T_UT for FLA_BACKWARD.
*/

#define LAPACK_ormtr(prefix, name)                                                      \
  int F77_ ## prefix ## name ## tr( char* side,                                         \
                                    char* uplo,                                         \
                                    char* trans,                                        \
                                    int*  m,                                            \
                                    int*  n,                                            \
                                    PREFIX2LAPACK_TYPEDEF(prefix) *buff_A, int* ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix) *buff_t,              \
                                    PREFIX2LAPACK_TYPEDEF(prefix) *buff_C, int* ldim_C, \
                                    PREFIX2LAPACK_TYPEDEF(prefix) *buff_w, int* lwork,  \
                                    int* info )

#define LAPACK_ormtr_body(prefix)                                       \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                \
  FLA_Side     side_fla;                                                \
  FLA_Uplo     uplo_fla;                                                \
  FLA_Trans    trans_fla;                                               \
  dim_t        m_d, m_e;                                                \
  FLA_Obj      A, C;                                                    \
  FLA_Error    init_result;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Param_map_netlib_to_flame_side( side, &side_fla );                \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );                \
  FLA_Param_map_netlib_to_flame_trans( trans, &trans_fla );             \
                                                                        \
  if      ( side_fla == FLA_LEFT )     m_d = *m;                        \
  else /* ( side_fla == FLA_RIGHT) */  m_d = *n;                        \
  m_e = ( m_d - 1 );                                                    \
                                                                        \
  FLASH_Obj_create_without_buffer( datatype,                            \
                                   *m,                                  \
                                   *n,                                  \
                                   FLASH_depth,                         \
                                   FLASH_blocksize,                     \
                                   &C );                                \
  FLASH_Obj_attach_buffer( buff_C, 1, *ldim_C, &C );                    \
                                                                        \
  FLASH_Obj_create_without_buffer( datatype,                            \
                                   m_d,                                 \
                                   m_d,                                 \
                                   FLASH_depth,                         \
                                   FLASH_blocksize,                     \
                                   &A );                                \
  FLASH_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                    \
                                                                        \
  if ( m_e > 0 ) {                                                      \
    FLA_Obj    ATL, ATR, ABL, ABR, T, W, t;                             \
    FLA_Direct direct;                                                  \
                                                                        \
    FLASH_Obj_create_without_buffer( datatype,                          \
                                     m_e,                               \
                                     1,                                 \
                                     FLASH_depth,                       \
                                     FLASH_blocksize,                   \
                                     &t );                              \
    FLASH_Obj_attach_buffer( buff_t, 1, m_e, &t );                      \
    PREFIX2FLAME_INVERT_TAU(prefix,t);                                  \
                                                                        \
    if ( uplo_fla == FLA_LOWER_TRIANGULAR ) {                           \
      FLA_Part_2x2( A, &ATL, &ATR,                                      \
                       &A,   &ABR, 1, 1, FLA_TR );                      \
      direct = FLA_FORWARD;                                             \
    } else {                                                            \
      FLA_Part_2x2( A, &ATL, &A,                                        \
                       &ABL, &ABR, 1, 1, FLA_BL );                      \
      direct = FLA_BACKWARD;                                            \
    }                                                                   \
    if ( side_fla == FLA_LEFT ) {                                       \
      FLA_Part_2x1( C, &W, &C, 1, FLA_TOP );                            \
    } else {                                                            \
      FLA_Part_1x2( C, &W, &C, 1, FLA_LEFT );                           \
    }                                                                   \
                                                                        \
    FLASH_QR_UT_create_T( A, &T );                                      \
    FLASH_Set( FLA_ZERO, T );                                           \
    FLASH_Apply_Q_UT_create_workspace_side( side_fla, T, C, &W);        \
    FLA_Accum_T_UT( direct, FLA_COLUMNWISE, A, t, T );                  \
                                                                        \
    if ( FLA_Obj_is_complex( A ) == TRUE ) {                            \
      FLA_Obj d2, e2, r;                                                \
                                                                        \
      /* Temporary vectors to store diagonal and subdiagonal */         \
      FLASH_Obj_create( datatype,                                       \
                        m_d,                                            \
                        1,                                              \
                        FLASH_depth,                                    \
                        FLASH_blocksize,                                \
                        &d2 );                                          \
      FLASH_Obj_create( datatype,                                       \
                        m_e,                                            \
                        1,                                              \
                        FLASH_depth,                                    \
                        FLASH_blocksize,                                \
                        &e2 );                                          \
                                                                        \
      /* Temporary vectors to store realifying transformation */        \
      FLASH_Obj_create( datatype,                                       \
                        m_d,                                            \
                        1,                                              \
                        FLASH_depth,                                    \
                        FLASH_blocksize,                                \
                        &r );                                           \
                                                                        \
      /* Extract diagonals (complex) and realify them. */               \
      FLA_Tridiag_UT_extract_diagonals( uplo_fla, A, d2, e2 );          \
      FLA_Tridiag_UT_realify_subdiagonal( e2, r );                      \
                                                                        \
      if      ( side_fla  == FLA_LEFT &&                                \
                trans_fla == FLA_NO_TRANSPOSE )                         \
        FLA_Apply_diag_matrix( FLA_LEFT, FLA_CONJUGATE, r, C );         \
      else if ( side_fla  == FLA_RIGHT &&                               \
                trans_fla == FLA_CONJ_TRANSPOSE )                       \
        FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, r, C );     \
                                                                        \
      FLASH_Apply_Q_UT( side_fla, trans_fla, direct, FLA_COLUMNWISE,    \
                      A, T, W, C );                                     \
                                                                        \
      if      ( side_fla  == FLA_LEFT &&                                \
                trans_fla == FLA_CONJ_TRANSPOSE )                       \
        FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, r, C );      \
      else if ( side_fla  == FLA_RIGHT &&                               \
                trans_fla == FLA_NO_TRANSPOSE )                         \
        FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, r, C );        \
                                                                        \
      FLASH_Obj_free( &r  );                                            \
      FLASH_Obj_free( &e2 );                                            \
      FLASH_Obj_free( &d2 );                                            \
    } else {                                                            \
      FLASH_Apply_Q_UT( side_fla, trans_fla, direct, FLA_COLUMNWISE,    \
                      A, T, W, C );                                     \
    }                                                                   \
                                                                        \
    FLASH_Obj_free( &W );                                               \
    FLASH_Obj_free( &T );                                               \
                                                                        \
    PREFIX2FLAME_INVERT_TAU(prefix,t);                                  \
    FLASH_Obj_free_without_buffer( &t );                                \
  }                                                                     \
                                                                        \
  FLASH_Obj_free_without_buffer( &A );                                  \
  FLASH_Obj_free_without_buffer( &C );                                  \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;

extern int sormtr_fla(char *side, char *uplo, char *trans, integer *m, integer *n, real *a,          integer *lda, real *tau,          real *c__,          integer *ldc, real *work,          integer *lwork, integer *info);
extern int dormtr_fla(char *side, char *uplo, char *trans, integer *m, integer *n, doublereal *a,    integer *lda, doublereal *tau,    doublereal *c__,    integer *ldc, doublereal *work,    integer *lwork, integer *info);
extern int cunmtr_fla(char *side, char *uplo, char *trans, integer *m, integer *n, complex *a,       integer *lda, complex *tau,       complex *c__,       integer *ldc, complex *work,       integer *lwork, integer *info);
extern int zunmtr_fla(char *side, char *uplo, char *trans, integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, integer *info);

LAPACK_ormtr(s, orm)
{
    {
        if ( *uplo == 'U' )
        {
            sormtr_fla( side, uplo, trans,
                        m, n,
                        buff_A, ldim_A,
                        buff_t,
                        buff_C, ldim_C,
                        buff_w, lwork,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( sormtr_check( side, uplo, trans,
                                           m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_C, ldim_C,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormtr_body(s)
    }
}
LAPACK_ormtr(d, orm)
{
    {
        if ( *uplo == 'U' )
        {
            dormtr_fla( side, uplo, trans,
                        m, n,
                        buff_A, ldim_A,
                        buff_t,
                        buff_C, ldim_C,
                        buff_w, lwork,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( dormtr_check(  side, uplo, trans,
                                            m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_C, ldim_C,
                                            buff_w, lwork,
                                            info ) )
    }
    {
        LAPACK_ormtr_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_ormtr(c, unm)
{
    {
        if ( *uplo == 'U' )
        {
            cunmtr_fla( side, uplo, trans,
                        m, n,
                        (complex*)buff_A, ldim_A,
                        (complex*)buff_t,
                        (complex*)buff_C, ldim_C,
                        (complex*)buff_w, lwork,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( cunmtr_check(  side, uplo, trans,
                                            m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_C, ldim_C,
                                            buff_w, lwork,
                                            info ) )
    }
    {
        LAPACK_ormtr_body(c)
    }
}
LAPACK_ormtr(z, unm)
{
    {
        if ( *uplo == 'U' )
        {
            zunmtr_fla( side, uplo, trans,
                        m, n,
                        (doublecomplex*)buff_A, ldim_A,
                        (doublecomplex*)buff_t,
                        (doublecomplex*)buff_C, ldim_C,
                        (doublecomplex*)buff_w, lwork,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( zunmtr_check(  side, uplo, trans,
                                            m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_C, ldim_C,
                                            buff_w, lwork,
                                            info ) )
    }
    {
        LAPACK_ormtr_body(z)
    }
}
#endif

#endif
