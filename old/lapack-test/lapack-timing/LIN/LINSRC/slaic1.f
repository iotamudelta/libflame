      SUBROUTINE SLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
*
*  -- LAPACK auxiliary routine (instrumented to count ops, version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            J, JOB
      REAL               C, GAMMA, S, SEST, SESTPR
*     ..
*     .. Array Arguments ..
      REAL               W( J ), X( J )
*     ..
*     .. Common block to return operation count ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      REAL               ITCNT, OPS
*     ..
*
*  Purpose
*  =======
*
*  SLAIC1 applies one step of incremental condition estimation in
*  its simplest version:
*
*  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j
*  lower triangular matrix L, such that
*           twonorm(L*x) = sest
*  Then SLAIC1 computes sestpr, s, c such that
*  the vector
*                  [ s*x ]
*           xhat = [  c  ]
*  is an approximate singular vector of
*                  [ L     0  ]
*           Lhat = [ w' gamma ]
*  in the sense that
*           twonorm(Lhat*xhat) = sestpr.
*
*  Depending on JOB, an estimate for the largest or smallest singular
*  value is computed.
*
*  Note that [s c]' and sestpr**2 is an eigenpair of the system
*
*      diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]
*                                            [ gamma ]
*
*  where  alpha =  x'*w.
*
*  Arguments
*  =========
*
*  JOB     (input) INTEGER
*          = 1: an estimate for the largest singular value is computed.
*          = 2: an estimate for the smallest singular value is computed.
*
*  J       (input) INTEGER
*          Length of X and W
*
*  X       (input) REAL array, dimension (J)
*          The j-vector x.
*
*  SEST    (input) REAL
*          Estimated singular value of j by j matrix L
*
*  W       (input) REAL array, dimension (J)
*          The j-vector w.
*
*  GAMMA   (input) REAL
*          The diagonal element gamma.
*
*  SESTPR  (output) REAL
*          Estimated singular value of (j+1) by (j+1) matrix Lhat.
*
*  S       (output) REAL
*          Sine needed in forming xhat.
*
*  C       (output) REAL
*          Cosine needed in forming xhat.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
      REAL               HALF, FOUR
      PARAMETER          ( HALF = 0.5E0, FOUR = 4.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               ABSALP, ABSEST, ABSGAM, ALPHA, B, COSINE, EPS,
     $                   NORMA, S1, S2, SINE, T, TEST, TMP, ZETA1, ZETA2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. External Functions ..
      REAL               SDOT, SLAMCH
      EXTERNAL           SDOT, SLAMCH
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      ALPHA = SDOT( J, X, 1, W, 1 )
*
      ABSALP = ABS( ALPHA )
      ABSGAM = ABS( GAMMA )
      ABSEST = ABS( SEST )
*
      IF( JOB.EQ.1 ) THEN
*
*        Estimating largest singular value
*
*        special cases
*
         IF( SEST.EQ.ZERO ) THEN
            S1 = MAX( ABSGAM, ABSALP )
            IF( S1.EQ.ZERO ) THEN
               S = ZERO
               C = ONE
               SESTPR = ZERO
            ELSE
               OPS = OPS + 9
               S = ALPHA / S1
               C = GAMMA / S1
               TMP = SQRT( S*S+C*C )
               S = S / TMP
               C = C / TMP
               SESTPR = S1*TMP
            END IF
            RETURN
         ELSE IF( ABSGAM.LE.EPS*ABSEST ) THEN
            OPS = OPS + 7
            S = ONE
            C = ZERO
            TMP = MAX( ABSEST, ABSALP )
            S1 = ABSEST / TMP
            S2 = ABSALP / TMP
            SESTPR = TMP*SQRT( S1*S1+S2*S2 )
            RETURN
         ELSE IF( ABSALP.LE.EPS*ABSEST ) THEN
            S1 = ABSGAM
            S2 = ABSEST
            IF( S1.LE.S2 ) THEN
               S = ONE
               C = ZERO
               SESTPR = S2
            ELSE
               S = ZERO
               C = ONE
               SESTPR = S1
            END IF
            RETURN
         ELSE IF( ABSEST.LE.EPS*ABSALP .OR. ABSEST.LE.EPS*ABSGAM ) THEN
            S1 = ABSGAM
            S2 = ABSALP
            IF( S1.LE.S2 ) THEN
               OPS = OPS + 8
               TMP = S1 / S2
               S = SQRT( ONE+TMP*TMP )
               SESTPR = S2*S
               C = ( GAMMA / S2 ) / S
               S = SIGN( ONE, ALPHA ) / S
            ELSE
               OPS = OPS + 8
               TMP = S2 / S1
               C = SQRT( ONE+TMP*TMP )
               SESTPR = S1*C
               S = ( ALPHA / S1 ) / C
               C = SIGN( ONE, GAMMA ) / C
            END IF
            RETURN
         ELSE
*
*           normal case
*
            OPS = OPS + 8
            ZETA1 = ALPHA / ABSEST
            ZETA2 = GAMMA / ABSEST
*
            B = ( ONE-ZETA1*ZETA1-ZETA2*ZETA2 )*HALF
            C = ZETA1*ZETA1
            IF( B.GT.ZERO ) THEN
               OPS = OPS + 5
               T = C / ( B+SQRT( B*B+C ) )
            ELSE
               OPS = OPS + 4
               T = SQRT( B*B+C ) - B
            END IF
*
            OPS = OPS + 12
            SINE = -ZETA1 / T
            COSINE = -ZETA2 / ( ONE+T )
            TMP = SQRT( SINE*SINE+COSINE*COSINE )
            S = SINE / TMP
            C = COSINE / TMP
            SESTPR = SQRT( T+ONE )*ABSEST
            RETURN
         END IF
*
      ELSE IF( JOB.EQ.2 ) THEN
*
*        Estimating smallest singular value
*
*        special cases
*
         IF( SEST.EQ.ZERO ) THEN
            SESTPR = ZERO
            IF( MAX( ABSGAM, ABSALP ).EQ.ZERO ) THEN
               SINE = ONE
               COSINE = ZERO
            ELSE
               SINE = -GAMMA
               COSINE = ALPHA
            END IF
            OPS = OPS + 7
            S1 = MAX( ABS( SINE ), ABS( COSINE ) )
            S = SINE / S1
            C = COSINE / S1
            TMP = SQRT( S*S+C*C )
            S = S / TMP
            C = C / TMP
            RETURN
         ELSE IF( ABSGAM.LE.EPS*ABSEST ) THEN
            S = ZERO
            C = ONE
            SESTPR = ABSGAM
            RETURN
         ELSE IF( ABSALP.LE.EPS*ABSEST ) THEN
            S1 = ABSGAM
            S2 = ABSEST
            IF( S1.LE.S2 ) THEN
               S = ZERO
               C = ONE
               SESTPR = S1
            ELSE
               S = ONE
               C = ZERO
               SESTPR = S2
            END IF
            RETURN
         ELSE IF( ABSEST.LE.EPS*ABSALP .OR. ABSEST.LE.EPS*ABSGAM ) THEN
            S1 = ABSGAM
            S2 = ABSALP
            IF( S1.LE.S2 ) THEN
               OPS = OPS + 9
               TMP = S1 / S2
               C = SQRT( ONE+TMP*TMP )
               SESTPR = ABSEST*( TMP / C )
               S = -( GAMMA / S2 ) / C
               C = SIGN( ONE, ALPHA ) / C
            ELSE
               OPS = OPS + 8
               TMP = S2 / S1
               S = SQRT( ONE+TMP*TMP )
               SESTPR = ABSEST / S
               C = ( ALPHA / S1 ) / S
               S = -SIGN( ONE, GAMMA ) / S
            END IF
            RETURN
         ELSE
*
*           normal case
*
            OPS = OPS + 14
            ZETA1 = ALPHA / ABSEST
            ZETA2 = GAMMA / ABSEST
*
            NORMA = MAX( ONE+ZETA1*ZETA1+ABS( ZETA1*ZETA2 ),
     $              ABS( ZETA1*ZETA2 )+ZETA2*ZETA2 )
*
*           See if root is closer to zero or to ONE
*
            TEST = ONE + TWO*( ZETA1-ZETA2 )*( ZETA1+ZETA2 )
            IF( TEST.GE.ZERO ) THEN
*
*              root is close to zero, compute directly
*
               OPS = OPS + 20
               B = ( ZETA1*ZETA1+ZETA2*ZETA2+ONE )*HALF
               C = ZETA2*ZETA2
               T = C / ( B+SQRT( ABS( B*B-C ) ) )
               SINE = ZETA1 / ( ONE-T )
               COSINE = -ZETA2 / T
               SESTPR = SQRT( T+FOUR*EPS*EPS*NORMA )*ABSEST
            ELSE
*
*              root is closer to ONE, shift by that amount
*
               OPS = OPS + 6
               B = ( ZETA2*ZETA2+ZETA1*ZETA1-ONE )*HALF
               C = ZETA1*ZETA1
               IF( B.GE.ZERO ) THEN
                  OPS = OPS + 5
                  T = -C / ( B+SQRT( B*B+C ) )
               ELSE
                  OPS = OPS + 4
                  T = B - SQRT( B*B+C )
               END IF
                  OPS = OPS + 10
               SINE = -ZETA1 / T
               COSINE = -ZETA2 / ( ONE+T )
               SESTPR = SQRT( ONE+T+FOUR*EPS*EPS*NORMA )*ABSEST
            END IF
            OPS = OPS + 6
            TMP = SQRT( SINE*SINE+COSINE*COSINE )
            S = SINE / TMP
            C = COSINE / TMP
            RETURN
*
         END IF
      END IF
      RETURN
*
*     End of SLAIC1
*
      END
