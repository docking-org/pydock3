      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
c------------------------------------------------------
      SUBROUTINE CROSS(U,V,W)
      dimension u(3),v(3),w(3)
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
      return
      end
c------------------------------------------------------
      SUBROUTINE DOT(U,V,RDOT)
      dimension u(3),v(3)
      rdot = 0.0
	rdot = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
      return
      end
c------------------------------------------------------
      SUBROUTINE NORM(U,SIZE)
      dimension u(3)
      call dot(u,u,size)
      if(size.ne.0.) then
        size = sqrt(size)
        u(1) = u(1)/size
        u(2) = u(2)/size
        u(3) = u(3)/size
      else
c        type *,'zero vector'
      end if
      return
      end
c-------------------------------------------------------
      SUBROUTINE PERP(U,V)
c	returns vector v normal to u
      dimension u(3),v(3)
      v(1) = u(2) - u(3)
      v(2) = u(3) - u(1)
      v(3) = u(1) - u(2)
	if(abs(v(1))+abs(v(2))+abs(v(3)).le.1.e-5)then
        v(1) = u(1)
        v(2) = -u(2)/2.
        v(3) = -u(3)/2.
	end if
      return
      end      
c-------------------------------------------------------
      function icomb(i,j)
      integer*4 icomb,i,j,imj,n

c     i!/j!(i-j)!
      icomb = 1
      do n = j+1,i
        icomb = icomb*n
      end do
      imj = i-j
      icomb = icomb/ifact(imj)
c     print *,'i,j,icomb: ',i,j,icomb
      return
      end

      function ifact(i)
      integer*4 ifact,i,n
      ifact = 1
      do n = 1,i
        ifact = ifact*n
      end do
c     print *,'i,ifact: ',i,ifact
      return
      end
c-------------------------------------------------------
      subroutine mtdet(rm1,det)
c calc determinant
      dimension rm1(3,3)
      d1 = rm1(1,1)*rm1(2,2)*rm1(3,3)
      d2 = rm1(1,1)*rm1(3,2)*rm1(2,3)
      d3 = rm1(2,1)*rm1(3,2)*rm1(1,3)
      d4 = rm1(2,1)*rm1(1,2)*rm1(3,3)
      d5 = rm1(3,1)*rm1(1,2)*rm1(2,3)
      d6 = rm1(3,1)*rm1(2,2)*rm1(1,3)
      det = d1 - d2 + d3 - d4 + d5 - d6
      return
      end
c--------------------------------------------------------

