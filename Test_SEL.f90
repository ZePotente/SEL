MODULE Test_SEL
    USE SEL_MET
    USE VYM_IO
    USE VYM_MANIP
    USE VYM_CALCULOS
    IMPLICIT NONE
CONTAINS
    SUBROUTINE PRUEBA(A, BMAT)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, BMAT
        REAL(8), DIMENSION(:), ALLOCATABLE :: B
        INTEGER :: N
        
        N = SIZE(A,1)
        !Procedimientos basicos
!        CALL SYSTEM("clear")
!        CALL PRUEBA_PROC_BASICOS(A, BMAT)
!        READ(*,*)

!        CALL SYSTEM("clear")        
!        CALL PRUEBA_NORMAS(A, B)
!        READ(*,*)

!        !Metodos directos
!        CALL SYSTEM("clear")
!        CALL PRUEBA_MET_GYGJ(A, BMAT)
!        READ(*,*)
        
        CALL SYSTEM("clear")
        CALL PRUEBA_THOMAS()
        READ(*,*)
        
!        ALLOCATE(B(SIZE(BMAT,1)))
!        B = BMAT(:,1)
!        CALL SYSTEM("clear")
!        CALL PRUEBA_CROUT(A, B)
!        READ(*,*)
!        CALL SYSTEM("clear")
!        CALL PRUEBA_RESIDUO(A, B)
!        READ(*,*)

!        CALL SYSTEM("clear")
!        CALL PRUEBA_REFINAMIENTO(A, B)
!        READ(*,*)
!        !Sensibilidad
!        CALL SYSTEM("clear")
!        CALL PRUEBA_ERR_REL()
!        READ(*,*)
        
!        CALL SYSTEM("clear")
!        CALL PRUEBA_COTA()
!        READ(*,*)
!        !Metodos Indirectos
!        CALL SYSTEM("clear")
!        CALL PRUEBA_MET_IND(A, B)
!        READ(*,*)
        IF (ALLOCATED(B)) DEALLOCATE(B)
    END SUBROUTINE
    
    SUBROUTINE PRUEBA_PROC_BASICOS(A, BMAT)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, BMAT
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: AUX
        
        PRINT *, 'PROCESOS BASICOS'
        CALL MOSTRARMATS(A, BMAT)
        !
        PRINT *, 'Matriz identidad'
        IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
        CALL MAT_IDENTIDAD(4, AUX)
        CALL MAT_MOSTRAR(AUX)
        !
        PRINT *, 'Matriz ampliada'
        CALL MAT_MOSTRAR(MATRIZAMPLIADA(A, BMAT))
        IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
    END SUBROUTINE
    
    SUBROUTINE PRUEBA_MET_GYGJ(A, BMAT)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, BMAT
        !
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: AUX, AUX2
        
        PRINT *, 'GAUSS Y GJ'
        CALL MOSTRARMATS(A, BMAT)
        
        !
        PRINT *, 'Gauss'
        IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
        CALL MET_GAUSS(A, BMAT, AUX)
        CALL MAT_MOSTRAR(AUX)
        PRINT *, 'Gauss con sustitucion regresiva'
        CALL SUST_REGRESIVA(AUX, AUX2)
        CALL MAT_MOSTRAR(AUX2)
        DEALLOCATE(AUX2)
        !
        PRINT *, 'Gauss-Jordan'
        IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
        CALL MET_GAUSSJORDAN(A, BMAT, AUX)
        CALL MAT_MOSTRAR(AUX)
        !
        PRINT *, 'Matriz inversa'
        IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
        CALL MAT_INVERSA(A, AUX)
        CALL MAT_MOSTRAR(AUX)
        PRINT *, 'Comprobacion de Matriz inversa'
        CALL MAT_MOSTRAR(MATMUL(A,AUX)) !Ver si A * A' = I
        IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
    END SUBROUTINE
    
    SUBROUTINE MOSTRARMATS(A, B)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A, B
        
        PRINT *, 'A'
        CALL MAT_MOSTRAR(A)
        PRINT *, 'B'
        CALL MAT_MOSTRAR(B)
    END SUBROUTINE
    
    SUBROUTINE PRUEBA_NORMAS(A, B)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8), DIMENSION(:), INTENT(IN) :: B
        !
        PRINT *, 'Normas vectoriales y matriciales'
        !
        CALL MOSTRARMYV(A, B)
        !
        PRINT *, 'Normas vectoriales'
        PRINT *, 'Norma M'
        PRINT *, VEC_NORMAM(B)
        PRINT *, 'Norma E'
        PRINT *, VEC_NORMAE(B)
        !
        PRINT *, 'Normas matriciales'
        PRINT *, 'Norma M'
        PRINT *, MAT_NORMAM(A)
        PRINT *, 'Norma L'
        PRINT *, MAT_NORMAL(A)
        PRINT *, 'Norma Frobenius'
        PRINT *, MAT_NORMAFROBENIUS(A)
    END SUBROUTINE
    
    
    SUBROUTINE PRUEBA_THOMAS()
        !Local
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: A, B
        REAL(8), DIMENSION(:), ALLOCATABLE :: L, D, U
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: AUX, AUX2
        INTEGER :: BANDERA
        !
        
        PRINT *, 'Metodos directos restantes'
        !Lectura de thomas matricial
        CALL LECTURA_THOMAS(A, B, U, D, L)
        
        PRINT *, 'Metodo de Thomas con matriz'
        CALL MET_THOMAS(A, B, AUX)
        CALL MAT_MOSTRAR(AUX)
        CALL MAT_GUARDAR(AUX, BANDERA, '3)Thomas.txt')
        PRINT *, 'Comparo con GJ'
        CALL MET_GAUSSJORDAN(A, B, AUX2)
        CALL MAT_MOSTRAR(AUX2(:,SIZE(A,1)+1:))
        CALL MAT_GUARDAR(AUX, BANDERA, '3)GJ.txt')
        DEALLOCATE(AUX, AUX2)
        
        PRINT *, 'Metodo de Thomas con vectores'
        CALL MET_THOMAS_VEC(U, D, L, B, AUX)
        CALL MAT_MOSTRAR(AUX)
        PRINT *, 'Comparo con GJ'
        CALL MET_GAUSSJORDAN(A, B, AUX2)
        PRINT*, 'ASD'
        CALL MAT_MOSTRAR(AUX2(:,SIZE(A,1)+1:))
        PRINT*, 'ASD'
        DEALLOCATE(AUX, AUX2)
        
        DEALLOCATE(A, B, U, D, L)
    END SUBROUTINE
    
    SUBROUTINE LECTURA_THOMAS(A_TH, B_TH, U, D, L)
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: A_TH, B_TH
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: U, D, L
        INTEGER :: BANDERA 
        PRINT *, 'Metodo de Thomas matricial'
        CALL MAT_LEER(A_TH, BANDERA, 'Prueba_A_Thomas_mat.txt')
        PRINT *, 'A Thomas'
!        CALL MAT_MOSTRAR(A_TH)
        
        CALL VEC_LEER(U, BANDERA, 'Prueba_U_Thomas.txt')
        CALL VEC_LEER(D, BANDERA, 'Prueba_D_Thomas.txt')
        CALL VEC_LEER(L, BANDERA, 'Prueba_L_Thomas.txt')
!        PRINT *, 'U, D, L de Thomas'
!        CALL VEC_MOSTRAR(U)
!        CALL VEC_MOSTRAR(D)
!        CALL VEC_MOSTRAR(L)
        
        CALL MAT_LEER(B_TH, BANDERA, 'Prueba_B_Thomas.txt')
!        PRINT *, 'B Thomas'
!        CALL MAT_MOSTRAR(B_TH)
    END SUBROUTINE
    
    SUBROUTINE PRUEBA_CROUT(A, B)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        !
        REAL(8), ALLOCATABLE :: CROUT(:), GJ(:,:), AUXB(:,:)
        
        PRINT *, 'Metodo de Crout'
        CALL MET_LU_CROUT(A, B, CROUT)
        CALL VEC_MOSTRAR(CROUT)
        PRINT *, 'Comparo con GJ'
        ALLOCATE(AUXB(SIZE(B),1)); AUXB(:,1) = B
        CALL MET_GAUSSJORDAN(A, AUXB, GJ)
        CALL VEC_MOSTRAR(GJ(:,SIZE(A,1)+1))
        DEALLOCATE(CROUT, GJ, AUXB)
    END SUBROUTINE
    
    SUBROUTINE PRUEBA_RESIDUO(A, B)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        !
        REAL(8), DIMENSION(:), ALLOCATABLE :: R
        REAL(8), ALLOCATABLE :: AUX(:)
        
        PRINT *, 'Prueba de Residuo'
        CALL MOSTRARMYV(A, B)
        
        CALL MET_LU_CROUT(A, B, AUX)
        R = VEC_RESIDUO(A, B, AUX)
        PRINT *, 'Residuo'
        CALL VEC_MOSTRAR(R)
    END SUBROUTINE
    
    SUBROUTINE PRUEBA_REFINAMIENTO(A, B)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        !
        REAL(8), ALLOCATABLE :: XORIG(:), XREF(:), XPERT(:)
        REAL(8) :: TOL
        PRINT *, 'Prueba de refinamiento iterativo'
        CALL MOSTRARMYV(A, B)
        
        CALL MET_LU_CROUT(A, B, XORIG)
        PRINT *, 'Solucion por Crout'
        CALL VEC_MOSTRAR(XORIG)
        
        PRINT *, 'Solucion para probar refinamiento y comparar'
        ALLOCATE(XPERT(SIZE(XORIG)))
!        XPERT(2) = 2.1; XPERT(3) = -1.21; XPERT(4) = 2.99;
        XPERT = 0.
!        XPERT(1) = 1.4;
        CALL VEC_MOSTRAR(XPERT)
        
        TOL = 1E-14
        CALL REF_ITER(A, B, XPERT, TOL, XREF)
        PRINT *, 'Resultado del refinamiento'
        CALL VEC_MOSTRAR(XREF)
        
        DEALLOCATE(XORIG, XREF, XPERT)
    END SUBROUTINE
    !----!
    SUBROUTINE PRUEBA_ERR_REL()
        REAL(8), DIMENSION(3,3) :: A, APERT
        REAL(8), DIMENSION(3) :: B
        REAL(8), ALLOCATABLE :: X(:), XREF(:), XPERT(:), XMAT(:,:)
        REAL(8), DIMENSION(3,1) :: BMAT
        REAL(8) :: TOL
        A(1,1) = 3.02; A(1,2) = -1.05; A(1,3) = 2.53;
        A(2,1) = 4.33; A(2,2) = 0.56; A(2,3) = -1.78;
        A(3,1) = -0.83; A(3,2) = -0.54; A(3,3) = 1.47;
        B(1) = -1.61; B(2) = 7.23; B(3) = -3.38;
        PRINT *, 'Prueba de Error relativo'
        CALL MOSTRARMYV(A, B)
        CALL MET_LU_CROUT(A, B, X)
        PRINT *, 'X'
        CALL VEC_MOSTRAR(X)
        BMAT(:,1) = B
        CALL MET_GAUSSJORDAN(A, BMAT, XMAT)
        TOL = 1E-8
        CALL REF_ITER(A, B, X, TOL, XREF)
        APERT(1,1) = 3.00; APERT(1,2) = -1.05; APERT(1,3) = 2.53;
        APERT(2,1) = 4.33; APERT(2,2) = 0.56; APERT(2,3) = -1.78;
        APERT(3,1) = -0.83; APERT(3,2) = -0.54; APERT(3,3) = 1.47;
        CALL MET_LU_CROUT(APERT, B, XPERT)
        
        PRINT *, 'APERT'
        CALL MAT_MOSTRAR(APERT)
        PRINT *, 'X GJ'
        CALL MAT_MOSTRAR(XMAT)
        PRINT *, 'XREF'
        CALL VEC_MOSTRAR(XREF)
        PRINT *, 'XPERT'
        CALL VEC_MOSTRAR(XPERT)
        
        WRITE(*,*) 'Error relativo en A', SENS_ERROR_REL_MAT(A, APERT)
        WRITE(*,*) 'Error relativo en X', SENS_ERROR_REL_VEC(X, XPERT)
    END SUBROUTINE
    
    SUBROUTINE PRUEBA_COTA()
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: A
        REAL(8), DIMENSION(:), ALLOCATABLE :: B, BPERT, X, XPERT
        INTEGER :: BANDERA, N
        N = SIZE(B)
        CALL MAT_LEER(A, BANDERA, '5)A.txt')
        CALL VEC_LEER(B, BANDERA, '5)B.txt')
        PRINT *, 'Prueba de sensibilidad'
        CALL MOSTRARMYV(A, B)
        PRINT *, 'APERT'
        CALL MAT_MOSTRAR(A)
        PRINT *, 'BPERT'
        CALL VEC_MOSTRAR(B)
        
        WRITE(*,'(A F20.10)') 'Nro de condicion: ', SENS_NROCOND(A)
        ALLOCATE(BPERT(N))
        BPERT(1) = 30.; BPERT(2) = 42.; BPERT(3) = 22.;
        
        CALL MET_LU_CROUT(A, B, X)
        CALL MET_LU_CROUT(A, BPERT, XPERT)
        PRINT *, 'X'
        CALL VEC_MOSTRAR(X)
        PRINT *, 'XPERT'
        CALL VEC_MOSTRAR(XPERT)
        
        WRITE(*,'(A F20.10)')'Error relativo de X: ', SENS_ERROR_REL_VEC(X, XPERT)
        WRITE(*,'(A F20.10)')'Cota teorica: ', SENS_COTAPERT_SOL(A, B = B, BPERT = BPERT)
        
        DEALLOCATE(BPERT, A, B)
    END SUBROUTINE
    !----!
    SUBROUTINE PRUEBA_MET_IND(AIN, BIN)
        REAL(8), INTENT(IN) :: AIN(:,:), BIN(:)
        !
        REAL(8), ALLOCATABLE :: AUX(:), AUX2(:,:), AUXB(:,:)
        REAL(8) :: A(3,3), B(3)
        REAL(8) :: TOLERANCIA, W
        TOLERANCIA = 0.001; W = 0.5
        A(1,1) = 3.02; A(1,2) = -1.05; A(1,3) = 2.53;
        A(2,1) = 4.33; A(2,2) = 0.56; A(2,3) = -1.78;
        A(3,1) = -0.83; A(3,2) = -0.54; A(3,3) = 1.47;
        B(1) = -1.61; B(2) = 7.23; B(3) = -3.38;
        !
        PRINT *, 'Metodo de Jacobi'
        CALL MET_JACOBI(A, B, TOLERANCIA, AUX)
        CALL VEC_MOSTRAR(AUX)
        PRINT *, 'Comparo con GJ'
        ALLOCATE(AUXB(SIZE(B),1)); AUXB(:,1) = B
        CALL MET_GAUSSJORDAN(A, AUXB, AUX2)
        CALL VEC_MOSTRAR(AUX2(:,SIZE(A,1)+1))
        DEALLOCATE(AUX, AUX2, AUXB)
        
        PRINT *, 'Metodo de Gauss-Seidel'
        CALL MET_GS(A, B, TOLERANCIA, AUX)
        CALL VEC_MOSTRAR(AUX)
        PRINT *, 'Comparo con GJ'
        ALLOCATE(AUXB(SIZE(B),1)); AUXB(:,1) = B
        CALL MET_GAUSSJORDAN(A, AUXB, AUX2)
        CALL VEC_MOSTRAR(AUX2(:,SIZE(A,1)+1))
        DEALLOCATE(AUX, AUX2, AUXB)
        
!        PRINT *, 'Metodo de SOR'
!        CALL MET_SOR(A, B, TOLERANCIA, W, AUX)
!        CALL VEC_MOSTRAR(AUX)
!        PRINT *, 'Comparo con GJ'
!        ALLOCATE(AUXB(SIZE(B),1)); AUXB(:,1) = B
!        CALL MET_GAUSSJORDAN(A, AUXB, AUX2)
!        CALL VEC_MOSTRAR(AUX2(:,SIZE(A,1)+1))
!        DEALLOCATE(AUX, AUX2, AUXB)
    END SUBROUTINE
    
    SUBROUTINE MOSTRARMYV(A, B)
        REAL(8), INTENT(IN) :: A(:,:), B(:)
        
        PRINT *, 'A'
        CALL MAT_MOSTRAR(A)
        PRINT *, 'B'
        CALL VEC_MOSTRAR(B) !No uso mostrarmats porque es vector
    END SUBROUTINE
END MODULE
