C $Header$
C**********************     UPDLDL_ADD       ***********************

      subroutine uldlad(ldl, nrdldl, n, x, z, nrdz, nz, y)

c     Esegue una modifica di rango unitario del problema dei minimi quadrati 
c     in ricorsione all'aggiunta di un vettore, ovvero aggiorna la
c     fattorizzazione LDL della matrice normale B = B + xx' e la matrice
c     dei termini noti Z = Z + xy'. 
c
c     Parametri:    ldl = matrice contenente L e D
c                nrdldl = numero di righe del dimension di LDL
c                     n = ordine di LDL
c                     x = vettore dei coefficienti da 'aggiungere'
c                     z = matrice dei termini noti
c                  nrdz = numero di righe del dimension di Z
c                    nz = numero di termini noti
c                     y = vettore dei termini noti da 'aggiungere'

      implicit none
      integer nrdldl, n, nrdz, nz
      double precision ldl(nrdldl,n), x(n), z(nrdz,nz), y(nz)
      
      double precision t_old, t_new, temp, beta
      integer i, j
      
      do i = 1, n
          do j = 1, nz
              z(i,j) = z(i,j) + x(i)*y(j)
          enddo
      enddo
      
      t_old = 1.D0
      do i = 1, n - 1
          temp = x(i)
          t_new = t_old + temp*temp*ldl(i,i)
          beta = ldl(i,i)*temp/t_new
          ldl(i,i) = ldl(i,i)*t_old/t_new
          t_old = t_new
          do j = i + 1, n
              x(j) = x(j) - temp*ldl(j,i)
              ldl(j,i) = ldl(j,i) + beta*x(j)
          enddo
      enddo
      ldl(n,n) = ldl(n,n)*t_old/(t_old + x(n)*x(n)*ldl(n,n))
      
      return
      end

C**********************       LDL_SOLVE       ***********************

      subroutine ldlsol(ldl, nrdldl, b, nrdb, n, nvet)

c     Solutore di un sistema fattorizzato LDL', LDL'x = b
c
c     Parametri:    ldl = matqrice contenente L e D
c                nrdldl = numero di righe del dimension di LDL
c                     b = matrice dei termini noti
c                  nrdb = numero di righe del dimension di B
c                     n = ordine di LDL, nonche' dei termini noti
c                  nvet = numero di termini noti (numero di colonne di B)
c
c     Restituisce: la soluzione sovrascritta in B. 
      
      implicit none
      integer nrdldl, nrdb, n, nvet
      double precision ldl(nrdldl,n), b(nrdb, nvet)
      
      double precision s
      integer i, j, k


c     risolve il sistema Lv = b con v = DL'x
      
      do k = 1, nvet
          do i = 2, n
              s = b(i,k)
              do j = 1, i - 1
                  s = s - ldl(i,j)*b(j,k)
              enddo
              b(i,k) = s
          enddo
      enddo

c      risolve il sistema L'x = inv(D)v
      
      do k = 1, nvet
          b(n,k) = b(n,k)*ldl(n,n)
          do i = n - 1, 1, -1
              s = b(i,k)*ldl(i,i)
              do j = i + 1, n
                  s = s - ldl(j,i)*b(j,k)
              enddo
              b(i,k) = s
          enddo
      enddo
      return
      end
