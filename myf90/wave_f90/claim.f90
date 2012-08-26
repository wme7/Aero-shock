      subroutine claim
      use mod_wave
      implicit none
       allocate(c(nx))   
      allocate(c11(nx-1))
      allocate(c12(nx-1))
      allocate(c22(nx-1))
      allocate(alpha(kt,0:kt-1))
       allocate(x(ks))
       allocate(xt(kt))
       allocate(w(ks))
       allocate(wt(kt))
       allocate(pl(ks,k+1))
       allocate(plt(kt,kt))
       allocate(u(nx,k+1,0:kt))
       allocate(v(nx,k+1,0:kt))
       allocate(fluxu(0:nx))
       allocate(fluxv(0:nx))
       allocate(f(kt))
       allocate(AA(kt))
      return
  
      end subroutine claim
