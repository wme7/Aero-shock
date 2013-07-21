!!!!!!! FLUXES2,  ROE FLUX !!!!!!!!!!!!!!!!
    subroutine FLUX_ROE
    use mainvar
    implicit none
    integer :: i
    real :: ro_L,u_L,p_L,ro_R,u_R,p_R, sr1,sr2,sr,rop,up,hp,a2p,  &
      HL,HR,ap,   &
      rom1,rom2,rom3,rv1(3),rv2(3),rv3(3), &
      d_ro,d_u,d_p,alp1,alp2,alp3, f1L,f2L,f3L,f1R,f2R,f3R,  &
      a_L,rosL,us,Er,ps,asL,rom1L,rom1R,rom3L,rom3R , &
	  z,aL,aR, asR

    do i=0,NN
      ro_L= roL(i)
      u_L = uL(i)
      p_L = pL(i)
      ro_R= roR(i)
      u_R = uR(i)
      p_R = pR(i)

      sr1=sqrt(ro_L)
      sr2=sqrt(ro_R)
      sr =sr1+sr2
      rop=sr1*sr2       !!!!!!!??????????????????????????????????????????, rop 可能不对????????
      up=(sr1*u_L+sr2*u_R)/sr
      HL=0.5*u_L*u_L+gama/(gama-1.)*p_L/ro_L
      HR=0.5*u_R*u_R+gama/(gama-1.)*p_R/ro_R
      Hp=(sr1*hL+sr2*hR)/(sr1+sr2)
      a2p=(gama-1.)*(Hp-0.5*up*up)
      ap=sqrt(a2p)

      rv1 = (/1.,up-ap,Hp-up*ap /)
      rv2 = (/1.,up,  0.5*up*up /)
      rv3 = (/1.,up+ap,Hp+up*ap /)

      d_ro= ro_R - ro_L
      d_u = u_R - u_L
      d_p = p_R - p_L
      alp1=(d_p-rop*ap*d_u)/2./a2p       !! need to deduce in one dimension
      alp2=(a2p*d_ro-d_p)/a2p
      alp3=(d_p+rop*ap*d_u)/2./a2p

	  rom1=(up-ap)  !!! ???? , abs?
      rom2=abs(up)
      rom3=(up+ap)
	  !! entropy fix 3:
	  z =(gama-1.)/(2.*gama)
	  aL=sqrt(gama*p_L/ro_L)
	  aR=sqrt(gama*p_R/ro_R)
	  ps=( (aL+aR-(gama-1.)/2.*(u_R-u_L))/(aL/p_L**z+aR/p_R**z) )**(1./z)
	  asL= aL*(ps/p_L)**z
	  us = u_L+2./(gama-1)*(aL-asL)
	  rom1L=u_L-aL
	  rom1R=us-asL
	  if( abs(rom1R-rom1L)>1.0e-16 ) rom1 =rom1L*(rom1R-rom1)/(rom1R-rom1L) !!
	  asR= aR*(ps/p_R)**z
	  us = u_R+2./(gama-1.)*(asR-aR)
	  rom3L=us+asR
	  rom3R=u_R+aR
	  if( abs(rom3R-rom3L)>1.0e-16 ) rom3 =rom3R*(rom3-rom3L)/(rom3R-rom3L)  !!
	  rom1=abs(rom1)
	  rom3=abs(rom3)

      F1L= ro_L*u_L
      F2L= ro_L*u_L*u_L + p_L
      F3L= u_L*p_L*gama/(gama-1.) + 0.5*ro_L*u_L*u_L*u_L
      F1R= ro_R*u_R
      F2R= ro_R*u_R*u_R + p_R
      F3R= u_R*p_R*gama/(gama-1.) + 0.5*ro_R*u_R*u_R*u_R

      F1(i)=0.5*(f1L+f1R - (alp1*rom1*rv1(1)+alp2*rom2*rv2(1)+alp3*rom3*rv3(1) ))
      F2(i)=0.5*(f2L+f2R - (alp1*rom1*rv1(2)+alp2*rom2*rv2(2)+alp3*rom3*rv3(2) ))
      F3(i)=0.5*(f3L+f3R - (alp1*rom1*rv1(3)+alp2*rom2*rv2(3)+alp3*rom3*rv3(3) ))
 
    enddo 

    end subroutine