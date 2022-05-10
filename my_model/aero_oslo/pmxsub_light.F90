
module pmxsub

  implicit none
  private
  public :: pmxsub_light


contains

  subroutine pmxsub_light(lchnk,pcols,pver,ncol,Nnatk,per_tau, per_tau_w, per_tau_w_g)


    use opttab
    use opttab_lw
    use optinterpol, only: interpol5to10
    use shr_kind_mod, only: r8 => shr_kind_r8
    use commondefinitions
    ! input arguments
    implicit none
    integer, intent(in) :: lchnk 
    integer, intent(in) :: pver
    integer, intent(in) :: pcols
    integer, intent(in) :: ncol


    ! in out arguments

    real(r8), intent(inout) :: Nnatk(pcols,pver,0:nmodes)! aerosol mode number concentration  
    
    ! Output arguments
    real(r8), intent(out) :: per_tau   (pcols,0:pver,nbands) ! aerosol extinction optical depth
    real(r8), intent(out) :: per_tau_w (pcols,0:pver,nbands) ! aerosol single scattering albedo * tau
    !real(r8), intent(out) :: per_lw_abs (pcols,pver,nlwbands) ! aerosol absorption optical depth (LW)
    real(r8), intent(out) :: per_tau_w_g(pcols,0:pver,nbands) ! aerosol assymetry parameter * w * tau



    ! Local variables
    integer  i, k, ib, icol, iloop

    logical  daylight(pcols)
    real(r8) deltah_km(pcols,pver) 
    real(r8) betot(pcols,pver,nbands)
    real(r8) ssatot(pcols,pver,nbands)
    real(r8) asymtot(pcols,pver,nbands) 
    real(r8) faqm(pcols,pver,nbmodes), fbcm(pcols,pver,nbmodes)
    real(r8) xrh(pcols,pver)
    real(r8) fnbc(pcols,pver), faitbc(pcols,pver), f_so4_cond(pcols,pver), &
    f_soa(pcols,pver),f_soana(pcols,pver), vnbc, vaitbc
    real(r8) xfbcbg(pcols,pver), xfbcbgn(pcols,pver), ifbcbgn1(pcols,pver)
    real(r8) xfombg(pcols,pver)
    real(r8) rhum(pcols,pver)       ! (trimmed) relative humidity for the aerosol calculations
    integer  irh1(pcols,pver)
    real(r8) focm(pcols,pver,4)
    real(r8) Cam(pcols,pver,nbmodes), fcm(pcols,pver,nbmodes)
    integer ict1(pcols,pver,nmodes)
    real(r8) xfac(pcols,pver,nbmodes)
    integer ifac1(pcols,pver,nbmodes)
    real(r8) xfbc(pcols,pver,nbmodes)
    integer ifbc1(pcols,pver,nbmodes)
    real(r8) xfaq(pcols,pver,nbmodes)
    integer ifaq1(pcols,pver,nbmodes)
    integer ifombg1(pcols,pver), ifbcbg1(pcols,pver)

    real(r8) xct(pcols,pver,nmodes)
    real(r8) ke(pcols,pver,0:nmodes,nbands)
    real(r8) kalw(pcols,pver,0:nmodes,nlwbands)
    real(r8) ssa(pcols,pver,0:nmodes,nbands), asym(pcols,pver,0:nmodes,nbands), & 
            be(pcols,pver,0:nmodes,nbands)

    logical, parameter :: lw_on=.false.  
    call initopt

    ! Set local variables that are not used to zero
    do k=1,pver
      do icol=1,ncol
        do ib=1,nbmodes
          Cam(icol,k,ib)=0.0_r8
          fcm(icol,k,ib)=0.00001_r8
          faqm(icol,k,ib)=0.00001_r8
          fbcm(icol,k,ib)=0.00001_r8
        end do
        rhum(icol, k)=0.00001_r8
        fnbc(icol, k)=0.00001_r8
        faitbc(icol,k)=0.00001_r8
        f_soana(icol,k)=0.00001_r8
        rhum(icol, k)=0.995_r8
        deltah_km(icol,k)=0.5_r8

      end do
    end do

    do icol=1,ncol
      daylight(icol) = .true. 
    end do

    Nnatk(:ncol,:,:) = Nnatk(:ncol,:,:)*1.e-6_r8
    cam(:ncol,:,:)=cam(:ncol,:,:)*1.e9_r8

    call inputForInterpol (lchnk, ncol, rhum, xrh, irh1,   &
    f_soana, xfombg, ifombg1, faitbc, xfbcbg, ifbcbg1,  &
    fnbc, xfbcbgn, ifbcbgn1, Nnatk, Cam, xct, ict1,     &
    focm, fcm, xfac, ifac1, fbcm, xfbc, ifbc1, faqm, xfaq, ifaq1)
    
    do iloop=1,1   
      !------------------------Get Optical Properties---------------------------------
        call interpol5to10 (lchnk, ncol, daylight, xrh, irh1, &
                            Nnatk, xct, ict1, xfac, ifac1, &
                            xfbc, ifbc1, xfaq, ifaq1, ssa, asym, be, ke, lw_on, kalw)
    end do

    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!     SW Optical properties of total aerosol:
    do ib=1,nbands
      do k=1,pver
        do icol=1,ncol
          betot(icol,k,ib)=0.0_r8
          ssatot(icol,k,ib)=0.0_r8
          asymtot(icol,k,ib)=0.0_r8
        end do
      end do  
    end do
    do ib=1,nbands
      do i=0,nmodes
      do k=1,pver
        do icol=1,ncol
          betot(icol,k,ib)=betot(icol,k,ib)+Nnatk(icol,k,i)*be(icol,k,i,ib)
          ssatot(icol,k,ib)=ssatot(icol,k,ib)+Nnatk(icol,k,i) &
                            *be(icol,k,i,ib)*ssa(icol,k,i,ib)           
          asymtot(icol,k,ib)=asymtot(icol,k,ib)+Nnatk(icol,k,i) &
                        *be(icol,k,i,ib)*ssa(icol,k,i,ib)*asym(icol,k,i,ib)
        end do
      end do
      end do
    end do

    do ib=1,nbands
      do k=1,pver
      do icol=1,ncol
        ssatot(icol,k,ib)=ssatot(icol,k,ib)/(betot(icol,k,ib)+eps)
        asymtot(icol,k,ib)=asymtot(icol,k,ib) &
                          /(betot(icol,k,ib)*ssatot(icol,k,ib)+eps)
      end do
      end do
    end do


    do i=1,ncol  ! zero aerosol in the top layer
      do ib=1,14 ! 1-nbands
          per_tau(i,0,ib)= 0._r8
          per_tau_w(i,0,ib)= 0.999_r8
          per_tau_w_g(i,0,ib)= 0.5_r8
      end do
    end do
    do i=1,ncol  ! zero aerosol in the top layer
      do ib=1,14  ! initialize also for the other layers
        do k=1,pver
          per_tau(i,k,ib)= 0._r8
          per_tau_w(i,k,ib)= 0.999_r8
          per_tau_w_g(i,k,ib)= 0.5_r8
        end do
      end do
    end do
  !      Remapping of SW wavelength bands from AeroTab to CAM5
      do i=1,ncol
      do ib=1,13
        do k=1,pver
          per_tau(i,k,ib)=deltah_km(i,k)*betot(i,k,14-ib)           
          per_tau_w(i,k,ib)=per_tau(i,k,ib)*max(min(ssatot(i,k,14-ib),0.999999_r8),1.e-6_r8)
          per_tau_w_g(i,k,ib)=per_tau_w(i,k,ib)*asymtot(i,k,14-ib)
        end do
      end do
        ib=14
        do k=1,pver
          per_tau(i,k,ib)=deltah_km(i,k)*betot(i,k,ib)
          per_tau_w(i,k,ib)=per_tau(i,k,ib)*max(min(ssatot(i,k,ib),0.999999_r8),1.e-6_r8)
          per_tau_w_g(i,k,ib)=per_tau_w(i,k,ib)*asymtot(i,k,ib)
        end do
    end do 
  end subroutine pmxsub_light
end module pmxsub