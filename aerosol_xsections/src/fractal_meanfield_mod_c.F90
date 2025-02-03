module fractal_meanfield_mod_c
  implicit none

contains

  subroutine fractal_meanfield_c(do_miess_in, xl_in, xk_in, xn_in, nb_in, alpha_in, &
    df_in, rmon,xv, ang, rcore, xksh_in, xnsh_in, Qext, Qsca, gfac, rc) bind(c,name='fractal_meanfield')
    use fractal_meanfield_mod, only: fractal_meanfield
    use iso_c_binding, only: f => c_double, c_bool, c_int
    logical(c_bool),value, intent(in)                :: do_miess_in   !! do core shell mie calculation for monomers
    real(kind=f),value,intent(in)            :: xl_in         !! Wavelength [microns]
    real(kind=f),value,intent(in)            :: xk_in         !! imaginary index of refraction ! if do_miess_in = .true, core imaginary index
    real(kind=f),value,intent(in)            :: xn_in         !! real index of refraction      ! if do_miess_in = .true, core real index
    real(kind=f),value,intent(in)            :: nb_in         !! number of monomers
    real(kind=f),value,intent(in)            :: alpha_in      !! Packing coefficient
    real(kind=f),value,intent(in)            :: df_in         !! Fractal dimension
    real(kind=f),value,intent(in)            :: rmon          !! monomer size [microns]
    real(kind=f),value,intent(in)            :: xv            !! set to 1
    real(kind=f),value,intent(in)            :: ang           !! angle set to zero
    real(kind=f),value,intent(in)            :: rcore         !! core radius [microns]               ! only used if do_miess_in = .true.  
    real(kind=f),value,intent(in)            :: xksh_in       !! shell imaginary index of refraction ! only used if do_miess_in = .true.  
    real(kind=f),value,intent(in)            :: xnsh_in       !! shell real index of refraction      ! only used if do_miess_in = .true.  
    real(kind=f),intent(out)           :: Qext          !! EFFICIENCY FACTOR FOR EXTINCTION
    real(kind=f),intent(out)           :: Qsca          !! EFFICIENCY FACTOR FOR SCATTERING
    real(kind=f),intent(out)           :: gfac          !! asymmetry factor
    integer(c_int),intent(inout)              :: rc            !! return code, negative indicates failure
    logical :: do_miess_in_f
    do_miess_in_f = do_miess_in
    call fractal_meanfield(do_miess_in_f, xl_in, xk_in, xn_in, nb_in, alpha_in, &
    df_in, rmon,xv, ang, rcore, xksh_in, xnsh_in, Qext, Qsca, gfac, rc)
  end subroutine

end module