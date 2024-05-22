module sulfuric_acid
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  
  public :: binary_saturation_pressure

contains

  !> So we can call from C
  subroutine binary_saturation_pressure_c(T, x_H2SO4, P_H2SO4, P_H2O) bind(c)
    use iso_c_binding
    real(c_double), value, intent(in) :: T
    real(c_double), value, intent(in) :: x_H2SO4
    real(c_double), intent(out) :: P_H2SO4
    real(c_double), intent(out) :: P_H2O
    call binary_saturation_pressure(T, x_H2SO4, P_H2SO4, P_H2O)
  end subroutine

  !> Computes the saturation vapor pressure of H2SO4 and H2O above a condensed sulfuric
  !> acid droplet. This a parameterization/fit to Equation 11 in Dai et al. (2023): 
  !> https://doi.org/10.1029/2021JE007060 , except ignoring their third term relevant to the
  !> curvature of the droplet's surface. Here, I assume the saturation vapor pressure of pure
  !> H2SO4 is given by Kulmala & Laaksonen (1990) (Equation 12 in Dai et al. (2023)), and
  !> the saturation vapor pressure of pure H2O is given by www.engineeringtoolbox.com .
  !> Furthermore, for the chemical potentials of Equation 11 of Dai et al. (2023), I use
  !> the thermodynamics laid out by Zeleznik (1991): https://doi.org/10.1063/1.555899 .
  !> Zeleznik (1991) is only valid from 150 K to 350 K, so here I have implemented a simple
  !> expontential extrapolation of their results to higher and lower temperatures.
  subroutine binary_saturation_pressure(T, x_H2SO4, P_H2SO4, P_H2O)
    real(dp), intent(in) :: T !! Temperature (K)
    real(dp), intent(in) :: x_H2SO4 !! Mole fraction of H2SO4 in a condensed droplet. Note that
                                    !! the mole fraction of H2O is 1 - x_H2SO4.
    real(dp), intent(out) :: P_H2SO4 !! Saturation vapor pressure of H2SO4 (dynes/cm^2)
    real(dp), intent(out) :: P_H2O !! Saturation vapor pressure of H2O (dynes/cm^2)

    ! Some constants for the parameterization
    real(dp), parameter :: mu = 18.01534_dp
    real(dp), parameter :: T_triple = 273.15_dp
    real(dp), parameter :: T_critical = 647.0_dp
    real(dp), parameter :: T_ref = 300.0_dp

    ! Grid of H2SO4 mole fractions that we will interpolate to.
    real(dp), parameter :: x_H2SO4_grid(*) = &
    [ &
      0.000000e+00_dp, 1.000000e-02_dp, 5.000000e-02_dp, 1.000000e-01_dp, 1.500000e-01_dp, &
      2.000000e-01_dp, 2.500000e-01_dp, 3.000000e-01_dp, 3.500000e-01_dp, 4.000000e-01_dp, &
      4.500000e-01_dp, 5.000000e-01_dp, 5.500000e-01_dp, 6.000000e-01_dp, 6.500000e-01_dp, &
      7.000000e-01_dp, 7.500000e-01_dp, 8.000000e-01_dp, 8.500000e-01_dp, 9.000000e-01_dp, &
      9.500000e-01_dp, 9.900000e-01_dp, 1.000000e+00_dp &
    ]

    ! Constants for H2SO4 saturation
    real(dp), parameter :: P_ref_H2SO4(*) = &
    [ &
      5.095729e-19_dp, 9.797323e-15_dp, 1.196130e-12_dp, 7.524852e-11_dp, 2.164168e-09_dp, &
      3.628208e-08_dp, 4.070169e-07_dp, 3.311650e-06_dp, 2.040525e-05_dp, 9.721374e-05_dp, &
      3.616182e-04_dp, 1.057906e-03_dp, 2.459897e-03_dp, 4.627574e-03_dp, 7.233160e-03_dp, &
      9.733141e-03_dp, 1.176242e-02_dp, 1.335216e-02_dp, 1.482752e-02_dp, 1.653790e-02_dp, &
      1.850000e-02_dp, 1.982962e-02_dp, 2.016085e-02_dp &
    ]
    real(dp), parameter :: A_v_H2SO4(*) = &
    [ &
      9.268866e+10_dp, 1.150058e+11_dp, 1.154358e+11_dp, 1.023884e+11_dp, 9.655652e+10_dp, &
      7.927856e+10_dp, 8.001911e+10_dp, 1.000506e+11_dp, 9.930541e+10_dp, 9.647878e+10_dp, &
      6.512191e+10_dp, 8.735501e+10_dp, 5.220553e+10_dp, 7.954502e+10_dp, 7.766621e+10_dp, &
      7.725585e+10_dp, 7.794300e+10_dp, 7.910873e+10_dp, 8.009908e+10_dp, 8.055334e+10_dp, &
      6.136681e+10_dp, 8.157246e+10_dp, 5.863745e+10_dp &
    ]
    real(dp), parameter :: B_v_H2SO4(*) = &
    [ &
      5.423812e+01_dp, -7.853300e+07_dp, -9.280866e+07_dp, -7.001541e+07_dp, -7.000210e+07_dp, &
      -3.277475e+07_dp, -4.999062e+07_dp, -1.253940e+08_dp, -1.334303e+08_dp, -1.328219e+08_dp, &
      -4.091893e+07_dp, -1.153205e+08_dp, -8.594643e+06_dp, -9.665820e+07_dp, -9.243321e+07_dp, &
      -9.233133e+07_dp, -9.540077e+07_dp, -9.985395e+07_dp, -1.037200e+08_dp, -1.058567e+08_dp, &
      -4.627851e+07_dp, -1.097433e+08_dp, -3.268146e+07_dp &
    ]
    real(dp), parameter :: A_s_H2SO4(*) = &
    [ &
      1.232585e+11_dp, 9.199536e+10_dp, 1.017863e+11_dp, 1.049495e+11_dp, 1.015327e+11_dp, &
      7.640003e+10_dp, 7.281877e+10_dp, 9.566832e+10_dp, 9.576816e+10_dp, 9.636312e+10_dp, &
      6.323842e+10_dp, 9.818473e+10_dp, 6.050828e+10_dp, 9.961581e+10_dp, 9.959902e+10_dp, &
      9.870468e+10_dp, 9.666506e+10_dp, 9.334849e+10_dp, 8.891694e+10_dp, 8.404251e+10_dp, &
      5.571668e+10_dp, 7.906592e+10_dp, 5.863745e+10_dp &
    ]
    real(dp), parameter :: B_s_H2SO4(*) = &
    [ &
      -6.085830e+01_dp, 0.000000e+00_dp, -4.099230e+07_dp, -7.863085e+07_dp, -8.816608e+07_dp, &
      0.000000e+00_dp, -1.016046e+02_dp, -1.103702e+08_dp, -1.216415e+08_dp, -1.336165e+08_dp, &
      -8.631619e+01_dp, -1.560730e+08_dp, 2.731753e+01_dp, -1.708874e+08_dp, -1.732486e+08_dp, &
      -1.711372e+08_dp, -1.640038e+08_dp, -1.518689e+08_dp, -1.357679e+08_dp, -1.183477e+08_dp, &
      -2.543527e+01_dp, -1.003070e+08_dp, -3.268146e+07_dp &
    ]
    real(dp), parameter :: A_c_H2SO4(*) = &
    [ &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp &
    ]
    real(dp), parameter :: B_c_H2SO4(*) = &
    [ &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp &
    ]

    ! Constants for H2O saturation
    real(dp), parameter :: P_ref_H2O(*) = &
    [ &
      3.518390e+04_dp, 3.446912e+04_dp, 3.003429e+04_dp, 2.154586e+04_dp, 1.336617e+04_dp, &
      7.362893e+03_dp, 3.654488e+03_dp, 1.651983e+03_dp, 6.891239e+02_dp, 2.704528e+02_dp, &
      1.025916e+02_dp, 3.892186e+01_dp, 1.535394e+01_dp, 6.549178e+00_dp, 3.121090e+00_dp, &
      1.690356e+00_dp, 1.028877e+00_dp, 6.655342e-01_dp, 4.046557e-01_dp, 1.857053e-01_dp, &
      4.427552e-02_dp, 3.347471e-03_dp, 6.456003e-07_dp &
    ]
    real(dp), parameter :: A_v_H2O(*) = &
    [ &
      2.841421e+10_dp, 4.086040e+10_dp, 4.081275e+10_dp, 4.185671e+10_dp, 4.264804e+10_dp, &
      4.259325e+10_dp, 4.198741e+10_dp, 4.151882e+10_dp, 4.189839e+10_dp, 4.361661e+10_dp, &
      4.680510e+10_dp, 5.118637e+10_dp, 5.610284e+10_dp, 6.062146e+10_dp, 6.371436e+10_dp, &
      6.451611e+10_dp, 6.265248e+10_dp, 5.123949e+10_dp, 5.395939e+10_dp, 5.087350e+10_dp, &
      -2.970220e+10_dp, 1.533113e+10_dp, -9.615150e+10_dp &
    ]
    real(dp), parameter :: B_v_H2O(*) = &
    [ &
      -1.399732e+07_dp, -5.570953e+07_dp, -5.515655e+07_dp, -5.694957e+07_dp, -5.681923e+07_dp, &
      -5.306463e+07_dp, -4.685036e+07_dp, -4.060079e+07_dp, -3.680970e+07_dp, -3.725898e+07_dp, &
      -4.256939e+07_dp, -5.202994e+07_dp, -6.367656e+07_dp, -7.460809e+07_dp, -8.154041e+07_dp, &
      -8.160080e+07_dp, -7.334625e+07_dp, -3.481727e+07_dp, -3.971756e+07_dp, -2.500159e+07_dp, &
      2.342411e+08_dp, 9.722722e+07_dp, 4.344690e+08_dp &
    ]
    real(dp), parameter :: A_s_H2O(*) = &
    [ &
      2.746884e+10_dp, 4.044188e+10_dp, 3.988834e+10_dp, 3.969948e+10_dp, 4.019651e+10_dp, &
      4.087980e+10_dp, 4.142348e+10_dp, 4.168605e+10_dp, 4.163211e+10_dp, 4.127069e+10_dp, &
      4.062853e+10_dp, 3.976110e+10_dp, 3.880194e+10_dp, 3.805062e+10_dp, 3.810030e+10_dp, &
      4.000389e+10_dp, 4.547076e+10_dp, 4.781122e+10_dp, 7.824908e+10_dp, 1.128670e+11_dp, &
      6.485085e+10_dp, 1.925886e+11_dp, 7.331063e+10_dp &
    ]
    real(dp), parameter :: B_s_H2O(*) = &
    [ &
      4.181527e+06_dp, -3.950446e+07_dp, -3.704729e+07_dp, -3.424137e+07_dp, -3.293380e+07_dp, &
      -3.177766e+07_dp, -2.967194e+07_dp, -2.601080e+07_dp, -2.056537e+07_dp, -1.336773e+07_dp, &
      -4.670571e+06_dp, 4.999876e+06_dp, 1.469284e+07_dp, 2.278076e+07_dp, 2.660376e+07_dp, &
      2.201467e+07_dp, 2.851450e+06_dp, 0.000000e+00_dp, -1.165092e+08_dp, -2.401535e+08_dp, &
      -1.688018e+01_dp, -5.392384e+08_dp, -2.328523e-03_dp &
    ]
    real(dp), parameter :: A_c_H2O(*) = &
    [ &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp, &
      1.793161e+12_dp, 1.793161e+12_dp, 1.793161e+12_dp &
    ]
    real(dp), parameter :: B_c_H2O(*) = &
    [ &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp, &
      0.000000e+00_dp, 0.000000e+00_dp, 0.000000e+00_dp &
    ]

    integer :: ind, ind1, ind2
    real(dp) :: x_H2SO4_, P_low, P_high

    ! Bound input H2SO4 mole fraction by numbers very close to 0 and 1
    x_H2SO4_ = min(max(x_H2SO4, tiny(1.0_dp)),1.0_dp - 1.0e-14_dp)

    ! Find the location in the grid that bounds input x_H2SO4
    ind1 = searchsorted(x_H2SO4_grid, x_H2SO4_) - 1 ! lower index
    ind2 = ind1 + 1 ! upper index

    !~~~ H2SO4 ~~~!
    ! SVP at lower x_H2SO4
    ind = ind1
    P_low = saturation_pressure_(mu, T_triple, T_critical, T_ref, P_ref_H2SO4(ind), &
                                 A_v_H2SO4(ind), B_v_H2SO4(ind), A_s_H2SO4(ind), &
                                 B_s_H2SO4(ind), A_c_H2SO4(ind), B_c_H2SO4(ind), T)
    ! SVP at upper x_H2SO4
    ind = ind2
    P_high = saturation_pressure_(mu, T_triple, T_critical, T_ref, P_ref_H2SO4(ind), &
                                 A_v_H2SO4(ind), B_v_H2SO4(ind), A_s_H2SO4(ind), &
                                 B_s_H2SO4(ind), A_c_H2SO4(ind), B_c_H2SO4(ind), T)

    ! log-linearly interpolate to estimate value
    P_H2SO4 = log(P_low) + (log(P_high) - log(P_low))*(x_H2SO4_ - x_H2SO4_grid(ind1))/(x_H2SO4_grid(ind2) - x_H2SO4_grid(ind1))
    P_H2SO4 = exp(P_H2SO4)

    !~~~ H2O ~~~!
    ! SVP at lower x_H2SO4
    ind = ind1
    P_low = saturation_pressure_(mu, T_triple, T_critical, T_ref, P_ref_H2O(ind), &
                                 A_v_H2O(ind), B_v_H2O(ind), A_s_H2O(ind), &
                                 B_s_H2O(ind), A_c_H2O(ind), B_c_H2O(ind), T)
    ! SVP at upper x_H2SO4
    ind = ind2
    P_high = saturation_pressure_(mu, T_triple, T_critical, T_ref, P_ref_H2O(ind), &
                                 A_v_H2O(ind), B_v_H2O(ind), A_s_H2O(ind), &
                                 B_s_H2O(ind), A_c_H2O(ind), B_c_H2O(ind), T)

    ! log-linearly interpolate to estimate value
    P_H2O = log(P_low) + (log(P_high) - log(P_low))*(x_H2SO4_ - x_H2SO4_grid(ind1))/(x_H2SO4_grid(ind2) - x_H2SO4_grid(ind1))
    P_H2O = exp(P_H2O)

  end subroutine

  function integral_(A, B, T) result(res)
    real(dp), intent(in) :: A, B, T
    real(dp) :: res
    res = -A/T + B*log(T)
  end function

  function saturation_pressure_(mu, T_triple, T_critical, T_ref, P_ref, A_v, B_v, A_s, B_s, A_c, B_c, T) result(Ps)
    real(dp), intent(in) :: mu, T_ref, P_ref, T_triple, T_critical, A_v, B_v, A_s, B_s, A_c, B_c, T
    real(dp) :: Ps !! dynes/cm^2
    real(dp), parameter :: R = 8.31446261815324e7_dp
    real(dp) :: tmp

    if (T > T_critical) then
      tmp = (integral_(A_v, B_v, T_critical) - integral_(A_v, B_v, T_ref)) + &
            (integral_(A_c, B_c, T) - integral_(A_c, B_c, T_critical))
    elseif (T <= T_critical .and. T > T_triple) then
      tmp = integral_(A_v, B_v, T) - integral_(A_v, B_v, T_ref)
    elseif (T <= T_triple) then
      tmp = (integral_(A_v, B_v, T_triple) - integral_(A_v, B_v, T_ref)) + &
              (integral_(A_s, B_s, T) - integral_(A_s, B_s, T_triple))
    endif

    Ps = P_ref*exp((mu/R)*(tmp))
  end function

  !> Mimics numpy.searchsorted
  pure function searchsorted(arr, val) result(ind)
    real(dp), intent(in) :: arr(:) !! Input sorted array
    real(dp), intent(in) :: val !! Value to compare to arr
    integer :: ind !! Index that satisfies arr(i-1) < val <= arr(i)

    integer :: low, high, mid

    if (val <= arr(1)) then
      ind = 1
      return
    endif

    if (val > arr(size(arr))) then
      ind = size(arr) + 1
      return
    endif

    low = 1
    high = size(arr)
    do
      mid = (low + high)/2
      if (val > arr(mid)) then
        low = mid
      else
        high = mid
      endif
      if (high-1 == low) exit
    enddo

    ind = high

  end function

end module
