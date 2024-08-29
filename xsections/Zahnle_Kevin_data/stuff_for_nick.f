      implicit real*8(A-H,O-Z)
      
c new ion SW photolysis     (not yet used)    added May 2023
      real*8 short_ion(177),long_ion(177),ion_flux(177)
      real*8 ion_sigma(16,177)
      integer ishort_ion(177),ilong_ion(177)
c new neutral SW photolysis   (not yet used)  added May 2023
      real*8 short(177),long(177),sw_flux(177),sw_sigma(18,177)
      integer ishort(177),ilong(177)
      real*8  swflux(177)  
c new bins 
      real*8 wave(177), uv_flux(177), sigma(7,177)

      open(54, file='SW_ions_for_stays.dat',status='OLD')  ! photo-ionization
      open(53, file='SW_neutrals_for_stays.dat',status='OLD')  ! short lambda photolysis 

 100  FORMAT(/)
        read(53,100)  ! This skips the header lines
        do L=1,177
          read(53,547) ishort(L),ilong(L),sw_flux(L),
     $     (sw_sigma(K,L),K=1,18)
        enddo
 547    FORMAT(3X,2I7,19E11.2)
c       do L=1,177,20
c         print 544, ishort(L), ilong(L), sw_flux(L), 
c    $    (sw_sigma(K,L), K=1,18)
c       enddo
 544    format(1x,2I5,1P19E12.3)

        read(54,100)   ! This skips the header lines
        do L=1,177
          read(54,547) ishort_ion(L),ilong_ion(L),ion_flux(L),
     $     (ion_sigma(K,L),K=1,16)
        enddo

c       do L=1,177,1
c         print 544, ishort_ion(L), ilong_ion(L), ion_flux(L), 
c    $    (ion_sigma(K,L), K=1,16)
c       enddo

c  put the new TOTAL cross sections into properly labeled bins
        do L=1,177
           wave(L)      = (ishort(L) + ilong(L))/2.
           uv_flux(L)   = sw_flux(L)
           sigma(LH,L)    = ion_sigma(1,L)
           sigma(LH2,L)   = sw_sigma(15,L) + ion_sigma(2,L)
           sigma(LN2,L)   = sw_sigma(17,L) + ion_sigma(10,L)
           sigma(LCO2,L)  = sw_sigma(2,L) + ion_sigma(6,L)
           sigma(LH2O,L)  = sw_sigma(3,L) + ion_sigma(15,L)

c for testing, only 1-5 are in use
           sigma(LCO,L)   = sw_sigma(18,L) + ion_sigma(5,L)
           sigma(LO,L)    = ion_sigma(3,L)
c          sigma(LHe,L)   = ion_sigma(16,L)
c          sigma(LN,L)    = ion_sigma(9,L)
c          sigma(LC,L)    = ion_sigma(7,L)
c          sigma(LO2,L)   = sw_sigma(1,L) + ion_sigma(4,L)
c          sigma(LCH4,L)  = sw_sigma(13,L) + ion_sigma(12,L)
c          sigma(LNO,L)   = sw_sigma(16,L) + ion_sigma(11,L)
c          sigma(LHCN,L)  = sw_sigma(6,L) + ion_sigma(13,L)
c          sigma(LC2H2,L) = sw_sigma(4,L) + ion_sigma(14,L)
c          sigma(LNH3,L)  = sw_sigma(7,L) 
        enddo

    