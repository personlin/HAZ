c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common001 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  1.2879,   1.2879,   1.1582,   -0.4289 / 
      data a2	  / -0.0102,  -0.0102,  -0.0673,    1.2224 / 
      data a3	  /  0.1627,   0.1627,   0.3759,    1.5170 / 
      data a4	  / -0.4696,  -0.4696,  -0.1079,   -0.0782 / 
      data a5	  / -1.4940,  -1.4940,  -1.2552,   -1.3684 / 
      data a6	  /  0.2493,   0.2493,   0.2181,    0.2369 / 
      data a7	  /  7.2372,   7.2372,   5.1799,    6.5817 / 
      data a8	  /  0.0473,   0.0473,   0.0766,    0.0162 / 
      data a9	  /  0.0064,   0.0064,   0.0247,    0.0096 / 
      data Mc1	/  5.6001,   5.6001,   5.7256,    5.5757 / 
      data Mc2	/  6.6637,   6.6637,   6.6473,    6.5889 / 
      data a10	/  0.0485,   0.0485,   0.0103,    0.1734 / 
      data a11	/ -0.1253,  -0.1253,  -0.1130,   -0.2023 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
         a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model001 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
	  
C     Compute the Hanging Wall term if needed. 
      a10h = 3.0 
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor 
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor 
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor 
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common002 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  0.5619,   0.5619,   1.7456,   -0.1204 / 
      data a2	  / -0.0808,  -0.0808,   0.1143,    0.3148 / 
      data a3	  /  0.0431,   0.0431,   0.2300,    1.8317 / 
      data a4	  / -0.2962,  -0.2962,  -0.4107,    0.3730 / 
      data a5	  / -1.4149,  -1.4149,  -1.5449,   -1.4235 / 
      data a6	  /  0.2683,   0.2683,   0.3057,    0.1607 / 
      data a7	  /  4.9661,   4.9661,   4.2839,    5.8676 / 
      data a8	  /  0.0609,   0.0609,   0.0698,    0.0116 / 
      data a9	  /  0.0177,   0.0177,   0.0334,    0.0147 / 
      data Mc1	/  5.3920,   5.3920,   5.1621,    6.1163 / 
      data Mc2	/  6.4915,   6.4915,   5.5450,    6.5876 / 
      data a10	/  0.0592,   0.0592,   0.0620,    0.1262 / 
      data a11	/ -0.1126,  -0.1126,  -0.1436,   -0.1086 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
         a11T      = a11(1)  

         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model002 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
      
C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common003 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )
	 
      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  1.0584,   1.0584,   1.4838,   -0.0535 / 
      data a2	  /  0.0359,   0.0359,  -0.0794,    1.2284 / 
      data a3	  /  0.2852,   0.2852,   0.4103,    1.6048 / 
      data a4	  / -0.2864,  -0.2864,  -0.2399,    0.1384 / 
      data a5	  / -1.4607,  -1.4607,  -1.3359,   -1.3429 / 
      data a6	  /  0.2387,   0.2387,   0.2168,    0.1991 / 
      data a7	  /  6.3262,   6.3262,   7.2651,    5.7542 / 
      data a8	  /  0.0593,   0.0593,   0.0838,    0.0449 / 
      data a9	  /  0.0012,   0.0012,   0.0141,   -0.0043 / 
      data Mc1	/  5.5229,   5.5229,   5.7782,    5.5103 / 
      data Mc2	/  6.5767,   6.5767,   6.6088,    6.5461 / 
      data a10	/  0.0932,   0.0932,   0.0457,    0.1645 / 
      data a11	/ -0.1192,  -0.1192,  -0.0700,   -0.1912 / 

C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model003 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
      
C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common004 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none
      
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  0.5050,   0.5050,   1.5348,    0.4774 / 
      data a2	  /  0.1514,   0.1514,   0.1776,    1.1181 / 
      data a3	  /  0.2210,   0.2210,   0.5088,    2.1317 / 
      data a4	  / -0.1429,  -0.1429,  -0.0490,    0.3187 / 
      data a5	  / -1.2796,  -1.2796,  -1.1567,   -1.2211 / 
      data a6	  /  0.2425,   0.2425,   0.1538,    0.0585 / 
      data a7	  /  3.6016,   3.6016,   6.2446,    7.8137 / 
      data a8	  /  0.0758,   0.0758,   0.0955,    0.0436 / 
      data a9	  /  0.0077,   0.0077,  -0.0003,    0.0025 / 
      data Mc1	/  5.6551,   5.6551,   5.5430,    5.7928 / 
      data Mc2	/  6.8449,   6.8449,   6.1635,    6.6020 / 
      data a10	/  0.0316,   0.0316,   0.1025,    0.1727 / 
      data a11	/ -0.1753,  -0.1753,  -0.0331,   -0.2930 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model004 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common005 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  1.7145,   1.7145,   1.4329,    0.0870 / 
      data a2	  / -0.1330,  -0.1330,  -0.0318,    1.0243 / 
      data a3	  / -0.1040,  -0.1040,   0.4805,    1.7607 / 
      data a4	  / -0.7006,  -0.7006,  -0.0462,    0.0683 / 
      data a5	  / -1.8911,  -1.8911,  -1.1945,   -1.2756 / 
      data a6	  /  0.3764,   0.3764,   0.2005,    0.1729 / 
      data a7	  /  7.2752,   7.2752,   5.4623,    5.9079 / 
      data a8	  /  0.0439,   0.0439,   0.1040,    0.0421 / 
      data a9	  /  0.0224,   0.0224,   0.0172,    0.0039 / 
      data Mc1	/  5.4114,   5.4114,   5.8048,    5.7941 / 
      data Mc2	/  6.9240,   6.9240,   7.3926,    6.7951 / 
      data a10	/  0.0586,   0.0586,   0.1171,    0.1868 / 
      data a11	/ -0.1329,  -0.1329,  -0.0691,   -0.0440 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model005 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common006 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  2.2338,   2.2338,   2.5560,   -0.2015 / 
      data a2	  /  0.3272,   0.3272,   0.1199,    0.5813 / 
      data a3	  /  0.1722,   0.1722,   0.2438,    1.6755 / 
      data a4	  / -0.2986,  -0.2986,  -0.4968,    0.0899 / 
      data a5	  / -1.7999,  -1.7999,  -1.6629,   -1.3164 / 
      data a6	  /  0.2576,   0.2576,   0.2686,    0.1831 / 
      data a7	  /  7.4611,   7.4611,   7.2618,    4.9050 / 
      data a8	  /  0.0145,   0.0145,   0.0480,    0.0132 / 
      data a9	  /  0.0175,   0.0175,   0.0128,    0.0159 / 
      data Mc1	/  5.6187,   5.6187,   5.5611,    5.8036 / 
      data Mc2	/  6.5426,   6.5426,   7.3868,    6.8294 / 
      data a10	/  0.1483,   0.1483,   0.0088,    0.0948 / 
      data a11	/ -0.0626,  -0.0626,  -0.1609,   -0.1001 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model006 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common007 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  1.2712,   1.2712,   1.7190,   -0.2125 / 
      data a2	  /  0.0299,   0.0299,   0.3245,    0.9319 / 
      data a3	  /  0.1272,   0.1272,   0.3021,    1.6181 / 
      data a4	  / -0.3903,  -0.3903,  -0.1940,    0.0114 / 
      data a5	  / -1.4163,  -1.4163,  -1.1634,   -1.2909 / 
      data a6	  /  0.2502,   0.2502,   0.1822,    0.2045 / 
      data a7	  /  6.6321,   6.6321,   5.6208,    7.8266 / 
      data a8	  /  0.0701,   0.0701,   0.0852,    0.0401 / 
      data a9	  /  0.0110,   0.0110,   0.0079,    0.0195 / 
      data Mc1	/  5.8451,   5.8451,   5.8578,    5.8669 / 
      data Mc2	/  7.0218,   7.0218,   7.1003,    6.7371 / 
      data a10	/  0.0471,   0.0471,   0.0964,    0.1633 / 
      data a11	/ -0.1436,  -0.1436,  -0.0575,   -0.3259 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model007 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common008 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none
 
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  0.6255,   0.6255,   1.1093,   -0.2083 / 
      data a2	  /  0.1363,   0.1363,   0.0328,    1.0637 / 
      data a3	  /  0.2042,   0.2042,   0.3703,    1.4599 / 
      data a4	  / -0.0380,  -0.0380,  -0.0014,   -0.1904 / 
      data a5	  / -1.1500,  -1.1500,  -1.1017,   -1.4072 / 
      data a6	  /  0.1896,   0.1896,   0.1940,    0.2663 / 
      data a7	  /  5.1050,   5.1050,   3.7748,    8.8245 / 
      data a8	  /  0.0783,   0.0783,   0.0883,    0.0179 / 
      data a9	  /  0.0054,   0.0054,   0.0178,    0.0081 / 
      data Mc1	/  5.4204,   5.4204,   5.7387,    5.7101 / 
      data Mc2	/  6.4630,   6.4630,   7.4409,    6.7753 / 
      data a10	/  0.1054,   0.1054,   0.0340,    0.1707 / 
      data a11	/ -0.0473,  -0.0473,  -0.1329,   -0.1207 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model008 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common009 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  0.7491,   0.7491,   1.0521,   -0.7153 / 
      data a2	  / -0.0910,  -0.0910,  -0.2409,    0.4855 / 
      data a3	  /  0.1419,   0.1419,   0.5491,    1.5375 / 
      data a4	  / -0.3045,  -0.3045,  -0.1489,    0.1429 / 
      data a5	  / -1.3451,  -1.3451,  -1.1965,   -1.2894 / 
      data a6	  /  0.2320,   0.2320,   0.1863,    0.1934 / 
      data a7	  /  5.7062,   5.7062,   4.4666,    7.0431 / 
      data a8	  /  0.0495,   0.0495,   0.0712,   -0.0491 / 
      data a9	  /  0.0129,   0.0129,   0.0197,    0.0171 / 
      data Mc1	/  5.7402,   5.7402,   5.7253,    5.9617 / 
      data Mc2	/  6.4223,   6.4223,   6.8326,    6.5867 / 
      data a10	/  0.0152,   0.0152,   0.0253,    0.1065 / 
      data a11	/ -0.1596,  -0.1596,  -0.1149,   -0.1333 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model009 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 

c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common010 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  1.2305,   1.2305,   0.5693,   -0.5469 / 
      data a2	  / -0.0556,  -0.0556,  -0.1030,    0.8337 / 
      data a3	  /  0.1208,   0.1208,   0.4818,    1.4414 / 
      data a4	  / -0.3485,  -0.3485,  -0.1026,   -0.2134 / 
      data a5	  / -1.6250,  -1.6250,  -1.1063,   -1.4756 / 
      data a6	  /  0.2401,   0.2401,   0.1864,    0.3035 / 
      data a7	  /  7.2356,   7.2356,   2.5148,    4.1190 / 
      data a8	  /  0.0269,   0.0269,   0.0795,    0.0428 / 
      data a9	  /  0.0208,   0.0208,   0.0215,    0.0226 / 
      data Mc1	/  5.5960,   5.5960,   5.7523,    5.8555 / 
      data Mc2	/  6.4851,   6.4851,   6.2045,    6.8400 / 
      data a10	/  0.0581,   0.0581,   0.0216,    0.1129 / 
      data a11	/ -0.0711,  -0.0711,  -0.1212,   -0.0647 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model010 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 	  
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common011 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  0.0111,   0.0111,   2.0817,   -0.3963 / 
      data a2	  / -0.2061,  -0.2061,   0.0647,    0.9235 / 
      data a3	  /  0.2553,   0.2553,   0.4920,    1.4640 / 
      data a4	  / -0.3341,  -0.3341,  -0.3568,   -0.0085 / 
      data a5	  / -1.1935,  -1.1935,  -1.5144,   -1.3363 / 
      data a6	  /  0.2469,   0.2469,   0.2195,    0.2533 / 
      data a7	  /  4.7454,   4.7454,   7.7493,    3.5875 / 
      data a8	  /  0.0959,   0.0959,   0.0742,    0.0780 / 
      data a9	  /  0.0171,   0.0171,   0.0134,    0.0220 / 
      data Mc1	/  5.6123,   5.6123,   5.7596,    5.8060 / 
      data Mc2	/  6.6120,   6.6120,   5.8510,    6.7200 / 
      data a10	/  0.0933,   0.0933,   0.0812,    0.1377 / 
      data a11	/ -0.0613,  -0.0613,  -0.0607,   -0.2187 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model011 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common012 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  2.1402,   2.1402,   0.9268,   -0.1636 / 
      data a2	  / -0.1807,  -0.1807,   0.2047,    1.6408 / 
      data a3	  /  0.0454,   0.0454,   0.3969,    1.6917 / 
      data a4	  / -0.4102,  -0.4102,  -0.0991,   -0.0713 / 
      data a5	  / -2.0904,  -2.0904,  -1.0360,   -1.0426 / 
      data a6	  /  0.3475,   0.3475,   0.1748,    0.1627 / 
      data a7	  /  7.7344,   7.7344,   3.1985,    6.1390 / 
      data a8	  / -0.0106,  -0.0106,   0.1058,    0.0834 / 
      data a9	  /  0.0309,   0.0309,   0.0076,    0.0003 / 
      data Mc1	/  5.9059,   5.9059,   5.6733,    5.7334 / 
      data Mc2	/  6.6330,   6.6330,   6.4234,    6.7503 / 
      data a10	/  0.0875,   0.0875,   0.0476,    0.2554 / 
      data a11	/ -0.1161,  -0.1161,  -0.1035,   -0.0004 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model012 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common013 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  1.1314,   1.1314,   1.0562,   -0.0580 / 
      data a2	  / -0.1224,  -0.1224,   0.2467,    0.3603 / 
      data a3	  /  0.2023,   0.2023,   0.3102,    1.8187 / 
      data a4	  / -0.4955,  -0.4955,  -0.0312,    0.4600 / 
      data a5	  / -1.6983,  -1.6983,  -0.8950,   -1.2158 / 
      data a6	  /  0.3373,   0.3373,   0.1682,    0.1157 / 
      data a7	  /  6.1112,   6.1112,   3.4947,    3.8245 / 
      data a8	  /  0.0667,   0.0667,   0.1193,    0.0428 / 
      data a9	  /  0.0369,   0.0369,   0.0086,    0.0144 / 
      data Mc1	/  6.2692,   6.2692,   5.3782,    6.1807 / 
      data Mc2	/  7.4743,   7.4743,   7.4795,    6.6854 / 
      data a10	/  0.1053,   0.1053,   0.1547,    0.1205 / 
      data a11	/ -0.0980,  -0.0980,  -0.0250,   -0.0428 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model013 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common014 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none
      
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  2.1282,   2.1282,   1.8784,   -0.5201 / 
      data a2	  / -0.3135,  -0.3135,   0.1079,    1.1397 / 
      data a3	  / -0.1293,  -0.1293,   0.1798,    1.7294 / 
      data a4	  / -0.6486,  -0.6486,  -0.2516,    0.1059 / 
      data a5	  / -2.1494,  -2.1494,  -1.3718,   -1.1108 / 
      data a6	  /  0.4193,   0.4193,   0.2680,    0.1402 / 
      data a7	  /  7.6013,   7.6013,   4.9614,    3.7622 / 
      data a8	  / -0.0217,  -0.0217,   0.0958,    0.0306 / 
      data a9	  /  0.0469,   0.0469,   0.0219,    0.0174 / 
      data Mc1	/  5.2166,   5.2166,   5.9196,    5.7281 / 
      data Mc2	/  7.0892,   7.0892,   7.0443,    6.5485 / 
      data a10	/  0.0564,   0.0564,   0.0821,    0.1436 / 
      data a11	/ -0.1733,  -0.1733,  -0.1046,   -0.1764 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model014 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common015 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  0.8856,   0.8856,   1.9212,   -0.4046 / 
      data a2	  /  0.3055,   0.3055,   0.1075,    0.8047 / 
      data a3	  /  0.1288,   0.1288,   0.2544,    1.5391 / 
      data a4	  / -0.3615,  -0.3615,  -0.2172,    0.0483 / 
      data a5	  / -1.2319,  -1.2319,  -1.4832,   -1.1845 / 
      data a6	  /  0.2481,   0.2481,   0.2691,    0.1686 / 
      data a7	  /  4.7634,   4.7634,   4.9037,    8.5936 / 
      data a8	  /  0.0830,   0.0830,   0.0672,   -0.0240 / 
      data a9	  /  0.0211,   0.0211,   0.0321,    0.0109 / 
      data Mc1	/  5.9347,   5.9347,   5.8247,    5.7817 / 
      data Mc2	/  6.9656,   6.9656,   7.4137,    6.6662 / 
      data a10	/  0.0301,   0.0301,   0.0912,    0.1245 / 
      data a11	/ -0.1534,  -0.1534,  -0.0989,   -0.2741 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model015 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common016 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  0.9286,   0.9286,   2.2734,   -0.7258 / 
      data a2	  /  0.3526,   0.3526,  -0.1283,    0.6455 / 
      data a3	  /  0.4498,   0.4498,   0.4025,    1.2554 / 
      data a4	  / -0.0640,  -0.0640,  -0.3103,   -0.0874 / 
      data a5	  / -1.2195,  -1.2195,  -1.5419,   -1.3902 / 
      data a6	  /  0.1574,   0.1574,   0.2117,    0.2820 / 
      data a7	  /  5.5203,   5.5203,   7.2034,    7.3332 / 
      data a8	  /  0.0586,   0.0586,   0.0366,    0.0036 / 
      data a9	  /  0.0223,   0.0223,   0.0182,    0.0182 / 
      data Mc1	/  6.1405,   6.1405,   6.1553,    5.8398 / 
      data Mc2	/  6.8394,   6.8394,   7.1726,    6.6093 / 
      data a10	/  0.1578,   0.1578,   0.0729,    0.1163 / 
      data a11	/ -0.0712,  -0.0712,  -0.0590,   -0.2427 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model016 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
                                                                                                          

C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_Crustal_Common017 ( m, Rrup, Rjb, ztor, ftype, dip, Width, Rx, HWFlag, 
     1                                   specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a10(MAXPER),  a11(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER)
      real Mc1(MAXPER), Mc2(MAXPER), Rrup, Ztor, HWfactor
      REAL M, specT, sigma, Rjb, dip, Width, Rx, a10h
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL a7T, a8T, a9T, Mc1T, Mc2T, a10T, a11T, Fn, Frv, fType
      integer iflag, count1, count2, nPer, i, HWFlag

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1	  /  1.3616,   1.3616,   2.2072 ,  -0.0297 / 
      data a2	  / -0.1853,  -0.1853,  -0.0613 ,   1.4595 / 
      data a3	  /  0.1332,   0.1332,   0.4127 ,   1.6152 / 
      data a4	  / -0.2241,  -0.2241,  -0.2848 ,  -0.0337 / 
      data a5	  / -1.6346,  -1.6346,  -1.5829 ,  -1.2537 / 
      data a6	  /  0.2659,   0.2659,   0.2231 ,   0.1882 / 
      data a7	  /  7.0751,   7.0751,   7.4694 ,   8.0630 / 
      data a8	  /  0.0543,   0.0543,   0.0526 ,   0.0306 / 
      data a9	  /  0.0195,   0.0195,   0.0145 ,  -0.0033 / 
      data Mc1	/  6.1016,   6.1016,   6.0370 ,   5.5830 / 
      data Mc2	/  6.9046,   6.9046,   7.5530 ,   6.6911 / 
      data a10	/  0.1299,   0.1299,   0.0544 ,   0.2095 / 
      data a11	/ -0.0323,  -0.0323,  -0.0735 ,  -0.1772 / 



C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T      = a1(1)
         a2T      = a2(1)
         a3T      = a3(1)
         a4T      = a4(1)
         a5T      = a5(1)
         a6T      = a6(1)
         a7T      = a7(1)
         a8T      = a8(1)
         a9T      = a9(1)
         Mc1T      = Mc1(1)
         Mc2T      = Mc2(1) 
         a10T      = a10(1) 
          a11T      = a11(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Crustal Common Model017 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +             specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +             specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +             specT,a3T,iflag)
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +             specT,a4T,iflag)
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +             specT,a5T,iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +             specT,a6T,iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +             specT,a7T,iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +             specT,a8T,iflag)
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +             specT,a9T,iflag)
      call S24_interp (period(count1),period(count2),Mc1(count1),Mc1(count2),
     +             specT,Mc1T,iflag)
      call S24_interp (period(count1),period(count2),Mc2(count1),Mc2(count2),
     +             specT,Mc2T,iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +             specT,a10T,iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +             specT,a11T,iflag)
   
 1011 period1 = specT                                                                                                              

      if ( fType .gt. 0.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .lt. 0.0 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif

	  
C     Set Constant Terms

C     Compute the Hanging Wall term if needed.
      a10h = 3.0
      if (HWFlag .eq. 1) then
         call S02_SWUS_HWModel (M, Rrup, Rjb, Rx, Dip, Width, a10h, ztor, specT, HWfactor)
      else 
         HWFactor = 0.0      
      endif
	  
C     Compute the ground motion for the given spectral period. 
      if (M .lt. Mc1T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(Mc1T-Mc2T) + a3T*(M-Mc1T)
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      elseif (M .ge. Mc2T) then
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a4T*(M-Mc2T)  
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      else
         lnY = a1T - a8T**2*Rrup + a9T*Ztor + a2T*(M-Mc2T) 
     1         + (a5T+a6T*(M-5))*alog(sqrt(a7T*a7T+Rrup*Rrup)) + a10T*Frv + a11T*Fn + HWFactor
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 

c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common001 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 7.6010, 7.6010, 6.5524, 1.0265 / 
      data a2 / -1.9563, -1.9563, -1.5544, -0.8219 / 
      data a9 / 0.6396, 0.6396, 0.4062, 0.3273 / 
      data a6 / -0.0129, -0.0129, -0.0384, -0.0474 / 
      data a4 / 1.1360, 1.1360, 0.4059, 0.9647 / 
      data a5 / 0.4968, 0.4968, -0.0164, 0.1406 / 
      data a13 / -0.0167, -0.0167, -0.0474, -0.0551 / 
      data Mref / 7.7734, 7.7734, 7.3962, 7.3556 / 
      data a11 / 0.0198, 0.0198, 0.0151, 0.0141 / 
      data a3 / 0.0625, 0.0625, 0.0972, 0.0868 / 
      data c4 / 17.3667, 17.3667, 12.5760, 9.4682 / 
      data dsi / 61.4544, 61.4544, 98.5737, 68.9010 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model001 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common002 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 6.6315, 6.6315, 4.9742, 0.3675 / 
      data a2 / -1.7536, -1.7536, -1.3114, -0.9139 / 
      data a9 / 0.5059, 0.5059, 0.5196, 0.5048 / 
      data a6 / -0.0270, -0.0270, -0.0414, -0.0348 / 
      data a4 / 1.0383, 1.0383, 0.3076, 1.7792 / 
      data a5 / 0.3197, 0.3197, 0.1108, 0.4393 / 
      data a13 / -0.0007, -0.0007, -0.0541, 0.0315 / 
      data Mref / 7.9432, 7.9432, 7.7129, 7.0885 / 
      data a11 / 0.0145, 0.0145, 0.0195, 0.0081 / 
      data a3 / 0.0748, 0.0748, 0.0810, 0.1024 / 
      data c4 / 14.1282, 14.1282, 7.3092, 6.5189 / 
      data dsi / 93.7570, 93.7570, 60.8882, 86.2021 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model002 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common003 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 3.7774, 3.7774, 7.0059, 1.4244 / 
      data a2 / -1.2463, -1.2463, -1.6520, -1.1243 / 
      data a9 / 0.6534, 0.6534, 0.3836, 0.5270 / 
      data a6 / -0.0499, -0.0499, -0.0487, -0.0255 / 
      data a4 / 0.3452, 0.3452, 0.2617, 1.8165 / 
      data a5 / 0.0039, 0.0039, 0.1626, 0.7126 / 
      data a13 / -0.0505, -0.0505, -0.0774, 0.0108 / 
      data Mref / 7.3320, 7.3320, 7.7770, 7.2068 / 
      data a11 / 0.0165, 0.0165, 0.0154, 0.0121 / 
      data a3 / 0.1120, 0.1120, 0.0755, 0.0516 / 
      data c4 / 1.3088, 1.3088, 12.9969, 9.3220 / 
      data dsi / 69.0272, 69.0272, 108.4136, 73.1003 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model003 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end  
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common004 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none
      
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 3.8283, 3.8283, 6.3262, 0.5708 / 
      data a2 / -1.2365, -1.2365, -1.5653, -0.6131 / 
      data a9 / 0.2588, 0.2588, 0.3817, 1.1698 / 
      data a6 / -0.0653, -0.0653, -0.0312, -0.0694 / 
      data a4 / 0.3462, 0.3462, 0.3981, 0.7018 / 
      data a5 / -0.1812, -0.1812, 0.0448, 0.1398 / 
      data a13 / -0.0160, -0.0160, -0.0587, -0.0928 / 
      data Mref / 7.1206, 7.1206, 7.9352, 7.6605 / 
      data a11 / 0.0068, 0.0068, 0.0152, 0.0108 / 
      data a3 / 0.1351, 0.1351, 0.0769, 0.1083 / 
      data c4 / 2.3068, 2.3068, 6.8757, 1.6528 / 
      data dsi / 140.0033, 140.0033, 94.0696, 80.3279 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model004 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common005 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 5.2566, 5.2566, 6.3325, 2.2480 / 
      data a2 / -1.4931, -1.4931, -1.6338, -1.0795 / 
      data a9 / 0.4911, 0.4911, 0.3316, 0.2266 / 
      data a6 / -0.0402, -0.0402, -0.0254, -0.0411 / 
      data a4 / 0.6733, 0.6733, 0.6955, 1.3152 / 
      data a5 / 0.0601, 0.0601, 0.0060, 0.5264 / 
      data a13 / -0.0253, -0.0253, -0.0205, -0.0367 / 
      data Mref / 7.4836, 7.4836, 7.6130, 7.4294 / 
      data a11 / 0.0160, 0.0160, 0.0170, 0.0088 / 
      data a3 / 0.0928, 0.0928, 0.0926, 0.0648 / 
      data c4 / 5.2613, 5.2613, 4.9944, 16.0688 / 
      data dsi / 87.1517, 87.1517, 81.8567, 83.9385 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model005 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common006 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.4922, 4.4922, 6.7169, 0.8300 / 
      data a2 / -1.3198, -1.3198, -1.6517, -0.7128 / 
      data a9 / 0.3089, 0.3089, 0.5400, 0.5829 / 
      data a6 / -0.0381, -0.0381, -0.0418, -0.0556 / 
      data a4 / 0.2321, 0.2321, 0.9988, 0.8864 / 
      data a5 / -0.4776, -0.4776, 0.3316, 0.2955 / 
      data a13 / -0.0319, -0.0319, -0.0205, -0.0656 / 
      data Mref / 7.2758, 7.2758, 7.3512, 7.3714 / 
      data a11 / 0.0128, 0.0128, 0.0163, 0.0046 / 
      data a3 / 0.1566, 0.1566, 0.0588, 0.1072 / 
      data c4 / 5.5454, 5.5454, 4.6080, 3.7810 / 
      data dsi / 102.6698, 102.6698, 93.1004, 85.4824 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model006 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common007 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.0993, 4.0993, 6.4873, 1.4771 / 
      data a2 / -1.2073, -1.2073, -1.5296, -0.8852 / 
      data a9 / 0.2077, 0.2077, 0.4319, 0.4509 / 
      data a6 / -0.0597, -0.0597, -0.0529, -0.0495 / 
      data a4 / -0.3052, -0.3052, 0.7997, 0.9847 / 
      data a5 / -0.3777, -0.3777, 0.1360, 0.3975 / 
      data a13 / -0.0811, -0.0811, -0.0307, -0.0606 / 
      data Mref / 6.5863, 6.5863, 7.3565, 7.1541 / 
      data a11 / 0.0074, 0.0074, 0.0139, 0.0041 / 
      data a3 / 0.1199, 0.1199, 0.0719, 0.0938 / 
      data c4 / 5.7668, 5.7668, 6.1596, 10.4666 / 
      data dsi / 126.3711, 126.3711, 100.5852, 88.7429 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model007 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common008 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none
 
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 5.1307, 5.1307, 5.9711, 1.5940 / 
      data a2 / -1.5487, -1.5487, -1.4710, -0.8534 / 
      data a9 / 0.3690, 0.3690, 0.2758, 0.2700 / 
      data a6 / -0.0331, -0.0331, -0.0319, -0.0504 / 
      data a4 / 0.4767, 0.4767, 0.5815, 0.6185 / 
      data a5 / 0.2524, 0.2524, 0.0830, 0.0163 / 
      data a13 / -0.0290, -0.0290, -0.0266, -0.0842 / 
      data Mref / 7.5607, 7.5607, 7.4862, 7.2653 / 
      data a11 / 0.0134, 0.0134, 0.0169, 0.0086 / 
      data a3 / 0.0787, 0.0787, 0.0915, 0.1143 / 
      data c4 / 15.5491, 15.5491, 8.1545, 13.8522 / 
      data dsi / 87.1045, 87.1045, 95.2334, 77.9879 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model008 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common009 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 5.3266, 5.3266, 6.6594, 0.4704 / 
      data a2 / -1.6325, -1.6325, -1.5621, -0.6435 / 
      data a9 / 0.6781, 0.6781, 0.5405, 0.8438 / 
      data a6 / -0.0177, -0.0177, -0.0412, -0.0600 / 
      data a4 / 1.1049, 1.1049, 0.4688, 0.7495 / 
      data a5 / 0.6079, 0.6079, 0.2201, 0.0232 / 
      data a13 / -0.0086, -0.0086, -0.0536, -0.0665 / 
      data Mref / 7.7873, 7.7873, 7.6401, 7.4228 / 
      data a11 / 0.0202, 0.0202, 0.0196, 0.0118 / 
      data a3 / 0.0393, 0.0393, 0.0768, 0.1085 / 
      data c4 / 11.6308, 11.6308, 13.3721, 4.8046 / 
      data dsi / 40.8043, 40.8043, 67.8297, 81.3104 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model009 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 

c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common010 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.1489, 4.1489, 4.7173, 0.1842 / 
      data a2 / -1.4127, -1.4127, -1.1855, -0.8482 / 
      data a9 / 0.4411, 0.4411, 0.6867, 0.7215 / 
      data a6 / -0.0332, -0.0332, -0.0703, -0.0333 / 
      data a4 / 0.4451, 0.4451, 0.1021, 1.5905 / 
      data a5 / -0.2411, -0.2411, 0.0132, 0.3161 / 
      data a13 / -0.0269, -0.0269, -0.0806, 0.0161 / 
      data Mref / 7.2384, 7.2384, 7.3130, 7.2780 / 
      data a11 / 0.0186, 0.0186, 0.0180, 0.0125 / 
      data a3 / 0.1433, 0.1433, 0.0927, 0.0915 / 
      data c4 / 3.7703, 3.7703, 9.9141, 4.6502 / 
      data dsi / 64.8118, 64.8118, 59.1310, 80.0262 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model010 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 	  
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common011 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 6.3214, 6.3214, 5.7330, 0.4810 / 
      data a2 / -1.7507, -1.7507, -1.4570, -0.7149 / 
      data a9 / 0.2199, 0.2199, 0.4098, 0.7081 / 
      data a6 / -0.0213, -0.0213, -0.0583, -0.0551 / 
      data a4 / 0.8063, 0.8063, 0.2276, 0.8726 / 
      data a5 / -0.0681, -0.0681, 0.2076, -0.0774 / 
      data a13 / -0.0052, -0.0052, -0.0747, -0.0566 / 
      data Mref / 7.7219, 7.7219, 7.2074, 7.4848 / 
      data a11 / 0.0182, 0.0182, 0.0151, 0.0135 / 
      data a3 / 0.1246, 0.1246, 0.0670, 0.1249 / 
      data c4 / 10.4318, 10.4318, 11.7409, 6.1431 / 
      data dsi / 84.5264, 84.5264, 104.5682, 77.1699 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model011 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common012 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 6.0139, 6.0139, 5.3662, 1.1045 / 
      data a2 / -1.6224, -1.6224, -1.4087, -0.8415 / 
      data a9 / 0.4369, 0.4369, 0.5209, 0.4863 / 
      data a6 / -0.0334, -0.0334, -0.0576, -0.0585 / 
      data a4 / 0.7709, 0.7709, 0.5013, 0.9889 / 
      data a5 / -0.1388, -0.1388, 0.3414, 0.2105 / 
      data a13 / -0.0255, -0.0255, -0.0544, -0.0642 / 
      data Mref / 7.6515, 7.6515, 7.4414, 7.2496 / 
      data a11 / 0.0151, 0.0151, 0.0148, 0.0089 / 
      data a3 / 0.1117, 0.1117, 0.0515, 0.1188 / 
      data c4 / 4.5970, 4.5970, 4.7594, 8.5847 / 
      data dsi / 88.4311, 88.4311, 101.8119, 80.4040 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model012 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common013 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 6.1544, 6.1544, 5.0165, 1.7735 / 
      data a2 / -1.5918, -1.5918, -1.3094, -0.9073 / 
      data a9 / 0.5370, 0.5370, 0.9003, 0.3817 / 
      data a6 / -0.0348, -0.0348, -0.0567, -0.0603 / 
      data a4 / 0.8008, 0.8008, 0.6633, 1.1443 / 
      data a5 / -0.0462, -0.0462, 0.2486, 0.5699 / 
      data a13 / -0.0247, -0.0247, -0.0476, -0.0672 / 
      data Mref / 7.6204, 7.6204, 7.4097, 7.4593 / 
      data a11 / 0.0148, 0.0148, 0.0177, 0.0086 / 
      data a3 / 0.0988, 0.0988, 0.0655, 0.0611 / 
      data c4 / 6.4457, 6.4457, 0.5961, 10.2027 / 
      data dsi / 94.1379, 94.1379, 65.7016, 81.5994 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model013 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common014 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none
      
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.6394, 4.6394, 6.4821, 1.8417 / 
      data a2 / -1.3217, -1.3217, -1.5732, -0.9433 / 
      data a9 / 0.1903, 0.1903, 0.4001, 0.1561 / 
      data a6 / -0.0397, -0.0397, -0.0407, -0.0518 / 
      data a4 / 0.0804, 0.0804, 0.7176, 1.1255 / 
      data a5 / -0.4410, -0.4410, -0.0111, 0.5776 / 
      data a13 / -0.0473, -0.0473, -0.0281, -0.0569 / 
      data Mref / 7.2315, 7.2315, 7.3453, 7.1780 / 
      data a11 / 0.0120, 0.0120, 0.0135, 0.0060 / 
      data a3 / 0.1369, 0.1369, 0.0893, 0.0663 / 
      data c4 / 5.9427, 5.9427, 2.4906, 13.8270 / 
      data dsi / 111.4712, 111.4712, 108.6521, 79.0768 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model014 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common015 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 5.3702, 5.3702, 6.5593, 3.0678 / 
      data a2 / -1.5065, -1.5065, -1.5348, -1.2869 / 
      data a9 / 0.6260, 0.6260, 0.4691, 0.5445 / 
      data a6 / -0.0368, -0.0368, -0.0300, -0.0363 / 
      data a4 / 0.8633, 0.8633, 0.6776, 1.4379 / 
      data a5 / 0.4420, 0.4420, -0.1784, 0.9867 / 
      data a13 / -0.0063, -0.0063, -0.0250, -0.0456 / 
      data Mref / 7.6397, 7.6397, 7.5294, 6.9958 / 
      data a11 / 0.0153, 0.0153, 0.0166, 0.0014 / 
      data a3 / 0.0650, 0.0650, 0.1268, 0.0573 / 
      data c4 / 14.1368, 14.1368, 3.6553, 18.4740 / 
      data dsi / 84.0762, 84.0762, 90.3494, 85.5721 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model015 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common016 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.8493, 4.8493, 6.1452, 1.4615 / 
      data a2 / -1.4707, -1.4707, -1.4313, -0.8889 / 
      data a9 / 0.6318, 0.6318, 0.2397, 0.3222 / 
      data a6 / -0.0518, -0.0518, -0.0592, -0.0435 / 
      data a4 / 0.3829, 0.3829, 0.3012, 0.9246 / 
      data a5 / 0.4741, 0.4741, -0.0550, 0.2915 / 
      data a13 / -0.0477, -0.0477, -0.0542, -0.0603 / 
      data Mref / 7.0911, 7.0911, 7.4798, 7.2608 / 
      data a11 / 0.0138, 0.0138, 0.0077, 0.0085 / 
      data a3 / 0.0710, 0.0710, 0.0771, 0.0840 / 
      data c4 / 16.0957, 16.0957, 10.1588, 11.4187 / 
      data dsi / 86.8582, 86.8582, 140.6040, 75.9709 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model016 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubIntra_Common017 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 5.9209, 5.9209, 6.2393, 2.4027 / 
      data a2 / -1.6806, -1.6806, -1.5484, -1.1370 / 
      data a9 / 0.3338, 0.3338, 0.3739, 0.4368 / 
      data a6 / -0.0453, -0.0453, -0.0449, -0.0410 / 
      data a4 / 0.3982, 0.3982, 0.6425, 1.0586 / 
      data a5 / 0.1463, 0.1463, 0.0951, 0.4804 / 
      data a13 / -0.0371, -0.0371, -0.0390, -0.0672 / 
      data Mref / 6.6640, 6.6640, 7.3289, 7.0792 / 
      data a11 / 0.0139, 0.0139, 0.0132, 0.0074 / 
      data a3 / 0.1106, 0.1106, 0.0695, 0.0793 / 
      data c4 / 12.6560, 12.6560, 8.6193, 16.0488 / 
      data dsi / 113.5005, 113.5005, 94.3497, 80.5859 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-intra Common Model017 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
      
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common001 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 3.7424, 3.7424, 3.8775, 1.0247 / 
      data a2 / -1.3644, -1.3644, -1.0503, -0.8087 / 
      data a9 / 0.2930, 0.2930, 0.4142, 0.3487 / 
      data a6 / 0.0443, 0.0443, 0.0626, 0.0616 / 
      data a4 / 0.3180, 0.3180, 0.0171, 0.9498 / 
      data a5 / 0.0771, 0.0771, -0.1192, 0.1362 / 
      data a13 / -0.0370, -0.0370, -0.0859, -0.0553 / 
      data Mref / 7.3602, 7.3602, 6.4274, 7.2721 / 
      data a11 / 0.0288, 0.0288, 0.0300, -0.0027 / 
      data a3 / 0.0946, 0.0946, 0.1029, 0.0975 / 
      data c4 / 10.4162, 10.4162, 7.3671, 11.5527 / 
      data dsi / 30.2891, 30.2891, 44.8829, 45.8685 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model001 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common002 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 5.3043, 5.3043, 6.2736, 1.4047 / 
      data a2 / -1.6592, -1.6592, -1.6213, -0.8688 / 
      data a9 / 0.4005, 0.4005, 0.3852, 0.3276 / 
      data a6 / 0.0429, 0.0429, 0.0422, 0.0546 / 
      data a4 / 0.3581, 0.3581, 0.2763, 0.7581 / 
      data a5 / 0.3253, 0.3253, 0.2038, 0.1679 / 
      data a13 / -0.0566, -0.0566, -0.0810, -0.0818 / 
      data Mref / 7.6478, 7.6478, 8.7756, 7.1762 / 
      data a11 / 0.0362, 0.0362, 0.0402, -0.0034 / 
      data a3 / 0.0839, 0.0839, 0.0649, 0.1098 / 
      data c4 / 13.7763, 13.7763, 15.7867, 11.6594 / 
      data dsi / 37.1124, 37.1124, 31.4039, 33.3315 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model002 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common003 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 3.8164, 3.8164, 5.7130, 2.0769 / 
      data a2 / -1.3597, -1.3597, -1.5499, -1.0501 / 
      data a9 / 0.3580, 0.3580, 0.2663, 0.3705 / 
      data a6 / 0.0441, 0.0441, 0.0464, 0.0528 / 
      data a4 / 0.1181, 0.1181, 0.1041, 1.0282 / 
      data a5 / -0.3053, -0.3053, 0.1188, 0.7027 / 
      data a13 / -0.0558, -0.0558, -0.0985, -0.0737 / 
      data Mref / 7.3755, 7.3755, 7.6501, 7.0961 / 
      data a11 / 0.0214, 0.0214, 0.0241, -0.0016 / 
      data a3 / 0.1440, 0.1440, 0.0591, 0.0704 / 
      data c4 / 8.7174, 8.7174, 12.4920, 14.5969 / 
      data dsi / 39.0604, 39.0604, 52.0196, 30.1867 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model003 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end  
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common004 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none
      
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 3.6331, 3.6331, 5.4627, 1.8227 / 
      data a2 / -1.2616, -1.2616, -1.5036, -0.8205 / 
      data a9 / 0.2588, 0.2588, 0.3936, 0.4540 / 
      data a6 / 0.0523, 0.0523, 0.0481, 0.0639 / 
      data a4 / -0.1672, -0.1672, 0.3739, 0.7088 / 
      data a5 / -0.5452, -0.5452, 0.0766, 0.2037 / 
      data a13 / -0.0764, -0.0764, -0.0613, -0.1046 / 
      data Mref / 7.5591, 7.5591, 7.9215, 7.7277 / 
      data a11 / 0.0227, 0.0227, 0.0259, -0.0046 / 
      data a3 / 0.1775, 0.1775, 0.0823, 0.0835 / 
      data c4 / 7.6829, 7.6829, 9.8460, 10.5834 / 
      data dsi / 40.7574, 40.7574, 44.8569, 94.5097 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model004 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common005 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.0315, 4.0315, 6.4900, 1.9982 / 
      data a2 / -1.4700, -1.4700, -1.6380, -0.9515 / 
      data a9 / 0.2128, 0.2128, 0.2840, 0.3020 / 
      data a6 / 0.0441, 0.0441, 0.0419, 0.0543 / 
      data a4 / 0.4509, 0.4509, 0.5884, 0.8981 / 
      data a5 / 0.2250, 0.2250, 0.3369, 0.4681 / 
      data a13 / -0.0333, -0.0333, -0.0617, -0.0733 / 
      data Mref / 7.7262, 7.7262, 8.1901, 7.2721 / 
      data a11 / 0.0224, 0.0224, 0.0177, -0.0013 / 
      data a3 / 0.0695, 0.0695, 0.0372, 0.0783 / 
      data c4 / 8.2992, 8.2992, 12.1733, 13.8116 / 
      data dsi / 28.7042, 28.7042, 56.4178, 16.7598 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model005 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common006 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.7161, 4.7161, 4.9221, 1.8387 / 
      data a2 / -1.5401, -1.5401, -1.3365, -0.9124 / 
      data a9 / 0.2282, 0.2282, 0.2841, 0.4761 / 
      data a6 / 0.0351, 0.0351, 0.0610, 0.0533 / 
      data a4 / 0.5886, 0.5886, 0.2559, 0.9936 / 
      data a5 / -0.0864, -0.0864, 0.2783, 0.3269 / 
      data a13 / -0.0160, -0.0160, -0.0985, -0.0533 / 
      data Mref / 7.4555, 7.4555, 6.5394, 7.5154 / 
      data a11 / 0.0223, 0.0223, 0.0125, -0.0027 / 
      data a3 / 0.1237, 0.1237, 0.0588, 0.0746 / 
      data c4 / 10.1710, 10.1710, 8.7643, 12.7814 / 
      data dsi / 44.8458, 44.8458, 44.2255, 26.3434 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model006 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common007 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.3371, 4.3371, 5.6258, 1.4876 / 
      data a2 / -1.4579, -1.4579, -1.3290, -0.9861 / 
      data a9 / 0.2150, 0.2150, 0.3967, 0.4154 / 
      data a6 / 0.0440, 0.0440, 0.0543, 0.0429 / 
      data a4 / 0.4709, 0.4709, 0.1021, 1.3009 / 
      data a5 / 0.1458, 0.1458, -0.1466, 0.5452 / 
      data a13 / -0.0216, -0.0216, -0.1026, -0.0232 / 
      data Mref / 7.3513, 7.3513, 7.5145, 7.1491 / 
      data a11 / 0.0168, 0.0168, 0.0204, -0.0022 / 
      data a3 / 0.1063, 0.1063, 0.0979, 0.0684 / 
      data c4 / 11.2785, 11.2785, 7.5206, 13.3581 / 
      data dsi / 45.9636, 45.9636, 70.1679, 30.6900 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model007 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common008 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none
 
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.2338, 4.2338, 4.5736, 2.2716 / 
      data a2 / -1.4621, -1.4621, -1.2444, -1.0188 / 
      data a9 / 0.3610, 0.3610, 0.2913, 0.6655 / 
      data a6 / 0.0459, 0.0459, 0.0496, 0.0600 / 
      data a4 / 0.6605, 0.6605, 0.3519, 0.9468 / 
      data a5 / 0.2848, 0.2848, 0.1728, 0.5838 / 
      data a13 / -0.0066, -0.0066, -0.0704, -0.0823 / 
      data Mref / 7.1221, 7.1221, 6.3425, 7.3565 / 
      data a11 / 0.0246, 0.0246, 0.0351, -0.0070 / 
      data a3 / 0.1004, 0.1004, 0.0715, 0.0556 / 
      data c4 / 11.1780, 11.1780, 8.3425, 16.8740 / 
      data dsi / 45.1281, 45.1281, 51.5860, 47.4536 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model008 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common009 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.6888, 4.6888, 3.6625, 1.5377 / 
      data a2 / -1.4726, -1.4726, -0.9014, -0.9677 / 
      data a9 / 0.3634, 0.3634, 1.3694, 0.3981 / 
      data a6 / 0.0257, 0.0257, 0.0533, 0.0557 / 
      data a4 / 0.4345, 0.4345, 0.3430, 1.0165 / 
      data a5 / 0.0385, 0.0385, -0.1962, 0.3546 / 
      data a13 / -0.0370, -0.0370, -0.0453, -0.0545 / 
      data Mref / 7.2518, 7.2518, 9.1572, 7.2854 / 
      data a11 / 0.0389, 0.0389, 0.0364, -0.0007 / 
      data a3 / 0.0976, 0.0976, 0.1430, 0.0859 / 
      data c4 / 12.3550, 12.3550, 3.2002, 13.0072 / 
      data dsi / 50.1469, 50.1469, 40.0120, 41.4477 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model009 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 

c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common010 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 3.9684, 3.9684, 6.5625, 1.9887 / 
      data a2 / -1.3964, -1.3964, -1.4112, -1.0311 / 
      data a9 / 0.5189, 0.5189, 0.9364, 0.3983 / 
      data a6 / 0.0415, 0.0415, 0.0381, 0.0568 / 
      data a4 / 0.3623, 0.3623, 0.7416, 1.1997 / 
      data a5 / 0.1068, 0.1068, 0.1502, 0.5326 / 
      data a13 / -0.0410, -0.0410, -0.0427, -0.0497 / 
      data Mref / 7.4692, 7.4692, 9.2243, 7.3419 / 
      data a11 / 0.0350, 0.0350, 0.0375, -0.0021 / 
      data a3 / 0.0947, 0.0947, 0.0820, 0.0802 / 
      data c4 / 11.4784, 11.4784, 10.1992, 14.5644 / 
      data dsi / 38.4534, 38.4534, 57.7496, 18.4030 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model010 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 	  
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common011 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.4088, 4.4088, 5.2053, 1.9462 / 
      data a2 / -1.6350, -1.6350, -1.5099, -0.9627 / 
      data a9 / 0.3923, 0.3923, 0.3651, 0.3418 / 
      data a6 / 0.0412, 0.0412, 0.0434, 0.0631 / 
      data a4 / 0.6072, 0.6072, -0.1079, 0.8899 / 
      data a5 / 0.4414, 0.4414, 0.0137, 0.4852 / 
      data a13 / -0.0193, -0.0193, -0.1055, -0.0844 / 
      data Mref / 7.2688, 7.2688, 7.2309, 7.2066 / 
      data a11 / 0.0187, 0.0187, 0.0352, -0.0046 / 
      data a3 / 0.0566, 0.0566, 0.0808, 0.1003 / 
      data c4 / 12.1529, 12.1529, 11.0067, 14.5843 / 
      data dsi / 41.6949, 41.6949, 36.6896, 78.2170 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model011 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common012 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.0061, 4.0061, 4.3774, 2.2394 / 
      data a2 / -1.4827, -1.4827, -1.4047, -0.8579 / 
      data a9 / 0.3026, 0.3026, 0.3299, 0.4018 / 
      data a6 / 0.0489, 0.0489, 0.0500, 0.0631 / 
      data a4 / 0.4282, 0.4282, 0.2378, 0.7610 / 
      data a5 / 0.0520, 0.0520, -0.0657, 0.0929 / 
      data a13 / -0.0428, -0.0428, -0.0519, -0.0832 / 
      data Mref / 7.5922, 7.5922, 7.8303, 7.8494 / 
      data a11 / 0.0223, 0.0223, 0.0177, -0.0048 / 
      data a3 / 0.0953, 0.0953, 0.0943, 0.1132 / 
      data c4 / 8.2679, 8.2679, 6.1358, 12.4935 / 
      data dsi / 33.0608, 33.0608, 32.6260, 26.0587 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model012 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common013 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 5.7752, 5.7752, 6.1411, 2.0259 / 
      data a2 / -1.8356, -1.8356, -1.5855, -0.8226 / 
      data a9 / 0.1917, 0.1917, 0.0399, 0.3818 / 
      data a6 / 0.0444, 0.0444, 0.0618, 0.0571 / 
      data a4 / 0.5433, 0.5433, 0.6053, 0.6445 / 
      data a5 / 0.3016, 0.3016, 0.2385, 0.0997 / 
      data a13 / -0.0339, -0.0339, -0.0455, -0.0976 / 
      data Mref / 7.2671, 7.2671, 8.3490, 7.7738 / 
      data a11 / 0.0245, 0.0245, 0.0135, -0.0040 / 
      data a3 / 0.1025, 0.1025, 0.0496, 0.0936 / 
      data c4 / 15.2827, 15.2827, 14.0772, 11.3501 / 
      data dsi / 29.5546, 29.5546, 46.9990, 58.7743 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model013 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common014 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none
      
      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 3.7669, 3.7669, 5.1466, 1.4016 / 
      data a2 / -1.1842, -1.1842, -1.2774, -0.9265 / 
      data a9 / 0.2310, 0.2310, 0.0930, 0.4161 / 
      data a6 / 0.0315, 0.0315, 0.0695, 0.0453 / 
      data a4 / 0.5934, 0.5934, 0.6047, 1.4587 / 
      data a5 / -0.2999, -0.2999, -0.0152, 0.6386 / 
      data a13 / -0.0203, -0.0203, -0.0378, -0.0147 / 
      data Mref / 7.8361, 7.8361, 7.7247, 7.2088 / 
      data a11 / 0.0158, 0.0158, 0.0130, -0.0040 / 
      data a3 / 0.1285, 0.1285, 0.0963, 0.0511 / 
      data c4 / 5.6574, 5.6574, 10.8209, 12.5254 / 
      data dsi / 82.5864, 82.5864, 54.4898, 37.5935 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model014 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common015 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.4670, 4.4670, 6.2871, 1.4544 / 
      data a2 / -1.3772, -1.3772, -1.5227, -0.9064 / 
      data a9 / 0.3084, 0.3084, 0.0780, 0.4053 / 
      data a6 / 0.0365, 0.0365, 0.0413, 0.0413 / 
      data a4 / 0.5852, 0.5852, 0.3567, 0.8821 / 
      data a5 / 0.1011, 0.1011, 0.1898, 0.1146 / 
      data a13 / -0.0301, -0.0301, -0.0944, -0.0546 / 
      data Mref / 7.6234, 7.6234, 6.4150, 7.0310 / 
      data a11 / 0.0136, 0.0136, 0.0238, -0.0040 / 
      data a3 / 0.0972, 0.0972, 0.0563, 0.1054 / 
      data c4 / 9.5482, 9.5482, 12.8497, 11.9022 / 
      data dsi / 60.1158, 60.1158, 64.9590, 35.2634 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model015 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common016 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.6382, 4.6382, 5.5387, 2.1083 / 
      data a2 / -1.4006, -1.4006, -1.4005, -1.2265 / 
      data a9 / 0.3332, 0.3332, 0.2745, 0.6133 / 
      data a6 / 0.0346, 0.0346, 0.0377, 0.0355 / 
      data a4 / 0.3453, 0.3453, 0.7070, 1.3955 / 
      data a5 / 0.0096, 0.0096, 0.1410, 0.9071 / 
      data a13 / -0.0428, -0.0428, -0.0154, -0.0396 / 
      data Mref / 7.5658, 7.5658, 8.0771, 6.9316 / 
      data a11 / 0.0236, 0.0236, 0.0425, -0.0029 / 
      data a3 / 0.1092, 0.1092, 0.0892, 0.0322 / 
      data c4 / 13.7206, 13.7206, 12.1381, 16.0940 / 
      data dsi / 48.1439, 48.1439, 33.7854, 11.8031 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model016 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 
                        
c----------------------------------------------------------------------                
      Subroutine S04_SubInter_Common017 ( m, Rrup, ztor, specT, lnY, sigma, iflag )

      implicit none

      integer MAXPER
      parameter (MAXPER=4)
      REAL Period(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER),  a11(MAXPER), dsi(MAXPER)
      REAL a4(MAXPER), a5(MAXPER), a6(MAXPER), c4(MAXPER), a9(MAXPER), a13(MAXPER)
      real Mref(MAXPER), Rrup, Ztor
      REAL M, specT, sigma
      REAL period2, lnY, period1, a1T, a2T, a3T, a4T, a5T, a6T
      REAL c4T, a9T, MrefT, dsiT, a11T, a13T
      integer iflag, count1, count2, nPer, i

      data period / 0.0,  0.01,  0.2,   2    / 
      data a1 / 4.7516, 4.7516, 6.4542, 1.8656 / 
      data a2 / -1.4691, -1.4691, -1.6687, -0.9552 / 
      data a9 / 0.3255, 0.3255, 0.3431, 0.3579 / 
      data a6 / 0.0375, 0.0375, 0.0240, 0.0536 / 
      data a4 / 0.2997, 0.2997, 0.7894, 0.9446 / 
      data a5 / -0.1630, -0.1630, 0.2290, 0.3337 / 
      data a13 / -0.0492, -0.0492, -0.0275, -0.0674 / 
      data Mref / 7.4556, 7.4556, 7.5839, 7.3076 / 
      data a11 / 0.0247, 0.0247, 0.0297, -0.0026 / 
      data a3 / 0.1396, 0.1396, 0.0499, 0.0817 / 
      data c4 / 11.0406, 11.0406, 11.3051, 13.2976 / 
      data dsi / 56.5213, 56.5213, 71.3392, 41.4432 / 


C Find the requested spectral period and corresponding coefficients
      nPer = 4
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         a1T   = a1(1)   
         a2T   = a2(1)   
         a3T   = a3(1)   
         a4T   = a4(1)   
         a5T   = a5(1)   
         a6T   = a6(1)   
         a9T   = a9(1)   
         a11T  = a11(1)  
         a13T  = a13(1)  
         c4T   = c4(1)  
         dsiT  = dsi(1)  
         MrefT = Mref(1) 
         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Sub-inter Common Model017 Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),a1(count1),a1(count2),  
     +                 specT, a1T  ,iflag) 
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),    
     +                 specT, a2T  ,iflag) 
      call S24_interp (period(count1),period(count2),a3(count1),a3(count2),    
     +                 specT, a3T  ,iflag) 
      call S24_interp (period(count1),period(count2),a4(count1),a4(count2),    
     +                 specT, a4T  ,iflag) 
      call S24_interp (period(count1),period(count2),a5(count1),a5(count2),    
     +                 specT, a5T  ,iflag) 
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),    
     +                 specT, a6T  ,iflag) 
      call S24_interp (period(count1),period(count2),a9(count1),a9(count2),    
     +                 specT, a9T  ,iflag) 
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),    
     +                 specT, a11T ,iflag) 
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),    
     +                 specT, a13T ,iflag) 
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),    
     +                 specT, c4T  ,iflag) 
      call S24_interp (period(count1),period(count2),dsi(count1),dsi(count2),    
     +                 specT, dsiT ,iflag) 
      call S24_interp (period(count1),period(count2),Mref(count1),Mref(count2),    
     +                 specT, MrefT,iflag) 
   
 1011 period1 = specT                                                                                                              

C     Set Constant Terms

C     Compute the ground motion for the given spectral period. 
      if (M .le. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a4T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      elseif (M .gt. MrefT) then
         lnY = a1T - a6T**2*Rrup + ( a2T + a3T *(M-7.8))*alog(Rrup + c4T*exp(a9T*(M-6)))
     1         + a5T*(M-MrefT) + a13T*(10-M)**2 + a11T*(min(Ztor,dsiT)-0.5*dsiT)
      endif

      sigma = 0.65

c     Convert from g to gal                                                     
      lnY = lnY + 6.89                                                          

      period2 = period1

      return
      end 