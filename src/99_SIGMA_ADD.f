C--------------------------------------------------------------------------------------  
C     Crustal SigmaSS Model - Taiwan SSHAC model , CV=0.2, 2018 revised 
C     August 2018
C     iBranch: 1 = Central, 2 = Upper, 3 = Lower

      subroutine S04_Taiwan_SigmaSS_Crustal ( mag, SigmaSS, iBranch)
      
      implicit none 
      real sigmaSS1(3), sigmaSS2(3), mag, sigmaSS
      integer iBranch

      data sigmaSS1 / 0.6062, 0.7703, 0.4560 /
      data sigmaSS2 / 0.5799, 0.7291, 0.4428 /


c     iBranch = 1 for central branch
c     iBranch = 2 for high branch
c     iBranch = 3 for low branch
      
c     Compute sigmaSS    
      if (mag .ge. 7.0 ) then 
        sigmaSS = sigmaSS2(iBranch)
      else
        sigmaSS = sigmaSS1(iBranch) + (mag-5.0)/2.0*(sigmaSS2(iBranch)-sigmaSS1(iBranch))
      endif

      return
      end

C--------------------------------------------------------------------------------------  
C     Subduction SigmaSS Model - Taiwan SSHAC model , CV=0.25, 2018 revised 
C     August 2018
C     iBranch: 1 = Central, 2 = Upper, 3 = Lower

      subroutine S04_Taiwan_SigmaSS_Sub ( mag, SigmaSS, iBranch)
      
      implicit none 
      real sigmaSS1(3), sigmaSS2(3), mag, sigmaSS, fM
      integer iBranch

      data sigmaSS1 / 0.5658, 0.7184, 0.4262 /
      data sigmaSS2 / 0.6101, 0.7606, 0.4713 /


c     iBranch = 1 for central branch
c     iBranch = 2 for high branch
c     iBranch = 3 for low branch
      
c     Compute sigmaSS    
      fM = min(7.0, max(6.0,mag))-6.0
      sigmaSS = sigmaSS1(iBranch) + (sigmaSS2(iBranch)-sigmaSS1(iBranch))*fM

      return
      end
	  
C--------------------------------------------------------------------------------------  
C     SigmaSS Model - Taiwan SSHAC model , CV=0.15
C     October 2017
C     iBranch: 1 = Central, 2 = Upper, 3 = Lower

      subroutine S04_Taiwan_SigmaSS__CV015 ( mag, SigmaSS, iBranch)
      
      implicit none 
      real sigmaSS1(3), sigmaSS2(3), mag, sigmaSS
      integer iBranch

      data sigmaSS1 /0.6098, 0.7406, 0.4878 /
      data sigmaSS2 /0.5819, 0.6954, 0.4753 /


c     iBranch = 1 for central branch
c     iBranch = 2 for high branch
c     iBranch = 3 for low branch
      
c     Compute sigmaSS    
      if (mag .ge. 7.0 ) then 
        sigmaSS = sigmaSS2(iBranch)
      else
        sigmaSS = sigmaSS1(iBranch) + (mag-5.0)/2.0*(sigmaSS2(iBranch)-sigmaSS1(iBranch))
      endif

      return
      end
	  
C--------------------------------------------------------------------------------------  
C     SigmaSS Model - Taiwan SSHAC model , CV=0.2
C     October 2017
C     iBranch: 1 = Central, 2 = Upper, 3 = Lower

      subroutine S04_Taiwan_SigmaSS_CV02 ( mag, SigmaSS, iBranch)
      
      implicit none 
      real sigmaSS1(3), sigmaSS2(3), mag, sigmaSS
      integer iBranch

      data sigmaSS1 / 0.6071, 0.7634, 0.4634 /
      data sigmaSS2 / 0.5788, 0.7232, 0.4457 /


c     iBranch = 1 for central branch
c     iBranch = 2 for high branch
c     iBranch = 3 for low branch
      
c     Compute sigmaSS    
      if (mag .ge. 7.0 ) then 
        sigmaSS = sigmaSS2(iBranch)
      else
        sigmaSS = sigmaSS1(iBranch) + (mag-5.0)/2.0*(sigmaSS2(iBranch)-sigmaSS1(iBranch))
      endif

      return
      end


C--------------------------------------------------------------------------------------  
C     PhiSS Model - Taiwan SSHAC model
C     June 2017
C     iBranch: 1 = Central, 3 = Lower, 2 = Upper

      subroutine S04_Taiwan_PhiSS_2017 ( specT, iBranch, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=24) 
      integer nPer, count1, count2, i, iflag, i1, iBranch
      real specT, period(MAXPER), PhiSSC(MAXPER), PhiSSH(MAXPER), PhiSSL(MAXPER)
      real PhiSS, PhiSSCT, PhiSSHT, PhiSSLT, period1
      
      data Period / 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 
     1              0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10 / 
      data PhissC / 4.71E-01, 4.71E-01, 4.77E-01, 4.84E-01, 4.87E-01, 4.85E-01, 4.79E-01,  
     1              4.73E-01, 4.71E-01, 4.68E-01, 4.73E-01, 4.79E-01, 4.85E-01, 4.96E-01,  
     1              5.03E-01, 4.99E-01, 4.88E-01, 4.72E-01, 4.54E-01, 4.38E-01, 4.29E-01,  
     1              4.18E-01, 4.00E-01, 3.81E-01 /
      data PhiSSH / 4.80E-01, 4.81E-01, 4.86E-01, 4.94E-01, 4.97E-01, 4.96E-01, 4.92E-01,  
     1              4.88E-01, 4.90E-01, 4.85E-01, 4.89E-01, 4.94E-01, 4.98E-01, 5.07E-01,  
     1              5.13E-01, 5.08E-01, 4.97E-01, 4.81E-01, 4.65E-01, 4.48E-01, 4.40E-01,  
     1              4.33E-01, 4.21E-01, 4.05E-01 /
      data PhiSSL / 4.62E-01, 4.61E-01, 4.68E-01, 4.74E-01, 4.77E-01, 4.74E-01, 4.66E-01,  
     1              4.58E-01, 4.53E-01, 4.51E-01, 4.57E-01, 4.64E-01, 4.72E-01, 4.85E-01,  
     1              4.93E-01, 4.90E-01, 4.79E-01, 4.63E-01, 4.43E-01, 4.28E-01, 4.18E-01,  
     1              4.03E-01, 3.79E-01, 3.57E-01 /    
	 
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSSCT = PhiSSC(i1)
        PhiSSHT = PhiSSH(i1)
        PhiSSLT = PhiSSL(i1)
        goto 5
      endif
      
      nPer = 24
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'Taiwan SSHAC PhiSS Model - June 2017'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSSC(count1), PhiSSC(count2),
     &    specT, PhiSSCT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSH(count1), PhiSSH(count2),
     &    specT, PhiSSHT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSL(count1), PhiSSL(count2),
     &    specT, PhiSSLT, iflag )

5     period1 = specT
      if (iBranch .eq. 1) then
         PhiSS = PhiSSCT
      elseif (iBranch .eq. 2) then 
         PhiSS = PhiSSHT
      elseif (iBranch .eq. 3) then 
         PhiSS = PhiSSLT
      endif
  
      return
      end 

C--------------------------------------------------------------------------------------  
C     PhiSS Model - Taiwan SSHAC model_Subduction (AGA16_C01)
C     June 2017
C     iBranch: 1 = Central, 2 = Upper, 3 = Lower

      subroutine S04_Taiwan_PhiSS_2017_Sub ( specT, iBranch, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=25) 
      integer nPer, count1, count2, i, iflag, i1, iBranch
      real specT, period(MAXPER), PhiSSC(MAXPER), PhiSSH(MAXPER), PhiSSL(MAXPER)
      real PhiSS, PhiSSCT, PhiSSHT, PhiSSLT, period1
      
      data Period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2,   
     1              0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10 / 
      data PhissC / 0.409117298, 0.408972331, 0.409568252, 0.413319126, 0.417848825, 0.423016386,   
     1               0.419957043, 0.417758955, 0.421538637, 0.421559765, 0.427086316, 0.424438357,   
     1               0.427894054, 0.432811208, 0.427497348, 0.426398424, 0.424738027, 0.419218648,   
     1               0.410078334, 0.413655599, 0.419596116, 0.38083348, 0.348258755, 0.32974697, 0.312483874 /
      data PhiSSH / 0.417855322, 0.41770523, 0.418314639, 0.422147358, 0.426774418, 0.432048719,    
     1              0.428929828, 0.426682762, 0.430542253, 0.430567451, 0.436212528, 0.433501769,    
     1              0.437021971, 0.442042054, 0.436617234, 0.43549417, 0.433811936, 0.42821915,    
     1              0.419032223, 0.423004254, 0.430251595, 0.391907903, 0.359846317, 0.345494802, 0.335111777 /
      data PhiSSL / 0.400188527, 0.400048842, 0.400630963, 0.404298168, 0.408728367, 0.413786937,    
     1              0.410788314, 0.408640318, 0.412338469, 0.412355357, 0.417760784, 0.415177135,    
     1              0.418567128, 0.423379153, 0.418178617, 0.417104376, 0.415465988, 0.410020621,    
     1              0.400924526, 0.404090721, 0.4086629, 0.369427226, 0.336272134, 0.313208349, 0.288084086 /    
	 
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSSCT = PhiSSC(i1)
        PhiSSHT = PhiSSH(i1)
        PhiSSLT = PhiSSL(i1)
        goto 5
      endif
      
      nPer = 24
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'Taiwan SSHAC PhiSS Model - June 2017'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSSC(count1), PhiSSC(count2),
     &    specT, PhiSSCT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSH(count1), PhiSSH(count2),
     &    specT, PhiSSHT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSL(count1), PhiSSL(count2),
     &    specT, PhiSSLT, iflag )

5     period1 = specT
      if (iBranch .eq. 1) then
         PhiSS = PhiSSCT
      elseif (iBranch .eq. 2) then 
         PhiSS = PhiSSHT
      elseif (iBranch .eq. 3) then 
         PhiSS = PhiSSLT
      endif      
  
      return
      end 
	  
C--------------------------------------------------------------------------------------  
C     Tau Model - NGA East Global model
C     Al-Atik 2015

      subroutine S04_Global_Tau_Linda ( mag, tau, iBranch)
      
      implicit none 
      real tau1(3), tau2(3), tau3(3), tau4(3), mag, tau
      integer iBranch
c     Tau values from Table 5.1 of NGA-East  Al-Atik 2015 report
      data tau1 /0.4436, 0.5706, 0.3280 /
      data tau2 /0.4169, 0.5551, 0.2928 /
      data tau3 /0.3736, 0.5214, 0.2439 /
      data tau4 /0.3415, 0.4618, 0.2343 /


c     iBranch = 1 for central branch
c     iBranch = 2 for high branch
c     iBranch = 3 for low branch
      
c     Compute tau      
      if (mag .gt. 6.5 ) then 
        tau = tau4(iBranch)
      elseif ( mag .le. 6.5 .and. mag .gt. 5.5  ) then
        tau = tau3(iBranch) + (tau4(iBranch)-tau3(iBranch))*(mag-5.5)
      elseif ( mag .le. 5.5 .and. mag .gt. 5.0 ) then
        tau = tau2(iBranch) + (tau3(iBranch)-tau2(iBranch))*(mag-5.0)/0.5
      elseif ( mag .le. 5.0 .and. mag .gt. 4.5 ) then
        tau = tau1(iBranch) + (tau2(iBranch)-tau1(iBranch))*(mag-4.5)/0.5
      else
        tau = tau1(iBranch)
      endif

      return
      end
	  
	  
C--------------------------------------------------------------------------------------  
C     Taiwan SSHAC Global Tau Model - Subduction
C     AGA16, ZH16

      subroutine S04_Global_Tau_Sub ( specT, tau, iBranch)
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      real tauC(MAXPER), tauH(MAXPER), tauL(MAXPER), tau
      integer iBranch, i1, nPer, i, count1, count2, iFlag
      real Period(MAXPER), specT, tauCT, tauHT, tauLT, period1
	  
      data Period/ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 
     1             0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5 /
      data tauC/ 0.422, 0.426, 0.436, 0.46, 0.482, 0.51, 0.492, 0.472, 0.44, 0.423,  
     1             0.406, 0.393, 0.381, 0.389, 0.399, 0.42, 0.424, 0.423, 0.412, 0.407, 0.396, 0.394 /
      data tauH/ 0.504, 0.508, 0.519, 0.566, 0.612, 0.679, 0.631, 0.582, 0.524, 0.479,  
     1             0.454, 0.46, 0.468, 0.466, 0.465, 0.464, 0.461, 0.481, 0.465, 0.452, 0.455, 0.456 /
      data tauL/ 0.339, 0.345, 0.352, 0.353, 0.352, 0.34, 0.353, 0.363, 0.356, 0.367,  
     1             0.358, 0.327, 0.294, 0.312, 0.333, 0.376, 0.387, 0.366, 0.36, 0.362, 0.337, 0.332 /


c     iBranch = 1 for central branch
c     iBranch = 2 for high branch
c     iBranch = 3 for low branch
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        tauCT = tauC(i1)
        tauHT = tauH(i1)
        tauLT = tauL(i1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'Taiwan SSHAC Global Tau Model - Subduction'
      write (*,*) 'Tau Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), tauC(count1), tauC(count2),
     &    specT, tauCT, iflag )
      call S24_interp ( period(count1), period(count2), tauH(count1), tauH(count2),
     &    specT, tauHT, iflag )
      call S24_interp ( period(count1), period(count2), tauL(count1), tauL(count2),
     &    specT, tauLT, iflag )
	 
5     period1 = specT
      if (iBranch .eq. 1) then
         tau = tauCT
      elseif (iBranch .eq. 2) then 
         tau = tauHT
      elseif (iBranch .eq. 3) then 
         tau = tauLT
      endif

      return
      end

