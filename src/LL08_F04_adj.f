c ------------------------------------------------------------------
C *** Adjusted Lin and Lee (2008) Horizontal for Subduction Zones, ****
c ------------------------------------------------------------------

      subroutine S04_LL08_F04 ( mag, rupdist, specT, period, lnY, sigma,
     1  iflag, Ztor, ftype, vs30)

      implicit none

      integer MAXPER, nPer, i
      parameter (MAXPER=23)
      real mag, rupDist, lnY, sigma, period
      real ftype, vs30, Ztor, fmag, fztor
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c6(MAXPER), c4(MAXPER), c7(MAXPER), 
     1     c8(MAXPER), c9(MAXPER), c11(MAXPER), c8a(MAXPER), c8b(MAXPER)
      real specT, c1T, c2T, c3T, c4T, c5, c6T, c7T, c8T, c9T, c10, c11T, sigT, c8aT, c8bT
      integer count1, count2, iflag
      real period1(MAXPER), sig(MAXPER)

      Data period1 / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 
     1      4, 5/
      Data c1 / 0.251917569769788, 0.251917569769788, 0.184928810489714, -0.00356699594929908, -0.116139847179393,
     1      -0.359205997670388, -0.618494316568272, 0.0248887792514275, 0.316179568518992, 0.355069600474626, 0.324182771381062,
     1      0.244285329050222, -0.850876179346134, -2.03266332645047, -3.55038622993569, -4.33158184045748, -6.1861295462799,
     1      -8.07773284436953, -9.59355609932776, -10.2410101192885, -11.8197430743799, -12.3899805660359, -11.2967044082474/
      Data c2 / 1.205, 1.205, 1.2, 1.155, 1.1, 1.09, 1.040234713, 1, 1.04, 1.045, 1.065, 1.085, 1.125123456, 1.215, 1.285, 1.365,
     1      1.465, 1.62, 1.705, 1.77, 1.83, 1.845, 1.805/
      Data c3 / -1.895, -1.895, -1.88, -1.875, -1.86, -1.855, -1.826241507, -1.795, -1.77, -1.73, -1.71, -1.675, -1.61902357,
     1      -1.57, -1.5, -1.465, -1.45, -1.45, -1.44, -1.43, -1.37, -1.26, -1.135/
      Data c6 / 0.0190085489909667, 0.0190085489909667, 0.0185931283145121, 0.0180754695677972, 0.0171797889011337,
     1      0.0167193132892576, 0.0164094510351137, 0.0156026600924293, 0.015593690561164, 0.0157954136494296, 0.015862378077345,
     1      0.0157281115774864, 0.0156001059686599, 0.0153175632854892, 0.017301105547853, 0.017252423599451, 0.0168567176022536,
     1      0.0156338216453473, 0.0163600126300559, 0.0152530865748235, 0.0101434708937324, 0.00763487493439363,
     1      0.00532464316318695/
      Data c8 / 1.29844688061264, 1.29844688061264, 1.2835194865868, 1.26799743889615, 1.23822332712607, 1.21958036706955,
     1      1.14922912205638, 1.10227413787241, 1.1108672551511, 1.11259460482704, 1.12053373611138, 1.13907016871582,
     1      1.13672583329478, 1.13371267626325, 1.12842094617372, 1.19344814904958, 1.37332618873098, 1.46401067034632,
     1      1.66265594663281, 1.72311049113089, 1.85023182441134, 1.91072173395565, 1.86912098319861/
      Data c9 / 0.0486333778494577, 0.0486333778494577, 0.0482154946217135, 0.0475826239873693, 0.0475694652133262,
     1      0.0474481959420581, 0.0479124556154867, 0.0501970201626954, 0.0506450856670888, 0.0494002298939118, 0.0502381974755613,
     1      0.0489867406047882, 0.049890773264954, 0.049003637866689, 0.0441952205455153, 0.041048694342018, 0.0307131546663227,
     1      0.0231171395193354, 0.0203770597524757, 0.0167666379486008, 0.00812986804206466, 0.000517797127688047,
     1      0.00385222909361174/
      Data c11 / -0.378053724277046, -0.378053724277046, -0.374068643626873, -0.361098846225156, -0.333731258276522,
     1      -0.303834200934664, -0.29609635351472, -0.317140653194989, -0.344785011141604, -0.39096980447971, -0.425857702744424,
     1      -0.444672585568288, -0.476617127319015, -0.504570787921981, -0.556131879476854, -0.550089918340281, -0.666683380683071,
     1      -0.693735359043865, -0.734307716172396, -0.708057413209543, -0.644740108862218, -0.603603881941717, -0.597188074861842/
      Data c4 / 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.352485413751673, 0.447854533027533,
     1      0.51552, 0.51552, 0.51552, 0.51552, 0.51552, 0.51552/
      Data c8a / 1.65355293187812, 1.65355293187812, 1.62736235457987, 1.56435667800801, 1.47706996359888, 1.37766571107588,
     1      1.27733907587263, 1.33837679454979, 1.34328098282119, 1.40599239759458, 1.47744479533333, 1.53084985739141,
     1      1.61153440302746, 1.75769420123297, 2.12165232630608, 2.35235234835414, 2.63083202918182, 2.78551811856462,
     1      3.08882658786425, 3.25257617832983, 3.08949547761579, 2.84396915436152, 2.69627621202097/
      Data c8b / 0.511622226208753, 0.511622226208753, 0.518809110964449, 0.564783626354318, 0.595527969814917, 0.664267061484588,
     1      0.74475103031188, 0.621026682547829, 0.557187053999304, 0.523473227825398, 0.513978132919416, 0.501853735087157,
     1      0.673294096714701, 0.850400204354929, 1.03499364271083, 1.11815076143448, 1.36785600868129, 1.65091170888597,
     1      1.77566937710392, 1.77458455729069, 1.89622633793847, 1.82000048123534, 1.41926944323687/
      Data sig / 0.5218, 0.5218, 0.5189, 0.5235, 0.5352, 0.537, 0.569479308, 0.5806, 0.5748, 0.5817, 0.5906, 0.6059, 0.637738271,
     1      0.6656, 0.7105, 0.7145, 0.7689, 0.7983, 0.8411, 0.8766, 0.859, 0.8055, 0.7654/
      Data c7 /-0.023024774, -0.023024774, -0.009555338, 0.023863931, 0.078064967, 0.123007013, 0.191175881, 0.238934319,
     1         0.214408068, 0.183558925, 0.127185951, 0.038416195, -0.073966903, -0.136315954, -0.243353466, -0.341654652,
     1         -0.48202467, -0.510018746, -0.559795963, -0.57197707, -0.441301253, -0.380235394, -0.321436112/
     
     
C Find the requested spectral period and corresponding coefficients
      nPer = 23

C First check for the PGA case (i.e., specT=0.0)
      if (specT .eq. 0.0) then
         period  = period1(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c6T     = c6(1)
         c8T     = c8(1)
         c9T     = c9(1)
C         c10T    = c10(1)
         c11T    = c11(1)
         sigT    = sig(1)
         c8aT     = c8a(1)
         c8bT     = c8b(1)
         c7T     = c7(1)
         c4T     = c4(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period1(i) .and. specT .le. period1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020
            endif
         enddo
      endif

C      write (*,*)
C      write (*,*) 'Lin and Lee (2008) Sub-Hor. Adjusted atttenuation model'
C      write (*,*) 'is not defined for a spectral period of: '
C      write (*,*)') ' Period = ',specT
C      write (*,*) 'This spectral period is outside the defined'
C      write (*,*) 'period range in the code or beyond the range'
C      write (*,*) 'of spectral periods for interpolation.'
C      write (*,*) 'Please check the input file.'
C      write (*,*)
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period1(count1),period1(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period1(count1),period1(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period1(count1),period1(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period1(count1),period1(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period1(count1),period1(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period1(count1),period1(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
c      call S24_interp (period1(count1),period1(count2),c10(count1),c10(count2),
c     +             specT,c10T,iflag)
      call S24_interp (period1(count1),period1(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period1(count1),period1(count2),sig(count1),sig(count2),
     +             specT,sigT,iflag)
      call S24_interp (period1(count1),period1(count2),c8a(count1),c8a(count2),
     +             specT,c8aT,iflag)
      call S24_interp (period1(count1),period1(count2),c8b(count1),c8b(count2),
     +             specT,c8bT,iflag)
      call S24_interp (period1(count1),period1(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period1(count1),period1(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)

 1011 period = specT

C     Compute the ground motions.
C      c4 = 0.51552
      c5 = 0.63255
C      c7 = 0.275
      c10 = 0   !This effectively removes the fattn term.
      
      if (mag .le. 5.0) then
       fmag = c8bT * mag 
      elseif (mag .le. 5.5) then
        fmag = c8aT*(mag) + (c8bT-c8aT)*5.0
      elseif (mag .le. 7.5) then
        fmag = c8T * mag + (c8bT-c8aT)*5.0 + (c8aT-c8T)*5.5
      else
        fmag = c2T*(mag) + (c8bT-c8aT)*5.0 + (c8aT-c8T)*5.5 + (c8T-c2T)*7.5
      endif
      
      if (ftype .eq. 0) then
       fztor = c9T*(min(Ztor, 30.0) - 20.0)
      else
        fztor = c6T*(min(Ztor, 80.0) - 20.0)
      endif
      
      lnY = c1T + fmag + c3T*alog(Rupdist+c4T*exp(c5*mag)) + 
     1      c7T*ftype + fztor + c10*Rupdist + c11T*alog(Vs30/760.0)


      sigma = sigT

C     Now convert to Ln Units in gals.
      lnY = lnY + 6.89

      return
      end
