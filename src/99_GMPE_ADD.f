c -------------------------------------------------------------------           
C **** Zhao et al. 2016 (BSSA, Vol 106, No.4) *************
c -------------------------------------------------------------------           

      subroutine S02_Zhaoetal2016_cru ( m, dist, ftype, lnY, sigma, sclass, specT,              
     1                   attenName, period1, iflag, sourcetype, depth, phiT, tauT )                                   

c     This  subroutine calculates the spectral acceleration from the 
c     Zhao et al. (2016) attenuation model which is defined for crustal
c     and upper mantle events. Four separate site classes are also modeled.
                                                                         
      parameter (MAXPER=37)                                                     
      real ftype, dist, m, lnY, sigma, specT, sclass, cm, cmax, r
      real period(MAXPER), Amax1(MAXPER), SRC1(MAXPER), Amax2(MAXPER)
      real SRC2(MAXPER), Amax3(MAXPER), SRC3(MAXPER), Amax4(MAXPER), SRC4(MAXPER)
      real AmSCI(MAXPER), fsrCR1(MAXPER), fsrCR2(MAXPER), fsrCR3(MAXPER)
      real fsrCR4(MAXPER), fsrUM1(MAXPER), fsrUM2(MAXPER), fsrUM3(MAXPER)
      real fsrUM4(MAXPER), c1(MAXPER), c2(MAXPER), ccr(MAXPER), dcr(MAXPER)
      real FCRN(MAXPER), FumRV(MAXPER), FumNS(MAXPER), bcr(MAXPER), gcr(MAXPER)
      real gUM(MAXPER), gcrN(MAXPER), gcrL(MAXPER), ecr(MAXPER), eum(MAXPER)
      real ecrv(MAXPER), gamma(MAXPER), S2(MAXPER), S3(MAXPER), S4(MAXPER)
      real phiS1(MAXPER), phiS2(MAXPER), phiS3(MAXPER), phiS4(MAXPER)
      real phi(MAXPER), tau(MAXPER)
      character*80 attenName                                                    
      integer nper, count1, count2, iflag, xcro
      real Amax1T, SRC1T, Amax2T, SRC2T, Amax3T, SRC3T, Amax4T, SRC4T
      real AmSCIT, fsrCR1T, fsrCR2T, fsrCR3T, sourcetype
      real fsrCR4T, fsrUM1T, fsrUM2T, fsrUM3T
      real fsrUM4T, c1T, c2T, ccrT, dcrT
      real FCRNT, FumRVT, FumNST, bcrT, gcrT
      real gUMT, gcrNT, gcrLT, ecrT, eumT
      real ecrvT, gammaT, S2T, S3T, S4T
      real phiS1T, phiS2T, phiS3T, phiS4T, phiT, tauT
      real mech, rockterm, ANmax, Amax, Sr, Src, fsr, Sreff, Sreffc, Snc, Sf, Smr
      real fm, gmterm, gmLterm, gN, eterm, logAn, Imf
                                                                                
      data period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,                                      
     1              0.1, 0.12, 0.14, 0.15, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 
     1              0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5,                                       
     1              3, 3.5, 4, 4.5, 5 /                     

      data Amax1 / 1.916, 1.919, 1.922, 1.925, 1.921, 1.959, 2.013, 2.049, 2.046,                                      
     1             2.066, 2.1, 2.143, 2.122, 2.092, 2.053, 1.923, 1.793, 1.793,                                       
     1             1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793,                                       
     1             1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793/ 
      data  SRC1 / 8.429, 8.09, 6.992, 6.35, 4.883, 5.043, 6.271, 7.667, 9.034,                                       
     1             11.251, 14.817, 14.817, 14.817, 14.817, 14.817, 14.817, 14.817,                                       
     1             14.817, 14.817, 14.817, 14.817, 14.817, 14.817, 14.817, 14.817,                                       
     1             14.817, 14.817, 14.817, 14.817, 14.817, 14.817, 14.817, 14.817,                                       
     1             14.817, 14.817, 14.817, 14.817/ 
      data Amax2 / 2.033, 2.027, 2.003, 1.989, 2.012, 2.017, 2.064, 2.103, 2.195,                                       
     1             2.219, 2.263, 2.329, 2.188, 2.214, 2.245, 2.324, 2.405, 2.554,                                       
     1             2.586, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718,                                       
     1             2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718,                                       
     1             2.718/ 
      data  SRC2 /1.91368,1.88256,1.77861,1.71781,
     &      2.05234,2.38713,2.83399,3.29447,3.99091,4.46576,5.04561,
     &      5.89960,5.05353,5.20490,5.38694,5.87165,6.57391,8.50000,
     &     10.67030,10.67030,10.67030,10.67030,10.67030,10.67030,
     &     10.67030,10.67030,10.67030,10.67030,10.67030,10.67030,
     &     10.67030,10.67030,10.67030,10.67030,10.67030,10.67030,
     &     10.67030/
      data Amax3 / 1.905, 1.908, 1.894, 1.886, 1.833, 1.854, 1.893, 1.924, 1.974,                                       
     1             2.032, 2.052, 2.066, 2.107, 2.14, 2.156, 2.132, 2.05, 1.925,                                       
     1             2.006, 2.179, 2.288, 2.402, 2.522, 2.648, 2.78, 2.919, 3.065,                                       
     1             3.217, 3.378, 3.547, 3.724, 3.909, 4.104, 4.309, 4.524, 4.75,                                       
     1             4.987/ 
      data  SRC3 /1.11714,1.11444,1.12437,1.13017,
     &     1.15080,1.23971,1.34819,1.45181,1.58315,1.73292,1.84134,
     &     2.03029,2.28133,2.44413,2.58017,2.74161,2.82587,2.71893,
     &     2.41759,2.30375,2.23625,2.21678,2.24338,2.80535,6.65839,
     &     30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,
     &     30.0/ 
      data Amax4 / 1.498, 1.498, 1.474, 1.46, 1.374, 1.363, 1.384, 1.425, 1.481,                                       
     1             1.525, 1.549, 1.603, 1.67, 1.706, 1.734, 1.773, 1.816, 1.843,                                       
     1             1.871, 1.878, 1.911, 1.899, 1.927, 1.987, 2.025, 2.043, 2.022,                                       
     1             1.97, 1.843, 1.729, 1.583, 1.504, 1.439, 1.391, 1.362, 1.34, 1.729/ 
      data  SRC4 /0.83644,0.83644,0.83000,0.82624,
     &     0.76758,0.78632,0.83775,0.92616,1.02228,1.11802,1.16578,
     &     1.28551,1.39808,1.44327,1.47177,1.54694,1.64401,1.79013,
     &     1.82345,1.79037,1.76844,1.67539,1.62539,1.52453,1.39724,
     &     1.32029,1.26637,1.22680,1.22065,1.31805,2.12485,14.38181,
     &     14.38181,14.38181,14.38181,14.38181,14.38181/
      data AmSCI / 1.381, 1.228, 1.087, 1.042, 1.035, 1.047, 1.071, 1.103, 1.141,                                       
     1             1.184, 1.231, 1.334, 1.448, 1.51, 1.573, 1.707, 1.833, 1.954,                                       
     1             2.034, 2.052, 2.025, 1.999, 1.975, 1.931, 1.891, 1.855, 1.822,                                       
     1             1.791, 1.724, 1.667, 1.574, 1.5, 1.439, 1.387, 1.341, 1.301, 1.265/ 
      data fsrCR1 /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     &     1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsrCR2/1.0,1.2,1.3,1.253,1.064,1.120,
     &     1.207,1.238,1.360,1.355,1.341,1.195,0.835,0.781,0.738,0.684,
     &     0.654,0.683,0.691,0.800,0.876,0.966,1.034,1.206,1.314,1.357,
     &     1.357,1.357,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsrCR3/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     &     1.0,1.080,1.093,0.948,0.908,0.862,0.745,0.623,0.436,0.400,
     &     0.433,0.460,0.496,0.529,0.578,0.578,0.578,0.578,0.000,0.0,
     &     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsrCR4 /1.0,1.0,0.949,0.550,0.477,0.492,
     &     0.531,0.613,0.693,0.780,0.816,0.998,0.954,0.942,0.927,0.927,
     &     0.949,0.970,0.961,0.948,0.965,0.958,0.984,1.055,1.114,1.175,
     &     1.216,1.230,1.192,0.942,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsrUM1 /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     &      1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsrUM2 /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     &      1.0,1.0,1.0,0.767,0.708,0.657,0.568,0.509,0.433,0.337,0.337,
     &      0.337,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0/
      data fsrUM3 /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.83,
     &      0.706,0.8,0.759,0.715,0.686,0.644,0.549,0.434,0.292,0.275,
     &      0.295,0.293,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0/
      data fsrUM4 /1.0,1.0,1.15,0.8,0.613,0.542,0.534,
     &      0.583,0.643,0.674,0.694,0.763,0.684,0.645,0.61,0.532,0.468,
     &      0.291,0.275,0.295,0.293,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &      0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data c1    /-3.224,-3.357,-3.552,-3.640,-3.758,-3.826,
     &        -3.890,-3.965,-4.055,-4.153,-4.255,-4.466,-4.677,-4.781,
     &        -4.883,-5.085,-5.233,-5.229,-5.226,-5.223,-5.221,-5.218,
     &        -5.216,-5.213,-5.210,-5.208,-5.206,-5.204,-5.200,-5.196,
     &        -5.191,-5.187,-5.183,-5.181,-5.178,-5.176,-5.174/          
      data c2     /0.900,0.909,0.927,0.937,0.944,0.948,0.956,
     &         0.967,0.980,0.995,1.009,1.040,1.070,1.085,1.100,1.129,
     &         1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,
     &         1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,
     &         1.151,1.151,1.151/
      data ccr    /1.07312,1.0785,1.07254,1.05736,1.03568,
     &    1.00578,0.98413,0.98059,0.98634,0.99121,1.00033,1.03440,
     &    1.08388,1.10602,1.12674,1.16455,1.19837,1.27000,1.32852,
     &    1.37801,1.42087,1.45868,1.49250,1.55102,1.60051,1.64337,
     &    1.68118,1.71500,1.78663,1.84515,1.93750,2.00913,2.06765,
     &    2.12357,2.13734,2.13734,2.13734/ 
      data dcr    /0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,
     &    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.19,0.178,0.162,
     &    0.148,0.136,0.125,0.101,0.083,0.053,0.02967,0.01078,0.0,0.0,
     &    0.0,0.0/
      data FcrN   /0.31282,0.31568,0.31854,0.32022,0.32141,
     &    0.32185,0.32471,0.32860,0.32984,0.33637,0.34096,0.34562,
     &    0.34596,0.34500,0.34347,0.33911,0.33355,0.31694,0.29912,
     &    0.28167,0.26521,0.24997,0.23598,0.21157,0.19133,0.17456,
     &    0.16062,0.14900,0.12761,0.11382,0.09918,0.09329,0.09106,
     &    0.09003,0.08897,0.08722,0.08444/
      data FumRV  /-0.20236,-0.21425,-0.22125,-0.22306,
     &   -0.22334,-0.22297,-0.22229,-0.22145,-0.22053,-0.21956,-0.21858,
     &   -0.21661,-0.21469,-0.21374,-0.21282,-0.21102,-0.20929,-0.20526,
     &   -0.20159,-0.19822,-0.19511,-0.19221,-0.18950,-0.18454,-0.18008,
     &   -0.17602,-0.17229,-0.16883,-0.16113,-0.15447,-0.14326,-0.13399,
     &   -0.12604,-0.11905,-0.11280,-0.10713,-0.10194/
      data FumNS  /0.25194,0.25847,0.26070,0.26025,0.25891,
     &    0.25756,0.25606,0.25456,0.25301,0.25560,0.24023,0.24775,
     &    0.24524,0.24410,0.24309,0.24098,0.23890,0.23447,0.23046,
     &    0.22674,0.22355,0.22041,0.21760,0.21250,0.20807,0.20398,
     &    0.20025,0.19671,0.18911,0.18269,0.17171,0.16289,0.15528,
     &    0.14853,0.14260,0.13724,0.13245/
      data bcr    /0.00907,0.00907,0.00907,0.00907,0.00907,
     &    0.00907,0.00907,0.00907,0.00907,0.00907,0.00907,0.00958,
     &    0.01055,0.01122,0.01170,0.01233,0.01346,0.01617,0.01831,
     &    0.01980,0.02078,0.02138,0.02168,0.02161,0.02094,0.01987,
     &    0.01854,0.01704,0.01292,0.00863,0.00042,-0.00687,-0.01315,
     &   -0.01849,-0.02298,-0.02672,-0.02980/
      data gcr    /-1.26034,-1.26949,-1.28775,-1.29619,-1.25147,
     &   -1.14724,-1.09126,-1.04676,-1.01364,-0.99232,-0.98315,-0.97311,
     &   -0.98324,-0.99261,-1.00423,-1.03299,-1.06493,-1.15136,-1.23723,
     &   -1.31736,-1.38989,-1.45461,-1.51219,-1.60834,-1.68374,-1.74306,
     &   -1.78967,-1.82603,-1.88672,-1.91705,-1.93178,-1.91787,-1.89632,
     &   -1.87695,-1.86168,-1.85415,-1.85292/
      data gUM    /-1.09985,-1.10720,-1.11766,-1.11392,-1.07971,
     &   -0.98606,-0.93901,-0.90036,-0.87064,-0.85080,-0.84175,-0.82925,
     &   -0.83484,-0.84153,-0.85033,-0.87328,-0.89945,-0.97276,-1.04809,
     &   -1.12034,-1.18735,-1.24859,-1.30439,-1.40097,-1.48074,-1.54711,
     &   -1.60263,-1.64914,-1.73810,-1.79802,-1.86912,-1.90393,-1.92241,
     &   -1.93474,-1.94360,-1.95343,-1.96352/
      data gcrN   /-0.49919,-0.48684,-0.44645,-0.42175,
     &   -0.37622,-0.43576,-0.46789,-0.50260,-0.53802,-0.57301,-0.60839,
     &   -0.67012,-0.72524,-0.74977,-0.77240,-0.81285,-0.84625,-0.90575,
     &   -0.93933,-0.95403,-0.95484,-0.94513,-0.92771,-0.87661,-0.81164,
     &   -0.73894,-0.66214,-0.58347,-0.40237,-0.24053,-0.03590,0.08743,
     &    0.17408,0.24501,0.31173,0.37949,0.45204/
      data gcrL   /1.26564,1.24149,1.19894,1.18654,1.14213,
     &  1.14137,1.15749,1.18679,1.22338,1.26404,1.30526,1.39284,1.47770,
     &  1.51881,1.55884,1.63476,1.70605,1.86247,1.99162,2.09823,2.18630,
     &  2.25934,2.31973,2.41072,2.47179,2.51100,2.53418,2.54538,2.54707,
     &  2.52581,2.49751,2.47889,2.46474,2.45267,2.44271,2.43391,2.42673/
      data ecr    /-0.00794,-0.00772,-0.00756,-0.00788,-0.00863,
     &     -0.00988,-0.01075,-0.01130,-0.01164,-0.01182,-0.01184,
     &     -0.01174,-0.01142,-0.01123,-0.01101,-0.01054,-0.01007,
     &     -0.00890,-0.00784,-0.00691,-0.00610,-0.00541,-0.00482,
     &     -0.00387,-0.00315,-0.00261,-0.00220,-0.00189,-0.00139,
     &     -0.00116,-0.00109,-0.00125,-0.00144,-0.00159,-0.00171,
     &     -0.00176,-0.00177/
      data eum    /-0.010830,-0.010580,-0.010505,-0.011005,
     &     -0.011488,-0.012292,-0.012739,-0.013080,-0.013313,-0.013455,
     &     -0.013479,-0.013488,-0.013335,-0.013231,-0.013115,-0.012840,
     &     -0.012558,-0.011817,-0.011100,-0.010436,-0.009831,-0.009285,
     &     -0.008788,-0.007922,-0.007193,-0.006568,-0.006025,-0.005551,
     &     -0.004580,-0.003842,-0.002838,-0.002226,-0.001844,-0.001604,
     &     -0.001488,-0.001452,-0.001499/
      data ecrv   /-0.00628,-0.00629,-0.00634,-0.00647,
     &     -0.00681,-0.00710,-0.00724,-0.00734,-0.00741,-0.00746,
     &     -0.00749,-0.00751,-0.00748,-0.00746,-0.00743,-0.00735,
     &     -0.00725,-0.00696,-0.00661,-0.00624,-0.00586,-0.00548,
     &     -0.00511,-0.00439,-0.00374,-0.00315,-0.00261,-0.00214,
     &     -0.00119,-0.00054,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data gamma  /-9.08724,-9.05206,-8.85775,-8.62794,
     &  -8.39822,-8.15307,-8.00327,-8.04783,-8.19691,-8.35079,-8.5243,
     &  -9.05478,-9.67640,-9.96595,-10.24160,-10.75200,-11.22492,
     &  -12.26275,-13.14181,-13.90254,-14.57355,-15.17384,-15.71567,
     &  -16.66219,-17.46728,-18.16453,-18.77798,-19.32562,-20.46705,
     &  -21.38822,-22.81170,-23.90512,-24.79448,-25.58872,-26.04635,
     &  -26.37008,-26.67797/
      data S2     /0.28877,0.29991,0.29778,0.22535,0.15865,
     &  0.08257,0.05542,0.04177,0.06040,0.09304,0.15210,0.26161,0.35467,
     &  0.39410,0.42523,0.47642,0.51169,0.55177,0.55325,0.53555,0.50855,
     &  0.47913,0.45094,0.39849,0.35254,0.31707,0.28896,0.26685,0.22889,
     &  0.20778,0.18675,0.17749,0.16807,0.15988,0.15433,0.15116,0.15275/
      data S3     /0.12210,0.10321,0.11109,0.08972,0.04037,
     &  -0.05934,-0.10956,-0.13262,-0.13537,-0.12513,-0.07090,0.00261,
     &   0.11202,0.15670,0.19678,0.27414,0.33785,0.47868,0.57391,
     &   0.63824,0.67945,0.70517,0.72065,0.72820,0.71582,0.69660,
     &   0.67366,0.64985,0.59394,0.55120,0.50247,0.48879,0.48245,
     &   0.47612,0.46978,0.46344,0.45711/
      data S4     /0.20813,0.21925,0.21713,0.15948,0.07029,
     &   -0.03529,-0.09017,-0.09583,-0.07321,-0.03629,0.01110,0.10298,
     &    0.21434,0.25574,0.29811,0.36360,0.42568,0.55643,0.65844,
     &    0.73536,0.79334,0.83853,0.87505,0.92621,0.95627,0.97619,
     &    0.98819,0.99485,0.99639,0.98814,0.96004,0.92903,0.89351,
     &    0.85603,0.81676,0.77379,0.72802/
      data phiS1  /0.4153,0.4142,0.4123,0.4089,0.4072,
     &      0.4059,0.4066,0.4171,0.4255,0.4315,0.4371,0.4515,0.4559,
     &      0.4594,0.4647,0.4721,0.4722,0.4780,0.4906,0.4918,0.4951,
     &      0.4927,0.4890,0.4792,0.4665,0.4600,0.4541,0.4475,0.4404,
     &      0.4284,0.4148,0.4003,0.4033,0.4026,0.3804,0.3736,0.3831/
      data phiS2  /0.4363,0.4356,0.4344,0.4307,0.4275,
     &      0.4294,0.4304,0.4329,0.4361,0.4417,0.4478,0.4476,0.4553,
     &      0.4547,0.4607,0.4697,0.4800,0.4882,0.4991,0.5043,0.5074,
     &      0.5074,0.5149,0.5177,0.5044,0.4954,0.4834,0.4720,0.4595,
     &      0.4499,0.4374,0.4343,0.4328,0.4214,0.4264,0.4290,0.4378/
      data phiS3  /0.4235,0.4230,0.4226,0.4207,0.4208,
     &      0.4168,0.4156,0.4215,0.4220,0.4161,0.4241,0.4332,0.4347,
     &      0.4380,0.4419,0.4514,0.4493,0.4626,0.4833,0.4783,0.4720,
     &      0.4944,0.4949,0.4848,0.4816,0.4751,0.4554,0.4415,0.4380,
     &      0.4219,0.4172,0.3955,0.3995,0.4131,0.4060,0.3921,0.3779/
      data phiS4  /0.4426,0.4419,0.4418,0.4415,0.4435,
     &      0.4446,0.4415,0.4407,0.4351,0.4374,0.4402,0.4594,0.4697,
     &      0.4707,0.4744,0.4784,0.4825,0.4843,0.4930,0.4914,0.4783,
     &      0.4787,0.4764,0.4854,0.4707,0.4837,0.4729,0.4709,0.4627,
     &      0.4564,0.4323,0.4128,0.4170,0.4265,0.4259,0.4256,0.4195/
      data phi /0.556,0.556,0.555,0.553,0.558,0.564,0.577,
     &     0.599,0.616,0.629,0.641,0.657,0.663,0.666,0.671,0.680,0.692,
     &     0.694,0.688,0.675,0.667,0.665,0.664,0.669,0.670,0.674,0.673,
     &     0.670,0.660,0.656,0.631,0.605,0.591,0.578,0.556,0.542,0.538/
      data tau /0.391,0.390,0.396,0.408,0.438,0.460,
     &     0.481,0.488,0.481,0.477,0.460,0.437,0.414,0.402,0.384,0.380,
     &     0.359,0.340,0.344,0.353,0.363,0.359,0.361,0.362,0.370,0.378,
     &     0.386,0.390,0.399,0.397,0.381,0.376,0.363,0.363,0.376,0.377,
     &     0.395/
c Set attenuation name                                                            
c     Sourcetype = 0 Crustal
c     Sourcetype = 1 Upper mantle
c     Sclass = 0 Hard Rock
c     Sclass = 1 SC I
c     Sclass = 2 SC II
c     Sclass = 3 SC III
c     Sclass = 4 SC IV
                                                                       
C Find the requested spectral period and corresponding coefficients
      nper = 37

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
       period1 = period(1)
       Amax1T =  Amax1(1) 
       SRC1T =   SRC1(1)  
       Amax2T =  Amax2(1) 
       SRC2T =   SRC2(1)  
       Amax3T =  Amax3(1) 
       SRC3T =   SRC3(1)  
       Amax4T =  Amax4(1) 
       SRC4T =   SRC4(1)  
       AmSCIT =  AmSCI(1) 
       fsrCR1T = fsrCR1(1)
       fsrCR2T = fsrCR2(1)
       fsrCR3T = fsrCR3(1)
       fsrCR4T = fsrCR4(1)
       fsrUM1T = fsrUM1(1)
       fsrUM2T = fsrUM2(1)
       fsrUM3T = fsrUM3(1)
       fsrUM4T = fsrUM4(1)
       c1T =     c1(1)    
       c2T =     c2(1)    
       ccrT =    ccr(1)   
       dcrT =    dcr(1)   
       FcrNT =   FcrN(1)  
       FumRVT =  FumRV(1) 
       FumNST =  FumNS(1) 
       bcrT =    bcr(1)   
       gcrT =    gcr(1)   
       gUMT =    gUM(1)   
       gcrNT =   gcrN(1)  
       gcrLT =   gcrL(1)  
       ecrT =    ecr(1)   
       eumT =    eum(1)   
       ecrvT =   ecrv(1)  
       gammaT =  gamma(1) 
       S2T =     S2(1)    
       S3T =     S3(1)    
       S4T =     S4(1)    
       phiS1T =  phiS1(1) 
       phiS2T =  phiS2(1) 
       phiS3T =  phiS3(1) 
       phiS4T =  phiS4(1) 
       phiT =  phi(1) 
       tauT =  tau(1) 
       goto 1011
      elseif (specT .ne. 0.0) then

C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1010 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Zhao et al. (2016) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1010    call S24_interp (period(count1),period(count2),Amax1(count1),Amax1(count2),
     +                specT,Amax1T,iflag)
         call S24_interp (period(count1),period(count2),SRC1(count1),SRC1(count2),
     +                specT,SRC1T,iflag)
         call S24_interp (period(count1),period(count2),Amax2(count1),Amax2(count2),
     +                specT,Amax2T,iflag)
         call S24_interp (period(count1),period(count2),SRC2(count1),SRC2(count2),
     +                specT,SRC2T,iflag)
         call S24_interp (period(count1),period(count2),Amax3(count1),Amax3(count2),
     +                specT,Amax3T,iflag)
         call S24_interp (period(count1),period(count2),SRC3(count1),SRC3(count2),
     +                specT,SRC3T,iflag)
         call S24_interp (period(count1),period(count2),Amax4(count1),Amax4(count2),
     +                specT,Amax4T,iflag)
         call S24_interp (period(count1),period(count2),SRC4(count1),SRC4(count2),
     +                specT,SRC4T,iflag)
         call S24_interp (period(count1),period(count2),AmSCI(count1),AmSCI(count2),
     +                specT,AmSCIT,iflag)
         call S24_interp (period(count1),period(count2),fsrCR1(count1),fsrCR1(count2),
     +                specT,fsrCR1T,iflag)
         call S24_interp (period(count1),period(count2),fsrCR2(count1),fsrCR2(count2),
     +                specT,fsrCR2T,iflag)
         call S24_interp (period(count1),period(count2),fsrCR3(count1),fsrCR3(count2),
     +                specT,fsrCR3T,iflag)
         call S24_interp (period(count1),period(count2),fsrCR4(count1),fsrCR4(count2),
     +                specT,fsrCR4T,iflag)
         call S24_interp (period(count1),period(count2),fsrUM1(count1),fsrUM1(count2),
     +                specT,fsrUM1T,iflag)
         call S24_interp (period(count1),period(count2),fsrUM2(count1),fsrUM2(count2),
     +                specT,fsrUM2T,iflag)
         call S24_interp (period(count1),period(count2),fsrUM3(count1),fsrUM3(count2),
     +                specT,fsrUM3T,iflag)
         call S24_interp (period(count1),period(count2),fsrUM4(count1),fsrUM4(count2),
     +                specT,fsrUM4T,iflag)
         call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                specT,c1T,iflag)
         call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +                specT,c2T,iflag)
         call S24_interp (period(count1),period(count2),ccr(count1),ccr(count2),
     +                specT,ccrT,iflag)
         call S24_interp (period(count1),period(count2),dcr(count1),dcr(count2),
     +                specT,dcrT,iflag)
         call S24_interp (period(count1),period(count2),FCRN(count1),FCRN(count2),
     +                specT,FcrNT,iflag)
         call S24_interp (period(count1),period(count2),FumRV(count1),FumRV(count2),
     +                specT,FumRVT,iflag)
         call S24_interp (period(count1),period(count2),FumNS(count1),FumNS(count2),
     +                specT,FumNST,iflag)
         call S24_interp (period(count1),period(count2),bcr(count1),bcr(count2),
     +                specT,bcrT,iflag)
         call S24_interp (period(count1),period(count2),gcr(count1),gcr(count2),
     +                specT,gcrT,iflag)
         call S24_interp (period(count1),period(count2),gUM(count1),gUM(count2),
     +                specT,gUMT,iflag)
         call S24_interp (period(count1),period(count2),gcrN(count1),gcrN(count2),
     +                specT,gcrNT,iflag)
         call S24_interp (period(count1),period(count2),gcrL(count1),gcrL(count2),
     +                specT,gcrLT,iflag)
         call S24_interp (period(count1),period(count2),ecr(count1),ecr(count2),
     +                specT,ecrT,iflag)
         call S24_interp (period(count1),period(count2),eum(count1),eum(count2),
     +                specT,eumT,iflag)
         call S24_interp (period(count1),period(count2),ecrv(count1),ecrv(count2),
     +                specT,ecrvT,iflag)
         call S24_interp (period(count1),period(count2),gamma(count1),gamma(count2),
     +                specT,gammaT,iflag)
         call S24_interp (period(count1),period(count2),S2(count1),S2(count2),
     +                specT,S2T,iflag)
         call S24_interp (period(count1),period(count2),S3(count1),S3(count2),
     +                specT,S3T,iflag)
         call S24_interp (period(count1),period(count2),S4(count1),S4(count2),
     +                specT,S4T,iflag)
         call S24_interp (period(count1),period(count2),phiS1(count1),phiS1(count2),
     +                specT,phiS1T,iflag)
         call S24_interp (period(count1),period(count2),phiS2(count1),phiS2(count2),
     +                specT,phiS2T,iflag)
         call S24_interp (period(count1),period(count2),phiS3(count1),phiS3(count2),
     +                specT,phiS3T,iflag)
         call S24_interp (period(count1),period(count2),phiS4(count1),phiS4(count2),
     +                specT,phiS4T,iflag)
         call S24_interp (period(count1),period(count2),tau(count1),tau(count2),
     +                specT,tauT,iflag)
         call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +                specT,phiT,iflag)

 1011 period1 = specT


C     Set mechanism term and source and near fault term.                

      xcro=2

      if (sourcetype .eq. 0.0 ) then
         if (ftype .lt. 0) then
            mech = FcrNT
         else
            mech = 0.0
         endif
         if (m .le. 7.1) then
           fm = bcrT*depth + mech + CcrT*m
         else
           fm = bcrT*depth + mech + ccrT*7.1 + dcrT*(m-7.1)
         endif
   
         if (dist .le. 30) then
           gN = gcrNT*alog(xcro+dist + exp(c1T + 6.5*c2T))
         else
           gN = gcrNT*alog(xcro+30 + exp(c1T + 6.5*c2T))
         endif

      elseif (sourcetype .eq. 1.0) then         
        if (ftype .gt. 0) then
          mech = FumRVT
         else
          mech = FumNST
        endif
         
         if (m .le. 7.1) then
           fm = mech + ccrT*m
         else
           fm = mech + ccrT*7.1 + dcrT*(m-7.1)
         endif
      endif
   
C    Set distance for geometric spreading term r
   
      Cmax = 7.1
      if (m .le. Cmax) then
       Cm = m
      else
       Cm = Cmax
      endif
      
      r = xcro + dist + exp(c1T + c2T*Cm)  
c     write(*,*) "r0 = ", r     

C    Set geometric attenuation rate term 
      if (sourcetype .eq. 0.0 ) then
      gmterm = gcrT*alog(r)
   
      elseif (sourcetype .eq. 1.0) then     
      gmterm = gumT*alog(r)
   
      endif

C    Set large distance geometric attenuation rate term 

      gmLterm = gcrLT*alog(dist+200)
   
C    Set anelastic attenuation rate term
C    Ignore volcanic path term
      if (sourcetype .eq. 0.0 ) then
      eterm = ecrT*dist
   
      elseif (sourcetype .eq. 1.0) then     
      eterm = eumT*dist
   
      endif   
   
C    Set Site class term 
   
      if (sourcetype .eq. 0.0 ) then

         if (sclass .eq. 1.0) then
            Amax = Amax1T
            Src = Src1T
         fsr = fsrCR1T
            Imf = 0.91   
            ANmax = AmSCIT 
         elseif (sclass .eq. 2.0) then           
            Amax = Amax2T   
            Src = Src2T   
         fsr = fsrCR2T
            Imf = 1.023   
            ANmax = AmSCIT * exp(S2T)
         elseif (sclass .eq. 3.0) then           
            Amax = Amax3T   
            Src = Src3T   
         fsr = fsrCR3T
            Imf = 1.034  
            ANmax = AmSCIT * exp(S3T)
         elseif (sclass .eq. 4.0) then           
            Amax = Amax4T   
            Src = Src4T   
         fsr = fsrCR4T
            Imf = 0.737  
            ANmax = AmSCIT * exp(S4T)
         endif

      elseif ( sourcetype .eq. 1.0 ) then

         if (sclass .eq. 1.0) then
            Amax = Amax1T
            Src = Src1T
         fsr = fsrUM1T
            Imf = 0.91   
            ANmax = AmSCIT 
         elseif (sclass .eq. 2.0) then           
            Amax = Amax2T   
            Src = Src2T   
         fsr = fsrUM2T
            Imf = 1.023   
            ANmax = AmSCIT * exp(S2T)
         elseif (sclass .eq. 3.0) then           
            Amax = Amax3T   
            Src = Src3T   
         fsr = fsrUM3T
            Imf = 1.034  
            ANmax = AmSCIT * exp(S3T)
         elseif (sclass .eq. 4.0) then           
            Amax = Amax4T   
            Src = Src4T   
         fsr = fsrUM4T
            Imf = 0.737  
            ANmax = AmSCIT * exp(S4T)
         endif
      endif
   
      rockterm = fm + gmterm + gmLterm + gN + eterm + gammaT - alog(AmSCIT)
      Sr = exp(rockterm) 
       
      Sreff = Sr *Imf
      Sreffc = Src *Imf
      Sf = ANmax/Amax
      
        if (ANmax .lt. 1.25) then
          ca=alog(Amax)/(alog(0.6)-alog(Sreffc**2 + 0.6))
          cb=-ca*alog(Sreffc**2+0.6)
          Snc=exp((ca*(2-1.0)*alog(0.6)*alog(10.0*0.6)-
     &        alog(10.0)*(cb+alog(Sf)))/(ca*(2*alog(10.0*0.6)-
     &        alog(0.6))))
        else
         Snc = (exp((alog(ANmax)*alog(Sreffc**2+0.6)-alog(Sf)* alog(0.6))/(alog(Amax)))-0.6)**(0.5)
        endif
           
        Smr = Sreff * (Snc/Sreffc) * fsr

        if (Smr .ne. 0.0) then
  
         logAn = alog(ANmax)-alog(Amax)*((alog(Smr**2+0.6)-alog(0.6))/(alog(Sreffc**2+0.6)-alog(0.6)))
        else
          logAn = alog(ANmax)  
        endif

      if (sclass .eq. 0.0) then
         lnY =  rockterm 
       else   
        lnY =  rockterm + logAn                                        
      endif

c      write(*,*) "--------------------------------- "
c      write(*,*) "dist = ", dist     
c      write(*,*) "sclass = ", sclass
c
c      write(*,*) "Sr = ", Sr
c      write(*,*) "Src = ", Src
c
c      write(*,*) "fm = ", fm
c      write(*,*) "gmterm = ", gmterm
c      write(*,*) "gmLterm = ", gmLterm
c      write(*,*) "gN = ", gN
c      write(*,*) "eterm = ", eterm
c      write(*,*) "gammaT = ", gammaT
c      write(*,*) "logAn = ", logAn
c      write(*,*) "Y(g) = ", exp(lnY)

C     Convert ground motion to units of gals.
      lnY = lnY + 6.89

c       write(*,*) "Y(gal) = ", exp(lnY)
                                                                               
c     Set standard error                                                       
      if (sclass .eq. 1.0) then
         phiT = phiS1T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 2.0) then
         phiT = phiS2T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 3.0) then
          phiT = phiS3T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 4.0) then
         phiT = phiS4T
         sigma = SQRT(phiT**2 + tauT**2)
      endif

      return                                                                    
      end       
          
c -------------------------------------------------------------------           
C **** Zhao et al. 2016 (BSSA, Vol 106, No.4) *************
c -------------------------------------------------------------------           

      subroutine S02_Zhaoetal2016_int ( m, dist, lnY, sigma, sclass, specT,              
     1                   attenName, period1, iflag, depth, phiT, tauT )                                   

c     This  subroutine calculates the spectral acceleration from the 
c     Zhao et al. (2016) attenuation model which is defined for subduction
c     interface events. Four separate site classes are also modeled.
                                                                         
      parameter (MAXPER=37)                                                     
      real dist, m, lnY, sigma, specT, sclass, cm, cmax, r
      real period(MAXPER), Amax1(MAXPER), SRC1(MAXPER), Amax2(MAXPER)
      real SRC2(MAXPER), Amax3(MAXPER), SRC3(MAXPER), Amax4(MAXPER), SRC4(MAXPER)
      real AmSCI(MAXPER), fsr1(MAXPER), fsr2(MAXPER), fsr3(MAXPER), fsr4(MAXPER)
      real c1(MAXPER), c2(MAXPER), gamma(MAXPER)
      real CintS(MAXPER), CintD(MAXPER), dint(MAXPER), bint(MAXPER), gint(MAXPER)
      real gintLS(MAXPER), gintLD(MAXPER), eintS(MAXPER), gammaS(MAXPER)
      real S2(MAXPER), S3(MAXPER), S4(MAXPER), S5(MAXPER), S6(MAXPER), S7(MAXPER)
      real phiS1(MAXPER), phiS2(MAXPER), phiS3(MAXPER), phiS4(MAXPER)
      real phi(MAXPER), tau(MAXPER)
      character*80 attenName                                                    
      integer nper, count1, count2, iflag, xinto
      real Amax1T, SRC1T, Amax2T, SRC2T, Amax3T, SRC3T, Amax4T, SRC4T
      real AmSCIT, fsr1T, fsr2T, fsr3T, fsr4T
      real c1T, c2T
      real CintST, CintDT, dintT, bintT, gintT, gintLST, gintLDT, eintST
      real gammaT, gammaST, S2T, S3T, S4T, S5T, S6T, S7T
      real phiS1T, phiS2T, phiS3T, phiS4T, phiT, tauT
      real mech, rockterm, ANmax, Amax, Sr, Src, fsr, Sreff, Sreffc, Snc, Sf, Smr
      real fm, gmterm, gmLterm, eterm, logAn, Imf
                                                                                
      data period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,                                      
     1              0.1, 0.12, 0.14, 0.15, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 
     1              0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5,                                       
     1              3, 3.5, 4, 4.5, 5 /                     

      data Amax1 / 1.916, 1.919, 1.922, 1.925, 1.921, 1.959, 2.013, 2.049, 2.046,                                      
     1             2.066, 2.1, 2.143, 2.122, 2.092, 2.053, 1.923, 1.793, 1.793,                                       
     1             1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793,                                       
     1             1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793/ 
      data  SRC1 /8.42900,8.09039,6.99211,6.34965,
     &   4.88339,5.04257,6.27132,7.66694,9.03436,11.25066,14.81739,
     &  14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,
     &  14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,
     &  14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,
     &  14.81739,14.81739,14.81739,14.81739,14.81739/ 
      data Amax2 / 2.033, 2.027, 2.003, 1.989, 2.012, 2.017, 2.064, 2.103, 2.195,                                       
     1             2.219, 2.263, 2.329, 2.188, 2.214, 2.245, 2.324, 2.405, 2.554,                                       
     1             2.586, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718,                                       
     1             2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718,                                       
     1             2.718/ 
      data  SRC2 /1.91368,1.88256,1.77861,1.71781,
     &  2.05234,2.38713,2.83399,3.29447,3.99091,4.46576,5.04561,5.89960,
     &  5.05353,5.20490,5.38694,5.87165,6.57391,8.50000,10.67030,
     & 10.67030,10.67030,10.67030,10.67030,10.67030,10.67030,10.67030,
     & 10.67030,10.67030,10.67030,10.67030,10.67030,10.67030,10.67030,
     & 10.67030,10.67030,10.67030,10.67030/
      data Amax3 / 1.905, 1.908, 1.894, 1.886, 1.833, 1.854, 1.893, 1.924, 1.974,                                       
     1             2.032, 2.052, 2.066, 2.107, 2.14, 2.156, 2.132, 2.05, 1.925,                                       
     1             2.006, 2.179, 2.288, 2.402, 2.522, 2.648, 2.78, 2.919, 3.065,                                       
     1             3.217, 3.378, 3.547, 3.724, 3.909, 4.104, 4.309, 4.524, 4.75,                                       
     1             4.987/ 
      data  SRC3 /1.11714,1.11444,1.12437,1.13017,
     &  1.15080,1.23971,1.34819,1.45181,1.58315,1.73292,1.84134,2.03029,
     %  2.28133,2.44413,2.58017,2.74161,2.82587,2.71893,2.41759,2.30375,
     &  2.23625,2.21678,2.24338,2.80535,6.65839,30.0,30.0,30.0,30.0,
     &  30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0/
      data Amax4 / 1.498, 1.498, 1.474, 1.46, 1.374, 1.363, 1.384, 1.425, 1.481,                                       
     1             1.525, 1.549, 1.603, 1.67, 1.706, 1.734, 1.773, 1.816, 1.843,                                       
     1             1.871, 1.878, 1.911, 1.899, 1.927, 1.987, 2.025, 2.043, 2.022,                                       
     1             1.97, 1.843, 1.729, 1.583, 1.504, 1.439, 1.391, 1.362, 1.34, 1.729/ 
      data  SRC4 /0.83644,0.83644,0.83000,0.82624,
     &  0.76758,0.78632,0.83775,0.92616,1.02228,1.11802,1.16578,1.28551,
     &  1.39808,1.44327,1.47177,1.54694,1.64401,1.79013,1.82345,1.79037,
     &  1.76844,1.67539,1.62539,1.52453,1.39724,1.32029,1.26637,1.22680,
     &  1.22065,1.31805,2.12485,14.38181,14.38181,14.38181,14.38181,
     &  14.38181,14.38181/
      data fsr1  /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     &        1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsr2  /1.0,1.0,1.0,1.0,0.843,0.663,0.841,
     &     1.029,1.235,1.144,1.092,0.945,0.624,0.577,0.545,0.527,0.546,
     &     0.596,0.623,0.701,0.708,0.737,0.748,0.728,0.634,0.0,0.0,0.0,
     &     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsr3   /1.165,0.944,1.012,1.1,0.959,0.889,
     &     0.946,1.006,1.065,1.093,1.077,1.036,0.895,0.860,0.822,0.743,
     &     0.663,0.487,0.447,0.473,0.487,0.511,0.536,0.540,0.477,0.0,
     &     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsr4   /1.0,1.0,1.0,1.0,0.557,0.543,0.574,
     &     0.648,0.721,0.809,1.015,0.972,0.967,0.963,0.953,0.967,1.005,
     &     1.045,1.035,1.008,1.007,0.981,0.990,1.016,1.022,1.023,0.997,
     &     0.948,0.802,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data c1     /-5.30119,-5.28844,-5.27568,-5.26822,
     &       -5.26293,-5.25882,-5.25547,-5.25263,-5.25017,-5.24801,
     &       -5.24607,-5.24271,-5.23988,-5.23861,-5.23742,-5.23525,
     &       -5.23331,-5.22921,-5.22585,-5.22302,-5.22056,-5.21839,
     &       -5.21645,-5.21310,-5.21026,-5.20781,-5.20564,-5.20370,
     &       -5.19959,-5.19624,-5.19095,-5.18684,-5.18349,-5.18065,
     &       -5.17819,-5.17602,-5.17409/
      data c2     /1.151,1.151,1.151,1.151,1.151,1.151,1.151,
     &        1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,
     &        1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,
     &        1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,1.151,
     &        1.151,1.151,1.151/
      data cintD  /1.09973,1.09848,1.09227,1.10690,1.11578,
     &        1.10234,1.08611,1.07291,1.06383,1.05863,1.05673,1.06051,
     &        1.07135,1.07856,1.08664,1.10470,1.12443,1.17689,1.22970,
     &        1.28058,1.32873,1.37394,1.41630,1.49310,1.56069,1.62055,
     &        1.67390,1.72171,1.82188,1.90081,2.01482,2.08892,2.13574,
     &        2.16254,2.17393,2.17301,2.16199/
      data cintS  /1.31479,1.31739,1.31919,1.34096,1.38051,
     &        1.43246,1.46238,1.47120,1.46432,1.44695,1.42323,1.36833,
     &        1.31558,1.29277,1.27319,1.24828,1.23715,1.22387,1.22846,
     &        1.24219,1.26077,1.28190,1.30432,1.35019,1.39523,1.43824,
     &        1.47881,1.51685,1.60148,1.67280,1.78372,1.86241,1.91711,
     &        1.95322,1.97450,1.98361,1.98255/
      data dint   /0.553,0.553,0.553,0.553,0.553,0.553,0.553,
     &        0.553,0.553,0.553,0.553,0.553,0.553,0.553,0.553,0.553,
     &        0.553,0.553,0.553,0.553,0.553,0.553,0.553,0.553,0.560,
     &        0.580,0.602,0.622,0.667,0.705,0.768,0.820,0.863,0.902,
     &        0.935,0.966,0.994/
      data gammaS /-3.89528,-3.89528,-3.89528,-3.89528,
     &         -3.89528,-3.89528,-3.89528,-3.89528,-3.89461,-3.90179,
     &         -3.90765,-3.91644,-3.92273,-3.92526,-3.92753,-3.93130,
     &         -3.93446,-3.94068,-3.94547,-3.94943,-3.95273,-3.95558,
     &         -3.95795,-3.96181,-3.96483,-3.96729,-3.96962,-3.97202,
     &         -3.97949,-3.99048,-4.02652,-4.08299,-4.15939,-4.25418,
     &         -4.36584,-4.49267,-4.63313/
      data bint  /0.01999,0.01999,0.01999,0.02066,0.02308,
     &       0.02709,0.02972,0.03207,0.03196,0.02972,0.02789,0.02470,
     &       0.02117,0.01951,0.01793,0.01505,0.01255,0.00769,0.00438,
     &       0.00215,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &       0.0,0.0,0.0,0.0,0.0/
      data gint   /-2.05587,-2.06565,-2.10231,-2.19232,
     &  -2.24640,-2.29341,-2.31172,-2.31101,-2.28778,-2.24680,-2.20410,
     &  -2.12011,-2.04337,-2.01088,-1.98301,-1.94607,-1.92702,-1.89876,
     &  -1.89141,-1.89299,-1.89531,-1.90578,-1.91467,-1.92743,-1.93452,
     &  -1.93739,-1.93725,-1.93505,-1.92471,-1.91189,-1.88859,-1.87252,
     &  -1.86345,-1.85969,-1.85953,-1.86151,-1.86451/
      data gintLD /0.54541,0.54975,0.56172,0.57887,
     &        0.49333,0.49100,0.50848,0.52750,0.54595,0.56307,0.57615,
     &        0.59258,0.60983,0.61961,0.63085,0.66196,0.69978,0.78447,
     &        0.85939,0.92338,0.98011,1.02220,1.05874,1.11803,1.16302,
     &        1.19735,1.22362,1.24367,1.27245,1.28541,1.28830,1.27727,
     &        1.26049,1.24107,1.22026,1.19859,1.17629/
      data gintLS /1.13364,1.13365,1.13365,1.13364,0.98811,
     &        0.90444,0.88766,0.90491,0.94206,0.98652,1.03553,1.13529,
     &        1.23424,1.28134,1.32659,1.41134,1.48854,1.65209,1.78127,
     &        1.88438,1.96763,2.03554,2.09142,2.17640,2.23600,2.27831,
     &        2.30847,2.32985,2.35853,2.36653,2.35536,2.33106,2.30409,
     &        2.27794,2.25371,2.23161,2.21149/
      data eintS  /-0.00628,-0.00625,-0.00616,-0.00572,
     &  -0.00532,-0.00503,-0.00528,-0.00569,-0.00619,-0.00673,-0.00718,
     &  -0.00793,-0.00853,-0.00879,-0.00902,-0.00927,-0.00942,-0.00959,
     %  -0.00952,-0.00933,-0.00911,-0.00888,-0.00866,-0.00824,-0.00787,
     &  -0.00755,-0.00726,-0.00700,-0.00644,-0.00597,-0.00518,-0.00451,
     &  -0.00393,-0.00344,-0.00302,-0.00267,-0.00240/
      data gamma  /-4.49858,-4.45894,-4.25807,-3.91800,
     &  -3.11423,-2.76035,-2.64088,-2.65622,-2.75266,-2.89920,-3.07698,
     &  -3.48283,-3.91612,-4.13481,-4.35243,-4.78030,-5.19439,-6.15803,
     &  -7.02003,-7.79153,-8.49551,-9.11351,-9.68515,-10.68946,
     & -11.54596,-12.28722,-12.93630,-13.51000,-14.69034,-15.60298,
     & -16.90010,-17.73658,-18.27135,-18.59264,-18.75474,-18.79351,
     & -18.73393/
      data AmSCI  / 1.358, 1.247, 1.149, 1.097, 1.065, 1.037, 1.038, 1.05, 1.103,  
     1              1.191, 1.278, 1.401, 1.525, 1.578, 1.626, 1.706, 1.768, 1.868,  
     1              1.917, 1.939, 1.944, 1.944, 1.943, 1.927, 1.912, 1.893, 1.874,  
     1              1.853, 1.799, 1.74, 1.619, 1.508, 1.416, 1.347, 1.303, 1.285, 1.267 / 
      data S2    /0.31288,0.30850,0.29297,0.22869,0.16316,
     &  0.12129,0.12345,0.13970,0.16385,0.20501,0.24449,0.32284,0.40120,
     &  0.43622,0.46736,0.51201,0.53926,0.58602,0.60465,0.60640,0.60277,
     &  0.58043,0.55692,0.50969,0.46501,0.42440,0.38841,0.35703,0.29673,
     &  0.25785,0.22262,0.21840,0.21595,0.21595,0.21595,0.21595,0.21595/
      data S3     /-0.00431,-0.01115,-0.02166,-0.11291,
     & -0.18870,-0.22831,-0.21920,-0.19009,-0.15499,-0.11250,-0.07508,
     &  0.01500,0.09699,0.14588,0.18793,0.25154,0.30298,0.42692,0.51620,
     &  0.56951,0.62366,0.65814,0.68671,0.71220,0.71239,0.69938,0.67996,
     &  0.65828,0.61928,0.58288,0.52620,0.48720,0.45698,0.42809,0.39524,
     &  0.35488,0.30468/
      data S4     /0.22838,0.22305,0.20892,0.13314,0.06959,
     &  0.02845,0.00070,-0.00949,-0.00485,0.03047,0.06080,0.14231,
     &  0.22696,0.25758,0.30748,0.35974,0.40313,0.50765,0.57779,0.63821,
     &  0.70323,0.75082,0.79378,0.84953,0.87979,0.89539,0.90256,0.90528,
     &  0.91793,0.92130,0.91709,0.90551,0.88668,0.85880,0.82025,0.76990,
     &  0.70705/
      data S5     /0.31288,0.30850,0.29297,0.22869,0.16316,
     &  0.12129,0.12345,0.13970,0.16385,0.20501,0.24449,0.32284,0.40120,
     &  0.43622,0.46736,0.51201,0.53926,0.58602,0.60465,0.60640,0.60277,
     &  0.58043,0.55692,0.50969,0.46501,0.42440,0.38841,0.35703,0.29673,
     &  0.25785,0.22262,0.21840,0.21595,0.21595,0.21595,0.21595,0.21595/
      data S6     /-0.00431,-0.01115,-0.02166,-0.08291,
     & -0.18870,-0.22831,-0.21920,-0.19009,-0.15499,-0.11250,-0.07508,
     &  0.01500,0.09699,0.13588,0.17293,0.24154,0.30298,0.42692,0.51620,
     &  0.57951,0.63366,0.65814,0.68671,0.71220,0.71239,0.69938,0.67996,
     &  0.65828,0.61928,0.58288,0.52620,0.48720,0.45698,0.42809,0.39524,
     &  0.35488,0.30468/
      data S7     /0.22838,0.22305,0.20892,0.16314,0.05959,
     &  0.00845,0.00070,0.00551,0.01415,0.04047,0.07080,0.14031,0.20196,
     &  0.23258,0.26248,0.31974,0.37313,0.48765,0.57779,0.64821,0.71323,
     &  0.75082,0.79378,0.84953,0.87979,0.89539,0.90256,0.90528,0.91793,
     &  0.92130,0.91709,0.90551,0.88668,0.85880,0.82025,0.76990,0.70705/
      data phiS1    / 0.388, 0.388, 0.386, 0.383, 0.386, 0.386, 0.386, 0.389, 0.397,  
     1              0.402, 0.406, 0.411, 0.413, 0.414, 0.416, 0.423, 0.426, 0.424,  
     1              0.422, 0.413, 0.413, 0.419, 0.409, 0.395, 0.395, 0.392, 0.395,  
     1              0.402, 0.408, 0.409, 0.415, 0.406, 0.397, 0.388, 0.39, 0.384,  
     1              0.394 / 
      data phiS2    / 0.421, 0.421, 0.42, 0.42, 0.42, 0.419, 0.419, 0.42, 0.424,  
     1              0.425, 0.424, 0.424, 0.439, 0.44, 0.446, 0.454, 0.469, 0.468,  
     1              0.473, 0.481, 0.47, 0.464, 0.457, 0.45, 0.442, 0.444, 0.434,  
     1              0.428, 0.425, 0.441, 0.437, 0.438, 0.44, 0.435, 0.442, 0.441,  
     1              0.461 / 
      data phiS3    / 0.393, 0.393, 0.394, 0.397, 0.394, 0.399, 0.398, 0.412, 0.406,  
     1              0.406, 0.403, 0.415, 0.407, 0.41, 0.404, 0.411, 0.414, 0.41,  
     1              0.398, 0.4, 0.393, 0.389, 0.383, 0.385, 0.385, 0.393, 0.392,  
     1              0.394, 0.402, 0.421, 0.412, 0.419, 0.399, 0.403, 0.398, 0.384,  
     1              0.375 / 
      data phiS4    / 0.458, 0.416, 0.415, 0.415, 0.421, 0.419, 0.42, 0.42, 0.419,  
     1              0.419, 0.423, 0.426, 0.438, 0.43, 0.437, 0.455, 0.456, 0.452,  
     1              0.439, 0.432, 0.423, 0.419, 0.42, 0.403, 0.412, 0.417, 0.407,  
     1              0.406, 0.419, 0.415, 0.409, 0.421, 0.414, 0.407, 0.396, 0.38,  
     1              0.371 / 
      data phi /0.553,0.554,0.553,0.555,0.565,0.570,0.583,
     &    0.602,0.614,0.625,0.637,0.646,0.654,0.659,0.663,0.672,0.678,
     &    0.659,0.640,0.634,0.627,0.620,0.612,0.613,0.625,0.628,0.628,
     &    0.633,0.636,0.644,0.635,0.619,0.599,0.581,0.568,0.551,0.563/
      data tau /0.669,0.670,0.673,0.682,0.709,0.734,
     &    0.760,0.784,0.793,0.797,0.796,0.789,0.773,0.773,0.774,0.776,
     &    0.778,0.753,0.729,0.729,0.720,0.719,0.712,0.720,0.739,0.743,
     &    0.743,0.750,0.753,0.754,0.741,0.734,0.712,0.692,0.682,0.667,
     &    0.675/
   
c Set attenuation name                                                            
c     Sourcetype = 0 Crustal
c     Sourcetype = 1 Subduction - Interface
c     Sourcetype = 2 Subduction - Slab
c     Sourcetype = 3 Upper mantle
c     Sclass = 0 Hard Rock
c     Sclass = 1 SC I
c     Sclass = 2 SC II
c     Sclass = 3 SC III
c     Sclass = 4 SC IV

                                                                 
C Find the requested spectral period and corresponding coefficients
      nper = 37

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
       period1 = period(1)
       Amax1T =  Amax1(1) 
       SRC1T =   SRC1(1)  
       Amax2T =  Amax2(1) 
       SRC2T =   SRC2(1)  
       Amax3T =  Amax3(1) 
       SRC3T =   SRC3(1)  
       Amax4T =  Amax4(1) 
       SRC4T =   SRC4(1)  
       AmSCIT =  AmSCI(1) 
       fsr1T = fsr1(1)
       fsr2T = fsr2(1)
       fsr3T = fsr3(1)
       fsr4T = fsr4(1)
       c1T =     c1(1)    
       c2T =     c2(1)    
       CintST =    CintS(1)   
       CintDT =    CintD(1)   
       dintT =    dint(1)   
       bintT =    bint(1)   
       gintT =    gint(1)   
       gintLST =   gintLS(1)  
       gintLDT =   gintLD(1)  
       eintST =    eintS(1)   
       gammaT =  gamma(1) 
       gammaST =  gammaS(1) 
       S2T =     S2(1)    
       S3T =     S3(1)    
       S4T =     S4(1)    
       S5T =     S5(1)    
       S6T =     S6(1)    
       S7T =     S7(1)    
       phiS1T =  phiS1(1) 
       phiS2T =  phiS2(1) 
       phiS3T =  phiS3(1) 
       phiS4T =  phiS4(1) 
       phiT =  phi(1) 
       tauT =  tau(1) 
       goto 1011
      elseif (specT .ne. 0.0) then

C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1010 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Zhao et al. (2016) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1010    call S24_interp (period(count1),period(count2),Amax1(count1),Amax1(count2),
     +                specT,Amax1T,iflag)
         call S24_interp (period(count1),period(count2),SRC1(count1),SRC1(count2),
     +                specT,SRC1T,iflag)
         call S24_interp (period(count1),period(count2),Amax2(count1),Amax2(count2),
     +                specT,Amax2T,iflag)
         call S24_interp (period(count1),period(count2),SRC2(count1),SRC2(count2),
     +                specT,SRC2T,iflag)
         call S24_interp (period(count1),period(count2),Amax3(count1),Amax3(count2),
     +                specT,Amax3T,iflag)
         call S24_interp (period(count1),period(count2),SRC3(count1),SRC3(count2),
     +                specT,SRC3T,iflag)
         call S24_interp (period(count1),period(count2),Amax4(count1),Amax4(count2),
     +                specT,Amax4T,iflag)
         call S24_interp (period(count1),period(count2),SRC4(count1),SRC4(count2),
     +                specT,SRC4T,iflag)
         call S24_interp (period(count1),period(count2),AmSCI(count1),AmSCI(count2),
     +                specT,AmSCIT,iflag)
         call S24_interp (period(count1),period(count2),fsr1(count1),fsr1(count2),
     +                specT,fsr1T,iflag)
         call S24_interp (period(count1),period(count2),fsr2(count1),fsr2(count2),
     +                specT,fsr2T,iflag)
         call S24_interp (period(count1),period(count2),fsr3(count1),fsr3(count2),
     +                specT,fsr3T,iflag)
         call S24_interp (period(count1),period(count2),fsr4(count1),fsr4(count2),
     +                specT,fsr4T,iflag)
         call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                specT,c1T,iflag)
         call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +                specT,c2T,iflag)
         call S24_interp (period(count1),period(count2),CintS(count1),CintS(count2), 
     +                specT,CintST,iflag)
         call S24_interp (period(count1),period(count2),CintD(count1),CintD(count2), 
     +                specT,CintDT,iflag)
         call S24_interp (period(count1),period(count2),dint(count1),dint(count2), 
     +                specT,dintT,iflag)
         call S24_interp (period(count1),period(count2),bint(count1),bint(count2), 
     +                specT,bintT,iflag)
         call S24_interp (period(count1),period(count2),gint(count1),gint(count2), 
     +                specT,gintT,iflag)
         call S24_interp (period(count1),period(count2),gintLS(count1),gintLS(count2), 
     +                specT,gintLST,iflag)
         call S24_interp (period(count1),period(count2),gintLD(count1),gintLD(count2), 
     +                specT,gintLDT,iflag)
         call S24_interp (period(count1),period(count2),eintS(count1),eintS(count2), 
     +                specT,eintST,iflag)
         call S24_interp (period(count1),period(count2),gammaS(count1),gammaS(count2), 
     +                specT,gammaST,iflag)
         call S24_interp (period(count1),period(count2),gamma(count1),gamma(count2),
     +                specT,gammaT,iflag)
         call S24_interp (period(count1),period(count2),S2(count1),S2(count2),
     +                specT,S2T,iflag)
         call S24_interp (period(count1),period(count2),S3(count1),S3(count2),
     +                specT,S3T,iflag)
         call S24_interp (period(count1),period(count2),S4(count1),S4(count2),
     +                specT,S4T,iflag)
         call S24_interp (period(count1),period(count2),S5(count1),S5(count2), 
     +                specT,S5T,iflag)
         call S24_interp (period(count1),period(count2),S6(count1),S6(count2), 
     +                specT,S6T,iflag)
         call S24_interp (period(count1),period(count2),S7(count1),S7(count2), 
     +                specT,S7T,iflag)
         call S24_interp (period(count1),period(count2),phiS1(count1),phiS1(count2),
     +                specT,phiS1T,iflag)
         call S24_interp (period(count1),period(count2),phiS2(count1),phiS2(count2),
     +                specT,phiS2T,iflag)
         call S24_interp (period(count1),period(count2),phiS3(count1),phiS3(count2),
     +                specT,phiS3T,iflag)
         call S24_interp (period(count1),period(count2),phiS4(count1),phiS4(count2),
     +                specT,phiS4T,iflag)
         call S24_interp (period(count1),period(count2),tau(count1),tau(count2),
     +                specT,tauT,iflag)
         call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +                specT,phiT,iflag)

 1011 period1 = specT


C     Set mechanism term and source and near fault term.                
      if (depth .le. 25 ) then

         if (m .le. 7.1) then
           fm = bintT*depth + gammaST + CintST*m
         else
           fm = bintT*depth + gammaST + CintST*7.1 + dintT*(m-7.1)
         endif
         
      else
         if (m .le. 7.1) then
           fm = bintT*depth + CintDT*m
         else
           fm = bintT*depth + CintDT*7.1 + dintT*(m-7.1)
         endif
          
      endif
   
C    Set distance for geometric spreading term r
   
      xinto=10
      Cmax = 7.1
      if (m .le. Cmax) then
       Cm = m
      else
       Cm = Cmax
      endif
      
      r = xinto + dist + exp(c1T + c2T*Cm)  
   
C    Set geometric attenuation rate term 
      gmterm = gintT*alog(r)

C    Set large distance geometric attenuation rate term 
      if (depth .le. 25 ) then

       gmLterm = gintLST*alog(dist+200)

      else
       
       gmLterm = gintLDT*alog(dist+200)
          
      endif
   
C    Set anelastic attenuation rate term
C    Ignore volcanic path term
      if (depth .le. 25  ) then
       eterm = eintST*dist
   
      else
       eterm = 0
   
      endif   
   
C    Set Site class term 
   
      if (depth .le. 25 ) then

         if (sclass .eq. 1) then
            Amax = Amax1T
            Src = Src1T
         fsr = fsr1T
            Imf = 0.91   
            ANmax = AmSCIT 
         elseif (sclass .eq. 2) then           
            Amax = Amax2T   
            Src = Src2T   
         fsr = fsr2T
            Imf = 1.023   
            ANmax = AmSCIT * exp(S2T)
         elseif (sclass .eq. 3) then           
            Amax = Amax3T   
            Src = Src3T   
         fsr = fsr3T
            Imf = 1.034  
            ANmax = AmSCIT * exp(S3T)
         elseif (sclass .eq. 4) then           
            Amax = Amax4T   
            Src = Src4T   
         fsr = fsr4T
            Imf = 0.737  
            ANmax = AmSCIT * exp(S4T)
         endif
      else
         if (sclass .eq. 1) then
            Amax = Amax1T
            Src = Src1T
         fsr = fsr1T
            Imf = 0.91   
            ANmax = AmSCIT 
         elseif (sclass .eq. 2) then           
            Amax = Amax2T   
            Src = Src2T   
         fsr = fsr2T
            Imf = 1.023   
            ANmax = AmSCIT * exp(S5T)
         elseif (sclass .eq. 3) then           
            Amax = Amax3T   
            Src = Src3T   
         fsr = fsr3T
            Imf = 1.034  
            ANmax = AmSCIT * exp(S6T)
         elseif (sclass .eq. 4) then           
            Amax = Amax4T   
            Src = Src4T   
         fsr = fsr4T
            Imf = 0.737  
            ANmax = AmSCIT * exp(S7T)
         endif
   
      endif   
   
      rockterm = fm + gmterm + gmLterm + eterm + gammaT - alog(AmSCIT)
      Sr = exp(rockterm) 
       
      Sreff = Sr *Imf
      Sreffc = Src *Imf
      Sf = ANmax/Amax

        if (ANmax .lt. 1.25) then
          ca=alog(Amax)/(alog(0.6)-alog(Sreffc**2 + 0.6))
          cb=-ca*alog(Sreffc**2+0.6)
          Snc=exp((ca*(2-1.0)*alog(0.6)*alog(10.0*0.6)-
     &        alog(10.0)*(cb+alog(Sf)))/(ca*(2*alog(10.0*0.6)-
     &        alog(0.6))))
        else
         Snc = (exp((alog(ANmax)*alog(Sreffc**2+0.6)-alog(Sf)* alog(0.6))/(alog(Amax)))-0.6)**(0.5)
        endif
           
      Smr = Sreff * (Snc/Sreffc) * fsr
 
        if (Smr .ne. 0.0) then
  
         logAn = alog(ANmax)-alog(Amax)*((alog(Smr**2+0.6)-alog(0.6))/(alog(Sreffc**2+0.6)-alog(0.6)))
        else
          logAn = alog(ANmax)  
        end if


      if (sclass .eq. 0) then
         lnY =  rockterm 
       else   
        lnY =  rockterm + logAn                                        
      endif                                       

c      write(*,*) "--------------------------------- "
c      write(*,*) "dist = ", dist     
c      write(*,*) "sclass = ", sclass
c
c      write(*,*) "Sr = ", Sr
c      write(*,*) "Src = ", Src
c
c      write(*,*) "fm = ", fm
c      write(*,*) "gmterm = ", gmterm
c      write(*,*) "gmLterm = ", gmLterm
c      write(*,*) "eterm = ", eterm
c      write(*,*) "gammaT = ", gammaT
c      write(*,*) "logAn = ", logAn
c      write(*,*) "Y(g) = ", exp(lnY)

C     Convert ground motion to units of gals.
      lnY = lnY + 6.89

c       write(*,*) "Y(gal) = ", exp(lnY)
                                                                               
c     Set standard error                                                       
      if (sclass .eq. 1) then
         phiT = phiS1T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 2) then
         phiT = phiS2T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 3) then
          phiT = phiS3T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 4) then
         phiT = phiS4T
         sigma = SQRT(phiT**2 + tauT**2)
      endif

      return                                                                    
      end       
          
c -------------------------------------------------------------------           
C **** Zhao et al. 2016 (BSSA, Vol 106, No.4) *************
c -------------------------------------------------------------------           

      subroutine S02_Zhaoetal2016_slab ( m, dist, lnY, sigma, sclass, specT,              
     1                   attenName, period1, iflag, depth, phiT, tauT )                                   

c     This  subroutine calculates the spectral acceleration from the 
c     Zhao et al. (2016) attenuation model which is defined for subduction
c     intraslab events. Four separate site classes are also modeled.
                                                                         
      parameter (MAXPER=37)                                                     
      real dist, m, lnY, sigma, specT, sclass, cm, cmax, r
      real period(MAXPER), Amax1(MAXPER), SRC1(MAXPER), Amax2(MAXPER)
      real SRC2(MAXPER), Amax3(MAXPER), SRC3(MAXPER), Amax4(MAXPER), SRC4(MAXPER)
      real AmSCI(MAXPER), fsr1(MAXPER), fsr2(MAXPER), fsr3(MAXPER), fsr4(MAXPER)
      real gamma(MAXPER), c1(MAXPER), cSL1(MAXPER), cSL2(MAXPER), dSL(MAXPER)
      real bSL(MAXPER), gSL(MAXPER), gSLL(MAXPER), eSL(MAXPER), eSLH(MAXPER)       
      real S2(MAXPER), S3(MAXPER), S4(MAXPER)
      real phiS1(MAXPER), phiS2(MAXPER),phiS3(MAXPER), phiS4(MAXPER)
      real phi(MAXPER), tau(MAXPER)
      character*80 attenName                                                    
      integer nper, count1, count2, iflag, xinto
      real Amax1T, SRC1T, Amax2T, SRC2T, Amax3T, SRC3T, Amax4T, SRC4T
      real AmSCIT, fsr1T, fsr2T, fsr3T, fsr4T
      real c1T, cSL1T, cSL2T, dSLT, bSLT, gSLT, gSLLT, eSLT, eSLHT
      real gammaT, S2T, S3T, S4T
      real phiS1T, phiS2T, phiS3T, phiS4T, phiT, tauT
      real mech, rockterm, ANmax, Amax, Sr, Src, fsr, Sreff, Sreffc, Snc, Sf, Smr
      real fm, gmterm, gmLterm, eterm, logAn, Imf, c2, qterm, qSLH
                                                                                
      data period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,                                      
     1              0.1, 0.12, 0.14, 0.15, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 
     1              0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5,                                       
     1              3, 3.5, 4, 4.5, 5 /                     

      data Amax1 / 1.916, 1.919, 1.922, 1.925, 1.921, 1.959, 2.013, 2.049, 2.046,                                      
     1             2.066, 2.1, 2.143, 2.122, 2.092, 2.053, 1.923, 1.793, 1.793,                                       
     1             1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793,                                       
     1             1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793, 1.793/ 
      data  SRC1 /8.42900,8.09039,6.99211,6.34965,
     &   4.88339,5.04257,6.27132,7.66694,9.03436,11.25066,14.81739,
     &  14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,
     &  14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,
     &  14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,
     &  14.81739,14.81739,14.81739,14.81739,14.81739,14.81739,
     &  14.81739,14.81739/
      data Amax2 / 2.033, 2.027, 2.003, 1.989, 2.012, 2.017, 2.064, 2.103, 2.195,                                       
     1             2.219, 2.263, 2.329, 2.188, 2.214, 2.245, 2.324, 2.405, 2.554,                                       
     1             2.586, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718,                                       
     1             2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718, 2.718,                                       
     1             2.718/ 
      data  SRC2 /1.91368,1.88256,1.77861,1.71781,
     &   2.05234,2.38713,2.83399,3.29447,3.99091,4.46576,5.04561,
     &   5.89960,5.05353,5.20490,5.38694,5.87165,6.57391,8.50000,
     &  10.6703,10.6703,10.6703,10.6703,10.6703,10.6703,10.6703,
     &  10.6703,10.6703,10.6703,10.6703,10.6703,10.6703,10.6703,
     &  10.6703,10.6703,10.6703,10.6703,10.6703/
      data Amax3 / 1.905, 1.908, 1.894, 1.886, 1.833, 1.854, 1.893, 1.924, 1.974,                                       
     1             2.032, 2.052, 2.066, 2.107, 2.14, 2.156, 2.132, 2.05, 1.925,                                       
     1             2.006, 2.179, 2.288, 2.402, 2.522, 2.648, 2.78, 2.919, 3.065,                                       
     1             3.217, 3.378, 3.547, 3.724, 3.909, 4.104, 4.309, 4.524, 4.75,                                       
     1             4.987/ 
      data  SRC3 /1.11714,1.11444,1.12437,1.13017,
     &   1.15080,1.23971,1.34819,1.45181,1.58315,1.73292,1.84134,
     &   2.03029,2.28133,2.44413,2.58017,2.74161,2.82587,2.71893,
     &   2.41759,2.30375,2.23625,2.21678,2.24338,2.80535,6.65839,
     &  30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0/
      data Amax4 / 1.498, 1.498, 1.474, 1.46, 1.374, 1.363, 1.384, 1.425, 1.481,                                       
     1             1.525, 1.549, 1.603, 1.67, 1.706, 1.734, 1.773, 1.816, 1.843,                                       
     1             1.871, 1.878, 1.911, 1.899, 1.927, 1.987, 2.025, 2.043, 2.022,                                       
     1             1.97, 1.843, 1.729, 1.583, 1.504, 1.439, 1.391, 1.362, 1.34, 1.729/ 
      data  SRC4 /0.83644,0.83644,0.83000,0.82624,
     &   0.76758,0.78632,0.83775,0.92616,1.02228,1.11802,1.16578,
     &   1.28551,1.39808,1.44327,1.47177,1.54694,1.64401,1.79013,
     &   1.82345,1.79037,1.76844,1.67539,1.62539,1.52453,1.39724,
     &   1.32029,1.26637,1.22680,1.22065,1.31805,2.12485,14.38181,
     &   14.38181,14.38181,14.38181,14.38181,14.38181/
      data AmSCI / 1.381, 1.228, 1.087, 1.042, 1.035, 1.047, 1.071, 1.103, 1.141,
     1             1.184, 1.231, 1.334, 1.448, 1.51, 1.573, 1.707, 1.833, 1.954, 
     1             2.034, 2.052, 2.025, 1.999, 1.975, 1.931, 1.891, 1.855, 1.822, 
     1             1.791, 1.724, 1.667, 1.574, 1.5, 1.439, 1.387, 1.341, 1.301, 1.265 / 
      data fsr1  /1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
     &   1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsr2  /1.0,1.0,1.0,1.0,1.006,0.851,0.803,
     &   0.918,1.062,1.106,1.071,0.9515,0.672,0.631,0.6,0.571,0.565,
     &   0.601,0.579,0.679,0.655,0.615,0.55,0.0,0.0,0.0,0.0,0.0,0.0,
     &   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data fsr3  /1.0,1.0,1.0,1.0,1.0,1.0,1.044,0.975,
     &   0.964,0.98,0.97,1.022,0.889,0.861,0.831,0.748,0.65,0.479,0.449,
     &   0.482,0.499,0.515,0.53,0.53,0.499,0.369,0.3,0.2,0.0,0.0,0.0,
     &   0.0,0.0,0.0,0.0,0.0,0.0/
      data fsr4  /1.0,1.0,1.05,0.58,0.482,0.472,0.506,
     &   0.587,0.683,0.782,0.823,1.029,0.991,0.983,0.973,0.979,1.006,
     &   1.027,1.021,1.003,1.01,0.985,0.99,1.006,1.0,1.0,0.96,0.904,
     &   0.738,0.535,0.358,0.0,0.0,0.0,0.0,0.0,0.0/
      data c1    /-5.30119,-5.28844,-5.27568,-5.26822,
     &   -5.26293,-5.25882,-5.25547,-5.25263,-5.25017,-5.24801,-5.24607,
     &   -5.24271,-5.23988,-5.23861,-5.23742,-5.23525,-5.23331,-5.22921,
     &   -5.22585,-5.22302,-5.22056,-5.21839,-5.21645,-5.21310,-5.21026,
     &   -5.20781,-5.20564,-5.20370,-5.19959,-5.19624,-5.19095,-5.18684,
     &   -5.18349,-5.18065,-5.17819,-5.17602,-5.17409/
      data cSL1  /1.44758,1.45400,1.46625,1.49246,1.50129,
     &    1.51051,1.51380,1.51111,1.50406,1.49423,1.48300,1.45559,
     &    1.44277,1.43314,1.43253,1.43710,1.44781,1.48260,1.51881,
     &    1.55291,1.58443,1.61360,1.64075,1.69020,1.73450,1.77474,
     &    1.81162,1.84561,1.92015,1.98274,2.08214,2.15841,2.22046,
     &    2.27406,2.32307,2.37009,2.37009/
      data cSL2  /0.37625,0.38099,0.39101,0.41976,0.45746,
     &    0.48601,0.50311,0.50704,0.50004,0.48071,0.45759,0.41355,
     &    0.37828,0.36308,0.34919,0.32464,0.30358,0.26174,0.23036,
     &    0.20580,0.18597,0.16960,0.15585,0.13405,0.11757,0.10476,
     &    0.09458,0.08636,0.07173,0.06258,0.05327,0.05036,0.04536,
     &    0.04536,0.04536,0.04536,0.04536/
      data dSL   /0.42646,0.42075,0.40055,0.36433,0.32072,
     &    0.30000,0.31147,0.32673,0.34289,0.35921,0.37000,0.40606,
     &    0.43450,0.45000,0.46055,0.48439,0.509,0.555,0.593,0.625,
     &    0.652,0.675,0.695,0.729,0.756,0.778,0.796,0.812,0.841,
     &    0.861,0.884,0.9,0.9,0.9,0.9,0.9,0.9/
      data bSL   /0.01826,0.01826,0.01826,0.01826,0.01826,
     &    0.01826,0.01826,0.01826,0.01826,0.01826,0.01826,0.01826,
     &    0.01826,0.01826,0.01826,0.01826,0.01826,0.01826,0.01826,
     &    0.01826,0.01826,0.01826,0.01826,0.01826,0.01826,0.01826,
     &    0.01826,0.01826,0.01808,0.01786,0.01718,0.01628,0.01549,
     &    0.01489,0.01458,0.01459,0.01459/
      data gSL   /-1.98471,-1.96360,-1.91839,-1.89271,
     &  -1.87260,-1.85351,-1.83395,-1.81345,-1.79189,-1.76931,-1.74581,
     &  -1.73746,-1.74463,-1.74972,-1.76259,-1.78989,-1.82110,-1.90412,
     &  -1.98439,-2.05756,-2.12282,-2.18047,-2.23118,-2.31475,-2.37885,
     &  -2.42769,-2.46450,-2.49170,-2.52758,-2.53359,-2.49565,-2.42623,
     &  -2.34726,-2.27002,-2.19947,-2.12528,-2.02646/
      data gSLL  /1.12071,1.03278,0.94715,0.93420,0.97168,
     &    1.01492,1.06854,1.13401,1.20364,1.25808,1.30112,1.39137,
     &    1.47084,1.50784,1.54326,1.60985,1.67146,1.80738,1.92242,
     &    2.02102,2.10642,2.18097,2.24651,2.35602,2.44331,2.51391,
     &    2.57166,2.61931,2.70638,2.76244,2.82205,2.84475,2.84988,
     &    2.84667,2.83992,2.82802,2.82521/
      data eSL   /-0.00340,-0.00331,-0.00345,-0.00391,
     &  -0.00454,-0.00510,-0.00552,-0.00588,-0.00615,-0.00635,-0.00652,
     &  -0.00660,-0.00652,-0.00647,-0.00636,-0.00614,-0.00590,-0.00526,
     &  -0.00468,-0.00415,-0.00369,-0.00327,-0.00290,-0.00227,-0.00178,
     &  -0.00139,-0.00109,-0.00086,-0.00052,-0.00043,-0.00070,-0.00127,
     &  -0.00198,-0.00271,-0.00341,-0.00421,-0.00500/
      data eSLH  /-0.00050,-0.00050,-0.00050,-0.00050,
     &  -0.00050,-0.00050,-0.00050,-0.00049,-0.00048,-0.00048,-0.00048,
     &  -0.00049,-0.00051,-0.00052,-0.00053,-0.00056,-0.00059,-0.00067,
     &  -0.00075,-0.00083,-0.00091,-0.00099,-0.00107,-0.00124,-0.00139,
     &  -0.00154,-0.00166,-0.00178,-0.00199,-0.00213,-0.00225,-0.00219,
     &  -0.00207,-0.00193,-0.00180,-0.00170,-0.00158/
      data gamma /-9.87956,-9.51269,-9.26626,-9.33150,
     &  -9.50798,-9.72858,-9.96628,-10.22583,-10.55111,-10.80721,
     & -11.02190,-11.36530,-11.73039,-11.88013,-12.05637,-12.42044,
     & -12.78542,-13.63537,-14.38086,-15.03511,-15.61599,-16.13830,
     & -16.61324,-17.45298,-18.18095,-18.82499,-19.40313,-19.92766,
     & -21.05818,-21.99633,-23.48839,-24.64741,-25.59713,-26.40997,
     & -27.13181,-27.79299,-28.31346/
      data S2    /0.23200,0.22886,0.21825,0.18737,
     &   0.12332,0.07207,0.02701,-0.00621,0.01565,0.05089,0.09560,
     &   0.20037,0.30372,0.34284,0.37400,0.42700,0.46297,0.50856,
     &   0.50776,0.49710,0.48065,0.46159,0.44224,0.40537,0.37342,
     &   0.34623,0.32364,0.30479,0.27026,0.24831,0.22529,0.21540,
     &   0.21154,0.20976,0.20875,0.20774,0.20672/
      data S3   /0.14371,0.13978,0.12600,0.06164,
     &  -0.01705,-0.06331,-0.10103,-0.14684,-0.14479,-0.12666,-0.09319,
     &  -0.00878,0.08926,0.13602,0.17751,0.25309,0.32005,0.45304,
     &   0.54875,0.61713,0.66634,0.70111,0.72558,0.75294,0.76245,
     &   0.76119,0.75384,0.74279,0.70833,0.67256,0.61067,0.56403,
     &   0.52612,0.49766,0.47685,0.46223,0.45267/
      data S4    /0.14704,0.13284,0.14431,0.06600,
     &  -0.01714,-0.07312,-0.11955,-0.16006,-0.12434,-0.07293,-0.01458,
     &   0.08252,0.17151,0.20932,0.24117,0.29896,0.34591,0.44231,
     &   0.51782,0.57596,0.62239,0.65976,0.69066,0.73796,0.77226,
     &   0.79736,0.81616,0.83009,0.85036,0.85732,0.84991,0.82757,
     &   0.79911,0.76782,0.73594,0.70407,0.67220/
      data phiS1 /  0.398,0.397,0.395,0.389,0.387,0.387,
     &   0.397,0.403,0.413,0.422,0.429,0.440,0.441,0.450,0.452,0.454,
     &   0.462,0.474,0.472,0.468,0.457,0.450,0.445,0.448,0.440,0.441,
     &   0.435,0.427,0.413,0.416,0.409,0.398,0.390,0.386,0.377,0.360,
     &   0.361/
      data phiS2   /0.417,0.417,0.417,0.418,0.420,0.422,
     &   0.416,0.412,0.412,0.412,0.418,0.429,0.435,0.440,0.444,0.447,
     &   0.450,0.470,0.475,0.478,0.484,0.477,0.470,0.460,0.460,0.460,
     &   0.456,0.448,0.442,0.448,0.437,0.429,0.423,0.410,0.410,0.412,
     &   0.447/
      data phiS3   /0.409,0.409,0.408,0.409,0.413,0.409,
     &   0.401,0.394,0.389,0.391,0.394,0.428,0.433,0.420,0.424,0.446,
     &   0.441,0.459,0.437,0.444,0.453,0.477,0.472,0.459,0.461,0.457,
     &   0.449,0.442,0.428,0.418,0.406,0.387,0.369,0.376,0.362,0.368,
     &   0.381/
      data phiS4   /0.415,0.415,0.416,0.417,0.420,0.422,
     &   0.423,0.420,0.423,0.427,0.426,0.445,0.442,0.440,0.437,0.436,
     &   0.432,0.432,0.432,0.432,0.416,0.408,0.409,0.407,0.403,0.407,
     &   0.408,0.408,0.412,0.418,0.413,0.414,0.412,0.402,0.394,0.385,
     &   0.380/
      Data phi /0.587,0.587,0.587,0.588,0.600,0.607,0.623,
     &   0.638,0.651,0.663,0.674,0.690,0.692,0.696,0.697,0.704,0.713,
     &   0.711,0.683,0.665,0.657,0.647,0.640,0.633,0.633,0.636,0.636,
     %   0.637,0.635,0.645,0.633,0.608,0.582,0.562,0.540,0.525,0.522/
      Data tau /0.457,0.459,0.465,0.480,0.520,
     &   0.555,0.584,0.598,0.598,0.585,0.567,0.534,0.504,0.486,0.465,
     &   0.430,0.406,0.385,0.365,0.373,0.383,0.391,0.403,0.412,0.432,
     &   0.436,0.437,0.436,0.444,0.448,0.424,0.413,0.407,0.394,0.381,
     &   0.365,0.378/
c Set attenuation name                                                            
c     Sourcetype = 0 Crustal
c     Sourcetype = 1 Subduction - Interface
c     Sourcetype = 2 Subduction - Slab
c     Sourcetype = 3 Upper mantle
c     Sclass = 0 Hard Rock
c     Sclass = 1 SC I
c     Sclass = 2 SC II
c     Sclass = 3 SC III
c     Sclass = 4 SC IV

                                                       
C Find the requested spectral period and corresponding coefficients
      nper = 37

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
       period1 = period(1)
       Amax1T =  Amax1(1) 
       SRC1T =   SRC1(1)  
       Amax2T =  Amax2(1) 
       SRC2T =   SRC2(1)  
       Amax3T =  Amax3(1) 
       SRC3T =   SRC3(1)  
       Amax4T =  Amax4(1) 
       SRC4T =   SRC4(1)  
       AmSCIT =  AmSCI(1) 
       fsr1T = fsr1(1)
       fsr2T = fsr2(1)
       fsr3T = fsr3(1)
       fsr4T = fsr4(1)
       c1T   =  c1(1)
       cSL1T =  cSL1(1)
       cSL2T =  cSL2(1)
       dSLT  =  dSL(1)
       bSLT  =  bSL(1)
       gSLT  =  gSL(1)
       gSLLT =  gSLL(1)
       eSLT  =  eSL(1)
       eSLHT =  eSLH(1)
       gammaT =  gamma(1) 
       S2T =     S2(1)    
       S3T =     S3(1)    
       S4T =     S4(1)    
       phiS1T =  phiS1(1) 
       phiS2T =  phiS2(1) 
       phiS3T =  phiS3(1) 
       phiS4T =  phiS4(1) 
       phiT =  phi(1) 
       tauT =  tau(1) 
       goto 1011
      elseif (specT .ne. 0.0) then

C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1010 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Zhao et al. (2016) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1010    call S24_interp (period(count1),period(count2),Amax1(count1),Amax1(count2),
     +                specT,Amax1T,iflag)
         call S24_interp (period(count1),period(count2),SRC1(count1),SRC1(count2),
     +                specT,SRC1T,iflag)
         call S24_interp (period(count1),period(count2),Amax2(count1),Amax2(count2),
     +                specT,Amax2T,iflag)
         call S24_interp (period(count1),period(count2),SRC2(count1),SRC2(count2),
     +                specT,SRC2T,iflag)
         call S24_interp (period(count1),period(count2),Amax3(count1),Amax3(count2),
     +                specT,Amax3T,iflag)
         call S24_interp (period(count1),period(count2),SRC3(count1),SRC3(count2),
     +                specT,SRC3T,iflag)
         call S24_interp (period(count1),period(count2),Amax4(count1),Amax4(count2),
     +                specT,Amax4T,iflag)
         call S24_interp (period(count1),period(count2),SRC4(count1),SRC4(count2),
     +                specT,SRC4T,iflag)
         call S24_interp (period(count1),period(count2),AmSCI(count1),AmSCI(count2),
     +                specT,AmSCIT,iflag)
         call S24_interp (period(count1),period(count2),fsr1(count1),fsr1(count2),
     +                specT,fsr1T,iflag)
         call S24_interp (period(count1),period(count2),fsr2(count1),fsr2(count2),
     +                specT,fsr2T,iflag)
         call S24_interp (period(count1),period(count2),fsr3(count1),fsr3(count2),
     +                specT,fsr3T,iflag)
         call S24_interp (period(count1),period(count2),fsr4(count1),fsr4(count2),
     +                specT,fsr4T,iflag)
         call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                specT,c1T,iflag)
         call S24_interp (period(count1),period(count2),cSL1(count1),cSL1(count2),
     +                specT,cSL1T,iflag)
         call S24_interp (period(count1),period(count2),cSL2(count1),cSL2(count2),
     +                specT,cSL2T,iflag)
         call S24_interp (period(count1),period(count2),dSL(count1),dSL(count2),
     +                specT,dSLT,iflag)
         call S24_interp (period(count1),period(count2),bSL(count1),bSL(count2),
     +                specT,bSLT,iflag)
         call S24_interp (period(count1),period(count2),gSL(count1),gSL(count2),
     +                specT,gSLT,iflag)
         call S24_interp (period(count1),period(count2),gSLL(count1),gSLL(count2),
     +                specT,gSLLT,iflag)
         call S24_interp (period(count1),period(count2),eSL(count1),eSL(count2),
     +                specT,eSLT,iflag)
         call S24_interp (period(count1),period(count2),eSLH(count1),eSLH(count2),
     +                specT,eSLHT,iflag)
         call S24_interp (period(count1),period(count2),gamma(count1),gamma(count2),
     +                specT,gammaT,iflag)
         call S24_interp (period(count1),period(count2),S2(count1),S2(count2),
     +                specT,S2T,iflag)
         call S24_interp (period(count1),period(count2),S3(count1),S3(count2),
     +                specT,S3T,iflag)
         call S24_interp (period(count1),period(count2),S4(count1),S4(count2),
     +                specT,S4T,iflag)
         call S24_interp (period(count1),period(count2),phiS1(count1),phiS1(count2),
     +                specT,phiS1T,iflag)
         call S24_interp (period(count1),period(count2),phiS2(count1),phiS2(count2),
     +                specT,phiS2T,iflag)
         call S24_interp (period(count1),period(count2),phiS3(count1),phiS3(count2),
     +                specT,phiS3T,iflag)
         call S24_interp (period(count1),period(count2),phiS4(count1),phiS4(count2),
     +                specT,phiS4T,iflag)
         call S24_interp (period(count1),period(count2),tau(count1),tau(count2),
     +                specT,tauT,iflag)
         call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +                specT,phiT,iflag)

 1011 period1 = specT

      if (depth .gt. 100.0) then 
          depth = 100.0
      endif  

C     Set mechanism term and source and near fault term.                
 
         if (m .le. 7.1) then
           fm = bSLT*depth + cSL1T*m + cSL2T*(m-6.3)**2
         else
           fm = bSLT*depth + cSL1T*7.1 + cSL2T*(7.1-6.3)**2 + dSLT*(m-7.1)
         endif
   
C    Set distance for geometric spreading term r
   
      Cmax = 7.1
      c2 = 1.151
      if (m .le. Cmax) then
       Cm = m
      else
       Cm = Cmax
      endif
      
      r = dist + exp(c1T + c2*Cm)  
   
C    Set geometric attenuation rate term 
      gmterm = gSLT*alog(r)

C    Set large distance geometric attenuation rate term 
  
      gmLterm = gSLLT*alog(dist+200)
          
   
C    Set anelastic attenuation rate term
C    Ignore volcanic path term
      if (depth .lt. 50  ) then
   
       qSLH = eSLHT*0
    
      else
       qSLH = eSLHT*(0.02*depth -1)
   
      endif   
      
      qterm = qSLH * dist

      eterm = eSLT * dist

   
C    Set Site class term 

         if (sclass .eq. 1) then
            Amax = Amax1T
            Src = Src1T
         fsr = fsr1T
            Imf = 0.91   
            ANmax = AmSCIT 
         elseif (sclass .eq. 2) then           
            Amax = Amax2T   
            Src = Src2T   
         fsr = fsr2T
            Imf = 1.023   
            ANmax = AmSCIT * exp(S2T)
         elseif (sclass .eq. 3) then           
            Amax = Amax3T   
            Src = Src3T   
         fsr = fsr3T
            Imf = 1.034  
            ANmax = AmSCIT * exp(S3T)
         elseif (sclass .eq. 4) then           
            Amax = Amax4T   
            Src = Src4T   
         fsr = fsr4T
            Imf = 0.737  
            ANmax = AmSCIT * exp(S4T)
         endif
   
      rockterm = fm + gmterm + gmLterm + eterm + gammaT - alog(AmSCIT)
      Sr = exp(rockterm) 
       
      Sreff = Sr *Imf
      Sreffc = Src *Imf
      Sf = ANmax/Amax
        if (ANmax .lt. 1.25) then
          ca=alog(Amax)/(alog(0.6)-alog(Sreffc**2 + 0.6))
          cb=-ca*alog(Sreffc**2+0.6)
          Snc=exp((ca*(2-1.0)*alog(0.6)*alog(10.0*0.6)-
     &        alog(10.0)*(cb+alog(Sf)))/(ca*(2*alog(10.0*0.6)-
     &        alog(0.6))))
        else
         Snc = (exp((alog(ANmax)*alog(Sreffc**2+0.6)-alog(Sf)* alog(0.6))/(alog(Amax)))-0.6)**(0.5)
        endif
           
      Smr = Sreff * (Snc/Sreffc) * fsr
  
        if (Smr .ne. 0.0) then
  
         logAn = alog(ANmax)-alog(Amax)*((alog(Smr**2+0.6)-alog(0.6))/(alog(Sreffc**2+0.6)-alog(0.6)))
        else
          logAn = alog(ANmax)  
        end if

      if (sclass .eq. 0) then
         lnY =  rockterm 
       else   
        lnY =  rockterm + logAn                                        
      endif
      
c      write(*,*) "--------------------------------- "
c      write(*,*) "dist = ", dist     
c      write(*,*) "sclass = ", sclass
c
c      write(*,*) "Sr = ", Sr
c      write(*,*) "Src = ", Src
c
c      write(*,*) "fm = ", fm
c      write(*,*) "gmterm = ", gmterm
c      write(*,*) "gmLterm = ", gmLterm
c      write(*,*) "qterm = ", qterm
c      write(*,*) "eterm = ", eterm
c      write(*,*) "gammaT = ", gammaT
c      write(*,*) "logAn = ", logAn
c      write(*,*) "Y(g) = ", exp(lnY)

C     Convert ground motion to units of gals.
      lnY = lnY + 6.89

c       write(*,*) "Y(gal) = ", exp(lnY)
                                                                               
c     Set standard error                                                       
      if (sclass .eq. 1) then
         phiT = phiS1T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 2) then
         phiT = phiS2T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 3) then
          phiT = phiS3T
         sigma = SQRT(phiT**2 + tauT**2)
      elseif (sclass .eq. 4) then
         phiT = phiS4T
         sigma = SQRT(phiT**2 + tauT**2)
      endif

      return                                                                    
      end       

c ------------------------------------------------------------------            
C *** Montalva2017 Subduction (2017 - Model) Horizontal ***********
c ------------------------------------------------------------------            
      subroutine S02_Montalva2017 ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth, disthypo )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT, depth1
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C

      
c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)

C     Compute Rock PGA
      call S02_Montalva2017_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, faba, depth, disthypo )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S02_Montalva2017_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
c ----------------------------------------------------------------------
      subroutine S02_Montalva2017_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo )

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=23)
      real a1(MAXPER), a2(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER), a16(MAXPER),
     1     a3(MAXPER), a4(MAXPER), a5(MAXPER), a9(MAXPER), dC1int(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo,
     1     mag
      real a1T, a2T, a6T, a7T, a8T, a3T, a4T, a5T, a9T, dC1intT
      real a10T, a11T, a12T, a13T, a14T, a15T, a16T, sigsT, sigtT
      real vLinT, b_soilT, sumgm, Ftype, tau, period1
      integer count1, count2, iflag
      real n, c, c4, c1, deltac1, faba, R, testmag, VsStar, depth, specT
      real base, fmag, fdepth, fsite, farc, dC1slab

      data period /  0.00, 0.02, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
     1               0.6,  0.75, 1.00,  1.5, 2.00, 2.5, 3.00, 4.0, 5.0, 
     2               6.0, 7.5, 10.00  /
  
      data vLin   / 865.1, 865.1, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 654.3, 587.1, 503,
     1             456.6, 430.3, 410.5, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400 / 
      data b_soil / -1.186, -1.186, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381, -2.518,
     1             -2.657, -2.669, -2.599, -2.401, -1.955, -1.025, -0.299, 0, 0, 0, 0, 0, 
     1            0, 0 / 
      data dC1int / 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.143682921, 0.1,  
     1            0.073696559, 0.04150375, 0, -0.05849625, -0.1, -0.155033971,  
     1            -0.2, -0.2, -0.2, -0.2, -0.2, -0.2 /
      data a1 / 5.87504394, 5.97631438, 7.45297044, 8.04759521, 7.76085108, 6.171919, 
     1            4.83403302, 4.42687615, 4.57008643, 3.98311294, 4.8603434, 4.67510367, 
     1            4.30862113, 3.57339281, 2.92216459, 2.39779653, 1.64147667, 1.66482796, 
     1            0.90564754, 0.6123444, 0.32672294, -0.24139803, -0.96313983 / 
      data a2 / -1.75359772, -1.77010766, -2.03336398, -2.10610081, -1.99370934, -1.58654201, 
     1            -1.2971103, -1.18774055, -1.24895678, -1.13377346, -1.38019755, -1.35362409, 
     1            -1.30799859, -1.23082022, -1.18750273, -1.16319283, -1.06543862, -1.12677535, 
     1            -1.07619985, -1.13079589, -1.1573438, -1.1407007, -1.09295336 / 
      data a3 / 0.13125248, 0.12246057, 0.08332151, 0.08012671, 0.0730312, 0.05481839, 0.05249728, 
     1            0.02995137, 0.03865827, 0.04682762, 0.03822425, 0.02523729, 0.00995253, 0.03605351, 
     1            0.02768934, 0.040113, 0.08310064, 0.09403648, 0.13838017, 0.15259121, 0.1242091, 
     1            0.10950824, 0.11343926 / 
      data a4 / 0.80276784, 0.84131709, 1.03131243, 1.03436999, 1.07565004, 1.17061492, 1.20531288, 
     1            1.37607187, 1.34990775, 1.3795388, 1.51949871, 1.66662746, 1.85625091, 1.81217177, 
     1            2.03469107, 2.04340485, 1.88987024, 1.9050392, 1.71178342, 1.59358719, 1.69183532, 
     1            1.71125604, 1.67160339 / 
      data a5 / -0.33486952, -0.28054559, -0.03954116, -0.01295063, 0.00758131, 0.10490549, 0.17968066, 
     1            0.22912175, 0.15592549, 0.11670946, 0.18347677, 0.21967977, 0.29782648, 0.24372341, 
     1            0.22521403, 0.27382886, 0.18739875, 0.13268085, 0.01379686, 0.06464958, 0.32368231, 
     1            0.60252124, 0.7762083 / 
      data a6 / -0.00039095, -0.00038903, 0, -0.00009638, -0.00078515, -0.00267532, -0.0033759, -0.00355237, 
     1            -0.00244847, -0.00207613, -0.00001896, 0, 0, 0, -0.00009996, -0.00033356, -0.00121364, 
     1            -0.00087595, -0.00061861, 0, 0, 0, 0 / 
      data a7 / 1.0988, 1.0988, 1.2536, 1.4175, 1.3997, 1.3582, 1.1648, 0.994, 0.8821, 0.7046, 0.5799, 
     1            0.5021, 0.3687, 0.1746, -0.082, -0.2821, -0.4108, -0.4466, -0.4344, -0.4368, -0.4586, 
     1            -0.4433, -0.4828 / 
      data a8 / -1.42, -1.42, -1.65, -1.8, -1.8, -1.69, -1.49, -1.3, -1.18, -0.98, -0.82, -0.7, -0.54, 
     1            -0.34, -0.05, 0.12, 0.25, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3 / 
      data a9 / 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
     1            0.4, 0.4, 0.4, 0.4, 0.4 / 
      data a10 / 4.53143081, 4.57416129, 4.56070915, 4.36639286, 3.90922953, 3.06236311, 3.50112817, 
     1            3.62815675, 3.87633808, 4.03388062, 4.31418239, 4.75196667, 4.70451938, 4.56020155, 
     1            4.83342978, 4.59028522, 4.13415056, 4.18978319, 4.50906779, 4.56385964, 4.55836575, 
     1            5.08281865, 5.49692364 / 
      data a11 / 0.0056735, 0.00565448, 0.00848068, 0.00921589, 0.00629627, 0.00558843, 0.00319554, 
     1            0.001817, 0.00212947, 0.00068979, 0.0006478, 0.0008707, -0.00031282, -0.00101097, 
     1            0.00009741, 0.00108512, 0.00035459, 0.0007295, 0.00084112, 0.00068188, 0.00137322, 
     1            0.00167053, -0.00070392 / 
      data a12 / 1.01494528, 1.03738201, 1.31034079, 1.48158019, 1.65618649, 1.93944484, 2.08901131, 
     1            2.25003086, 2.283387, 2.3140873, 2.33333479, 2.23421777, 2.05217228, 1.63506217, 
     1            0.69338467, -0.09761879, -0.34931995, -0.33269783, -0.41320697, -0.42395126, -0.38759507, 
     1            -0.32638288, -0.25811162 / 
      data a13 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 / 
      data a14 / -0.73080261, -0.73868917, -0.69848828, -0.65335577, -0.5505116, -0.42997222, -0.53087673, 
     1            -0.58085678, -0.66280655, -0.72244113, -0.79644275, -0.90120145, -0.89829099, -0.87330858, 
     1            -0.94685865, -0.90845421, -0.80518214, -0.81689247, -0.87331394, -0.87800447, -0.88436295, 
     1            -0.98803311, -1.05008478 / 
      data a15 / 0.9969, 0.9969, 1.103, 1.2732, 1.3042, 1.26, 1.223, 1.16, 1.05, 0.8, 0.662, 0.58, 0.48, 
     1            0.33, 0.31, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3 / 
      data a16 / -1, -1, -1.18, -1.36, -1.36, -1.3, -1.25, -1.17, -1.06, -0.78, -0.62, -0.5, -0.34, -0.14, 
     1            0, 0, 0, 0, 0, 0, 0, 0, 0 / 
      data sigs   / 0.6911808, 0.69938258, 0.70173433, 0.71412373, 0.741128, 0.74606525, 0.7451527, 
     1            0.72855743, 0.72093248, 0.71005053, 0.66934213, 0.66733247, 0.66329494, 0.63504015, 
     1            0.60012607, 0.56961713, 0.55384735, 0.53658882, 0.51345287, 0.51417184, 0.49080507, 
     1            0.4706381, 0.46023151 / 
      data sigt   / 0.47462209, 0.47631913, 0.53776165, 0.56188074, 0.52707475, 0.50642417, 0.44618739, 
     1            0.45040229, 0.42549471, 0.42945015, 0.43333698, 0.44599448, 0.46723155, 0.50143305, 
     1            0.51633193, 0.50688464, 0.51465398, 0.50365207, 0.45311429, 0.43900131, 0.4208419, 
     1            0.41701232, 0.38872242 / 

C Constant parameters            
      n = 1.18
      c = 1.88
      dC1slab = -0.3
      c4 = 10.0
      c1 = 7.2
 
C Find the requested spectral period and corresponding coefficients
      nPer = 23

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a3T = a3(i1)
         a4T = a4(i1)
         a5T = a5(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a9T = a9(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         a15T = a15(i1)
         a16T = a16(i1)
         dC1intT = dC1int(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Montalva Sub (2017 Model) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +                   specT,a5T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +                   specT,a9T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,a16T,iflag)
            call S24_interp (period(count1),period(count2),dC1int(count1),dC1int(count2),
     +                   specT,dC1intT,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)

 1011 period1 = specT                                                                                                              

C     Compute the R term and base model based on either Rupture Distance 
c         (Interface events) of Hypocentral distance (Intraslab events). 
      if (ftype .eq. 0) then
      deltaC1 = dC1intT
         R = rRup + c4*exp( (mag-6.0)*a9T ) 
         base = a1T + a4T*deltaC1 + (a2T + a14T*ftype + a3T*(mag - 7.2))*alog(R) + a6T*rRup + a10T*ftype
      elseif (ftype .eq. 1) then
      deltac1 = dC1slab
         R = disthypo + c4*exp( (mag-6.0)*a9T ) 
         base = a1T + a4T*deltaC1 + (a2T + a14T*ftype + a3T*(mag - 7.2))*alog(R) + a6T*disthypo + a10T*ftype
      else
         write (*,*) 'BC Hydro V3 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif
      
C     Base model for Magnitude scaling.      
      testmag = (c1 + deltaC1)
      if (mag .le. testmag ) then
         fmag = a4T*(mag-testmag) + a13T*(10.0-mag)**2.0
      else
         fmag = a5T*(mag-testmag) + a13T*(10.0-mag)**2.0
      endif      
      
C     Depth Scaling
        fdepth =  a11T*(min(depth, 120.0) -60.0 )*ftype

C     Forearc/Backarc scaling      
      if (ftype .eq. 1) then
         farc =  (a7T +a8T*alog(max(disthypo,85.0)/40.0))*faba
      elseif (ftype .eq. 0) then   
         farc =  (a15T +a16T*alog(max(rRup,100.0)/40.0))*faba
      endif 

C     Site Response 
      if (vs30 .gt. 1000.0) then
          VsStar = 1000.0
      else
          VsStar = vs30
      endif
       
      if (vs30 .ge. VlinT) then
         fsite = a12T*alog(VsStar/vLinT) + b_soilT*n*alog(VsStar/vLinT)
      else
         fsite = a12T*alog(VsStar/vLinT) - b_soilT*alog(pgarock + c) +
     1          b_soilT*alog(pgarock + c*(VsStar/vlinT)**n)     
      endif

      sumgm = base + fmag + fdepth + farc + fsite

c      write(*,*) "deltaC1 = ", deltaC1
c      write(*,*) "testmag = ", testmag
c      write(*,*) "R = ", R
c      write(*,*) "base = ", base
c      write(*,*) "fmag = ", fmag
c      write(*,*) "fdepth = ", fdepth
c      write(*,*) "f_faba = ", farc
c      write(*,*) "fsite = ", fsite
c      write(*,*) "lnYSa = ", sumgm
c      write(*,*) "Sa = ", exp(sumgm)
   

C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = sumgm

      return
      end

C====================================================================================================
      subroutine S02_Arroyo2010 ( mag, rRup, lnSa, specT, iflag, phi1, tau1, sigma)

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=57)
      real a1(MAXPER), a2(MAXPER), a3(MAXPER), a4(MAXPER)
      real period(MAXPER), sig(MAXPER), tau(MAXPER), phi(MAXPER)
      real lnSa, rRup, sigma, mag
      real a1T, a2T, a3T, a4T
      real sigT, tauT, phiT, tau1, phi1, a4r,a4Rr0
      real period1, E1a, E1b, r0square
      integer count1, count2, iflag
      real specT

      data period / 0, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 
     1           0.09, 0.095, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 
     1           0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
     1           0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 
     1           2.5, 3, 3.5, 4, 4.5, 5 / 
      data a1 / 2.4862, 3.8123, 4.044, 4.1429, 4.3092, 4.377, 4.5185, 4.4591, 4.5939, 4.4832, 4.5062, 
     1          4.4648, 4.394, 4.3391, 4.0505, 3.5599, 3.1311, 2.7012, 2.5485, 2.2699, 1.913, 1.7181, 
     1          1.4039, 1.108, 1.0652, 0.8319, 0.4965, 0.3173, 0.2735, 0.099, -0.0379, -0.3512, -0.6897,  
     1          -0.6673, -0.7154, -0.7015, -0.8581, -0.9712, -1.097, -1.2346, -1.26, -1.7687, -2.1339,  
     1          -2.4122, -2.5442, -2.8509, -3.0887, -3.4884, -3.7195, -4.0141, -4.1908, -5.1104, -5.5926,  
     1          -6.1202, -6.5318, -6.9744, -7.1389 / 
      data a2 / 0.9392, 0.8636, 0.8489, 0.858, 0.8424, 0.8458, 0.8273, 0.8394, 0.8313, 0.8541, 0.8481,  
     1          0.8536, 0.858, 0.862, 0.8933, 0.9379, 0.9736, 1.003, 0.9988, 1.0125, 1.045, 1.0418,  
     1          1.0782, 1.1038, 1.0868, 1.1088, 1.1408, 1.1388, 1.1533, 1.1662, 1.2206, 1.2445, 1.2522,  
     1          1.2995, 1.3263, 1.2994, 1.3205, 1.3375, 1.3532, 1.3687, 1.3652, 1.4146, 1.4417, 1.4577,  
     1          1.4618, 1.492, 1.5157, 1.575, 1.5966, 1.6162, 1.6314, 1.7269, 1.7515, 1.8077, 1.8353,  
     1          1.8685, 1.8721 / 
      data a3 / 0.5061, 0.5578, 0.5645, 0.5725, 0.5765, 0.5798, 0.5796, 0.5762, 0.5804, 0.5792, 0.5771,  
     1          0.5742, 0.5712, 0.5666, 0.5546, 0.535, 0.5175, 0.4985, 0.485, 0.471, 0.4591, 0.445, 0.4391,  
     1          0.4287, 0.4208, 0.4142, 0.4044, 0.393, 0.4067, 0.4127, 0.4523, 0.4493, 0.4421, 0.4785,  
     1          0.5068, 0.5056, 0.5103, 0.5201, 0.5278, 0.5345, 0.5426, 0.5342, 0.5263, 0.5201, 0.5242,  
     1          0.522, 0.5215, 0.5261, 0.5255, 0.5187, 0.5199, 0.5277, 0.5298, 0.5402, 0.5394, 0.5328, 0.5376 / 
      data a4 / 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,  
     1          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,  
     1          0.015, 0.015, 0.0134, 0.0117, 0.0084, 0.0076, 0.0067, 0.0051, 0.0034, 0.0029, 0.0023, 0.0018,  
     1          0.0012, 0.0007, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,  
     1          0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001 / 
      data sig / 0.75, 0.8228, 0.8429, 0.8512, 0.8583, 0.8591, 0.8452, 0.8423, 0.8473, 0.8421, 0.8344, 0.8304,  
     1          0.8294, 0.8254, 0.796, 0.7828, 0.7845, 0.7717, 0.7551, 0.7431, 0.7369, 0.7264, 0.7209, 0.7198,  
     1          0.7206, 0.7264, 0.7255, 0.7292, 0.7272, 0.7216, 0.7189, 0.7095, 0.7084, 0.7065, 0.707, 0.7092,  
     1          0.6974, 0.6906, 0.6923, 0.6863, 0.6798, 0.6701, 0.6697, 0.6801, 0.6763, 0.6765, 0.6674, 0.648,  
     1          0.6327, 0.6231, 0.6078, 0.6001, 0.6029, 0.6137, 0.6201, 0.6419, 0.6701 / 
      data tau / 0.4654, 0.5179, 0.5246, 0.5199, 0.5253, 0.5563, 0.527, 0.5241, 0.5205, 0.5148, 0.5115, 0.5273,  
     1          0.5309, 0.5116, 0.4768, 0.465, 0.4523, 0.4427, 0.4428, 0.4229, 0.4223, 0.4356, 0.4191, 0.4281,  
     1          0.4384, 0.425, 0.4348, 0.4419, 0.4574, 0.4249, 0.4265, 0.4215, 0.4304, 0.4096, 0.3999, 0.4113,  
     1          0.3923, 0.4047, 0.398, 0.3921, 0.3842, 0.3871, 0.3931, 0.3939, 0.4146, 0.4159, 0.4187, 0.4164,  
     1          0.3985, 0.4062, 0.3828, 0.3936, 0.4148, 0.4273, 0.4393, 0.4577, 0.5011 / 
      data phi / 0.5882, 0.6394, 0.6597, 0.674, 0.6788, 0.6547, 0.6607, 0.6594, 0.6685, 0.6664, 0.6593, 0.6415,  
     1          0.6373, 0.6477, 0.6374, 0.6298, 0.6409, 0.6321, 0.6116, 0.6109, 0.6039, 0.5814, 0.5865, 0.5787,  
     1          0.5719, 0.5891, 0.5808, 0.58, 0.5653, 0.5833, 0.5788, 0.5707, 0.5627, 0.5756, 0.583, 0.5778,  
     1          0.5766, 0.5596, 0.5665, 0.5632, 0.5608, 0.5471, 0.5422, 0.5544, 0.5343, 0.5335, 0.5197, 0.4965,  
     1          0.4914, 0.4726, 0.4721, 0.453, 0.4375, 0.4405, 0.4376, 0.45, 0.4449 / 

 
C Find the requested spectral period and corresponding coefficients
      nPer = 57

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a3T = a3(i1)
         a4T = a4(i1)
         sigT = sig(i1)
         phiT = phi(i1)
         tauT = tau(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Arroyo Sub-Interface (2010 Model) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),sig(count1),sig(count2),
     +                   specT,sigT,iflag)
            call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +                   specT,phiT,iflag)
            call S24_interp (period(count1),period(count2),tau(count1),tau(count2),
     +                   specT,tauT,iflag)

 1011 period1 = specT                                                                                                              
 
   
      r0square = 1.4447*0.00001*exp(2.3026*mag)
      
      a4R = a4T*rRup
      a4Rr0 = a4T*sqrt(rRup**2+r0square)
   
      call S02_fEnx(a4R,E1a,1)
      call S02_fEnx(a4Rr0,E1b,1)
   
      lnSa = a1T + a2T*mag + a3T*alog((E1a-E1b)/r0square) 

c   write(*,*) "r0square = ", r0square
c   write(*,*) "a4R = ", a4R
c  write(*,*) "a4Rr0 = ", a4Rr0
c   write(*,*) "E1a = ", E1a
c   write(*,*) "E1b = ", E1b
c   write(*,*) "lnYSa = ", lnSa
c   write(*,*) "Sa = ", exp(lnSa)
c   write(*,*) "Sa = ", exp(lnSa)/980

C     Set sigma values to return
      sigma = sigT
      tau1 = tauT
      phi1 = phiT

      return
      end

      subroutine S02_fEnx( x, Enx, n)
      
      Integer MaxIt, i, ii, nm1 As 
      real expint, Eps, FpMin, Euler, Enx
      real a, b, c, d, del, Fact, H, psi 
      
c      'PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)
       MaxIt = 100
       Eps = 0.0000001
       FpMin = 1E-30
       Euler = 0.5772156649
      
       nm1 = n - 1
      
      If (n .le. 0 .or. x .le. 0 ) Then
      
          write(*,*) "fEnx Input Error"
      goto 99
    
      elseif ( x .eq. 0 .And. n .eq.0 ) then
          write(*,*) "fEnx Input Error"
      goto 99

      elseif ( x .eq. 0 .And.  n .eq. 1) then
            write(*,*) "fEnx Input Error"
      goto 99
     
      End If
      
      If ( n .eq. 0) Then
          expint = Exp(-x) / x
      Else
          
          If (x .eq. 0) Then
              expint = 1.0/ nm1
          Else
      
              If (x .ge. 1) Then
                  b = x + n
                  c = 1.0 / FpMin
                  d = 1.0 / b
                  H = d
              
                  Do i = 1 , MaxIt
                      a = -i * (nm1 + i)
                      b = b + 2.0
                      d = 1.0 / (a * d + b)
                      c = b + a / c
                      del = c * d
                      H = H * del
                
                      If (Abs(del - 1.0) <= Eps) Then
                          expint = H * Exp(-x)
                          Enx = expint
                          goto 99
                      End If
                   enddo
              
                  write(*,*) "fEnx Input Error"
                  goto 99
          
              Else
                  
                  If (nm1 .ne. 0) Then
                      expint = 1.0 / nm1
                  Else
                      expint = -Log(x) - Euler
                  End If
                      
                  Fact = 1.0
              
                  Do i = 1, MaxIt
                      Fact = -Fact * x / i
                          
                      If (i .ne. nm1) Then
                          del = -Fact / (i - nm1)
                      Else
                          psi = -Euler
                              
                          DO ii = 1, nm1
                              psi = psi + 1.0 / ii
                          enddo
                              
                          del = Fact * (-Log(x) + psi)
                      End If
                          
                      expint = expint + del
                          
                      If (Abs(del) .le. Abs(expint) * Eps) Then
                          Enx = expint
                          goto 99
                      End If
                  enddo
              
                  write(*,*)  "fEnx Input Error"
                  goto 99
                  
              End If
              
          End If
      
      End If
      
      
      
      Enx = expint
      
      
99    continue
      return
      End 

c ------------------------------------------------------------------            
C *** Chao2017 (Crustal and Subduction - Model) Horizontal ***********
c ------------------------------------------------------------------            
      subroutine S04_Chao2017 ( mag, dist, ftype, lnY, sigma, specT, vs, Ztor, Z10,           
     1            vs30_class, attenName, period2, iflag, sourcetype, phi, tau, msasflag )         

      implicit none

      real mag, dip, fType, dist, vs, SA1100,
     1      Z10,  ZTOR, fltWidth, lnSa, sigma, lnY, vs30_rock, sourcetype
      real Fn, Frv, specT, period2, CRjb, phi, tau, z10_rock, SA_rock
      integer hwflag, iflag, vs30_class, regionflag, msasflag
      character*80 attenName                                                    

C     Set the reference spectrum.                
c     sourcetype = 0 for crustal
c                  1 for Subduction 
c     Vs30_class = 0 for estimated
c     Vs30_class = 1 for measured 
C     Mainshock and Aftershocks included based on MSASFlag
C         0 = Mainshocks
C         1 = Aftershocks

c     Compute SA1100
      vs30_rock = 1100.
      z10_rock = 0.0058959
      SA_rock = 0.
      
         call S04_Chaoetal2017 ( mag, dist, ftype, sigma, specT, vs30_rock, Ztor, z10_rock,
     1             SA_rock, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag)
      Sa1100 = exp(lnSa)

c     Compute Sa at spectral period for given Vs30

         call S04_Chaoetal2017 ( mag, dist, ftype, sigma, specT, vs, Ztor, Z10,
     1             sa1100, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag )

C     Convert ground motion to units of gals.
      lnY = lnSa + 6.89

      period2 = specT

      return
      end
c -------------------------------------------------------------------           
C **** Chao et al. 2017 (SSHAC model) *************
c -------------------------------------------------------------------           

      subroutine S04_Chaoetal2017 ( mag, dist, ftype, sigma, specT, vs, Ztor, Z10,           
     1            sa1100, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag )                                   

      implicit none
      
      integer MAXPER                                                                            
      parameter (MAXPER=21)                                                     
      real ftype, dist, mag, lnSa, sigma, specT, lnYref, vs, Ztor, Z10, period1
      real period(MAXPER), c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER), c11(MAXPER)
      real c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER), c16(MAXPER), c17(MAXPER)
      real c18(MAXPER), c19(MAXPER), c20(MAXPER), c21(MAXPER), c22(MAXPER), c23(MAXPER)
      real c24(MAXPER), c25(MAXPER), c26(MAXPER), c27(MAXPER), taucr1(MAXPER), taucr2(MAXPER)
      real tausb1(MAXPER), tausb2(MAXPER), phisscr1(MAXPER), phisscr2(MAXPER), phisssb1(MAXPER)
      real phisssb2(MAXPER), arfacr(MAXPER), arfasb(MAXPER), phis2s(MAXPER)
      character*80 attenName                                                    
      integer nper, count1, count2, C11flag, C20flag, C26flag, iflag, C10flag
      integer vs30_class, h, n, i, msasflag
      integer Fcr, Fsb, Fcrss, Fcrno, Fcrro, Fsbintra, Fsbinter, Fas, Fkuo17, Fks17, Frf, Fmanila 
      real Mc, Mref, Mmax, Rrupref, Vs30ref, Zref, sourcetype
      real c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, c10T, c11T, c12T, c13T, c14T, c15T
      real c16T, c17T, c18T, c19T, c20T, c21T, c22T, c23T, c24T, c25T, c26T, c27T
      real taucr1T, taucr2T, tausb1T, tausb2T, phisscr1T, phisscr2T, phisssb1T, phisssb2T
      real arfacrT, arfasbT, phis2sT, phi, tau, fm, SA1100, Z10ref
      real Ssource, Spath, Ssite, Ssitelin, Ssitenon, Sztor, Smag, Sgeom, Sanel
      real taucr, tausb, phisscr, phisssb, phiss, sigmass
                                                                                
      data period  / 0, -2.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,  
     &          0.5, 0.75, 1, 1.5, 2, 3, 4, 5 / 
      data c1 / -0.42682, -0.32770488, -0.42617, -0.38267, -0.30514, -0.17682, -0.06398, -0.01395,  
     &          0.00575, -0.04985, -0.12963, -0.21231, -0.38179, -0.54385, -0.90047,  
     &          -1.20935, -1.72357, -2.15733, -2.89398, -3.49691, -4.05677 / 
      data c2 / -0.57576, -0.478751, -0.57399, -0.52904, -0.44828, -0.29911, -0.15800, -0.09117,  
     &          -0.07933, -0.16741, -0.28157, -0.39457, -0.60337, -0.78423, -1.15145,  
     &          -1.44846, -1.93170, -2.32776, -2.98541, -3.51400, -4.01754 / 
      data c3 / -0.56614, -0.464954, -0.56367, -0.51227, -0.42253, -0.24552, -0.08661, -0.02380,  
     &          -0.04936, -0.17938, -0.33569, -0.48978, -0.75722, -0.98086, -1.39035,  
     &          -1.68527, -2.13907, -2.50355, -3.12984, -3.64054, -4.22973 / 
      data c4 / -0.59679, -0.501209, -0.59948, -0.54896, -0.45758, -0.23559, -0.02770, 0.07032,  
     &          0.06111, -0.05191, -0.19320, -0.33347, -0.56778, -0.76385, -1.20988,  
     &          -1.56650, -2.15102, -2.60304, -3.33785, -3.93845, -4.56577 / 
      data c5 / -0.27336, -0.214600, -0.28004, -0.23702, -0.14609, 0.09781, 0.34210, 0.47686,  
     &          0.49113, 0.32613, 0.12375, -0.07993, -0.43822, -0.72362, -1.30331,  
     &          -1.75247, -2.42949, -2.88125, -3.49756, -3.93043, -4.31969 / 
      data c6 / -0.15066, -0.143591, -0.14965, -0.14920, -0.14601, -0.13802, -0.13783, -0.14249,  
     &          -0.15125, -0.15390, -0.15240, -0.14882, -0.14188, -0.13795, -0.13962,  
     &          -0.14688, -0.14972, -0.13286, -0.09571, -0.06254, -0.04672 / 
      data c7 / 0.18374, 0.1740786, 0.17900, 0.18516, 0.17620, 0.17060, 0.17877, 0.19418, 0.24573,  
     &          0.27648, 0.28582, 0.28807, 0.24809, 0.18660, 0.07403, 0.00114, -0.06290,  
     &          -0.09683, -0.19885, -0.31814, -0.42574 / 
      data c8 / 0.69358, 0.6707077, 0.69599, 0.66643, 0.61448, 0.53997, 0.54911, 0.61070, 0.75072,  
     &          0.87959, 0.98255, 1.06785, 1.21051, 1.32141, 1.51422, 1.64692, 1.82708,  
     &          1.94116, 2.06926, 2.13724, 2.18800 / 
      data c9 / 0.65185, 0.5581807, 0.64948, 0.60930, 0.56920, 0.55149, 0.59056, 0.65929, 0.84132,  
     &          1.00469, 1.14291, 1.24782, 1.37894, 1.46404, 1.56895, 1.63417, 1.65603,  
     &          1.62994, 1.54933, 1.48562, 1.40692 / 
      data c10 / -0.13872, -0.134137, -0.13920, -0.13329, -0.12290, -0.10799, -0.10982, -0.12214,  
     &           -0.15014, -0.17592, -0.19651, -0.21357, -0.24210, -0.26428, -0.30284,   
     &           -0.32938, -0.36542, -0.38823, -0.41385, -0.42745, -0.43760 / 
      data c11 / -0.04953, -0.039918, -0.04495, -0.04775, -0.05401, -0.04126, -0.01951, -0.01249,   
     &           -0.08485, -0.18925, -0.29110, -0.38255, -0.49316, -0.52701, -0.46679,   
     &           -0.35296, -0.16278, -0.07360, -0.00361, 0.00000, 0.00000 / 
      data c12 / 0.02846, 0.0288010, 0.02853, 0.02907, 0.03050, 0.03437, 0.03732, 0.03774, 0.03435,   
     &           0.02885, 0.02380, 0.01946, 0.01337, 0.00967, 0.00559, 0.00419, 0.00260,   
     &           0.00116, -0.00263, -0.00665, -0.01553 / 
      data c13 / 0.00919, 0.0101329, 0.00923, 0.00988, 0.01075, 0.01238, 0.01338, 0.01332, 0.01179,   
     &           0.00982, 0.00785, 0.00624, 0.00411, 0.00268, 0.00117, 0.00087, 0.00059,   
     &           0.00024, -0.00096, -0.00228, -0.00431 / 
      data c14 / -1.90982, -1.935253, -1.90965, -1.93011, -1.97278, -1.95803, -1.87814, -1.77781,   
     &           -1.65632, -1.57576, -1.52224, -1.49528, -1.46281, -1.44274, -1.41615,   
     &           -1.38993, -1.35565, -1.33932, -1.32171, -1.31274, -1.28503 / 
      data c15 / -1.63583, -1.634563, -1.63218, -1.64190, -1.67487, -1.74968, -1.79673, -1.79306,   
     &           -1.72580, -1.64287, -1.57054, -1.51230, -1.43027, -1.37373, -1.25896,   
     &           -1.18777, -1.10246, -1.06588, -1.02734, -0.99037, -0.94811 / 
      data c16 / 0.34806, 0.3584079, 0.34729, 0.35611, 0.37509, 0.39998, 0.38893, 0.36491, 0.31819,   
     &           0.27725, 0.24717, 0.22297, 0.18919, 0.17368, 0.17194, 0.18538, 0.20432,   
     &           0.21646, 0.23573, 0.25375, 0.26914 / 
      data c17 / 0.20176, 0.2166604, 0.20014, 0.21494, 0.22245, 0.20085, 0.15835, 0.12373, 0.08516,   
     &           0.06822, 0.05435, 0.04964, 0.05306, 0.05612, 0.08094, 0.10162, 0.14942,   
     &           0.19708, 0.26366, 0.30340, 0.31352 / 
      data c18 / 0.00000, 0.0003798, 0.00000, 0.00000, -0.00003, -0.00161, -0.00363, -0.00514, -0.00570,   
     &           -0.00512, -0.00429, -0.00334, -0.00191, -0.00106, -0.00012, 0.00000, 0.00000,   
     &           0.00000, -0.00003, -0.00006, -0.00024 / 
      data c19 / -0.00315, -0.003266, -0.00316, -0.00334, -0.00337, -0.00324, -0.00312, -0.00315, -0.00320,   
     &           -0.00306, -0.00283, -0.00259, -0.00222, -0.00198, -0.00196, -0.00196, -0.00200,   
     &           -0.00211, -0.00255, -0.00318, -0.00362 / 
      data c20 / -3.37443, -3.129942, -3.37957, -3.26276, -3.12039, -2.77699, -2.45023, -2.19165, -1.85215,   
     &           -1.62798, -1.49980, -1.44580, -1.45380, -1.52110, -1.59305, -1.46595, -0.89634,   
     &           -0.48495, -0.02762, 0.00000, 0.00000 / 
      data c21 / -0.53674, -0.534043, -0.53652, -0.52911, -0.51328, -0.48078, -0.46572, -0.46949, -0.49343,   
     &           -0.51786, -0.54454, -0.57046, -0.62522, -0.68213, -0.78909, -0.85442, -0.91005,   
     &           -0.92335, -0.91464, -0.89283, -0.85483 / 
      data c22 / 0.04523, 0.0458655, 0.04544, 0.04662, 0.04976, 0.05873, 0.06324, 0.06127, 0.05315, 0.04813,   
     &           0.04907, 0.05354, 0.06521, 0.07589, 0.09600, 0.11018, 0.12800, 0.13712, 0.13748,   
     &           0.12787, 0.11250 / 
      data c23 / -0.58641, -0.544826, -0.58670, -0.55816, -0.50549, -0.40157, -0.32642, -0.30593, -0.34304,   
     &           -0.42439, -0.51573, -0.60347, -0.76233, -0.89923, -1.15933, -1.34357, -1.58317,   
     &           -1.71310, -1.79927, -1.79038, -1.70330 / 
      data c24 / -0.66005, -0.614151, -0.66006, -0.63650, -0.59453, -0.50578, -0.43321, -0.40897, -0.43638,   
     &           -0.50625, -0.58308, -0.65623, -0.79171, -0.91483, -1.16552, -1.35415, -1.59529,   
     &           -1.72062, -1.79429, -1.77800, -1.68840 / 
      data c25 / -0.57611, -0.521623, -0.57563, -0.54296, -0.47735, -0.34285, -0.25612, -0.24798, -0.33656,   
     &           -0.46536, -0.58483, -0.68520, -0.84288, -0.96167, -1.17506, -1.33319, -1.54295,   
     &           -1.65505, -1.71968, -1.70446, -1.64548 / 
      data c26 / -0.09885, -0.098846, -0.09648, -0.05630, -0.01620, 0.00151, -0.03756, -0.10629, -0.28832,   
     &           -0.45169, -0.58991, -0.69482, -0.82594, -0.91104, -0.99861, -1.01217, -0.95103,   
     &           -0.86194, -0.68633, -0.55062, -0.41292 / 
      data c27 / -0.22539, -0.225386, -0.22873, -0.20875, -0.20487, -0.25149, -0.25548, -0.28929, -0.39132,   
     &           -0.49569, -0.58791, -0.65482, -0.72694, -0.76904, -0.80158, -0.82217, -0.79503,   
     &           -0.74594, -0.64933, -0.58562, -0.50692 / 
      data taucr1 / 0.32270, 0.333360, 0.32290, 0.33291, 0.34636, 0.36504, 0.36988, 0.36764, 0.36766, 
     &              0.38435, 0.40767, 0.42422, 0.44443, 0.45178, 0.44390, 0.43348, 0.41861,  
     &              0.39601, 0.37594, 0.36949, 0.42843 / 
      data taucr2 / 0.32595, 0.323792, 0.32442, 0.32587, 0.32941, 0.34135, 0.35775, 0.36540, 0.35112,  
     &              0.32432, 0.30171, 0.28995, 0.29053, 0.31136, 0.38256, 0.43299, 0.48885,  
     &              0.51573, 0.53354, 0.54041, 0.49538 / 
      data tausb1 / 0.25683, 0.223224, 0.25386, 0.25123, 0.24759, 0.24218, 0.24447, 0.25710, 0.31560,  
     &              0.38679, 0.44282, 0.48555, 0.53440, 0.54554, 0.52292, 0.48796, 0.43807,  
     &              0.39735, 0.34954, 0.30732, 0.33465 / 
      data tausb2 / 0.59159, 0.627475, 0.59451, 0.61853, 0.64800, 0.69230, 0.71020, 0.69445, 0.60944,  
     &              0.53922, 0.49368, 0.46672, 0.46083, 0.48019, 0.52884, 0.58232, 0.65015,  
     &              0.68099, 0.66195, 0.60955, 0.52706 / 
      data phisscr1 / 0.53254, 0.547197, 0.53259, 0.52733, 0.52005, 0.50843, 0.50827, 0.52071, 0.55228,  
     &              0.57885, 0.59736, 0.60769, 0.61086, 0.60090, 0.55589, 0.51336, 0.45771,  
     &              0.42721, 0.39838, 0.38461, 0.37698 / 
      data phisscr2 / 0.42278, 0.434184, 0.42262, 0.42746, 0.43800, 0.45841, 0.46418, 0.45535, 0.43212,  
     &              0.41800, 0.41410, 0.41806, 0.43197, 0.44586, 0.47073, 0.48275, 0.48885,  
     &              0.48432, 0.46852, 0.45047, 0.44101 / 
      data phisssb1 / 0.46609, 0.466874, 0.46818, 0.46163, 0.44776, 0.41437, 0.39915, 0.40914, 0.45847,  
     &              0.49875, 0.51974, 0.53135, 0.53489, 0.52593, 0.50290, 0.48257, 0.45684,  
     &              0.44057, 0.42007, 0.40229, 0.37888 / 
      data phisssb2 / 0.46213, 0.472038, 0.46184, 0.46785, 0.47969, 0.50154, 0.50900, 0.50192, 0.47897,  
     &              0.45854, 0.44536, 0.43689, 0.43234, 0.43663, 0.45184, 0.46434, 0.46862,  
     &              0.46433, 0.43642, 0.40445, 0.34992 / 
      data arfacr / 0.25318, 0.229696, 0.25414, 0.26085, 0.26648, 0.25936, 0.23528, 0.21591, 0.21253,  
     &              0.21733, 0.21368, 0.20487, 0.18226, 0.16450, 0.16298, 0.17668, 0.18885,  
     &              0.18576, 0.18142, 0.19052, 0.19394 / 
      data arfasb / 0.19692, 0.172664, 0.20383, 0.20481, 0.20131, 0.19567, 0.19501, 0.19415, 0.18761,  
     &              0.18374, 0.18200, 0.17591, 0.16516, 0.15715, 0.15460, 0.17358, 0.22851,  
     &              0.26757, 0.31216, 0.33304, 0.35954 / 
      data phis2s / 0.31600, 0.312931, 0.31613, 0.32192, 0.33755, 0.38533, 0.42476, 0.43829, 0.42157,  
     &              0.39303, 0.36995, 0.35307, 0.33529, 0.32910, 0.33174, 0.34038, 0.35318,  
     &              0.36003, 0.36784, 0.37191, 0.38065 /

c Set attenuation name                                                            
c     Sourcetype = 0 Crustal
c     Sourcetype = 1 Subduction 
                                                                       
C Find the requested spectral period and corresponding coefficients
      nper = 21

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
        c1T = c1(1)
        c2T = c2(1)
        c3T = c3(1)
        c4T = c4(1)
        c5T = c5(1)
        c6T = c6(1)
        c7T = c7(1)
        c8T = c8(1)
        c9T = c9(1)
        c10T = c10(1)
        c11T = c11(1)
        c12T = c12(1)
        c13T = c13(1)
        c14T = c14(1)
        c15T = c15(1)
        c16T = c16(1)
        c17T = c17(1)
        c18T = c18(1)
        c19T = c19(1)
        c20T = c20(1)
        c21T = c21(1)
        c22T = c22(1)
        c23T = c23(1)
        c24T = c24(1)
        c25T = c25(1)
        c26T = c26(1)
        c27T = c27(1)
        taucr1T = taucr1(1)
        taucr2T = taucr2(1)
        tausb1T = tausb1(1)
        tausb2T = tausb2(1)
        phisscr1T = phisscr1(1)
        phisscr2T = phisscr2(1)
        phisssb1T = phisssb1(1)
        phisssb2T = phisssb2(1)
        arfacrT = arfacr(1)
        arfasbT = arfasb(1)
        phis2sT = phis2s(1)
       goto 1011
C   Function Form for PGAraw, max Regression     
       elseif (specT .eq. -2.0 .or. specT .eq. -1.0) then
         period1 = period(2)
         c1T = c1(2)
         c2T = c2(2)
         c3T = c3(2)
         c4T = c4(2)
         c5T = c5(2)
         c6T = c6(2)
         c7T = c7(2)
         c8T = c8(2)
         c9T = c9(2)
         c10T = c10(2)
         c11T = c11(2)
         c12T = c12(2)
         c13T = c13(2)
         c14T = c14(2)
         c15T = c15(2)
         c16T = c16(2)
         c17T = c17(2)
         c18T = c18(2)
         c19T = c19(2)
         c20T = c20(2)
         c21T = c21(2)
         c22T = c22(2)
         c23T = c23(2)
         c24T = c24(2)
         c25T = c25(2)
         c26T = c26(2)
         c27T = c27(2)
         taucr1T = taucr1(2)
         taucr2T = taucr2(2)
         tausb1T = tausb1(2)
         tausb2T = tausb2(2)
         phisscr1T = phisscr1(2)
         phisscr2T = phisscr2(2)
         phisssb1T = phisssb1(2)
         phisssb2T = phisssb2(2)
         arfacrT = arfacr(2)
         arfasbT = arfasb(2)
         phis2sT = phis2s(2)
         goto 1011      
       endif
C Now loop over the spectral period range of the attenuation relationship.
         do i=3,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1010 
            endif
         enddo
        
      write (*,*) 
      write (*,*) 'Chao et al. (2017) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1010    call S24_interp (period(count1),period(count2),c1(count1),c1(count2), 
     +                specT,c1T,iflag)
         call S24_interp (period(count1),period(count2),c2(count1),c2(count2), 
     +                specT,c2T,iflag)
         call S24_interp (period(count1),period(count2),c3(count1),c3(count2), 
     +                specT,c3T,iflag)
         call S24_interp (period(count1),period(count2),c4(count1),c4(count2), 
     +                specT,c4T,iflag)
         call S24_interp (period(count1),period(count2),c5(count1),c5(count2), 
     +                specT,c5T,iflag)
         call S24_interp (period(count1),period(count2),c6(count1),c6(count2), 
     +                specT,c6T,iflag)
         call S24_interp (period(count1),period(count2),c7(count1),c7(count2), 
     +                 specT,c7T,iflag)
         call S24_interp (period(count1),period(count2),c8(count1),c8(count2), 
     +                 specT,c8T,iflag)
         call S24_interp (period(count1),period(count2),c9(count1),c9(count2), 
     +                 specT,c9T,iflag)
         call S24_interp (period(count1),period(count2),c10(count1),c10(count2), 
     +                 specT,c10T,iflag)
         call S24_interp (period(count1),period(count2),c11(count1),c11(count2), 
     +                 specT,c11T,iflag)
         call S24_interp (period(count1),period(count2),c12(count1),c12(count2), 
     +                 specT,c12T,iflag)
         call S24_interp (period(count1),period(count2),c13(count1),c13(count2), 
     +                 specT,c13T,iflag)
         call S24_interp (period(count1),period(count2),c14(count1),c14(count2), 
     +                 specT,c14T,iflag)
         call S24_interp (period(count1),period(count2),c15(count1),c15(count2), 
     +                 specT,c15T,iflag)
         call S24_interp (period(count1),period(count2),c16(count1),c16(count2), 
     +                 specT,c16T,iflag)
         call S24_interp (period(count1),period(count2),c17(count1),c17(count2), 
     +                 specT,c17T,iflag)
         call S24_interp (period(count1),period(count2),c18(count1),c18(count2), 
     +                 specT,c18T,iflag)
         call S24_interp (period(count1),period(count2),c19(count1),c19(count2), 
     +                 specT,c19T,iflag)
         call S24_interp (period(count1),period(count2),c20(count1),c20(count2), 
     +                 specT,c20T,iflag)
         call S24_interp (period(count1),period(count2),c21(count1),c21(count2), 
     +                 specT,c21T,iflag)
         call S24_interp (period(count1),period(count2),c22(count1),c22(count2), 
     +                 specT,c22T,iflag)
         call S24_interp (period(count1),period(count2),c23(count1),c23(count2), 
     +                 specT,c23T,iflag)
         call S24_interp (period(count1),period(count2),c24(count1),c24(count2), 
     +                 specT,c24T,iflag)
         call S24_interp (period(count1),period(count2),c25(count1),c25(count2), 
     +                 specT,c25T,iflag)
         call S24_interp (period(count1),period(count2),c26(count1),c26(count2), 
     +                 specT,c26T,iflag)
         call S24_interp (period(count1),period(count2),c27(count1),c27(count2), 
     +                 specT,c27T,iflag)
         call S24_interp (period(count1),period(count2),taucr1(count1),taucr1(count2), 
     +                 specT,taucr1T,iflag)
         call S24_interp (period(count1),period(count2),taucr2(count1),taucr2(count2), 
     +                 specT,taucr2T,iflag)
         call S24_interp (period(count1),period(count2),tausb1(count1),tausb1(count2), 
     +                 specT,tausb1T,iflag)
         call S24_interp (period(count1),period(count2),tausb2(count1),tausb2(count2), 
     +                 specT,tausb2T,iflag)
         call S24_interp (period(count1),period(count2),phisscr1(count1),phisscr1(count2), 
     +                 specT,phisscr1T,iflag)
         call S24_interp (period(count1),period(count2),phisscr2(count1),phisscr2(count2), 
     +                 specT,phisscr2T,iflag)
         call S24_interp (period(count1),period(count2),phisssb1(count1),phisssb1(count2), 
     +                 specT,phisssb1T,iflag)
         call S24_interp (period(count1),period(count2),phisssb2(count1),phisssb2(count2), 
     +                 specT,phisssb2T,iflag)
         call S24_interp (period(count1),period(count2),arfacr(count1),arfacr(count2), 
     +                 specT,arfacrT,iflag)
         call S24_interp (period(count1),period(count2),arfasb(count1),arfasb(count2), 
     +                 specT,arfasbT,iflag)
         call S24_interp (period(count1),period(count2),phis2s(count1),phis2s(count2), 
     +                specT,phis2sT,iflag)
   
 1011 period1 = specT

      h = 10.0
      n = 2.0
      Mc = 7.1
      Mref = 5.5
      Mmax = 8
      Rrupref = 0.0
      Vs30ref = 760.0
    
C     Set the reference spectrum.                
c     sourcetype = 0 for crustal
c                  1 for Subduction 
c     Vs30_class = 0 for estimated
c     Vs30_class = 1 for measured 

      Fcr=0
      Fsb=0
      Fcrss = 0
      Fcrno = 0
      Fcrro = 0
      Fsbintra = 0
      Fsbinter = 0
      Fas = 0
      Fkuo17 = 0
      Fks17 = 0
      Frf = 0
      Fmanila = 0
      C11flag = 0
      C20flag = 0
      C26flag = 0
   
      if (sourcetype .eq. 0.0 ) then
       Fcr = 1
       Zref = 15
         if(ftype .gt. 0) then
              Fcrro = 1
           elseif(ftype .lt. 0) then
              Fcrno = 1
           else
              Fcrss = 1
         endif
      elseif (sourcetype .eq. 1.0 ) then
        Fsb = 1
          Zref = 50
         if(ftype .eq. 0) then
              Fsbinter = 1
           elseif(ftype .eq. 1) then
              Fsbintra = 1
         endif
      endif
   
C     Add aftershock factor 
      if (msasflag .eq. 1) then
           Fas = 1
      endif 

C     choose Site ref by Vs30 class
        if (vs30_class .eq. 0 ) then
         Fks17 = 1
        elseif (vs30_class .eq. 1) then
         Fkuo17 = 1
        endif

      lnYref = c1T*Fcrro + c2T*Fcrss + c3T*Fcrno + c4T*Fsbinter + c5T*Fsbintra +
     &         c6T*Fas + c7T*Fmanila + c23T*Fkuo17 + c24T*Fks17 + c25T*Frf

C     Set Source scaling term 
     
      if(mag .LE. 5 ) then  
       C11flag=1
      endif
      if(mag .GE. Mc ) then  
       C26flag=1
      endif
      if(mag .GE. 7.6 ) then  
       C10flag=1
      endif   
   
      if (sourcetype .eq. 0.0 ) then
        Smag = c8T*(mag - Mref) + c10T*(mag - Mref)**2 
     1         - c10T*(mag-7.6)**2*C10flag + c11T*(5-mag)*C11flag  
      elseif (sourcetype .eq. 1.0 ) then
        Smag = c9T*(mag - Mref) + c26T*Fsbinter*(Mag-Mc)*c26flag + c27T*Fsbintra*(Mag-Mc)*c26flag 
      endif

      Sztor = c12T * Fcr *(Ztor-Zref) + c13T * Fsb * ((Ztor-Zref))   
      Ssource = Smag + Sztor

C     Set Path scaling term

      if (sourcetype .eq. 0.0 ) then
          Sgeom = (c14T + c16T*(min(mag,Mmax)- Mref )) * alog(SQRT(dist**2 + h**2)/SQRT(Rrupref**2 + h**2))
      elseif (sourcetype .eq. 1.0 ) then
          Sgeom = (c15T + c17T*(min(mag,Mc)- Mref )) * alog(SQRT(dist**2 + h**2)/SQRT(Rrupref**2 + h**2))
      endif

      Sanel = c18T*Fcr*(dist-Rrupref) + c19T*Fsb*(dist-Rrupref)
      Spath = Sgeom + Sanel 
    
C     Set Site scaling term 
    
      Z10ref = exp((-4.08/2)*alog((vs**2+355.4**2)/(1750**2+355.4**2)))
      Ssitelin = c21T * alog(vs/vs30ref) + c22T*alog(Z10*1000/Z10ref)

      if(vs .LT. vs30ref ) then  
           C20flag=1
      endif
     
      Ssitenon = c20T * C20flag * (-1.5*alog(vs/vs30ref)-alog(SA1100+2.4)+alog(SA1100+2.4*(vs/vs30ref)**1.5))  
      Ssite = Ssitenon + Ssitelin

      lnSa =  lnYref + Ssource + Spath + Ssite                                        
   
c    write(*,*) "lnYref = ", lnYref
c    write(*,*) "Ssource = ", Ssource
c    write(*,*) "Spath = ", Spath
c    write(*,*) "Ssitelin = ", Ssitelin
c    write(*,*) "Ssitenon = ", Ssitenon
c    write(*,*) "lnSa = ", lnSa
c    write(*,*) "Sa = ", exp(lnSa)

   
C     Set the event-specific residual term
 
      fm = 0.5*(min(6.5, max(4.5, mag))-4.5)
      
      taucr = taucr1T + (taucr2T - taucr1T)*fm
      tausb = tausb1T + (tausb2T - tausb1T)*fm 
      
      tau = taucr*Fcr + tausb*Fsb   
   
C     Set Site-specific residual term

      

C     Set Recoed-specific residual term

      phisscr = phisscr1T + (phisscr2T -phisscr1T)*fm
      phisssb = phisssb1T + (phisssb2T -phisssb1T)*fm

      phiss = phisscr*Fcr + phisssb*Fsb
      
      phi=(phis2sT**2+phiss**2)**0.5
      sigma=(tau**2+phi**2)**0.5
      sigmass=(tau**2+phiss**2)**0.5



c       write(*,*) "Y(gal) = ", exp(lnSa)

      return                                                                    
      end       
          
c ------------------------------------------------------------------            
C *** Phung2017 Crust and Subduction Model- Horizontal ***********
c ------------------------------------------------------------------            

      Subroutine S04_PhungCrust2017 ( m, Rrup, specT, period2, lnY, sigma, iflag, 
     1                     vs, Delta, DTor, Ftype, depthvs10, vs30_class,
     2                       regionflag, msasflag, phi, tau )


      implicit none
      
      integer MAXPER, i, nPer
      parameter (MAXPER=23)
      REAL Period(MAXPER), C1(MAXPER), C1a(MAXPER), C1b(MAXPER), C1c(MAXPER), C1d(MAXPER)
      REAL cn(MAXPER), cm(MAXPER), c3(MAXPER), c5(MAXPER), c6(MAXPER)
      REAL c7(MAXPER), C7b(MAXPER), C11(MAXPER), C11b(MAXPER), CHM(MAXPER)
      REAL phi1(MAXPER), phi2(MAXPER), phi3(MAXPER), phi4(MAXPER), phi5(MAXPER)
      REAL sigma1inf(MAXPER), sigma2inf(MAXPER)
      REAL tau1(MAXPER), tau2(MAXPER), sigma1, sigma2
      REAL sigma3(MAXPER), c9a(MAXPER), deltac5(MAXPER)
      REAL cgglb1(MAXPER), cgglb2(MAXPER), cgglb3(MAXPER), cge(MAXPER)
      Real tauT1(MAXPER), phiT1(MAXPER), cgtw1(MAXPER), cgtw2(MAXPER), cgtw3(MAXPER)
      real phiss(MAXPER), phis2s(MAXPER), sigma1mea(MAXPER), sigma2mea(MAXPER)
      real phiss1M(MAXPER), phiss2M(MAXPER)
      real phiglb1(MAXPER),phiglb2(MAXPER),phiglb3(MAXPER),phitw1(MAXPER),phitw2(MAXPER),phitw3(MAXPER)
      real vs, phi6
      real Finferred, Fmeasured
      REAL c1T, c1aT, c1bT, c1cT, c1dT,cnT, cmT, c5T, c6T, c3T
      REAL phi1T, phi2T, phi3T, phi4T, cgeT, deltac5T, sigma3T
      REAL phi5T, tau1T, tau2T, phiT, tauT
      real c7T, c7bT, c11T, c11bT, cHMT, cgglb1T, cgglb2T, cgglb3T
      real cgtw1T, cgtw2T, cgtw3T, sigma1infT, sigma2infT
      real phiglb1T,phiglb2T,phiglb3T,phitw1T,phitw2T,phitw3T
      REAL c2, c4, c4a, cRB, pi, d2r, term14, term15, term16, NL0
      REAL term1, term2, term3, term5, term4, term6, term8, term9, term10 
      REAL phissT, phis2sT, sigma1meaT, sigma2meaT, phiss1MT, phiss2MT
      real CNS, cosdelta, psa_ref, psa, cg1T, cg2T, cg3T
      integer iflag, count1, count2, vs30_class, regionflag, msasflag
      REAL M, RRUP, DTOR, Delta, specT, sigma, Ftype
      REAL period2, lnY, F_RV, F_NM, tau, phi, rkdepth
      real c8, c8a, c8bT, fd, lnpsa_ref, lnpsa, sa
      real sigmaNL0, F_Measured, F_Inferred, mz_TOR, deltaZ_TOR, coshM
      real period1 ,Ez1, term7, deltaZ1, depthvs10

C     Mainshock and Aftershocks included based on MSASFlag
C         0 = Mainshocks
C         1 = Aftershocks
C
C     regionflag  
C           = 1 for Taiwan
C           = 0 for global
C
C     vs30_class     Note
C     -------------------------
C      0         estimated
C      1         measured
C

      data period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 
     1            0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5.0 /
      
      
      data c1  / -1.5095, -1.5137, -1.4795, -1.3729, -1.2484, -1.0969, -0.7407, -0.5699, -0.5059, 
     1           -0.5061, -0.5538, -0.5932, -0.7551, -0.9194, -1.2439, -1.4716, -2.0247, -2.3813,  
     1           -2.9856, -3.3075, -3.7822, -4.0489, -4.2928  / 
      data c3  / 1.2468, 1.253, 1.2231, 1.1917, 1.1374, 1.0586, 1.1135, 1.18, 1.2719, 1.3754, 1.4184,  
     1           1.523, 1.6667, 1.8271, 2.029, 2.2367, 2.5124, 2.6398, 2.7256, 2.7831, 2.8165, 2.8036,  
     1           2.9213  / 
      data cn  / 12.3145, 12.3349, 12.0592, 10.9614, 12.7233, 16.2717, 16.493, 16.1585, 17.1615,  
     1           13.5472, 14.929, 15.8416, 12.7457, 9.6295, 6.6784, 4.5576, 3.4064, 3.1612, 2.8078,  
     1           2.4631, 2.2111, 1.9468, 1.8671  / 
      data cm  / 5.0417, 5.0419, 5.045, 5.0712, 5.0761, 5.0774, 5.0706, 5.0378, 5.031, 5.0801,  
     1           5.0752, 5.0915, 5.133, 5.1763, 5.22, 5.2948, 5.4189, 5.5442, 5.7193, 5.8612, 6.0824,  
     1           6.2674, 6.4197  / 
      data c5  / 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.8305, 7.1333, 7.3621,  
     1           7.4365, 7.4972, 7.54E+00, 7.56, 7.5735, 7.5778, 7.5808, 7.5814, 7.5817, 7.5818,  
     1           7.5818, 7.5818, 7.5818  / 
      data deltac5  / -1.8316, -1.8293, -1.8424, -1.8292, -1.9728, -1.9608, -1.867, -1.7551, -1.7174,  
     1           -1.6555, -1.5751, -1.5154, -1.48E+00, -1.6291, -0.9786, -1.1988, -1.6463, -1.8809,  
     1           -2.3783, -2.5278, -2.1459, -2.3076, -2.3007  / 
      data c6  / 0.4908, 0.4908, 0.4925, 0.4992, 0.5037, 0.5048, 0.5048, 0.5048, 0.5048, 0.5045,  
     1           0.5036, 0.5016, 4.97E-01, 0.4919, 0.4807, 0.4707, 0.4575, 0.4522, 0.4501, 0.45,  
     1           0.45, 0.45, 0.45  / 
      data cHM  / 3.0956, 3.0956, 3.0963, 3.0974, 3.0988, 3.1011, 3.1094, 3.2381, 3.3407, 3.43,  
     1           3.4688, 3.5146, 3.57E+00, 3.6232, 3.6945, 3.7401, 3.7941, 3.8144, 3.8284, 3.833,  
     1           3.8361, 3.8369, 3.8376  / 
      data c7  / 0.0141, 0.0142, 0.0142, 0.0151, 0.015, 0.0157, 0.0158, 0.0169, 0.018, 0.0187,  
     1           0.0194, 0.0187, 1.72E-02, 0.0171, 0.0164, 0.0153, 0.0133, 0.0129, 0.0104, 0.0078,  
     1           0, -0.0149, -0.0208  / 
      data c7b  / 0.0127, 0.0127, 0.0131, 0.0152, 0.0149, 0.0154, 0.0166, 0.0134, 0.0098, 0.0053,  
     1           0.0054, 0.0022, 3.00E-03, 0.0018, -0.0013, 0.0016, 0.0029, 0.004, 0.0034, 0.0067,  
     1           0.0068, 0.0243, 0.037  / 
      data c1a  / 0.1502, 0.1466, 0.1392, 0.1312, 0.1196, 0.0981, 0.1021, 0.1027, 0.1109, 0.1251,  
     1           0.1256, 0.1561, 1.59E-01, 0.1673, 0.1645, 0.2297, 0.2002, 0.215, 0.1334, 0.07,  
     1           0.0373, 0, 0  / 
      data c1c  / -0.0243, -0.0199, -0.0122, 0.0072, 0.0379, 0.0649, 0.0714, 0.0168, -0.0157, -0.0643,  
     1           -0.0658, -0.0806, -5.11E-02, -0.0819, -0.0389, -0.1289, -0.0511, -0.0223, 0.0563,  
     1           0.108, 0.1459, 0.1726, 0.1868  / 
      data c1b  / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.1549, -0.1815, -0.1998, -0.1209,  
     1           -0.081, -0.0867, -0.0867  / 
      data c1d  / -0.1712, -0.1712, -0.1712, -0.1712, -0.1712, -0.1712, -0.1712, -0.1712, -0.1712,  
     1           -0.1712, -0.1712, -0.1712, -0.1712, -0.1712, -0.1712, -0.1712, -0.1712, -0.1712,  
     1           -0.1712, -0.1712, -0.1712, -0.1712, -0.1712  / 
      data c11  / -0.1234, -0.1072, -0.0854, -0.0852, -0.0787, -0.0522, -0.1591, -0.1788, -0.1599,  
     1           -0.1231, -0.0645, -0.1604, -0.1152, -0.1613, -0.0896, -0.2766, -0.2254, -0.2945,  
     1           -0.0468, 0.0748, 0.10033, 0.19, 0.2357  / 
      data c11b  / -0.1677, -0.1522, -0.1808, -0.2026, -0.2264, -0.2688, -0.0894, 0.0172, -0.0037,  
     1           -0.0467, -0.1003, -0.0431, -0.1682, -0.0271, -0.1761, 0.1431, 0.0775, 0.1553, -0.1704,  
     1           -0.2993, -0.3915, -0.489, -0.2708  / 
      data cgglb1  / -0.00777287, -0.00777287, -7.87E-03, -0.008513485, -0.0090007, -0.0094918,  
     1           -0.010315, -0.010597, -0.0108897, -0.010839, -0.010645, -1.03E-02, -8.79E-03,  
     1           -8.31E-03, -7.19E-03, -0.006207, -4.47E-03, -3.85E-03, -3.11E-03, -3.05E-03,  
     1           -0.0025055, -0.00210659, -0.00194175  / 
      data cgglb2  / -0.008307401, -0.008307401, -0.008224304, -0.008319595, -0.008383, -0.0084453,  
     1           -0.0075936, -0.0061128, -0.0041812, -0.002678, -0.001964, -0.0015011, -0.0013538,  
     1           -0.0010318, -0.0004922, -0.001285, -0.0011589, -0.001824, -0.0014489, -0.0034758,  
     1           -0.024638, -0.019669, -0.0141327  / 
      data cgglb3  / 4.052466395, 4.052466395, 4.043855472, 4.04820785, 4.04402, 4.1283, 4.3754,  
     1           4.6228, 4.8756, 4.9326, 4.95E+00, 4.7776, 5.5509, 5.2382, 5.3261, 4.38, 4.6225,  
     1           4.636, 4.18648, 4.0464, 1.21957, 1.79675, 1.52486  / 
      data cgtw1  / -0.007289154, -0.007289154, -0.007391413, -0.008045598, -0.0084688, -0.00899689,  
     1           -0.01003, -0.010551, -0.011115, -0.011158, -0.0104477, -0.009978, -0.0090725,  
     1           -0.0082395, -0.00699, -0.005702, -0.0044667, -0.003847, -0.003111, -0.002558,  
     1           -0.0018576, -0.0015354, -0.0014361  / 
      data cgtw2  / -0.007103831, -0.007103831, -0.007027146, -0.007258757, -0.0072455, -0.0075512,  
     1           -0.007104, -0.005758, -0.0038376, -0.002273, -0.002101, -0.0015013, -0.00069154,  
     1           -0.0006181, -0.00021158, -0.0008504, -0.0011589, -0.001824, -0.0014489, -0.0027134,  
     1           -0.002776, -0.0087011, -0.007748  / 
      data cgtw3  / 4.22959614, 4.22959614, 4.222091708, 4.218021689, 4.2343, 4.27912, 4.4406, 4.6068,  
     1           4.7495, 4.6934, 4.9878, 5.1215, 5.419, 5.3499, 5.4856, 5.1327, 4.6225, 4.636, 4.18648,  
     1           4.443, 3.94986, 2.4836, 2.1575  / 
      data cge  / -0.002154124, -0.002154962, -0.002161527, -0.002149373, -0.002236978, -0.002258625,  
     1           -0.002171027, -0.002242017, -0.002190855, -0.00216606, -0.002071088, -0.002076124,  
     1           -0.002032334, -0.002000855, -0.001908347, -0.001877458, -0.0017069, -0.001368707,  
     1           -0.001190811, -0.000899898, -0.001337119, -0.00138601, -0.001628408  / 
      data phiglb1  / -0.5088, -0.5088, -0.4995, -0.4863, -0.4686, -4.47E-01, -0.428, -0.4577, -0.4793, 
     1            -0.5078, -0.5298, -0.5513, -0.5801, -0.6311, -0.6803, -0.7113, -0.8294, -0.9068,  
     1           -0.9733, -0.9776, -0.9775, -0.9659, -0.9258  / 
      data phiglb2  / -0.1417, -0.1417, -0.1364, -0.1403, -0.1591, -0.1862, -0.2538, -0.2943, -0.3077, 
     1            -0.3113, -0.3062, -0.2927, -0.2662, -0.2405, -0.1975, -0.1633, -0.1028, -0.0699, 
     1            -0.0425, -0.0302, -0.0129, -0.0016, 0  / 
      data phiglb3  / -0.00701, -0.00701, -0.007279, -0.007354, -0.006977, -0.006467, -0.0057, -0.0056, 
     1            -0.0057, -0.0058, -0.006, -0.0061, -0.0064, -0.0067, -0.0071, -0.0074, -0.0081, 
     1            -0.0084, -0.0077, -0.0048, -0.0018, -0.0015, -0.0014  / 
      data phitw1  / -0.529702044, -0.529656391, -0.521077213, -0.512073039, -0.498112888, -4.88E-01, 
     1            -0.453835736, -0.476664485, -4.96E-01, -0.507200973, -0.519826139, -0.534817785, 
     1            -0.559543859, -0.605127453, -0.67522556, -0.706879609, -0.853566618, -0.883170749, 
     1            -0.964946569, -0.979694815, -0.994511901, -0.997658667, -0.968404707  / 
      data phitw2  / -0.2226, -0.2226, -0.2178, -0.2221, -0.2513, -0.2958, -0.3235, -0.3378, -0.3723, 
     1            -0.3992, -0.3634, -0.4163, -0.3133, -0.3093, -0.2983, -0.2167, -0.2075, -0.0699, 
     1            -0.0425, -0.0302, -0.0129, -0.0016, 0  / 
      data phitw3  / -0.0076, -0.0076, -0.0078, -0.0079, -0.0077, -0.0074, -0.006, -0.0063, -0.0061, 
     1            -0.0056, -0.0061, -0.0055, -0.0075, -0.0077, -0.0079, -0.0092, -0.0095, -0.0084, 
     1            -0.0077, -0.0048, -0.0018, -0.0015, -0.0014  / 
      data phi4   / 0.1022, 0.1022, 0.1084, 0.1199, 0.1336, 0.1489, 0.1906, 0.2307, 0.2532, 0.2665, 
     1            0.2651, 0.2553, 0.2315, 0.2073, 0.1655, 0.1338, 0.0852, 0.0586, 0.0318, 0.0197, 
     1            0.0096, 0.0054, 0.0032  / 
      data phi5  / 0.040582663, 0.040603911, 0.040793653, 0.048275069, 0.057526547, 0.064938488,  
     1           0.082071634, 0.107265026, 0.093297608, 0.070415091, 0.060452567, 0.050253034,  
     1           0.042040485, 0.034861903, 0.0293967, 0.040453687, 0.070852531, 0.086614739,  
     1           0.160460769, 0.22535196, 0.272648022, 0.248210653, 0.221285413  / 
      data tau1  / 0.328705918, 0.32824584, 0.331061862, 0.358256174, 0.381906279, 0.404559567,  
     1           0.388100354, 0.366211703, 0.330648493, 0.329813563, 0.322463538, 0.315621419,  
     1           0.30519909, 0.311612501, 0.359831951, 0.397825282, 0.377486688, 0.39578383,  
     1           0.394229339, 0.405498239, 0.402609845, 0.381774468, 0.342688475  / 
      data tau2  / 0.372296285, 0.37547105, 0.372011554, 0.375978951, 0.390221346, 0.418424682,  
     1           0.474682772, 0.493852233, 0.491612285, 0.468071578, 0.444187615, 0.389796361,  
     1           0.355660861, 0.343225478, 0.310855105, 0.330660711, 0.396500599, 0.420526746,  
     1           0.490908629, 0.469648447, 0.510213646, 0.461885197, 0.465570287  / 
      data sigma3  / 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7999, 0.7997,  
     1           0.7988, 0.7966, 0.7792, 0.7504, 0.7136, 0.7035, 0.7006, 0.7001, 0.7  / 
      data sigma1mea  / 0.440462645, 0.440237517, 0.440552924, 0.454924612, 0.471237638, 0.500836323,  
     1           0.493188218, 0.472404434, 0.471527522, 0.450726803, 0.447200284, 0.433579051,  
     1           0.434137568, 0.435776383, 0.425775914, 0.413369858, 0.384583422, 0.359153292,  
     1           0.369920274, 0.361673931, 0.362299163, 0.359960063, 0.363572219  / 
      data sigma2mea  / 0.383905228, 0.383988636, 0.385173862, 0.398286282, 0.406936953, 0.419087157,  
     1           0.423852303, 0.432312248, 0.424703436, 0.41528495, 0.41504729, 0.417764253,  
     1           0.420837422, 0.429538727, 0.418804302, 0.431150448, 0.429025695, 0.45285511,  
     1           0.466343119, 0.477325584, 0.461878034, 0.444550934, 0.470269388  / 
      data sigma1inf  / 0.398123649, 0.398039907, 0.395199187, 0.405412008, 0.420115937, 0.421842238,  
     1           0.443501473, 0.455007817, 0.393576914, 0.398285761, 0.401960163, 0.478368061, 
     1           0.501817351, 0.441522586, 0.460174194, 0.422873766, 0.38843105, 0.394761518,  
     1           0.366757709, 0.398227355, 0.437440799, 0.375321005, 0.354401467  / 
      data sigma2inf  / 0.340567681, 0.340684182, 0.347366065, 0.306296924, 0.308395153, 0.29729216,  
     1           0.337069527, 0.320906139, 0.368035681, 0.341855716, 0.393292587, 0.370778536, 
     1            0.382037432, 0.371645616, 0.412904011, 0.385657898, 0.455887057, 0.448865694,  
     1           0.427656647, 0.475169411, 0.440811156, 0.459470633, 0.459536833  / 
      data phiss  / 0.447212907, 0.447102969, 0.446397587, 0.453095524, 0.46011904, 0.464069345,  
     1           0.468531341, 0.463696461, 0.457623166, 0.452211849, 0.52114735, 0.518086792,  
     1           0.524082259, 0.529838949, 0.537774644, 0.546980781, 0.547514088, 0.549963272,  
     1           0.55739441, 0.557527617, 0.545825976, 0.538054814, 0.552545987  / 
      data phis2s  / 0.289092069, 0.289175184, 0.287996213, 0.296191534, 0.310294276, 0.330992478,  
     1           0.36510545, 0.371308208, 0.367758691, 0.347494163, 0.226689915, 0.20097083,  
     1           0.164392064, 0.149065552, 0.138765534, 0.13996281, 0.126256451, 0.120197996,  
     1           0.12257396, 0.129125085, 0.13379128, 0.135468401, 0.133353227  / 
      data phiss1M  / 0.506329173, 0.506540456, 0.504282764, 0.522559148, 0.543975792, 0.564926662,  
     1           0.556630958, 0.537009107, 0.535309967, 0.53305874, 0.537071855, 0.533736576,  
     1           0.55317255, 0.554129774, 0.521646021, 0.488907503, 0.471282128, 0.449064043,  
     1           0.447660372, 0.414099867, 0.442978244, 0.401052695, 0.409456295  / 
      data phiss2M  / 0.415029056, 0.415220494, 0.415874458, 0.418944451, 0.421676471, 0.426776548, 
     1            0.427088307, 0.444029842, 0.432955207, 0.429122746, 0.440129094, 0.444863496,  
     1           0.461534752, 0.467060561, 0.443935541, 0.451123589, 0.464532391, 0.469747607, 
     1            0.487492461, 0.545482623, 0.4770667, 0.555766814, 0.567308177  / 


C Find the requested spectral period and corresponding coefficients
      nPer = 23
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T  =   c1(1)
         c3T  =   c3(1)
         cnT  =   cn(1)
         cmT  =   cm(1)
         c5T  =   c5(1)
         deltac5T  =   deltac5(1)
         c6T  =   c6(1)
         cHMT  =   cHM(1)
         c7T  =   c7(1)
         c7bT  =   c7b(1)
         c1aT  =   c1a(1)
         c1cT  =   c1c(1)
         c1bT  =   c1b(1)
         c1dT  =   c1d(1)
         c11T  =   c11(1)
         c11bT  =   c11b(1)
         cgglb1T  =   cgglb1(1)
         cgglb2T  =   cgglb2(1)
         cgglb3T  =   cgglb3(1)
         cgtw1T  =   cgtw1(1)
         cgtw2T  =   cgtw2(1)
         cgtw3T  =   cgtw3(1)

         phiglb1T  =   phiglb1(1)
         phiglb2T  =   phiglb2(1)
         phiglb3T  =   phiglb3(1)
         phitw1T  =   phitw1(1)
         phitw2T  =   phitw2(1)
         phitw3T  =   phitw3(1)
         phi4T  =   phi4(1)
         phi5T  =   phi5(1)

         tau1T  =   tau1(1)
         tau2T  =   tau2(1)
         sigma3T  =   sigma3(1)
         sigma1meaT  =   sigma1mea(1)
         sigma2meaT  =   sigma2mea(1)
         sigma1infT  =   sigma1inf(1)
         sigma2infT  =   sigma2inf(1)
         phissT  =   phiss(1)
         phis2sT  =   phis2s(1)
         phiss1MT  =   phiss1M(1)
         phiss2MT  =   phiss2M(1)

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
      write (*,*) 'Phung et al. 2017 horizontal'
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
 1020       call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +                   specT,c3T,iflag)
            call S24_interp (period(count1),period(count2),cn(count1),cn(count2),
     +                   specT,cnT,iflag)
            call S24_interp (period(count1),period(count2),cm(count1),cm(count2),
     +                   specT,cmT,iflag)
            call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +                   specT,c5T,iflag)
            call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +                   specT,c6T,iflag)
            call S24_interp (period(count1),period(count2),cHM(count1),cHM(count2),
     +                   specT,cHMT,iflag)
            call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +                   specT,c7T,iflag)
            call S24_interp (period(count1),period(count2),c7b(count1),c7b(count2),
     +                   specT,c7bT,iflag)

            call S24_interp (period(count1),period(count2),c1a(count1),c1a(count2),
     +                   specT,c1aT,iflag)
            call S24_interp (period(count1),period(count2),c1b(count1),c1b(count2),
     +                   specT,c1bT,iflag)
            call S24_interp (period(count1),period(count2),c1c(count1),c1c(count2),
     +                   specT,c1cT,iflag)
            call S24_interp (period(count1),period(count2),c1d(count1),c1d(count2),
     +                   specT,c1dT,iflag)
            call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +                   specT,c11T,iflag)
            call S24_interp (period(count1),period(count2),c11b(count1),c11b(count2),
     +                   specT,c11bT,iflag)

            call S24_interp (period(count1),period(count2),cgglb1(count1),cgglb1(count2),
     +                   specT,cgglb1T,iflag)
            call S24_interp (period(count1),period(count2),cgglb2(count1),cgglb2(count2),
     +                   specT,cgglb2T,iflag)
            call S24_interp (period(count1),period(count2),cgglb3(count1),cgglb3(count2),
     +                   specT,cgglb3T,iflag)
            call S24_interp (period(count1),period(count2),cgtw1(count1),cgtw1(count2),
     +                   specT,cgtw1T,iflag)
            call S24_interp (period(count1),period(count2),cgtw2(count1),cgtw2(count2),
     +                   specT,cgtw2T,iflag)
            call S24_interp (period(count1),period(count2),cgtw3(count1),cgtw3(count2),
     +                   specT,cgtw3T,iflag)
            call S24_interp (period(count1),period(count2),cge(count1),cge(count2),
     +                   specT,cgeT,iflag)
     
            call S24_interp (period(count1),period(count2),phiglb1(count1),phiglb1(count2),
     +                   specT,phiglb1T,iflag)
            call S24_interp (period(count1),period(count2),phiglb2(count1),phiglb2(count2),
     +                   specT,phiglb2T,iflag)
            call S24_interp (period(count1),period(count2),phiglb3(count1),phiglb3(count2),
     +                   specT,phiglb3T,iflag)
            call S24_interp (period(count1),period(count2),phitw1(count1),phitw1(count2),
     +                   specT,phitw1T,iflag)
            call S24_interp (period(count1),period(count2),phitw2(count1),phitw2(count2),
     +                   specT,phitw2T,iflag)
            call S24_interp (period(count1),period(count2),phitw3(count1),phitw3(count2),
     +                   specT,phitw3T,iflag)

            call S24_interp (period(count1),period(count2),phi4(count1),phi4(count2),
     +                   specT,phi4T,iflag)
            call S24_interp (period(count1),period(count2),phi5(count1),phi5(count2),
     +                   specT,phi5T,iflag)

            call S24_interp (period(count1),period(count2),phiss(count1),phiss(count2),
     +                   specT,phissT,iflag)
            call S24_interp (period(count1),period(count2),phis2s(count1),phis2s(count2),
     +                   specT,phis2sT,iflag)
            call S24_interp (period(count1),period(count2),tau1(count1),tau1(count2),
     +                   specT,tau1T,iflag)
            call S24_interp (period(count1),period(count2),tau2(count1),tau2(count2),
     +                   specT,tau2T,iflag)

            call S24_interp (period(count1),period(count2),sigma1mea(count1),sigma1mea(count2),
     +                   specT,sigma1meaT,iflag)
            call S24_interp (period(count1),period(count2),sigma2mea(count1),sigma2mea(count2),
     +                   specT,sigma2meaT,iflag)
            call S24_interp (period(count1),period(count2),sigma1inf(count1),sigma1inf(count2),
     +                   specT,sigma1infT,iflag)
            call S24_interp (period(count1),period(count2),sigma2inf(count1),sigma2inf(count2),
     +                   specT,sigma2infT,iflag)
          
            call S24_interp (period(count1),period(count2),phiss1M(count1),phiss1M(count2),
     +                   specT,phiss1MT,iflag)
            call S24_interp (period(count1),period(count2),phiss2M(count1),phiss2M(count2),
     +                   specT,phiss2MT,iflag)

 1011 period1 = specT                                                                                                              

c     Set the fault mechanism term.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
         if (ftype .eq. -1) then
            F_RV = 0.0
            F_NM = 1.0
         elseif (ftype .ge. 0.5) then
            F_RV = 1.0
            F_NM = 0.0
         else
            F_RV = 0.0
            F_NM = 0.0
         endif

C     Constant terms
        c2 = 1.06
        c4 = -2.1
        c4a = -0.5
        cRB = 50.0
        phi6 = 300

  
      if(regionflag .eq. 0) then
       
        cg1T = cgglb1T
        cg2T = cgglb2T
        cg3T = cgglb3T
        phi1T = phiglb1T
        phi2T = phiglb2T
        phi3T = phiglb3T

      elseif(regionflag .eq. 1) then
       
        cg1T = cgtw1T
        cg2T = cgtw2T
        cg3T = cgtw3T
        phi1T = phitw1T
        phi2T = phitw2T
        phi3T = phitw3T
       
      endif

C     Current code set for Measured Vs30 values (i.e., Vs30class=1)
      if (vs30_class .eq. 0) then
         Fmeasured = 0.0
         FInferred = 1.0
         sigma1 = sigma1infT
         sigma2 = sigma2infT

      elseif (vs30_class .eq. 1) then      
         Fmeasured = 1.0
         FInferred = 0.0
         sigma1 = sigma1meaT
         sigma2 = sigma2meaT

      endif       
  
c Center Z_TOR on the Z_TOR-M relation
        if (F_RV.EQ.1) then

            mZ_TOR = max(4.069-1.992*max(M-5.776,0.0),0.0)
            mZ_TOR = mZ_TOR * mZ_TOR

        else
            mZ_TOR = max(3.419-1.961*max(M-5.403,0.0),0.0)
            mZ_TOR = mZ_TOR * mZ_TOR
        endif
c        if (Z_TOR .EQ. -999) Z_TOR = mZ_TOR
        deltaZ_TOR = Dtor - mZ_TOR

c Reference motion  
  
        pi = atan(1.0)*4.0
        d2r = pi/180.0
        term1 = c1T
  
c Magnitude scaling
        term6 = c2 * (M-6.0) 
        term7 = (c2-c3T)/cnT * alog(1.0 + exp(cnT*(cMT-M)))  

c Near-field magnitude and distance scaling
        if (Dtor <= 20) then
          CNS = c5T* cosh(c6T * max((M-cHMT),0.0))   
        else
         CNS = (c5T + deltac5T) * cosh(c6T * max((M-cHMT),0.0))  
        endif
       
        term8 = c4 * alog(Rrup + CNS)                 

c Distance scaling at large distance
        term9 = (c4a-c4) * alog( sqrt(Rrup*Rrup+cRB*cRB) )  
        term10 = (cg1T + cg2T/cosh(max((M-cg3T),0.0)))*Rrup  


c Scaling with other source variables (F_RV, F_NM, deltaZ_TOR, and Dip)
        coshM = cosh(2*max(M-4.5,0.0))
        cosDELTA = cos(DELTA*d2r)
        term2 = (c1aT+c1cT/coshM) * F_RV 
        term3 = (c1bT+c1dT/coshM) * F_NM 
        term4 = (c7T +c7bT/coshM) * deltaZ_TOR 
        term5 = (c11T+c11bT/coshM)* cosDELTA**2   
        
c Predicted median Sa on reference condition (Vs=1130 m/sec)
        lnpsa_ref = term1+term2+term3+term5+term4+term6+term7+term8+term9+term10
        psa_ref = exp(lnpsa_ref)
  
c Linear soil amplification
        term14 = phi1T * min(alog(Vs/1130.0), 0.0)   

c Nonlinear soil amplification
        term15 = phi2T *
     1      (exp(phi3T*(min(Vs,1130.0)-360.0)) - exp(phi3T*(1130.0-360.0)))*
     1      alog((psa_ref+phi4T)/phi4T)

C Deviation from ln(Vs30) scaling: bedrock depth (Z1) effect.
        Ez1 = exp(-2.63/4.0 * alog((VS**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))
        deltaZ1 = depthvs10*1000.0 - Ez1
C     1    exp(-2.63/4.0 * alog((VS**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))
 
        if (regionflag .eq. 0) then
            term16 = 0.0
         elseif (regionflag .eq. 1) then
            term16 = phi5T*( 1.0 -exp(-deltaZ1/phi6))
        endif
  
c Sa on soil condition
        lnpsa = lnpsa_ref + term14 + term15 + term16
        sa = exp(lnpsa_ref + term14 + term15 + term16)
        psa = psa_ref * exp(term14 + term15 + term16)

c        write(*,*) "term1 = " , term1
c        write(*,*) "term2 = " , term2
c        write(*,*) "term3 = " , term3
c        write(*,*) "term4 = " , term4
c        write(*,*) "term5 = " , term5
c        write(*,*) "term6 = " , term6
c        write(*,*) "term7 = " , term7
c        write(*,*) "term8 = " , term8
c        write(*,*) "term9 = " , term9
c        write(*,*) "term10 = ", term10
c        write(*,*) "term14 = ", term14
c        write(*,*) "term15 = ", term15
c        write(*,*) "term16 = ", term16
c        write(*,*) "Ez1 = " , Ez1
c        write(*,*) "deltaZ1 = " , deltaZ1
c        write(*,*) "lnpsa_ref = ", lnpsa_ref
c        write(*,*) "lnpsa = ", lnpsa
c        write(*,*) "psa = ", psa
  
C Compute the sigma term

        tau = tau1T +(tau2T-tau1T)/1.5*(min(max(M-0.0,5.0),6.5)-5.0)

       NL0=phi2T*(exp(phi3T*(min(Vs,1130.0)-360.0))-exp(phi3T*(1130.0-360.0)))
     1    *(psa_ref/(psa_ref+phi4T))
  
       sigmaNL0 = (sigma1+(sigma2 - sigma1)/1.5*(min(max(M-0.0,5.0),6.5)-5.0))*
     1           sqrt((sigma3T*Finferred + 0.7* Fmeasured) + (1.0+NL0)**2.0)

     
C     Current code set for Measured Vs30 values (i.e., Vs30class=1)

      if (M .lt. 5.0) then
          phiss = phiss1MT
      elseif (M .le. 6.5) then
          phiss = phiss1MT + (phiss2MT-phiss1MT)*((M-5.0)/(6.5-5.5))
      else
          phiss = phiss2MT
      endif

        sigma = sqrt((1+NL0)**2.0*(tau)**2.0+sigmaNL0**2.0)

      phi = sigmaNL0

C     Convert ground motion to units of gals.
      lnY = lnpsa + 6.89
      period2 = period1

      return
      end 
 
 
c ------------------------------------------------------------------            
C *** Adjusted BCHydro model by Phung and Loh ***********
c ------------------------------------------------------------------            

      subroutine S04_PhungSub2017 ( mag, rRup, vs30, Z10, ZTor, lnY, sigma,  
     2                     specT, period2, iflag, regionflag, ftype )

      implicit none
     
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=21)
      real period(MAXPER), a5(MAXPER), a13(MAXPER), Mref(MAXPER), a2(MAXPER), a6tw(MAXPER), a12tw(MAXPER), 
     1     a6jp(MAXPER), a12jp(MAXPER), a1jp(MAXPER), a4jp(MAXPER), a1tw(MAXPER), a4tw(MAXPER),  
     1     a11sitw(MAXPER), a11sijp(MAXPER), a11sstw(MAXPER), a11ssjp(MAXPER), si(MAXPER),  
     1     ss(MAXPER), dsitw(MAXPER), dsstw(MAXPER), dsijp(MAXPER), dssjp(MAXPER), b12tw(MAXPER),  
     1     b12jp(MAXPER), tau1(MAXPER), phi1(MAXPER), phisstw(MAXPER), phis2stw(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo, mag 

      real a5T, a13T, MrefT, a2T, a6twT, a12twT, a6jpT, a12jpT, a1jpT, a4jpT, a1twT, a4twT,  
     1     a11sitwT, a11sijpT, a11sstwT, a11ssjpT, siT, ssT, dsitwT, dsstwT, dsijpT, dssjpT,  
     1     b12twT, b12jpT, tau1T, phi1T, phisstwT, phis2stwT    
      real Ez1, fz10, fmag, frup, fsite, fztor, fevt
      real period1, a3, Z10, ZTor, a9, d, b12, a1, a4, a6, a12, lnY
      integer count1, count2, iflag, regionflag
      real n, c, c4, c1, faba, R, depth, specT, tau, phi, ftype, period2


      data period  /0, 0.01, 0.02, 0.05, 0.075, 0.1, 0.15, 0.20, 0.25, 0.3, 0.4,
     1           0.5, 0.6, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5  /
      data a1jp  / 5.968739524, 5.968739524, 6.065937902, 6.883190759, 7.227465993, 7.497559544, 7.233835008, 
     1          6.935196245, 6.521314733, 6.156716362, 5.514367143, 4.98146131, 4.658399353, 4.008772718, 
     1          2.9524273, 1.828678126, 0.746782153, -0.00663011, -0.650990819, -1.431118975, -1.755209284 / 
      data a1tw  / 5.436921725, 5.481983362, 5.538686092, 6.089139138, 6.242050014, 6.500741577, 6.504662887, 
     1          6.389354121, 6.078155863, 5.827385493, 5.42827832, 5.006275132, 4.765700364, 4.170084662, 
     1          3.159305932, 2.174570827, 1.117040118, 0.315923466, -0.420579077, -1.301652283, -1.501990409 /       
      data a2  / -1.6170844, -1.6189049, -1.6263518, -1.7079854, -1.6933903, -1.7187589, -1.6686977, 
     1          -1.6457193, -1.6008945, -1.5604226, -1.5007535, -1.4451253, -1.4488396, -1.3851978, 
     1          -1.2532309, -1.1369522, -0.9895625, -0.8912422, -0.8027285, -0.7196396, -0.7258237 / 
      data a4jp  / 0.794352858, 0.792505723, 0.783529563, 0.849352353, 0.900836388, 0.872871415, 0.829902097, 
     1          0.783290162, 0.756119852, 0.730950329, 0.722697329, 0.759119408, 0.804186335, 0.969649141, 
     1          1.066170023, 1.227507791, 1.313655898, 1.354466562, 1.412705471, 1.401011805, 1.437155747 /    
      data a4tw  / 0.566548102, 0.567853502, 0.564490754, 0.51715887, 0.486404981, 0.488707369, 0.527569432, 
     1          0.553252818, 0.566033333, 0.589679636, 0.682072341, 0.767327775, 0.882682225, 1.106610265, 
     1          1.167032801, 1.394743299, 1.474103504, 1.511297382, 1.518900018, 1.460828925, 1.514666676 /    
      data a5  / 0.03849929, 0.04033665, 0.04190178, 0.04509359, 0.04623140, 0.04819708, 
     1          0.04325090, 0.03692059, 0.06597319, 0.06197944, 0.06979644, 0.04783791, 
     1          0.10612877, 0.14511870, 0.20744484, 0.24136621, 0.25767232, 0.27153377, 
     1          0.28822087, 0.32589322, 0.26383949 / 
      data a6jp  / -0.006852048, -0.006876357, -0.006874997, -0.007336676, -0.007786774, -0.007760518, 
     1          -0.007938507, -0.007685433, -0.007463007, -0.007083461, -0.006233276, -0.005626905, 
     1          -0.004826234, -0.0043472, -0.004003156, -0.003371602, -0.003456088, -0.003403299, 
     1          -0.003638611, -0.003796507, -0.003298056 /
      data a6tw  / -0.00099225, -0.00093636, -0.00088804, -0.00075076, -0.00064516, -0.00054289, 
     1          -0.00036864, -0.00022201, -0.00010404, -1.32E-05, -7.96E-05, -2.36E-05, -9.86E-06, 
     1          -6.74E-05, -0.00015376, -0.00032041, -0.00063504, -0.00097344, -0.00119025, -0.001296, -0.00100489 / 
      data a12jp  / -0.7516020, -0.7500948, -0.7307185, -0.4831132, -0.3413025, -0.4948081, -0.8669192, 
     1          -1.0634892, -1.1789740, -1.2253631, -1.2073943, -1.1299835, -1.0859780, -1.0233555, 
     1          -0.9766258, -0.9437327, -0.8880212, -0.8546554, -0.7803988, -0.6937169, -0.6499105 /    
      data a12tw  / -0.4528715, -0.4516550, -0.4403449, -0.2766783, -0.2833841, -0.3205012, -0.4471684, 
     1          -0.5552021, -0.6466667, -0.7124316, -0.7599690, -0.7702118, -0.8037457, -0.8730668, 
     1          -0.9821700, -1.0045641, -0.9337591, -0.9174852, -0.9334706, -0.8808471, -0.9343411 /    
      data a13  / -0.0256568, -0.0259617, -0.0262528, -0.0270426, -0.0276048, -0.0280794, 
     1          -0.0287650, -0.0291017, -0.0290970, -0.0287552, -0.0269993, -0.0235859, 
     1          -0.0180673, -0.0059789, -0.0031849, 0.0016328, 0.0038157, 0.0049829, 
     1          0.0063214, 0.0055941, 0.0007480 /   
      data Mref  / 7.42, 7.48, 7.48, 7.48, 7.48, 7.49, 7.50, 7.46, 7.42, 7.42, 7.34, 7.29, 
     1          7.24, 7.12, 7.12, 7.12, 7.12, 7.12, 7.12, 7.12, 7.12 /       
      data a11sitw  / 0.042901300, 0.043750861, 0.044802304, 0.049068196, 0.057175742, 0.058933025, 0.055861367, 
     1          0.050972069, 0.046435903, 0.04197234, 0.032788712, 0.020606033, 0.016652911, 0.012881203, 
     1          0.005555861, 0.001024621, -0.001859106, -0.004604082, -0.006390933, -0.006873052, -0.004145292 /        
      data a11sijp  / 0.013435789, 0.014125461, 0.014739968, 0.016119648, 0.016739089, 0.016921614, 1.55E-02, 
     1          1.38E-02, 1.31E-02, 1.28E-02, 1.21E-02, 1.08E-02, 9.87E-03, 9.18E-03, 9.08E-03, 0.008774591, 
     1          0.008850491, 0.008710436, 0.008597954, 0.008187343, 0.007571521 /   
      data a11sstw  / 0.017644058, 0.017707758, 0.017831492, 0.018544001, 0.019599185, 0.020556507, 0.020780166, 
     1          0.019426551, 0.017169719, 0.015414566, 0.01343689, 0.012308457, 0.011834766, 0.011954438, 0.01211556, 
     1          0.011352728, 0.011039715, 0.011055828, 0.011036725, 0.011031382, 0.011076041 /        
      data a11ssjp  / 0.01360108, 0.014106421, 0.014584367, 0.015798974, 0.016504644, 0.017134713, 0.01538687, 
     1          0.013705495, 0.012894641, 0.012707692, 0.012153928, 0.010798309, 0.009795145, 0.009095072, 
     1          0.009079087, 0.008626541, 0.008866849, 0.008843388, 0.008657879, 0.008229409, 0.0075585 /  
      data si  / -0.374408902, -3.75E-01, -0.376960941, -0.393806881, -0.417431719, -0.461212658, -0.47904196, 
     1          -0.436329045, -0.375236021, -0.329114488, -2.84E-01, -0.257984301, -2.22E-01, -0.16365332, 
     1          -0.095853886, -0.013344152, 0.006417125, 0.00428927, 0.004621522, -0.003842767, -0.012188147 /    
      data ss  / 0.447307171, 4.62E-01, 0.475632001, 0.512832173, 0.538728435, 0.561037907, 0.567754698, 0.531740717, 
     1          0.471888484, 0.411758139, 0.315569795, 0.252665953, 0.223480274, 0.200898826, 0.114288024, 0.01711266, 
     1          0.000179932, -0.008553978, -0.007969413, -0.004929665, 0.020215779 /  
      data dsitw  / 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25 /      
      data dsstw  / 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60 /      
      data dsijp  / 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 /      
      data dssjp  / 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 /      
      data b12tw  / -0.072698986, -0.072700471, -0.07332583, -0.10357132, -0.116533595, -0.116567166, -0.107681879, 
     1          -0.092269089, -0.072089752, -0.056308961, -0.006357987, 0.040510944, 0.071593315, 0.090540171, 
     1          0.126396911, 0.157972926, 0.160804716, 0.15178163, 0.131622296, 0.105217306, 0.103258831 /  
      data b12jp  / 0.002650526, 0.002810000, 0.002376243, 0.000548389, 0.003562477, 0.006709533, 0.007165476, 
     1          0.004341874, 0.004771685, 0.005747192, 0.005294191, 0.002624524, 0.00028237, -0.001824156, 
     1          -0.003649632, -0.005143864, -0.005623872, -0.004503918, -0.005225851, -0.006465579, -0.004203115 /    
      data tau1  / 0.359598302, 0.364722689, 0.361825448, 0.387821608, 0.410744393, 0.399061125, 0.393877013, 
     1          0.364957377, 0.388594983, 0.398189617, 0.418613789, 0.419158414, 0.395628973, 0.426568299, 0.444125779, 
     1          0.4450853, 0.439742589, 0.448589081, 0.469825447, 0.476038164, 0.490620146 /     
      data phi1  / 0.593221182, 0.594151268, 0.594258617, 0.638810416, 0.68778566, 0.704335066, 0.680147269, 
     1          0.653868401, 0.633684301, 0.613119026, 0.594846014, 0.587374669, 0.592607356, 0.589920833, 0.594625406, 
     1          0.606739824, 0.606858839, 0.604371225, 0.590755145, 0.563957325, 0.528088844 /    
      data phisstw  / 0.431332331, 0.437722489, 0.438659734, 0.443126663, 0.431998764, 0.423756788, 0.436286652, 
     1          0.44955763, 0.449672545, 0.449868674, 0.449033202, 0.44443679, 0.43483679, 0.433100211, 0.433210407, 
     1          0.412389349, 0.428153057, 0.401153057, 0.396086007, 0.335863656, 0.325693232 /        
      data phis2stw  / 0.350967141, 0.351110041, 0.350750998, 0.404608769, 0.452049131, 0.45379124, 0.416881992, 
     1          0.39160663, 0.3807666, 0.366885001, 0.360524884, 0.352567905, 0.363479045, 0.38126146, 0.369861574, 
     1          0.386518485, 0.378274223, 0.368173523, 0.3632874, 0.30811597, 0.302949265 /  
  
C Constant parameters            

      c4 = 10
      a3 = 0.1
      a9 = 0.268

C     regionflag     Note
C     -------------------------
C      0         for Japan
C      1         for Taiwan
C

C Find the requested spectral period and corresponding coefficients
      nPer = 21

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a5T = a5(i1)
         a13T = a13(i1)
         MrefT = Mref(i1)
         a2T = a2(i1)
         a6twT = a6tw(i1)
         a12twT = a12tw(i1)
         a6jpT = a6jp(i1)
         a12jpT = a12jp(i1)
         a1jpT = a1jp(i1)
         a4jpT = a4jp(i1)
         a1twT = a1tw(i1)
         a4twT = a4tw(i1)
         a11sitwT = a11sitw(i1)
         a11sijpT = a11sijp(i1)
         a11sstwT = a11sstw(i1)
         a11ssjpT = a11ssjp(i1)
         siT = si(i1)
         ssT = ss(i1)
         dsitwT = dsitw(i1)
         dsstwT = dsstw(i1)
         dsijpT = dsijp(i1)
         dssjpT = dssjp(i1)
         b12twT = b12tw(i1)
         b12jpT = b12jp(i1)
         tau1T = tau1(i1)
         phi1T = phi1(i1)
         phisstwT = phisstw(i1)
         phis2stwT = phis2stw(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Phung et al. interface (2017 Model) Horizontal'
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
 1020 call S24_interp (period(count1),period(count2),a5(count1), a5(count2),
     +                 specT, a5T, iflag)
      call S24_interp (period(count1),period(count2),a13(count1), a13(count2),
     +                 specT, a13T, iflag)
      call S24_interp (period(count1),period(count2),Mref(count1), Mref(count2),
     +                 specT, MrefT, iflag)
      call S24_interp (period(count1),period(count2),a2(count1), a2(count2),
     +                 specT, a2T, iflag)
      call S24_interp (period(count1),period(count2),a6tw(count1), a6tw(count2),
     +                 specT, a6twT, iflag)
      call S24_interp (period(count1),period(count2),a12tw(count1), a12tw(count2),
     +                 specT, a12twT, iflag)
      call S24_interp (period(count1),period(count2),a6jp(count1), a6jp(count2),
     +                 specT, a6jpT, iflag)
      call S24_interp (period(count1),period(count2),a12jp(count1), a12jp(count2),
     +                 specT, a12jpT, iflag)
      call S24_interp (period(count1),period(count2),a1jp(count1), a1jp(count2),
     +                 specT, a1jpT, iflag)
      call S24_interp (period(count1),period(count2),a4jp(count1), a4jp(count2),
     +                 specT, a4jpT, iflag)
      call S24_interp (period(count1),period(count2),a1tw(count1), a1tw(count2),
     +                 specT, a1twT, iflag)
      call S24_interp (period(count1),period(count2),a4tw(count1), a4tw(count2),
     +                 specT, a4twT, iflag)
      call S24_interp (period(count1),period(count2),a11sitw(count1), a11sitw(count2),
     +                 specT, a11sitwT, iflag)
      call S24_interp (period(count1),period(count2),a11sijp(count1), a11sijp(count2),
     +                 specT, a11sijpT, iflag)
      call S24_interp (period(count1),period(count2),a11sstw(count1), a11sstw(count2),
     +                 specT, a11sstwT, iflag)
      call S24_interp (period(count1),period(count2),a11ssjp(count1), a11ssjp(count2),
     +                 specT, a11ssjpT, iflag)
      call S24_interp (period(count1),period(count2),si(count1), si(count2),
     +                 specT, siT, iflag)
      call S24_interp (period(count1),period(count2),ss(count1), ss(count2),
     +                 specT, ssT, iflag)
      call S24_interp (period(count1),period(count2),dsitw(count1), dsitw(count2),
     +                 specT, dsitwT, iflag)
      call S24_interp (period(count1),period(count2),dsstw(count1), dsstw(count2),
     +                 specT, dsstwT, iflag)
      call S24_interp (period(count1),period(count2),dsijp(count1), dsijp(count2),
     +                 specT, dsijpT, iflag)
      call S24_interp (period(count1),period(count2),dssjp(count1), dssjp(count2),
     +                 specT, dssjpT, iflag)
      call S24_interp (period(count1),period(count2),b12tw(count1), b12tw(count2),
     +                 specT, b12twT, iflag)
      call S24_interp (period(count1),period(count2),b12jp(count1), b12jp(count2),
     +                 specT, b12jpT, iflag)
      call S24_interp (period(count1),period(count2),tau1(count1), tau1(count2),
     +                 specT, tau1T, iflag)
      call S24_interp (period(count1),period(count2),phi1(count1), phi1(count2),
     +                 specT, phi1T, iflag)
      call S24_interp (period(count1),period(count2),phisstw(count1), phisstw(count2),
     +                 specT, phisstwT, iflag)
      call S24_interp (period(count1),period(count2),phis2stw(count1), phis2stw(count2),
     +                 specT, phis2stwT, iflag)       


 1011 period1 = specT                                                                                                              

C     Regional term
      if(ftype .eq. 0.0) then 
        fevt = siT
      elseif(ftype .eq. 1.0) then 
       fevt = ssT
      endif
     
C     Ztor Scaling        
      if  (ftype .eq. 0.0 .and. regionflag .eq. 1 ) then
         d = dsitwT
         fztor = a11sitwT *(min(Ztor,40.0)-d)
      elseif (ftype .eq. 0.0 .and. regionflag .eq. 0 ) then
         d = dsijpT
         fztor = a11sijpT *(min(Ztor,40.0)-d)
      elseif (ftype .eq. 1.0 .and. regionflag .eq. 1 ) then
         d = dsstwT
         fztor = a11sstwT *(min(Ztor,80.0)-d)
      elseif (ftype .eq. 1.0 .and. regionflag .eq. 0 ) then
         d = dssjpT
         fztor = a11ssjpT*(min(Ztor,80.0)-d)
      endif
 
C  Regional term and  Basin Depth term
      if(regionflag .eq. 0) then
       
        a1 = a1jpT
        a4 = a4jpT
        a6 = a6jpT
        a12 = a12jpT

        Ez1 = exp(-5.23/2.0 * alog((vs30**2.0 + 412.39**2.0)/(1360.0**2.0 + 412.39**2.0)))
        b12 = b12jpT

      elseif(regionflag .eq. 1) then
       
        a1 = a1twT
        a4 = a4twT
        a6 = a6twT
        a12 = a12twT
    
        Ez1 = exp(-2.63/4.0 * alog((vs30**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))
        b12 = b12twT
  
      endif
      
      fz10 = b12*(min(alog(Z10*1000.0/Ez1),0.0))
      
C     Magnitude Scaling
      if (mag .le. MrefT ) then
        fmag = a4*(mag-MrefT) + a13T*(10.0-mag)**2.0
      else
        fmag = a5T*(mag-MrefT) + a13T*(10.0-mag)**2.0
      endif   
      
C     Path Scaling
       R = rRup + c4*exp( (mag-6.0)*a9 ) 
       frup = a1 + (a2T + a3*(mag - 7.8))*alog(R) + a6*rRup 
      
C     Site Effect
       fsite = a12*min(alog(vs30/760),0.0)
       

       lnSa = fmag + frup + fztor + fsite + fz10 + fevt
       
C     Set sigma values to return
       tau = tau1T
       phi = phi1T
       sigma = sqrt(tau**2+phi**2)
    
c     write(*,*) "fz10 = ", fz10
c      write(*,*) "fmag = ", fmag
c     write(*,*) "X = ", frup
c     write(*,*) "fsite = ", fsite
c     write(*,*) "fztor = ", fztor
c     write(*,*) "fevt = ", fevt
c     write(*,*) "lnSa = ", lnSa
c     write(*,*) "Sa = ", exp(lnSa)
 
C     Convert ground motion to units of gals.
      lnY = lnSa + 6.89
      period2 = period1
      return
      END

C--------------------------------------------------------------------------------------  
C     Style of Faulting Model - Taiwan SSHAC model
C     June 2017
C     megnitude dependent

      subroutine S04_Taiwan_Sof_Mag ( mag, Sof, ftype, iBranch )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=3) 
      integer nPer, count1, count2, i, iflag, iBranch
      real period(MAXPER), NM_a1(MAXPER), NM_a2(MAXPER), NM_a3(MAXPER)
      real RV_a1(MAXPER), RV_a2(MAXPER), RV_a3(MAXPER), mag, ftype
      real Sof, period1
      
      data NM_a1 / 0.03237403, 0.47049283, -0.40574477 /
      data NM_a2 /-0.1312677, 0.1482281, -0.4107635 /
      data NM_a3 / -0.1502311, 0.1306667, -0.4311289 /
      data RV_a1 / 0.08524778, 0.35955318, -0.18905762 /
      data RV_a2 / 0.1512955, 0.3856995, -0.0831085 /
      data RV_a3 / 0.1531903, 0.3856243, -0.0792437  /

c     iBranch = 1 for central branch
c     iBranch = 2 for high branch
c     iBranch = 3 for low branch

c     Compute SoF    
      If  ( ftype .lt. 0. ) then
        if (mag .gt. 6.5 ) then 
          SoF = NM_a3(iBranch)
        elseif ( mag .le. 6.5 .and. mag .gt. 5.5  ) then
          SoF = NM_a2(iBranch) + (NM_a3(iBranch)-NM_a2(iBranch))*(mag-5.5)
        elseif ( mag .le. 5.5 .and. mag .gt. 4.5 ) then
          SoF = NM_a1(iBranch) + (NM_a2(iBranch)-NM_a1(iBranch))*(mag-4.5)
        else
           SoF = NM_a1(iBranch)
       endif
      elseif ( ftype .gt. 0. ) then
        if (mag .gt. 6.5 ) then 
          SoF = RV_a3(iBranch)
        elseif ( mag .le. 6.5 .and. mag .gt. 5.5  ) then
          SoF = RV_a2(iBranch)+ (RV_a3(iBranch)-RV_a2(iBranch))*(mag-5.5)
        elseif ( mag .le. 5.5 .and. mag .gt. 4.5 ) then
          SoF = RV_a1(iBranch) + (RV_a2(iBranch)-RV_a1(iBranch))*(mag-4.5)
        else
           SoF = RV_a1(iBranch)
        endif
      Endif
   
      return
      end 

C     Magnitude Independent
   
      subroutine S04_Taiwan_Sof_MagI ( iBranch, Sof, ftype )
      
      implicit none 

      integer nPer, iflag, i1, iBranch
      real NM(3), RV(3)
      real mag, ftype
      real Sof

      
      data NM / -1.05E-01, 6.40E-02, -2.74E-01 /
      data RV /  1.48E-01, 2.94E-01,  2.74E-03 /

c     iBranch = 1 for central branch
c     iBranch = 2 for high branch
c     iBranch = 3 for low branch

c     Compute SoF    
      If  ( ftype .lt. 0. ) then
          SoF = NM(iBranch)
        elseif ( ftype .gt. 0. ) then
          SoF = RV(iBranch)
      Endif
   
      return
      end 

c-------------------- Adjusted in Taiwan SSHAC Project--------------------------------------------------  
c-------------------- C01  --------------------------------------------------  
      
c ------------------------------------------------------------------            
C *** BCHydro Subduction (06/2010 - adjustefd Model Version C ) Horizontal ***********
c ------------------------------------------------------------------            
      subroutine S04_AGA16_TW_C01 ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, Ztor, depth, disthypo )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT, Ztor
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance

c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.

C     Compute Rock PGA
      call S04_AGA16_TW_C01_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, Ztor, depth, disthypo )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S04_AGA16_TW_C01_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, Ztor, depth, disthypo )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
c ----------------------------------------------------------------------
      subroutine S04_AGA16_TW_C01_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, Ztor, depth, disthypo )

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=25)
      real a1(MAXPER), a2(MAXPER),dC1_itf(MAXPER) ,dC1_itb(MAXPER),
     1     a6(MAXPER), a10(MAXPER), a11(MAXPER), a11a(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo,
     1     mag, a3, a4, a5, a9 ,a7 , a8 , a15, a16, a6a
      real a1T, a2T, a6T, a11aT, Ztor
      real a10T, a11T, a12T, a13T, a14T, sigsT, sigtT, dC1_itfT, dC1_itbT
      real vLinT, b_soilT, sumgm, Ftype, tau, period1
      integer count1, count2, iflag
      real n, c, c4, c1, deltac1, faba, R, testmag, VsStar, depth, specT
      real base, fmag, fdepth, fsite

      Data Period(1:25) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 
     1            0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10/
       data a1 / 4.442324272, 4.454995106, 4.476989613, 4.570375455, 4.68085322, 4.795283362, 5.187708654, 
     1           5.366885037, 5.493314388, 5.606487923, 5.588152124, 5.56106773, 5.458627336, 5.256968986,  
     1           4.863221356, 4.406005361, 3.632386593, 3.163230307, 2.454565259, 1.873974482, 1.161774809,  
     1           0.622785806, 0.063732135, -1.21868348, -2.213522559 /
      data a2  / -1.35, -1.35, -1.35, -1.372125352, -1.38782354, -1.4, -1.45, -1.45, -1.45, -1.45, -1.428246273,  
     1           -1.4, -1.35, -1.28, -1.18, -1.08, -0.91, -0.85, -0.77, -0.71, -0.64, -0.58, -0.54, -0.46, -0.4 /
      data a6  / -0.000279112, -0.000236461, -0.000254929, -0.000277444, -0.000309597, -0.000328968, -0.000147104,  
     1           -0.000143273, -0.000101224, -0.000228865, -0.000398518, -0.000509299, -0.00056124, -0.00099964,  
     1           -0.001455436, -0.001931196, -0.003474407, -0.003950625, -0.004406306, -0.004095739, -0.003534331,  
     1           -0.003889444, -0.003579808, -0.002904953, -0.004265509 /
      data a10 /  2.614780773, 2.594860589, 2.601252405, 2.624840044, 2.66846963, 2.719062888, 2.82423238,  
     1           2.870663541, 2.873690705, 2.854306237, 2.690740299, 2.465778348, 2.135003365, 1.878180458,  
     1           1.517225292, 1.226658613, 0.719712006, 0.37882339, -0.02664084, -0.106219895, -0.17035158,  
     1           -0.079613443, 0.070314474, -0.07113619, 0.117820808 /
      data a11 /  0.019824826, 0.018076965, 0.01785939, 0.018011975, 0.017752615, 0.017565997, 0.017961591,  
     1           0.018186982, 0.018242659, 0.018510535, 0.018229181, 0.017858183, 0.017389833, 0.017450025,  
     1           0.018642358, 0.018912508, 0.019302661, 0.018150184, 0.015597208, 0.014257969, 0.011651582,  
     1           0.009288783, 0.008337491, 0.005121553, 0.002275984 /
      data a11a /  0.04081182, 0.036578254, 0.03708527, 0.038337152, 0.040004028, 0.041961023, 0.045225749,  
     1           0.048061124, 0.048918252, 0.047045382, 0.046243618, 0.044247831, 0.04247908, 0.038177175,  
     1           0.029427572, 0.022530308, 0.011216329, 0.004626401, -0.000948396, -0.004863106, -0.005063613, 
     1            -0.007774634, -0.004944396, -0.007533656, -0.017715245 /
      data a12  /  0.917324645, 0.918587417, 0.922959838, 1.009910498, 1.087815359, 1.150363857, 1.300591025,  
     1           1.450179555, 1.601462189, 1.776459199, 1.87491032, 2.031540841, 2.229875976, 2.379374019,  
     1           2.492260723, 2.537023323, 2.124665878, 1.57816376, 0.450025289, -0.391751696, -0.716939762,  
     1           -0.66351539, -0.670875563, -0.587002599, -0.56025643 /
      data a13  /  -0.0135, -0.0135, -0.0135, -0.013632752, -0.013726941, -0.0138, -0.0142, -0.0145, -0.014859728,  
     1           -0.0153, -0.015691567, -0.0162, -0.0172, -0.0183, -0.0206, -0.0231, -0.0296, -0.0363, -0.0493,  
     1           -0.061, -0.0798, -0.0935, -0.098, -0.098, -0.098 /
      data a14  /  -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.378246273, -0.35, -0.31, -0.28,  
     1           -0.23, -0.19, -0.12, -0.07, 0, 0, 0, 0, 0, 0, 0 /
      data Vlin /  865.1, 865.1, 865.1, 948.4683281, 1007.619098, 1053.5, 1085.7, 1032.5, 962.8476216, 877.6,  
     1           821.3013556, 748.2, 654.3, 587.1, 503, 456.6, 410.5, 400, 400, 400, 400, 400, 400, 400, 400 /
      data b_soil / -1.186, -1.186, -1.186, -1.256801128, -1.307035328, -1.346, -1.471, -1.624, -1.762045708,  
     1           -1.931, -2.042814155, -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299,  
     1           0, 0, 0, 0, 0 /
      data dC1_itf / 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.14, 0.1,  
     1           0.04, 0, -0.06, -0.1, -0.2, -0.2, -0.2, -0.2, -0.2 /
      data dC1_itb / -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3,  
     1           -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3 /
      data sigs /  0.60,0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 
     1            0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60 /
      data sigt /  0.43,0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 
     1            0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43 /


C Constant parameters            
      n = 1.18
      c = 1.88
      a3 = 0.1
      a4 = 0.9
      a5 = 0.0
      a9 = 0.4
      c4 = 10.0
      c1 = 7.8
      a7 = 0
      a8 = 0
      a15 = 0
      a16 = 0
      a6a = 0
      
C Find the requested spectral period and corresponding coefficients
      nPer = 25

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a6T = a6(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a11aT = a11a(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         dC1_itfT = dC1_itf(i1)
         dC1_itbT = dC1_itb(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
 
C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'AGA16_TW_C01 Subduction Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a11a(count1),a11a(count2),
     +                   specT,a11aT,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),dC1_itf(count1),dC1_itf(count2),
     +                   specT,dC1_itfT,iflag)
            call S24_interp (period(count1),period(count2),dC1_itb(count1),dC1_itb(count2),
     +                   specT,dC1_itbT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)

 1011 period1 = specT                                                                                                              

C     Compute the R term and base model based on either Rupture Distance 
c         (Interface events) of Hypocentral distance (Intraslab events). 
      if (ftype .eq. 0.0) then
         deltaC1 = dC1_itfT
         R = rRup + c4*exp( (mag-6.0)*a9 ) 
         base = a1T + a4*deltaC1 + (a2T + a14T*ftype + a3*(mag - 7.8))*alog(R) + a6T*rRup + a10T*ftype
      elseif (ftype .eq. 1.0) then
         deltaC1 = dC1_itbT
         R = disthypo + c4*exp( (mag-6.0)*a9 ) 
         base = a1T + a4*deltaC1 + (a2T + a14T*ftype + a3*(mag - 7.8))*alog(R) + a6T*disthypo + a10T*ftype
      else
         write (*,*) 'AGA16_TW_C01 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif
      
C     Base model for Magnitude scaling.      
      testmag = (7.8 + deltaC1)
      if (mag .le. testmag ) then
         fmag = a4*(mag-testmag) + a13T*(10.0-mag)**2.0
      else
         fmag = a5*(mag-testmag) + a13T*(10.0-mag)**2.0
      endif      
      
C     Depth Scaling
      if (ftype .eq. 0.0) then
        fdepth = a11aT*(Ztor - 20.0)
      elseif (ftype .eq. 1.0) then
C        fdepth = a11T*(depth - 60.0 )
        fdepth =  a11T*(min(depth, 80.0) -60.0 )
      else
         write (*,*) 'AGA16_TW_C01 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif

C     Site Response 
      if (vs30 .ge. 1000.0) then
          VsStar = 1000.0
      else
          VsStar = vs30
      endif
       
      if (vs30 .ge. VlinT) then
         fsite = a12T*alog(VsStar/vLinT) + b_soilT*n*alog(VsStar/vLinT)
      else
         fsite = a12T*alog(VsStar/vLinT) - b_soilT*alog(pgarock + c) +
     1          b_soilT*alog(pgarock + c*(VsStar/vlinT)**n)     
      endif

      sumgm = base + fmag + fdepth + fsite
c      write(*,*) "deltaC1 = ", deltaC1
c      write(*,*) "testmag = ", testmag
c      write(*,*) "base = ", base
c      write(*,*) "fmag = ", fmag
c      write(*,*) "fdepth = ", fdepth
c      write(*,*) "fsite = ", fsite
c      write(*,*) "lnYSa = ", sumgm
c      write(*,*) "Sa = ", exp(sumgm)
   
C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = sumgm

      return
      end

c ------------------------------------------------------------------
C *** Adjusted Lin and Lee (2008) Horizontal for Subduction Zones, ****
c ------------------------------------------------------------------

      subroutine S04_LL08_C02 ( mag, rupdist, specT, period, lnY, sigma,
     1  iflag, Ztor, ftype, vs30)

      implicit none

      integer MAXPER, nPer, i
      parameter (MAXPER=23)
      real mag, rupDist, lnY, sigma, period
      real ftype, vs30, Ztor, eq1, eq2
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c6(MAXPER), 
     1     c8(MAXPER), c9(MAXPER), c10(MAXPER), c11(MAXPER)
      real specT, c1T, c2T, c3T, c4, c5, c6T, c7, c8T, c9T, c10T, c11T, sigT
      integer count1, count2, iflag
      real period1(MAXPER), sig(MAXPER), tau(MAXPER)

      Data period1 / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 
     1      4, 5/
      Data c1 / -2.67215141251388, -2.71177472200959, -2.72221832187676, -2.39784018709008, -2.03538207060753, -1.9088731533554,
     1      -1.51924693067264, -1.22860451814747, -1.46548753025297, -1.60073705472033, -1.73022355739054, -1.92135557019542,
     1      -2.39746675094267, -3.1626558324922, -3.92684084497968, -4.66300148921766, -5.70953049123902, -7.1524588205603,
     1      -8.48237255072763, -9.52408882517922, -11.1234179760534, -12.2807829729212, -13.0556468090273/
      Data c2 / 1.205, 1.205, 1.2, 1.155, 1.1, 1.09, 1.04023471290541, 1, 1.04, 1.045, 1.065, 1.085, 1.12512345582376, 1.215,
     1      1.285, 1.365, 1.465, 1.62, 1.705, 1.77, 1.83, 1.845, 1.805/
      Data c3 / -1.905, -1.895, -1.88, -1.875, -1.86, -1.855, -1.82624150716967, -1.795, -1.77, -1.73, -1.71, -1.675,
     1      -1.61902356955068, -1.57, -1.5, -1.465, -1.45, -1.45, -1.44, -1.43, -1.37, -1.26, -1.135/
      Data c6 / 0.0156102889995883, 0.0155937815787548, 0.0153748969231095, 0.0151405607461744, 0.0148129562446649,
     1      0.0146316205990692, 0.0154910951819304, 0.0168002641290514, 0.0170361616812206, 0.017589045794533, 0.0173841662544598,
     1      0.0169136478032877, 0.0146367531026353, 0.0133391252630137, 0.0130416388541367, 0.0132968472446285, 0.0110158907384788,
     1      0.00840352664801909, 0.0052779340533018, 0.00391178925933505, -0.000827499686994609, -0.00149591256410688,
     1      0.00017541935878931/
      Data c8 / -0.0307330537090917, -0.0298343520996223, -0.0274126868068662, -0.0364602430745523, -0.0407807117320402,
     1      -0.031390816808513, -0.0287137244543388, -0.0351687921689559, -0.0323399626797574, -0.0437727054144336,
     1      -0.0498562436473804, -0.0639911797350464, -0.0669593104826402, -0.0635239987476825, -0.0947684446353348,
     1      -0.114035873731828, -0.148789496848836, -0.133337325622367, -0.139809507698123, -0.129201814425729,
     1      -0.0715007654732053, -0.0372707415321892, -0.000597669002690569/
      Data c9 / 0.0266593626202029, 0.0266334030802883, 0.0269669701408332, 0.0280072960151164, 0.0289904354524744,
     1      0.030604317861146, 0.0332015363157863, 0.0354619001161029, 0.0372753772673137, 0.0357195049604629, 0.0351815190940847,
     1      0.0338209794345063, 0.0339195582609039, 0.0320745613914819, 0.0267583477942808, 0.0213670091589494, 0.0111478707446052,
     1      0.00553780059760793, -0.0016119522124433, -0.00691562108364793, -0.0176232269444349, -0.0263408475728589,
     1      -0.0281361511687443/
      Data c10 / -0.000949647771928277, -0.00101200055609113, -0.00110255232367069, -0.0011326533367276, -0.00109413518804539,
     1      -0.00100199972311216, -0.00090996785184189, -0.00131476815827683, -0.00193723592589072, -0.00227333008442422,
     1      -0.00287937008272209, -0.00340654996752039, -0.00332677600001939, -0.00365410072006342, -0.00477028778490143,
     1      -0.00558825116399655, -0.00552557581642036, -0.00489580988712914, -0.00302787737429278, -0.00260229160594944,
     1      -0.000675577853054268, -0.00107307865445684, -0.00190354645418688/
      Data c11 / -0.359741188020741, -0.358470777504636, -0.34775073846162, -0.315351550224917, -0.254330134095103,
     1      -0.193789527625835, -0.123352653335684, -0.148748667363803, -0.211818594094982, -0.295503232941766, -0.332260362282402,
     1      -0.362579617763326, -0.406482874390383, -0.460061886261643, -0.589573384248183, -0.668328275187303, -0.828947940168923,
     1      -0.929268961596905, -0.915372798229036, -0.901918483623071, -0.851893512362305, -0.771271120200326, -0.707883845911145/
      Data tau / 0.362709384160633, 0.362228781621561, 0.360765115235053, 0.364305472555706, 0.369410732844507, 0.371386260224396,
     1      0.375711844802132, 0.378373296413057, 0.375256347952233, 0.370246674428955, 0.371237804386007, 0.365428318168295,
     1      0.380040605508848, 0.389791614996634, 0.402000747210632, 0.41366610274721, 0.42388131060612, 0.451342780155734,
     1      0.454659536359888, 0.473606264026091, 0.470378444115232, 0.470153829776495, 0.364032908959872/
      Data sig / 0.488925125190048, 0.488832360888311, 0.488325414244525, 0.491514552262097, 0.504342137392063,
     1      0.524922820418443, 0.54658776104583, 0.539020452290274, 0.532707465059637, 0.520338086029458, 0.528640131814135,
     1      0.52886057587926, 0.526673481453196, 0.527190629182218, 0.537849220431881, 0.542923431321146, 0.533182979055384,
     1      0.520721584299849, 0.532634760118896, 0.54532659397375, 0.53310130534746, 0.49937482519724, 0.453578088334651/

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
         c10T    = c10(1)
         c11T    = c11(1)
         sigT    = sig(1)
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
      call S24_interp (period1(count1),period1(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period1(count1),period1(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period1(count1),period1(count2),sig(count1),sig(count2),
     +             specT,sigT,iflag)

 1011 period = specT

C     Compute the ground motions.
      c4 = 0.51552
      c5 = 0.63255
      c7 = 0.275
      
      if (mag .lt. 7.2) then
       eq1 = c8T*(mag - 7.2)**2.0
      else
        eq1 = 0
      endif
      
      if (ftype .eq. 0) then
       eq2 = c9T*(min(Ztor, 30.0) - 20.0)
      else
        eq2 = c6T*(min(Ztor, 80.0) - 20.0)
      endif
      
      lnY = c1T + c2T*mag + eq1 + c3T*alog(Rupdist+c4*exp(c5*mag)) + 
     1      c7*ftype + eq2 + c10T*Rupdist + c11T*alog(Vs30/760.0)


      sigma = sigT

C     Now convert to Ln Units in gals.
      lnY = lnY + 6.89

      return
      end



   
c ---------------------------------------------------------------------            
C     *** Akkar, Sandikkaya, and Bommer (2013) *** Adjusted in Taiwan SSHAC Project
c ---------------------------------------------------------------------            
      subroutine S04_ASB14_TW_C01 ( mag, Rbjf, specT, 
     1                     period2, lnY, sigma, iflag, ftype, Vs, phiT, tauT ) 
      implicit none
      integer MAXPER
      parameter (MAXPER=22)
      REAL Period(MAXPER), a1(MAXPER), a3(MAXPER), a4(MAXPER), a8(MAXPER)
      Real a9(MAXPER), b1(MAXPER), b2(MAXPER), phi(MAXPER), tau(MAXPER)
      real specT, a1T, a3T, a4T, a8T, a9T, b1T, b2T, phiT, tauT, period2
      real mag, Rbjf, Ftype, Fn, Fr, Vs, lnY
      real a5, a6, a7, c1, c, n, sigma, period1, pgaref
      INTEGER iFlag, count1, count2, nPer, i
      real a2(MAXPER), a2T
 
 
      Data Period(1:22) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2,
     1      3, 4/
      Data a1(1:22) / 1.95754376344839, 2.00817307775539, 2.10206659629286, 2.20961153930911, 2.35078005202385, 2.49873114534198,
     1      2.81381004271636, 3.12694681898327, 3.21378881748549, 3.32779032291873, 3.24345450186723, 3.12402349228009,
     1      2.91284344812683, 2.72736664472706, 2.24534448743599, 1.99211376932964, 1.1530939925732, 0.641079591580465,
     1      -0.0254787065133857, -0.511129081837452, -1.18463752251569, -1.73520928065602/
      Data a3(1:22) / -0.02807, -0.0274, -0.02715, -0.02403, -0.01797, -0.01248, -0.00532, -0.00925, -0.0149516924364565, -0.02193,
     1      -0.0274510958074301, -0.03462, -0.0467825076620119, -0.05672, -0.07684, -0.0949, -0.12347, -0.14345, -0.17187,
     1      -0.19029, -0.21392, -0.23848/
      Data a4(1:22) / -1.23452, -1.23698, -1.25363, -1.27525, -1.30123, -1.32632, -1.35722, -1.38182, -1.37407684986153, -1.3646,
     1      -1.33160829826025, -1.28877, -1.22380239685518, -1.17072, -1.0653, -1.01909, -0.88393, -0.81838, -0.75751, -0.72033,
     1      -0.69085, -0.66482/
      Data a8(1:22) / -0.103636044256603, -0.132928920885736, -0.12986863151264, -0.129851892594279, -0.116141160848401,
     1      -0.119177416631283, -0.139260487428174, -0.171750171226412, -0.157099766506712, -0.120064436335966, -0.10406684558815,
     1      -0.104366241763941, -0.0899951148151869, -0.0994122752510732, -0.144947920374219, -0.196997416408109,
     1      -0.222073676123046, -0.25282351912814, -0.234590527565048, -0.232998579238557, -0.0848380615438779, 0.0861377135345235/
      Data a9(1:22) / 0.126400204982005, 0.0932411706137936, 0.0919391699827972, 0.0864920444615963, 0.0846063286281077,
     1      0.0774891547676303, 0.0493136670345936, 0.0343857396319973, 0.0339561453478299, 0.0418349456816038, 0.0612488290723383,
     1      0.0952480389354383, 0.122571526364686, 0.122639270021399, 0.13891186952902, 0.144080973849871, 0.164385885840223,
     1      0.17991213158502, 0.163888527122665, 0.144017834272244, 0.108612885472365, 0.0432613427525789/
      Data b1(1:22) / -0.511626065987955, -0.53046724276345, -0.52254494554547, -0.50864646002686, -0.485614368483321,
     1      -0.453434312321338, -0.407260400776915, -0.434402182634647, -0.46871799480745, -0.50684048713477, -0.53389023952211,
     1      -0.553869172268956, -0.561144973177327, -0.608010203560465, -0.666054146807456, -0.700286156154561, -0.855477635291906,
     1      -0.966259736674558, -1.04096576213125, -1.02116484088795, -1.02432695483163, -1.01420693530889/
      Data b2(1:22) / -0.28846, -0.28685, -0.28241, -0.26842, -0.24759, -0.22385, -0.17525, -0.29293, -0.339056152218589, -0.39551,
     1      -0.417668345900112, -0.44644, -0.452416689285495, -0.4573, -0.43008, -0.37408, -0.28957, -0.28702, -0.24695, -0.17336,
     1      -0.13336, -0.07749/
      Data phi(1:22) / 0.332271753816561, 0.331853645328173, 0.332226638023841, 0.340352054536531, 0.356032340893303,
     1      0.375882001454717, 0.414401540225185, 0.425576834343117, 0.421822057126099, 0.401352020442047, 0.390089235488791,
     1      0.370563407684413, 0.347128326879158, 0.339900164897318, 0.32815538601622, 0.317687018977947, 0.325315598182983,
     1      0.327742581374142, 0.354898012625689, 0.371465570115995, 0.390612539946501, 0.395777102651996/
      Data tau(1:22) / 0.350418686172617, 0.347748165151219, 0.350377709334341, 0.358987239286868, 0.369064296128902,
     1      0.378518325684224, 0.396400441287509, 0.396733147119714, 0.390535011398694, 0.386521582565263, 0.384648579313126,
     1      0.375129445685108, 0.362144468664536, 0.366581082849823, 0.386276478084319, 0.406760171807368, 0.436916198216006,
     1      0.454080668424532, 0.4699258063812, 0.48322880880426, 0.506656166693695, 0.521040518558988/
       Data a2(1:22) / 0.050335468, 0.106913336, 0.106685064, 0.106873132, 0.117318196, 0.111561424, 0.10736431 , 0.103185029,
     1       0.095335574, 0.131059036, 0.127769367, 0.131100025, 0.131821032, 0.153387307, 0.129558757, 0.125270741,
     1       0.116455398, 0.143364126, 0.146078527, 0.131832393, 0.128727888, 0.007683513 /

 
 
C First check for the PGA case (i.e., specT=0.0) 
      nPer = 22
      if (specT .eq. 0.0) then
         period1 = period(1)
         a1T = a1(1)
         a2T = a2(1)
         a3T = a3(1)
         a4T = a4(1)
         a8T = a8(1)
         a9T = a9(1)
         b1T = b1(1)
         b2T = b2(1)
         phiT = phi(1)
         tauT = tau(1)
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
      write (*,*) 'Akkar,Sandikkaya&Bommer (2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +                   specT,a9T,iflag)
            call S24_interp (period(count1),period(count2),b1(count1),b1(count2),
     +                   specT,b1T,iflag)
            call S24_interp (period(count1),period(count2),b2(count1),b2(count2),
     +                   specT,b2T,iflag)
            call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +                   specT,phiT,iflag)
            call S24_interp (period(count1),period(count2),tau(count1),tau(count2),
     +                   specT,tauT,iflag)
 1011 period1 = specT                                                                                                              
C.....Set the mechanism terms based on ftype............
C     Set mechanism term and corresponding Frv and Fnm values.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                    -120 < Rake <  -60
C     -0.5      Normal/Oblique            -150 < Rake < -120
C                                          -60 < Rake <  -30
C       0       Strike-Slip               -180 < Rake < -150
C                                          -30 < Rake <   30
C                                          150 < Rake <  180
C      0.5      Reverse/Oblique             30 < Rake <   60
C                                          120 < Rake <  150
C       1       Reverse                     60 < Rake <  120   
      if (ftype .eq. -1.0) then
         Fr = 0.0
         Fn = 1.0
      elseif (ftype .eq. -0.5) then 
         Fr = 0.0
         Fn = 1.0
      elseif (ftype .eq. 0.0) then 
         Fr = 0.0
         Fn = 0.0
      elseif (ftype .eq. 0.5) then
         Fr = 1.0
         Fn = 0.0
      elseif (ftype .eq. 1.0) then
         Fr = 1.0
         Fn = 0.0
      endif 
C     Set frequency independent terms
c      a2 = 0.0029
      a5 = 0.2529
      a6 = 7.5
      a7 = -0.5096
      c1 = 6.75
      c = 2.5
      n = 3.2
C     Compute the PGA for reference Vs=750m/s.
      if (mag .lt. c1 ) then
         pgaref = a1(1) + a2(1)*(mag-c1) + a3(1)*(8.5-mag)**2.0 + 
     1                 (a4(1)+a5*(mag-c1))*alog(sqrt(Rbjf*Rbjf+a6*a6)) + 
     2                  a8(1)*Fn + a9(1)*Fr
      else
         pgaref = a1(1) + a7*(mag-c1) + a3(1)*(8.5-mag)**2.0 + 
     1                 (a4(1)+a5*(mag-c1))*alog(sqrt(Rbjf*Rbjf+a6*a6)) + 
     2                  a8(1)*Fn + a9(1)*Fr      
      endif
      pgaref = exp(pgaref) 
C.....Now compute the ground motion value........
      if (mag .lt. c1 ) then
         lnY = a1T + a2T*(mag-c1) + a3T*(8.5-mag)**2.0 + 
     1                 (a4T+a5*(mag-c1))*alog(sqrt(Rbjf*Rbjf+a6*a6)) + 
     2                  a8T*Fn + a9T*Fr
      else
         lnY = a1T + a7*(mag-c1) + a3T*(8.5-mag)**2.0 + 
     1                 (a4T+a5*(mag-c1))*alog(sqrt(Rbjf*Rbjf+a6*a6)) + 
     2                  a8T*Fn + a9T*Fr      
      endif
C.....Now apply site amplification term......
      if (vs .le. 750.0) then
         lnY = lnY + b1T*alog(Vs/750.0) + 
     1         b2T*alog( (pgaref + c*(Vs/750.0)**n) / ((pgaref+c)*(Vs/750.0)**n) )
      else
         lnY = lnY + b1T*alog( min(Vs,1000.0)/750.0)
      endif
C.....Set Sigma value.........
      sigma = sqrt (phiT*phiT + tauT*tauT)
C     Convert ground motion to units of gals.
      lnY = lnY + 6.89
      period2 = period1
      return
      END
      
c ------------------------------------------------------------------            
C *** Abrahamson, Silva, and Kamai (NGA-West2 2013) Horizontal ****
C     Earthquake Spectra Paper:
C        Summary of the Abrahamson, Silva, and Kamai NGA-West2 
C            Ground-Motion Relations for Active Crustal REgions
C         N. A. Abrahamson, S. J. Silva, and R. Kamai
C     Notes:
C        Applicable Range (see Abstract):  
C           3 <= M <= 8.5
C           Rrup <= 300 km
C        Regional attenuation included based on Regionflag
C             0 = Global
C             1 = Taiwan
C             2 = China
C             3 = Japan
C         Mainshock and Aftershocks included based on MSASFlag
C             0 = Mainshocks
C             1 = Aftershocks
C         Sigma dependent on estimated or measured Vs30m based on 
C             Vs30_Class
C             0 = Estimated Vs30m
C             1 = Measured Vs30m
c ------------------------------------------------------------------            
      subroutine S04_ASK14_TW_C01 ( mag, dip, fType, fltWidth, rRup, Rjb,  
     1                     vs30, hwflag, lnY, sigma, specT, period2, ztor,
     2                     iflag, vs30_class, z10, Rx, Ry0, regionflag, msasflag,
     1                     phi, tau )
C     Last Updated: 8/1/13
      implicit none
 
      real mag, dip, fType, rRup, rjb, Rx, Ry0, vs30, SA1180,
     1      Z10,  ZTOR, fltWidth, lnSa, sigma, lnY, vs30_rock
      real Fn, Frv, specT, period2, CRjb, phi, tau, z10_rock, SA_rock
      integer hwflag, iflag, vs30_class, regionflag, msasflag
 
c     Vs30 class is to distinguish between the sigma if the Vs30 is measured
c     vs the VS30 being estimated from surface geology.
c         Vs30_class = 0 for estimated
c         Vs30_class = 1 for measured 
 
C     Current version is not programmed for Aftershock cases. 
C       For implementation of Aftershock a new distance metric, CRjb
C       will need to be computed and passed along to this subroutine. 
 
      CRjb = 999.9 
 
C     Set mechanism term and corresponding Frv and Fnm values.     
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
C 
      if ( fType .eq. 1.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .eq. 0.5 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .eq. -1.0 ) then
        Frv = 0.0
        Fn = 1.0
      elseif ( fType .eq. -0.5 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
 
c     Compute SA1180
      vs30_rock = 1180.
      z10_rock = 0.005
      SA_rock = 0.
      
      call S04_ASK14_TW_C01_model ( mag, dip, fltWidth, ZTOR, Frv, Fn, rRup, rjb, rx, Ry0, 
     1                     vs30_rock, SA_rock, Z10_rock, hwflag, vs30_class,
     2                     specT, lnSa, phi, tau, iflag, regionflag, msasflag, CRjb )
      Sa1180 = exp(lnSa)
 
c     Compute Sa at spectral period for given Vs30
 
      call S04_ASK14_TW_C01_model ( mag, dip, fltWidth, ZTOR, Frv, Fn, rRup, rjb, rx, Ry0, 
     1                     vs30, SA1180, Z10, hwflag, vs30_class,
     2                     specT, lnSa, phi, tau, iflag, regionflag, msasflag, CRjb )
 
c     compute Sa (given the PGA rock value)
      sigma = sqrt( phi**2 + tau**2 )
 
      lnY = lnSa + 6.89
 
      period2 = specT
 
      return
      end
 
c ----------------------------------------------------------------------
      subroutine S04_ASK14_TW_C01_model ( mag, dip, FltWidth, ZTOR, Frv, Fn, rRup, rjb, Rx, Ry0, 
     1                     vs30, Sa1180, Z1, hwflag, vs30_class,
     3                     specT, lnSa, phi, tau, iflag, regionflag, msasflag, CRjb)
 
      implicit none
      
      integer MAXPER     
      parameter (MAXPER=25)
      real Vlin(MAXPER), b(MAXPER), c4(MAXPER), M1(MAXPER), a1(MAXPER)
      real a2(MAXPER), a3(MAXPER), a6(MAXPER), a8(MAXPER), a10(MAXPER)
      real a11(MAXPER), a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER)
      real a17(MAXPER), a43(MAXPER), a44(MAXPER), a45(MAXPER), a46(MAXPER)
      real a25(MAXPER), a28(MAXPER), a29(MAXPER), a31(MAXPER), a36(MAXPER)
      real a37(MAXPER), a38(MAXPER), a39(MAXPER), a40(MAXPER), a41(MAXPER), a42(MAXPER)
      real s1est(MAXPER), s2est(MAXPER), s1msr(MAXPER), s2msr(MAXPER)
      real s3(MAXPER), s4(MAXPER), s5(MAXPER), s6(MAXPER), period(MAXPER)
      real a4(MAXPER), a5, a7
      real VlinT, bT, c4T, M1T, a1T
      real a2T, a3T, a6T, a8T, a10T, a4T
      real a11T, a12T, a13T, a14T, a15T
      real a17T, a43T, a44T, a45T, a46T
      real a25T, a28T, a29T, a31T, a36T
      real a37T, a38T, a39T, a40T, a41T, a42T
      real s1estT, s2estT, s1msrT, s2msrT
      real s3T, s4T, s5T, s6T, c4_mag
      real phiA_est, phiA_msr, period1
      
      real M2
      real lnSa, SA1180, rjb, rRup, Rx, Ry0, dip, mag, vs30
      real HW_taper1, HW_taper2, HW_taper3, HW_taper4, HW_taper5
      real damp_dSA1180, sigAmp, fltWidth
      real f1, f4, f5, f6, f7, f8, f10, f11, fReg, f12, f13
      real Ry1, ZTOR, Frv, Fn, SpecT
      real phiA, phiB, tauA, tauB, phi, tau
      integer vs30_class, hwflag, iflag, nPer, regionflag, msasflag
      real n, c, z1, z1_ref
      real R, V1, Vs30Star, hw_a2, h1, h2, h3, R1, R2, CRjb
      integer count1, count2, i
      real y1, y2, x1, x2, y1z, y2z, x1z, x2z
 
 
      Data Period(1:25) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2,
     1      3, 4, 5, 7.5, 10/
      Data Vlin(1:25) / 660, 660, 680, 770, 851.659765220817, 915, 960, 910, 833.557751246246, 740, 674.738820243143, 590, 495,
     1      430, 360, 340, 330, 330, 330, 330, 330, 330, 330, 330, 330/
      Data b(1:25) / -1.47, -1.47, -1.459, -1.39, -1.2936977941189, -1.219, -1.152, -1.23, -1.39052872238288, -1.587,
     1      -1.77190667597776, -2.012, -2.411, -2.757, -3.278, -3.599, -3.8, -3.5, -2.4, -1, 0, 0, 0, 0, 0/
      Data c4(1:25) / 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     1      4.5, 4.5, 4.5/
      Data M1(1:25) / 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75,
     1      6.75, 6.75, 6.82, 6.92, 7, 7.15, 7.25/
      Data a1(1:25) / 0.593866860165754, 0.593866860165754, 0.607350205626504, 0.623596721976998, 0.652048054162018,
     1      0.706900122976899, 0.977092520718666, 1.19353371909298, 1.32002501982816, 1.41943704493427, 1.45760409041355,
     1      1.48593774481291, 1.5482692399218, 1.59530501168755, 1.52577976900992, 1.42330982305614, 1.28852744087696,
     1      1.21763682762334, 0.926314398721161, 0.752132501511978, 0.306567902675289, -0.0776771579300297, -0.336929967399733,
     1      -1.3124274646718, -1.96852544897435/
      Data a2(1:25) / -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79,
     1      -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.765, -0.634, -0.529/
      Data a3(1:25) / 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275,
     1      0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275/
      Data a4(1:25) / 0.028559622, 0.028559622, 0.024114472, 0.029331155, 0.019548333, 0.004946714, 0.005303378, -0.005794728,
     1       -0.013326944, -0.023008825, -0.033586465, -0.066345205, -0.057650858, -0.029948392, -0.002088881, 0.0203011,
     1       0.031432073, 0.135117619, 0.263006151, 0.363764714, 0.299931624, 0.171388832, 0.125319644, 0.000172165, 0.001611885/
      Data a6(1:25) / 1.0971052478214, 1.0971052478214, 1.07034252625442, 1.01191011946962, 0.960151214531469, 0.895279929309867,
     1      0.92654506977541, 1.00221368734431, 1.08991983886901, 1.22240081684137, 1.29468215932786, 1.40361090204586,
     1      1.51069901081663, 1.63441282368804, 1.74891808838416, 1.86175559265212, 2.18640072401506, 2.42168821321506,
     1      2.68100338021949, 2.76037541026386, 2.74674409095451, 2.81072195031604, 2.96079416196568, 2.54253432152637,
     1      1.22631618381355/
      Data a8(1:25) / -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.0181476220075075, -0.022,
     1      -0.0254805962536991, -0.03, -0.038, -0.045, -0.055, -0.065, -0.095, -0.11, -0.124, -0.138, -0.172, -0.197, -0.218,
     1      -0.255, -0.285/
      Data a10(1:25) / 1.88550206478838, 1.88550206478838, 1.87041425559732, 1.75625542618031, 1.60356687034528, 1.49217324726696,
     1      1.39247214525391, 1.4779307664696, 1.62803769217999, 1.82902654690301, 2.07844213755755, 2.42387695711654,
     1      3.05112199129081, 3.36681275006325, 4.14386296817676, 4.64423228581923, 4.84231546795258, 4.30714161619116,
     1      2.63511186125103, 0.617331058415303, -0.758430542317734, -0.737742033861981, -0.713264155594585, -0.777366444748309,
     1      -0.79623295261254/
      Data a11(1:25) / 0.0865290017628497, 0.0865290017628497, 0.0863920308121283, 0.0833475966285011, 0.0805940005440329,
     1      0.0739970397632909, 0.0317420642709333, 0.00950566126829064, 0.0278573208474259, 0.0617341652642051,
     1      0.0765627168872937, 0.100927646749515, 0.133863284021228, 0.158502746654088, 0.174755835638996, 0.189203321064511,
     1      0.196322660498526, 0.194565992382258, 0.163916257091526, 0.138191034588059, 0.119417922134445, 0.0518858982566291,
     1      -0.000299209835122006, -0.126816643183004, -0.0952237645759153/
      Data a12(1:25) / -0.0851813014732252, -0.0851813014732252, -0.0811946125443235, -0.0736733480467424, -0.0447695469303866,
     1      -0.0326016592728057, -0.057299086049061, -0.109198803849594, -0.0892438233690082, -0.0570349718601036,
     1      -0.0425580668821486, -0.0512274722133986, -0.0561371570080197, -0.0617139982746329, -0.091740938717726,
     1      -0.151914306681818, -0.178531415315115, -0.210716153796932, -0.216476851197754, -0.22661557721623, -0.0950274132257672,
     1      0.129608160368009, -0.154555139232628, 0.0114922876378795, -0.0640609795492863/
      Data a13(1:25) / 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.58, 0.56, 0.53, 0.5, 0.42, 0.35,
     1      0.2, 0, 0, 0, 0/
      Data a14(1:25) / -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.24, -0.19, -0.11, -0.04, 0.07,
     1      0.15, 0.27, 0.35, 0.46, 0.54, 0.61, 0.72, 0.8/
      Data a15(1:25) / 1.35951090458584, 1.35951090458584, 1.36654581947039, 1.40869323733483, 1.48532549730173, 1.57986239181509,
     1      1.739195456708, 1.79908761029304, 1.77441289606063, 1.69212932619259, 1.672990391914, 1.55489999718983,
     1      1.39852777161035, 1.3162179619873, 1.2259862575015, 1.11569501263915, 0.991158857090293, 0.958925297374564,
     1      0.810318882076326, 0.769601371751858, 0.378047007927471, 0.227924124763207, -0.0330462634934479, -0.0979718229415627,
     1      -0.181661927158766/
      Data a17(1:25) / -0.00781955144257678, -0.00781955144257678, -0.00794388960832398, -0.00861792272478135,
     1      -0.00918743653244077, -0.0098216355840628, -0.010534120645979, -0.0103669729589076, -0.0100175235859223,
     1      -0.00894850816177607, -0.008255299022226, -0.00709893304178311, -0.00581931156279194, -0.00505808157463984,
     1      -0.00368321099920911, -0.0028063952303506, -0.00289623937801703, -0.00387904788105345, -0.00384632217298883,
     1      -0.00493839857385291, -0.00347153307789942, -0.00226353012210327, -0.00179670583741962, -0.00252889897970746,
     1      -0.00310008394168946/
      Data a25(1:25) / -0.0015, -0.0015, -0.0015, -0.0016, -0.00182526831785053, -0.002, -0.0027, -0.0033, -0.00338993205735736,
     1      -0.0035, -0.00341298509365752, -0.0033, -0.0029, -0.0027, -0.0023, -0.002, -0.001, -5e-04, -4e-04, -2e-04, 0, 0, 0, 0,
     1      0/
      Data a28(1:25) / 0.0025, 0.0025, 0.0024, 0.0023, 0.00252526831785053, 0.0027, 0.0032, 0.0036, 0.00346510191396396, 0.0033,
     1      0.00303895528097257, 0.0027, 0.0024, 0.002, 0.001, 8e-04, 7e-04, 7e-04, 6e-04, 3e-04, 0, 0, 0, 0, 0/
      Data a29(1:25) / -0.0034, -0.0034, -0.0033, -0.0034, -0.00334368292053737, -0.0033, -0.0029, -0.0025, -0.0025, -0.0025,
     1      -0.00276104471902743, -0.0031, -0.0036, -0.0039, -0.0048, -0.005, -0.0041, -0.0032, -0.002, -0.0017, -0.002, -0.002,
     1      -0.002, -0.002, -0.002/
      Data a31(1:25) / -0.1503, -0.1503, -0.1479, -0.1447, -0.137885633385021, -0.1326, -0.1353, -0.1128, -0.0448563306665158,
     1      0.0383, 0.0553549216431254, 0.0775, 0.0741, 0.2548, 0.2136, 0.1542, 0.0787, 0.0476, -0.0163, -0.1203, -0.2719, -0.2958,
     1      -0.2718, -0.14, -0.0216/
      Data a36(1:25) / 0.265, 0.265, 0.255, 0.249, 0.222530972652563, 0.202, 0.126, 0.022, -0.049046325312313, -0.136,
     1      -0.110765677160682, -0.078, 0.037, -0.091, 0.129, 0.31, 0.505, 0.358, 0.131, 0.123, 0.109, 0.135, 0.189, 0.15, 0.092/
      Data a37(1:25) / 0.337, 0.337, 0.328, 0.32, 0.302541705366584, 0.289, 0.275, 0.256, 0.213731933042042, 0.162,
     1      0.188974620966168, 0.224, 0.248, 0.203, 0.232, 0.252, 0.208, 0.208, 0.108, 0.068, -0.023, 0.028, 0.031, -0.07, -0.159/
      Data a38(1:25) / 0.188, 0.188, 0.184, 0.18, 0.172678779669858, 0.167, 0.173, 0.189, 0.15257751677027, 0.108,
     1      0.111045521721987, 0.115, 0.122, 0.096, 0.123, 0.134, 0.129, 0.152, 0.118, 0.119, 0.093, 0.084, 0.058, 0, -0.05/
      Data a39(1:25) / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      Data a40(1:25) / 0.088, 0.088, 0.088, 0.093, 0.115526831785053, 0.133, 0.186, 0.16, 0.118631253615615, 0.068,
     1      0.0592985093657524, 0.048, 0.055, 0.073, 0.143, 0.16, 0.158, 0.145, 0.131, 0.083, 0.07, 0.101, 0.095, 0.151, 0.124/
      Data a41(1:25) / -0.196, -0.196, -0.194, -0.175, -0.127130482456762, -0.09, 0.09, 0.006, -0.0668449664594602, -0.156,
     1      -0.207338794742061, -0.274, -0.248, -0.203, -0.154, -0.159, -0.141, -0.144, -0.126, -0.075, -0.021, 0.072, 0.205,
     1      0.329, 0.301/
      Data a42(1:25) / 0.044, 0.044, 0.061, 0.162, 0.324756359647008, 0.451, 0.506, 0.335, 0.146592339836334, -0.084,
     1      -0.124897005980964, -0.178, -0.187, -0.159, -0.023, -0.029, 0.061, 0.062, 0.037, -0.143, -0.028, -0.097, 0.015, 0.299,
     1      0.243/
      Data a43(1:25) / -0.122962023498924, -0.122962023498924, -0.120092998613746, -0.126886260792928, -0.126035664024335,
     1      -0.0904729506991718, -0.0149049072273536, 0.0103808659653804, -0.0154811955589853, -0.0841622352036735,
     1      -0.119282978541206, -0.12421065727408, -0.107040732978698, -0.0960315277316025, -0.157890828883787, -0.19410442396108,
     1      -0.211688638835866, -0.182842209269601, 0.0357796738204077, 0.0978715683661098, 0.180188769627975, 0.267208253179032,
     1      0.288512616632419, 0.336857781708233, 0.402190602744037/
      Data a44(1:25) / 0.172338591249922, 0.172338591249922, 0.170022426697282, 0.17246825329291, 0.177160008120209,
     1      0.183647781047614, 0.188806700354148, 0.243698862114109, 0.222677566880558, 0.237785176023185, 0.237613022368251,
     1      0.218308740069599, 0.2032746665851, 0.181666749213636, 0.161804828013144, 0.122191096979769, 0.0662625613004855,
     1      0.0475907780093519, 0.069351375345784, 0.129000718045731, 0.190504182454582, 0.170342502732787, 0.157276582209913,
     1      0.137031837000311, 0.0589963511005088/
      Data a45(1:25) / 0.0614738216205322, 0.0614738216205322, 0.0583202285935086, 0.0619929094404484, 0.0716204505668039,
     1      0.0735921825736064, 0.0738034667488054, 0.064417927749121, 0.0598856791348016, 0.0357966722878199, 0.0346058994769634,
     1      0.0318550184537276, 0.0222552236580935, 0.023693014881416, 0.0490897682923249, 0.105075768222844, 0.180265627533596,
     1      0.222035554590515, 0.260100426376672, 0.254045901077365, 0.236586467752453, 0.214938951003774, 0.178100653589404,
     1      0.12514948661591, 0.116411079489916/
      Data a46(1:25) / 0.0626052993621598, 0.0626052993621598, 0.058744653649071, 0.0544910136649466, 0.0302800958064618,
     1      0.00511890031399667, -0.0298335598698967, -0.0480238700565291, -0.0326498810238118, 0.0142571550229835,
     1      0.0402025828982244, 0.0551568774847881, 0.0972265261342537, 0.151626103452716, 0.233058972791216, 0.265871303338564,
     1      0.297582627738386, 0.289181720276543, 0.235785512188508, 0.27699299739719, 0.249507530382689, 0.219942760963378,
     1      0.254640765867001, 0.162535435459097, 0.0900289349548069/
      Data s1est(1:25) / 0.754, 0.754, 0.76, 0.781, 0.797331953044163, 0.81, 0.81, 0.81, 0.805953057418919, 0.801,
     1      0.795779105619451, 0.789, 0.77, 0.74, 0.699, 0.676, 0.631, 0.609, 0.578, 0.555, 0.548, 0.527, 0.505, 0.457, 0.429/
      Data s2est(1:25) / 0.52, 0.52, 0.52, 0.52, 0.525631707946263, 0.53, 0.54, 0.55, 0.554496602867868, 0.56, 0.562175372658562,
     1      0.565, 0.57, 0.58, 0.59, 0.6, 0.615, 0.63, 0.64, 0.65, 0.64, 0.63, 0.63, 0.63, 0.63/
      Data s1msr(1:25) / 0.741, 0.741, 0.747, 0.769, 0.785331953044163, 0.798, 0.798, 0.795, 0.785107473690691, 0.773,
     1      0.764298509365752, 0.753, 0.729, 0.693, 0.644, 0.616, 0.566, 0.541, 0.506, 0.48, 0.472, 0.447, 0.425, 0.378, 0.359/
      Data s2msr(1:25) / 0.501, 0.501, 0.501, 0.501, 0.50719487874089, 0.512, 0.522, 0.527, 0.523402717705706, 0.519,
     1      0.516824627341438, 0.514, 0.513, 0.519, 0.524, 0.532, 0.548, 0.565, 0.576, 0.587, 0.576, 0.565, 0.568, 0.575, 0.585/
      Data s3(1:25) / 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47,
     1      0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47/
      Data s4(1:25) / 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36,
     1      0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36/
      Data s5(1:25) / 0.54, 0.54, 0.54, 0.55, 0.555631707946263, 0.56, 0.57, 0.57, 0.574496602867868, 0.58, 0.584350745317124,
     1      0.59, 0.61, 0.63, 0.66, 0.69, 0.73, 0.77, 0.8, 0.8, 0.8, 0.76, 0.72, 0.67, 0.64/
      Data s6(1:25) / 0.63, 0.63, 0.63, 0.63, 0.641263415892527, 0.65, 0.69, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.69,
     1      0.68, 0.66, 0.62, 0.55, 0.52, 0.5, 0.5, 0.5/
 
 
C Find the requested spectral period and corresponding coefficients
      nPer = 25
 
C First check for the PGA, PGV, PGD cases 
      if (specT .eq. 0.0) then
         period1 = period(1)
         a1T = a1(1)
         a2T = a2(1)
         a3T = a3(1)
         a4T = a4(1)
         a6T = a6(1)
         a8T = a8(1)
         M1T = M1(1)
         a10T = a10(1)
         a11T = a11(1)
         a12T = a12(1)
         a13T = a13(1)
         a14T = a14(1)
         a15T = a15(1)
         a17T = a17(1)
         a43T = a43(1)
         a44T = a44(1)
         a45T = a45(1)
         a46T = a46(1)
         a25T = a25(1)
         a28T = a28(1)
         a29T = a29(1)
         a31T = a31(1)
         a36T = a36(1)
         a37T = a37(1)
         a38T = a38(1)
         a39T = a39(1)
         a40T = a40(1)
         a41T = a41(1)
         a42T = a42(1)
         VlinT = Vlin(1)
         bT = b(1)
         c4T = c4(1)
         s1estT = s1est(1)
         s2estT = s2est(1)
         s1msrT = s1msr(1)
         s2msrT = s2msr(1)
         s3T = s3(1)
         s4T = s4(1)
         s5T = s5(1)
         s6T = s6(1)
         goto 1011
      endif
C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Abrahamson, Silva, and Kamai (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),M1(count1),M1(count2),
     +                   specT,M1T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a17(count1),a17(count2),
     +                   specT,a17T,iflag)
            call S24_interp (period(count1),period(count2),a43(count1),a43(count2),
     +                   specT,a43T,iflag)
            call S24_interp (period(count1),period(count2),a44(count1),a44(count2),
     +                   specT,a44T,iflag)
            call S24_interp (period(count1),period(count2),a45(count1),a45(count2),
     +                   specT,a45T,iflag)
            call S24_interp (period(count1),period(count2),a46(count1),a46(count2),
     +                   specT,a46T,iflag)
            call S24_interp (period(count1),period(count2),a25(count1),a25(count2),
     +                   specT,a25T,iflag)
            call S24_interp (period(count1),period(count2),a28(count1),a28(count2),
     +                   specT,a28T,iflag)
            call S24_interp (period(count1),period(count2),a29(count1),a29(count2),
     +                   specT,a29T,iflag)
            call S24_interp (period(count1),period(count2),a31(count1),a31(count2),
     +                   specT,a31T,iflag)
            call S24_interp (period(count1),period(count2),a36(count1),a36(count2),
     +                   specT,a36T,iflag)
            call S24_interp (period(count1),period(count2),a37(count1),a37(count2),
     +                   specT,a37T,iflag)
            call S24_interp (period(count1),period(count2),a38(count1),a38(count2),
     +                   specT,a38T,iflag)
            call S24_interp (period(count1),period(count2),a39(count1),a39(count2),
     +                   specT,a39T,iflag)
            call S24_interp (period(count1),period(count2),a40(count1),a40(count2),
     +                   specT,a40T,iflag)
            call S24_interp (period(count1),period(count2),a41(count1),a41(count2),
     +                   specT,a41T,iflag)
            call S24_interp (period(count1),period(count2),a42(count1),a42(count2),
     +                   specT,a42T,iflag)
            call S24_interp (period(count1),period(count2),Vlin(count1),Vlin(count2),
     +                   specT,VlinT,iflag)
            call S24_interp (period(count1),period(count2),b(count1),b(count2),
     +                   specT,bT,iflag)
            call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +                   specT,c4T,iflag)
            call S24_interp (period(count1),period(count2),s1est(count1),s1est(count2),
     +                   specT,s1estT,iflag)
            call S24_interp (period(count1),period(count2),s2est(count1),s2est(count2),
     +                   specT,s2estT,iflag)
            call S24_interp (period(count1),period(count2),s1msr(count1),s1msr(count2),
     +                   specT,s1msrT,iflag)
            call S24_interp (period(count1),period(count2),s2msr(count1),s2msr(count2),
     +                   specT,s2msrT,iflag)
            call S24_interp (period(count1),period(count2),s3(count1),s3(count2),
     +                   specT,s3T,iflag)
            call S24_interp (period(count1),period(count2),s4(count1),s4(count2),
     +                   specT,s4T,iflag)
            call S24_interp (period(count1),period(count2),s5(count1),s5(count2),
     +                   specT,s5T,iflag)
            call S24_interp (period(count1),period(count2),s6(count1),s6(count2),
     +                   specT,s6T,iflag)
 1011 period1 = specT                                                                                                              
C     Constant values
      n = 1.5
      M2 = 5.0
c      a4 = -0.1
      a5 = -0.41
      a7 = 0.0
C     Set C term
      if (period1 .eq. -1.0) then
         c = 2400.0
      else
         c = 2.4
      endif
C     Magnitude dependent taper for C4 (eq. 4.4)
      if (mag .ge. 5.0) then
         c4_mag = c4T
      elseif (mag .ge. 4.0) then
         c4_mag = c4T - (c4T-1.0) * (5.0-mag)
      else
         c4_mag = 1.0
      endif 
     
c     Set distance (eq 4.3)
      R = sqrt(rRup**2 + c4_mag**2)
          
C     Base Model (eq 4.2)
      if ( mag .lt. M2 ) then
        f1 = a1T + a6T*(Mag-M2) + a7*(Mag-M2)**2 + a4T*(M2-M1T) + a8T*(8.5-M2)**2 +
     1                (a2T + a3T*(M2-M1T)) * alog(R) + a17T*Rrup
      elseif ( mag .le. M1T ) then
        f1 = a1T + a4T*(Mag-M1T) + a8T*(8.5-Mag)**2 + (a2T + a3T*(Mag-M1T)) * alog(R) + a17T*Rrup
      else
        f1 = a1T + a5*(Mag-M1T) + a8T*(8.5-Mag)**2 + (a2T + a3T*(Mag-M1T)) * alog(R) + a17T*Rrup
      endif
   
c     style of faulting (eq 4.5 and 4.6) 
      if ( mag .gt. 5. ) then
        f7 = Frv * a11T
        f8 = Fn * a12T
      elseif ( mag .ge. 4. ) then
        f7 = Frv * a11T * (mag-4.)
        f8 = Fn * a12T * (mag-4.)
      else 
        f7 = 0
        f8 = 0
      endif
c     ZTOR (eq 14) 
c     form modified"Extend the upper bound ZTOR to 50km:"
      if (ZTOR .lt. 50.) then 
        f6 = a15T * ZTOR/50.0
      else
        f6 = a15T
      endif    
c     Set VS30_star (eq 4.8 and 4.9)
      if ( specT .gt. 3.0 ) then
        V1 = 800.
      elseif ( specT .gt. 0.5 ) then
        V1 = exp( -0.35 * alog(specT/0.5)  + alog(1500.) )
      else
        V1=1500.
      endif      
      if ( vs30 .lt. v1 ) then 
         vs30Star = vs30
      else
      vs30Star = v1
      endif  
c     Compute site amplification (Eq. 4.7)  
      if (vs30 .lt. vLinT) then
        f5 = a10T*alog(vs30Star/vLinT) - bT*alog(c+Sa1180) 
     1              + bT*alog(Sa1180+c*((vs30Star/vLinT)**(n)) )
      else
      f5 = (a10T + bT*n) * alog(vs30Star/vLinT)
      endif
      if (vs30 .eq. 1180.) then
      f5 = (a10T + bT*n) * alog(1180/vLinT)
      endif
   
c     Set Regional z1 reference (eq 4.18)
      if (regionflag .eq. 1) then  
         z1_ref =  exp(-2.629 / 4.0 * alog((vs30**4.0 + 253.299**4.0)/(2491.945**4.0 + 253.299**4.0)))/ 1000.
      elseif (regionflag .eq. 3) then
         z1_ref = exp ( -5.23/2. * alog( (Vs30**2.0 + 412.**2.0)/(1360.**2.0+412.**2.0) ) ) / 1000.
      else
         z1_ref = exp ( -7.67/4. * alog( (Vs30**4.0 + 610.**4.0)/(1360.**4.0+610.**4.0) ) ) / 1000.
      endif 
      
C     Soil Depth Model (eq 4.17)
C     Updated 8/1/13
      if ( vs30 .lt. 150.0 ) then
         y1z = a43T
         y2z = a43T 
         x1z = 50.0
         x2z = 150.0
      elseif ( vs30 .lt. 250.0 ) then
         y1z = a43T
         y2z = a44T 
         x1z = 150.0
         x2z = 250.0
      elseif ( vs30 .lt. 400.0 ) then
         y1z = a44T
         y2z = a45T 
         x1z = 250.0
         x2z = 400.0
      elseif ( vs30 .lt. 700.0 ) then
         y1z = a45T
         y2z = a46T 
         x1z = 400.0
         x2z = 700.0
      else
         y1z = a46T
         y2z = a46T 
         x1z = 700.0
         x2z = 1000.0
      endif
C     Calculation f10 term and set it equal to zero for Vs=1180m/s (i.e., reference condition)
      if (vs30 .eq. 1180.0) then
          f10 = 0.0
      else
         f10 = ( y1z + (Vs30-x1z)*(y2z-y1z)/(x2z-x1z))*alog( (z1 + 0.01) / (z1_ref+0.01) )
      endif
c     Compute HW taper1 (eq 4.11) 
      if ( dip .le. 30. ) then
        HW_taper1 = 60./ 45.
      else
        HW_taper1 = (90.-dip)/45.
      endif
c     Compute HW taper2 (eq. 4.12)
      hw_a2 = 0.2
      if( mag .ge. 6.5 ) then
        HW_taper2 = 1. + hw_a2 * (mag-6.5) 
      elseif ( mag .gt. 5.5 ) then
        HW_taper2 = 1. + HW_a2 * (mag-6.5) - (1.0 - HW_a2)*(mag-6.5)**2
      else
        HW_taper2 = 0.
      endif
c     Compute HW taper 3 (eq. 4.13)
C   April 11, correction by ronnie for HW_taper3 when Rx.gt.R2
      h1 = 0.25
      h2 = 1.5
      h3 = -0.75
      R1 = fltWidth * cos(dip*3.1415926/180.)
      R2 = 3.*R1
      if ( Rx .lt. R1 ) then
        HW_taper3 = h1 + h2*(Rx/R1) + h3*(Rx/R1)**2
      elseif ( Rx . le. R2 ) then
        HW_taper3 = 1. - (Rx-R1)/(R2-R1)
      else
        HW_taper3 = 0.
      endif 
c     Compute HW taper 4 (eq 4.14)
      if ( ZTOR .le. 10. ) then
        HW_taper4 = 1. - (ZTOR**2) / 100.
      else
        HW_taper4 = 0.
      endif
      
c     Compute HW taper 5 (eq. 13)  **** Ry0 version ***
      Ry1 = Rx * tan(20.*3.1415926/180.)
      if ( Ry0 .lt. Ry1 ) then
        HW_taper5 = 1.
      elseif ( Ry0-Ry1 .lt. 5. ) then
        HW_taper5 = 1. - (Ry0-Ry1) / 5.
      else
        HW_taper5 = 0.
      endif
c     Compute HW taper 5 (eq. 4.15b)  **** No Ry0 version ***     
c      if (Rjb .eq. 0. ) then
c        HW_taper5 = 1. 
c      elseif ( Rjb .lt. 30. ) then
c        HW_taper5 = 1 - Rjb/30.
c      else
c        HW_taper5 = 0.
c      endif
c     Hanging wall Model (eq 4.10)
      if ( HWFlag .eq. 1 ) then
        f4 = a13T * HW_taper1 * HW_taper2 * HW_taper3 * HW_taper4 * HW_taper5
      else
        f4 = 0.
      endif
C     Add aftershock factor (eq 4.21)
      if (msasflag .eq. 1) then
         if (CRjb .ge. 15.0) then
             f11 = 0.0
         elseif (CRjb .le. 5.0) then
             f11 = a14T
         else
             f11 = a14T * ( 1.0 - (CRjb - 5.0) /10.0) 
         endif
      elseif (msasflag .eq. 0) then
         f11 = 0.0
      endif 
C     Now apply the regional attenuation differences (eq 4.22)
C     Global No Change
      if (regionflag .eq. 0) then
         freg = 0.0
C     Taiwan
      elseif (regionflag .eq. 1) then
         f12 = a31T * alog(vs30Star/VlinT) 
         freg = f12 + a25T*Rrup 
      elseif (regionflag .eq. 1 .and. vs30 .eq. 1180.0) then
         f12 = a31T * alog(1180/VlinT) 
         freg = f12 + a25T*Rrup 
C     China
      elseif (regionflag .eq. 2) then
         freg = a28T*Rrup 
C     Japan
      elseif (regionflag .eq. 3) then
         freg = a42T + a29T*Rrup 
      endif
C     Set the Sigma Values
      if (regionflag .ne. 3)  then
c     Compute within-event term, phiA, at the surface for linear site response (eq 7.1)
        if (mag .lt. 4.0) then
           phiA_est = s1estT
        elseif (mag .le. 6.0) then
           phiA_est = s1estT + ((s2estT-s1estT)/2.0)*(mag-4.0)
        else
           phiA_est = s2estT
        endif
c     Compute within-event term, phiA, for known Vs30
        if (mag .lt. 4.0) then
           phiA_msr = s1msrT
        elseif (mag .le. 6.0) then
           phiA_msr = s1msrT + ((s2msrT-s1msrT)/2.0)*(mag-4.0)
        else
           phiA_msr = s2msrT
        endif
C     choose phiA by Vs30 class
        if (vs30_class .eq. 0 ) then
      phiA = phiA_est
        elseif (vs30_class .eq. 1) then
      phiA = phiA_msr
        else
      stop 99
        endif
C     Set Sigma values for Japan Region
      else
C calculate phi_A for Japan (eq. 7.3)
        if (Rrup .lt. 30) then
           phiA = s5T        
        elseif (Rrup .le. 80) then
           phiA = s5T + (s6T-s5T)/50*(Rrup-30)
        else
           phiA = s6T
        endif
      endif
   
c     Compute between-event term, tau (eq. 7.2)
      if (mag .lt. 5.0) then
         tauA = s3T
      elseif (mag .le. 7.0) then
         tauA = s3T + ((s4T-s3T)/2.0)*(mag-5.0)
      else
         tauA = s4T
      endif
      tauB = tauA
c     Compute phiB, within-event term with site amp variability removed (eq. 7.7)
c     with fix to model for small mag at long periods - limit sigAmp to be less than phiA
      sigAmp = 0.4
      if (phiA .le. sigAmp) then
        sigAmp = phiA*0.99
      endif
      phiB = sqrt( phiA**2 - sigAmp**2)
      
c     Compute partial derivative of alog(soil amp) w.r.t. alog(Sa1180) (eq. 7.10)
      if ( vs30 .ge. vLinT) then
        dAmp_dSa1180 = 0.
      else
        dAmp_dSa1180 = bT*Sa1180 * ( -1. / (Sa1180+c) 
     1              + 1./ (Sa1180 + c*(vs30/vLinT)**(n)) )
      endif
C     Compute phi, with non-linear effects (eq. 7.8)
      phi = sqrt( phiB**2 * (1. + dAmp_dSa1180)**2 + sigAmp**2 )
C     Compute tau, with non-linear effects (eq. 7.9)
      tau = tauB * (1. + dAmp_dSa1180)
      
c     Compute median ground motion (eq. 1)
      lnSa = f1 + f4 + f5 + f6 + f7 + f8 + f10 + f11 + freg 
      return
      end
 
 
c-------------------- Adjusted in Taiwan SSHAC Project--------------------------------------------------  
              
      Subroutine S04_Bindi14_TW_C01 ( m, jbDist, ftype, specT,
     1                     period2, lnY, sigma, iflag, vs, phiT, tauT )
      implicit none
      integer MAXPER
      parameter (MAXPER=21)
      REAL Period(MAXPER), e1(MAXPER), c1(MAXPER), c2(MAXPER), h(MAXPER)
      REAL c3(MAXPER), b1(MAXPER), b2(MAXPER), b3(MAXPER), gamma(MAXPER)
      REAL sofN(MAXPER), sofR(MAXPER), sofS(MAXPER), phi(MAXPER), tau(MAXPER), sig(MAXPER), sigs2s(MAXPER)
      real e1T, c1T, c2T, hT, c3T, b1T, b2T, b3T, gammaT, sofNT, sofRT, sofST, sigs2sT
      real phiT, tauT, sigT, period1
      real Rref, Mref, Mh, R, Vref, vs
 
      REAL M, jbDist, specT, sigma, termsof
      REAL period2, lnY, ftype
      integer iflag, count1, count2, nPer, i
      real f_D, f_M, f_S
  
 
      Data period(1:21) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2,
     1      3/
      Data e1(1:21) / 0.659508496312107, 0.659508496312107, 0.684257395856917, 0.792185875253012, 0.910245329949904,
     1      1.04124346176327, 1.36695980240164, 1.58175260431507, 1.55527075935291, 1.48876800771044, 1.47418275968898,
     1      1.43106965774598, 1.32383819300195, 1.3382650551676, 1.08674468845801, 0.744502966314904, 0.339631620365221,
     1      0.267909172341515, -0.363413009167707, -0.8677036818281, -1.77042537043066/
      Data c1(1:21) / -1.26358, -1.26358, -1.26358, -1.29088019990866, -1.31025, -1.30237878943392, -1.28882132582112, -1.28178,
     1      -1.23465110534188, -1.17697, -1.14479188763455, -1.10301, -1.08792194443592, -1.10591, -1.09538, -1.05767,
     1      -1.04831378793046, -1.0527, -0.983388, -0.979215, -0.940373/
      Data c2(1:21) / 0.220527, 0.220527, 0.220527, 0.234653259429915, 0.244676, 0.239572871690745, 0.229465497782482, 0.219406,
     1      0.202883682422306, 0.182662, 0.161122330083983, 0.133154, 0.118226690912901, 0.108276, 0.101111, 0.112197,
     1      0.122345650356548, 0.103471, 0.109072, 0.163344, 0.227241/
      Data h(1:21) / 5.20082, 5.20082, 5.20082, 5.0346146046701, 4.91669, 5.09314627212803, 5.50666272693079, 6.12146,
     1      5.95062506384396, 5.74154, 5.55812998041133, 5.31998, 5.16226984536375, 5.12846, 4.95386, 4.43205, 4.18619685740845,
     1      4.41613, 4.56697, 4.58186, 5.74173/
      Data c3(1:21) / -0.00162234204269646, -0.00162234204269646, -0.00169729704805731, -0.00182364633843366, -0.00211014463275481,
     1      -0.00287513335217997, -0.004184784769582, -0.00466828237187935, -0.00500654836196918, -0.00499803696576266,
     1      -0.00473066619716198, -0.00423022629441665, -0.00300999356000803, -0.00189118393920752, -0.000589441966065091,
     1      -0.000241978980141242, 0.000285085906072067, -1.74038270368973e-05, -0.00134660884745955, -0.00283883648666463,
     1      -0.00339965384388641/
      Data b1(1:21) / -0.270897727235864, -0.270897727235864, -0.281154290911621, -0.291943091967611, -0.312942644598605,
     1      -0.279796632295181, -0.224987839767471, -0.275839067753525, -0.276546876529609, -0.247413951462292, -0.15365006490395,
     1      -0.0299424838422354, -0.0227660513080655, 0.0145880936350126, 0.0160023913376982, 0.0708600769194794,
     1      0.304811546059666, 0.646443935653183, 0.889724237533995, 0.819623500158111, 0.957007301722072/
      Data b2(1:21) / -0.210951793811668, -0.210951793811668, -0.210951793811668, -0.193123478069862, -0.180474086745271,
     1      -0.160842272902352, -0.139406648181843, -0.173459491517974, -0.201518836294038, -0.235860698830659, -0.238616639987212,
     1      -0.242195110421486, -0.285257767878105, -0.319087637016929, -0.373226017723405, -0.377872634441067, -0.365550258824207,
     1      -0.339764851152015, -0.322055669201798, -0.341431922759343, -0.256611595688721/
      Data b3(1:21) / 0, 0, 0, 0, 0, 0, 0, 0, 0.0765190099149165, 0.170170709229651, 0.246865015276161, 0.346449255676977,
     1      0.402141458721742, 0.436758945609297, 0.517683298702572, 0.356325043140829, 0.186190244946956, 0.213765783053816,
     1      0.227350804843972, 0, 0/
      Data gamma(1:21) / -0.422676821032021, -0.422676821032021, -0.414978878595817, -0.406699004557271, -0.3914970371187,
     1      -0.369326047622375, -0.343869768741454, -0.333289584723655, -0.352050474022221, -0.367605646663696, -0.387393323332123,
     1      -0.396818066941227, -0.39952338981378, -0.445165123281344, -0.512409231991832, -0.561508199422803, -0.752231535522345,
     1      -0.873980613150435, -0.968710126360196, -0.981012736821664, -0.992114119974754/
      Data sofN(1:21) / -0.117251075170762, -0.117251075170762, -0.115338543884633, -0.119347338965019, -0.108122399642431,
     1      -0.110836575215843, -0.127904216218309, -0.151584972241392, -0.131927184713145, -0.091333399271005,
     1      -0.0734541746220918, -0.0720743866300066, -0.053504945030001, -0.0577148145517801, -0.0967772518776173,
     1      -0.153197353204361, -0.176638749464974, -0.202427135498213, -0.184012887842084, -0.180956841720691,
     1      -0.0658311164233819/
      Data sofR(1:21) / 0.0950331165879268, 0.0950331165879268, 0.0942147492408604, 0.0886194576576785, 0.0865931352500485,
     1      0.0798032805953013, 0.0518078943711761, 0.0377736339116048, 0.0374540870711511, 0.0461070357350431, 0.0650643320300671,
     1      0.0987719958171821, 0.125118256865787, 0.123740513623839, 0.137614216748372, 0.143277552748404, 0.160418527309156,
     1      0.172380927444367, 0.156550741034289, 0.135252483671118, 0.10353955455201/
      Data sofS(1:21) / -0.096283977, -0.096283977, -0.096283977, -0.0966103375235148, -0.096841894, -0.0963783235595202,
     1      -0.0979532926457108, -0.107435167, -0.101695221063626, -0.094670095, -0.0968825621131357, -0.099755355,
     1      -0.0954761145312324, -0.090850567, -0.108074134, -0.098696165, -0.0876704883298241, -0.103411629, -0.096162861,
     1      -0.058467241, -0.005358461/
      Data tau(1:21) / 0.317898754427414, 0.317898754427414, 0.319323264507343, 0.332912367257089, 0.345138786542856,
     1      0.356539794838403, 0.373304099976226, 0.365831192336457, 0.352942682997377, 0.334646212108993, 0.328499912385758,
     1      0.313551883565158, 0.297330303430244, 0.299006855270763, 0.33361918283012, 0.364066847460064, 0.40039358200669,
     1      0.420084326301466, 0.436932806015766, 0.45401623183321, 0.490721750019259/
      Data phi(1:21) / 0.32571333838296, 0.32571333838296, 0.325174907932113, 0.334158210279457, 0.35048253853852,
     1      0.370382687881984, 0.409261273800818, 0.41946450409783, 0.415563072394708, 0.393759932519684, 0.382737515801976,
     1      0.363483599894857, 0.341015171374459, 0.335643324687049, 0.327306984921592, 0.31669506610895, 0.328111960636172,
     1      0.334986782672972, 0.363636845023503, 0.388384142446741, 0.404771608028183/
      Data sigs2s(1:21) / 0.474767396333149, 0.474767396333149, 0.474134438654724, 0.480278729542088, 0.486514622964184,
     1      0.488556279560537, 0.487077378450021, 0.481594309301403, 0.477211581548162, 0.47594992226976, 0.472293215881846,
     1      0.475976807032338, 0.48288546933477, 0.488981047072412, 0.499582923250037, 0.506548849370516, 0.50222265741314,
     1      0.491260159723756, 0.472398005717859, 0.454546252692053, 0.434481072336867/
      Data sig(1:21) / 0.745772773, 0.745772773, 0.745772773, 0.753543188479092, 0.759056386, 0.767508805717544, 0.783602925138285,
     1      0.797567122, 0.792929650007148, 0.787253843, 0.780874405181787, 0.772590981, 0.777647499576309, 0.775374807,
     1      0.77430871, 0.78683247, 0.795506751916842, 0.819874566, 0.835458462, 0.84932463, 0.829789498/
 
C Find the requested spectral period and corresponding coefficients
      nPer = 21
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         e1T      = e1(1)
         c1T      = c1(1)
         c2T      = c2(1)
         hT       = h(1)
         c3T      = c3(1)
         b1T      = b1(1)
         b2T      = b2(1)
         b3T      = b3(1)
         gammaT   = gamma(1)
         sofNT    = sofN(1)
         sofRT    = sofR(1)
         sofST    = sofS(1)
         sigs2sT  = sigs2s(1)
         phiT     = phi(1)
         tauT     = tau(1)
         sigT     = sig(1)
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
      write (*,*) 'Bindi et al. (2013) Horizontal'
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
 1020 call S24_interp (period(count1),period(count2),e1(count1),e1(count2),
     +             specT,e1T,iflag)
      call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +             specT,hT,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),b1(count1),b1(count2),
     +             specT,b1T,iflag)
      call S24_interp (period(count1),period(count2),b2(count1),b2(count2),
     +             specT,b2T,iflag)
      call S24_interp (period(count1),period(count2),b3(count1),b3(count2),
     +             specT,b3T,iflag)
      call S24_interp (period(count1),period(count2),gamma(count1),gamma(count2),
     +             specT,gammaT,iflag)
      call S24_interp (period(count1),period(count2),sofN(count1),sofN(count2),
     +             specT,sofNT,iflag)
      call S24_interp (period(count1),period(count2),sofR(count1),sofR(count2),
     +             specT,sofRT,iflag)
      call S24_interp (period(count1),period(count2),sofS(count1),sofS(count2),
     +             specT,sofST,iflag)
      call S24_interp (period(count1),period(count2),sigs2s(count1),sigs2s(count2),
     +             specT,sigs2sT,iflag)
      call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +             specT,phiT,iflag)
      call S24_interp (period(count1),period(count2),tau(count1),tau(count2),
     +             specT,tauT,iflag)
      call S24_interp (period(count1),period(count2),sig(count1),sig(count2),
     +             specT,sigT,iflag)
   
 1011 period1 = specT                                                                                                              
C     Set Constant Terms
      Mref = 5.5
      Rref = 1.0
      Mh = 6.75
      Vref = 800.0
C     Set the mechanism term. 
      if (ftype .eq. 0 ) then
         termsof = 0
      elseif (ftype .ge. 0.5) then
         termsof = sofRT
      elseif (ftype .le. -0.5) then
         termsof = sofNT
      endif
      R = sqrt (jbdist**2 + hT**2)
      f_D = (c1T+c2T*(M-Mref))*alog(R/Rref) + c3T*(R-Rref)
C     Compute the ground motion for the given spectral period. 
      if (M .le. Mh) then
         f_M = b1T*(M-Mh) + b2T*(M-Mh)**2.0 
      else
         f_M = b3T*(M-Mh)  
      endif
   
      f_S = gammaT*alog(vs/vref)
   
      lnY = e1T + f_D + f_M + f_S + termsof
   
C     Set the sigma value and convert from log10 to Ln units
c      phiT = phiT*alog(10.0)
c      tauT = tauT*alog(10.0)
c      sigma = sigT*alog(10.0)
c      sigs2sT = sigs2sT*alog(10.0)
C     Convert ground motion to units of gals in natural log units.
c      lnY = lnY*alog(10.0)
      lnY = lnY + 6.89
      period2 = period1
      return
      end
 
 
c ---------------------------------------------------------------------------            
C     *** Boore, Stewart, Seyhan and Atkinson NGA West 2 (NGA West2-2013) ***
C         Earthquake Spectra Report:
C            NGA-West2 Equations for Predicting PGA, PGV, and 5%-Damped
C                PSA for Shallow Crustal Earthquakes.
C             D. M. Boore, J. P. Stewart, E. Seyhan, and G. M. Atkinson
C     Notes:
C        Applicable Range:
C            3.0 < M < 8.5 (Strike-Slip)
C            3.0 < M < 7.0 (Normal)
C            Distance < 300km
C            150 < Vs < 1500 m/s
C            0.0 < Z1 < 3.0 km
C            Region Flag:
C               0 = Global 
C               1 = China-Turkey
C               2 = Italy-Japan
c ---------------------------------------------------------------------------            
      subroutine S04_BSSA14_TW_C01 ( mag, Rbjf, specT, 
     1        period2, lnY, sigma, iflag, vs, ftype, pga4nl, z10, regionflag, basinflag,
     1        phi, tau ) 
C     Last Updated: 9/16/13
      parameter (MAXPER=25)
      REAL Period(MAXPER), c1(MAXPER), c2(MAXPER), c3(MAXPER)
      real h(MAXPER), DC3ChinaTrk(MAXPER), DC3ItalyJapan(MAXPER), e0(MAXPER)
      real e1(MAXPER), e2(MAXPER), e3(MAXPER), e4(MAXPER)
      real e5(MAXPER), e6(MAXPER), mh(MAXPER), c(MAXPER), Vc(MAXPER)
      real phi2(MAXPER), phi3(MAXPER), f4(MAXPER)
      real l1(MAXPER), l2(MAXPER), t1(MAXPER), t2(MAXPER)
      real f5(MAXPER), rjbbar(MAXPER), Dfr(MAXPER), Dfv(MAXPER)
      real R1(MAXPER), R2(MAXPER), DC3Global(MAXPER), f6(MAXPER), f7(MAXPER)
      real Mref, Rref, Vref, f1, f3, specT
      REAL MAG, RBJF, VS, z10
      real ftype, Rp, R
      INTEGER iFlag, count1, count2, regionflag, basinflag
      real lnY, mechS, mechN, mechR, pga4nl
      real f2, flin, fBasin, phi, tau
      real c1T, c2T, c3T, hT, e0T, e1T, e2T, e3T, e4T
      real e5T, e6T, mhT, cT, VcT, phi2T, phi3T, f4T, l1T, l2T, t1T, t2T
      real deltaz1, f5T, rjbbarT, DfrT, DfvT, R1T, R2T, DC3GlobalT
      real DC3ChinaTrkT, DC3ItalyJapanT, f6T, f7T
 
 
      Data Period(1:25) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2,
     1      3, 4, 5, 7.5, 10/
      Data e0(1:25) / 0.4534, 0.48598, 0.56916, 0.673459231164795, 0.75436, 0.96447, 1.1268, 1.20895293439595, 1.3095,
     1      1.3164611925074, 1.3255, 1.2766, 1.2217, 1.1046, 0.96991, 0.66903, 0.3932, -0.14954, -0.58669, -1.1898, -1.6388,
     1      -1.966, -2.5865, -3.0702, 0.4473/
      Data e1(1:25) / -0.0192358695063516, 0.0567227727912004, 0.103830570426141, 0.198897708101522, 0.271900769282252,
     1      0.367398657969924, 0.606318878890423, 0.823883469673358, 0.913809415580094, 0.972718537695808, 0.972741127823899,
     1      0.952673927272273, 0.929588626686887, 0.906677969774733, 0.843243142498517, 0.70851892231752, 0.369210918626143,
     1      0.162950871861978, -0.505519010009077, -0.791847279739693, -1.52940409045199, -2.08919352516209, -2.29332932326379,
     1      -2.99734840119961, -3.09836742174658/
      Data e2(1:25) / -0.0466127974661881, -0.0899522328065243, -0.0889941216873682, -0.0881837044024162, -0.0738384916140176,
     1      -0.0774606068618486, -0.0933753870789509, -0.119805547559339, -0.102166706267363, -0.0619031495690598,
     1      -0.0438630892907241, -0.0424937532989739, -0.0281214761287716, -0.0309767673715999, -0.081060821379759,
     1      -0.135927788790708, -0.141847724215483, -0.157893131542867, -0.119300477644201, -0.103636530614135, 0.0302804751375015,
     1      0.197391958820006, -0.105591250006451, 0.09472882175673, 0.0180272297249192/
      Data e3(1:25) / 0.113150415770813, 0.100661193118287, 0.0982992168992328, 0.0926856182254271, 0.0910574057573448,
     1      0.0841869161547781, 0.0561265104794706, 0.0417968407844044, 0.0410646271818258, 0.0492102889262502, 0.0678307028722051,
     1      0.101248016855977, 0.129758071123678, 0.128286336016744, 0.141561631803489, 0.145442944817284, 0.163752332665633,
     1      0.176717020297631, 0.158522053758111, 0.13897165041878, 0.10401554209148, 0.0505868993201825, -0.0453840386998217,
     1      -0.150207517493787, -0.13838956231916/
      Data e4(1:25) / 0.634956463990215, 0.619860972185183, 0.633560416066013, 0.626773923260574, 0.605074753181471,
     1      0.566320478533912, 0.565076405914992, 0.596866780802841, 0.589028550987566, 0.607990040954619, 0.624159967438116,
     1      0.623636047860314, 0.588632504600515, 0.678668378713168, 0.67378138305232, 0.762124819833728, 1.11558730063412,
     1      1.40281398367416, 1.861987715669, 2.03842856513544, 2.2278862997874, 2.24750560056704, 2.27392203189865,
     1      2.02295395240319, 1.47472914989199/
      Data e5(1:25) / 0.05053, 0.04932, 0.05339, 0.06144, 0.0647739711041879, 0.06736, 0.07355, 0.05523, 0.0114780540956452,
     1      -0.04207, -0.072042284489666, -0.11096, -0.16213, -0.1959, -0.22608, -0.23522, -0.21591, -0.18983, -0.1467, -0.11237,
     1      -0.04332, -0.01464, -0.01486, -0.08161, -0.15096/
      Data e6(1:25) / -0.1662, -0.1659, -0.16561, -0.1669, -0.174739337461198, -0.18082, -0.19665, -0.19838, -0.19116744899994,
     1      -0.18234, -0.171976524654611, -0.15852, -0.12784, -0.09286, -0.02319, 0.02912, 0.10829, 0.17895, 0.33896, 0.44788,
     1      0.62694, 0.76303, 0.87314, 1.0121, 1.0651/
      Data Mh(1:25) / 5.75939535845878, 5.93161026644103, 5.89125554011655, 5.90885411844775, 5.92355835882231, 5.93448092637676,
     1      5.91330519700299, 5.90874187926377, 5.92092742673115, 5.94829040352112, 5.93622578077096, 5.96033962320899,
     1      6.06055764485135, 6.03454555393634, 6.20028144493362, 6.22545222018765, 6.17573383439527, 6.18463084988356,
     1      6.01489943068217, 6.0514358625268, 6.07860067368602, 6.1072608317045, 6.18885674149109, 6.23781972475555,
     1      6.57470966466657/
      Data c1(1:25) / -1.134, -1.134, -1.1394, -1.1421, -1.12734492518079, -1.1159, -1.0831, -1.0652, -1.05980407655856, -1.0532,
     1      -1.05646305898784, -1.0607, -1.0773, -1.0948, -1.1243, -1.1459, -1.1777, -1.193, -1.2063, -1.2159, -1.2179, -1.2162,
     1      -1.2189, -1.2543, -1.3253/
      Data c2(1:25) / 0.1917, 0.1916, 0.18962, 0.18842, 0.187670982843147, 0.18709, 0.18225, 0.17203, 0.163927121632102, 0.15401,
     1      0.150042120270783, 0.14489, 0.13925, 0.13388, 0.12512, 0.12015, 0.11054, 0.10248, 0.09645, 0.09636, 0.09764, 0.10218,
     1      0.10353, 0.12507, 0.15183/
      Data c3(1:25) / -0.00788875587408985, -0.0078818509689674, -0.00780428543929034, -0.00832664074171499, -0.0092200593557657,
     1      -0.0100375197535989, -0.0115480726222495, -0.0118959621513459, -0.0115081781258884, -0.0105994476660068,
     1      -0.00984494177643466, -0.00870503597848327, -0.00695258543010702, -0.00578513843981538, -0.00361675198164066,
     1      -0.00202594754252705, -0.000393963454749957, -0.000268145047578852, 0.000338647265255224, -0.000274615471785213,
     1      0.000881504173671615, 0.00121342431999415, 0.00142297452256839, 0.00116987226508603, 0.00197956130485759/
      Data h(1:25) / 4.5, 4.5, 4.5, 4.49, 4.32668046955837, 4.2, 4.04, 4.13, 4.24691167456457, 4.39, 4.48571639697672, 4.61, 4.78,
     1      4.93, 5.16, 5.34, 5.6, 5.74, 6.18, 6.54, 6.93, 7.32, 7.78, 9.48, 9.66/
      Data Dc3Global(1:25) / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      Data Dc3ChinaTrk(1:25) / 0.002816, 0.00278, 0.002765, 0.00287312879256825, 0.002957, 0.002957, 0.002879, 0.00283718159332883,
     1      0.002786, 0.00271029703148205, 0.002612, 0.002444, 0.002196, 0.002107, 0.002348, 0.00269, 0.002921, 0.003039, 0.002923,
     1      0.002616, 0.002605, 0.002604, 0.0026, 0.00303, 0.002858/
      Data Dc3ItalyJapan(1:25) / -0.002437, -0.00234, -0.002168, -0.00206831876935114, -0.001991, -0.002159, -0.002439,
     1      -0.00255905929657207, -0.002706, -0.00282085967637207, -0.00297, -0.00314, -0.003297, -0.003212, -0.002907, -0.002527,
     1      -0.002089, -0.001518, -0.00117, -0.001188, -0.001083, -0.000571, 0.000385, 0.00149, -0.00255/
      Data c(1:25) / -0.494639525126874, -0.492386790829841, -0.486959917611639, -0.483734070461185, -0.473331126329196,
     1      -0.453061661114652, -0.428535271473178, -0.422096875367947, -0.442766131734041, -0.460824990697695, -0.480490619799534,
     1      -0.48968579619744, -0.492658432074322, -0.536122783765441, -0.599481092255171, -0.641720588625172, -0.829810827332613,
     1      -0.948492569720881, -1.03094032845303, -1.03669695507734, -1.05585320583303, -1.07883768198775, -1.0420295546124,
     1      -0.981155106166575, -0.88092955963527/
      Data Vc(1:25) / 1500, 1500.2, 1500.36, 1502.95, 1502.08834868422, 1501.42, 1494, 1479.12, 1462.81082139824, 1442.85,
     1      1420.99185552677, 1392.61, 1356.21, 1308.47, 1252.66, 1203.91, 1147.59, 1109.95, 1072.39, 1009.49, 922.43, 844.48,
     1      793.13, 771.01, 775/
      Data f4(1:25) / -0.15, -0.14833, -0.1471, -0.15485, -0.1782103245611, -0.19633, -0.22866, -0.24916, -0.252743792485691,
     1      -0.25713, -0.252539963690434, -0.24658, -0.23574, -0.21912, -0.19582, -0.17041, -0.13866, -0.10521, -0.06794, -0.03614,
     1      -0.01358, -0.00321, -0.00025, -5e-05, 0/
      Data f5(1:25) / -0.00701, -0.00701, -0.00728, -0.00735, -0.00685440970072883, -0.00647, -0.00573, -0.0056,
     1      -0.0057124150716967, -0.00585, -0.00597617161419659, -0.00614, -0.00644, -0.0067, -0.00713, -0.00744, -0.00812,
     1      -0.00844, -0.00771, -0.00479, -0.00183, -0.00152, -0.00144, -0.00137, -0.00136/
      Data f6(1:25) / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.092, 0.367, 0.638, 0.871, 1.135, 1.271, 1.329, 1.329,
     1      1.183/
      Data f7(1:25) / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.059, 0.208, 0.309, 0.382, 0.516, 0.629, 0.738, 0.809,
     1      0.703/
      Data R1(1:25) / 111.667, 113.105, 112.133, 104.132595691538, 97.927, 85.989, 79.587, 80.3689592387222, 81.326,
     1      85.4944490883363, 90.907, 97.039, 103.152, 106.018, 105.536, 108.388, 116.388, 125.38, 130.369, 130.365, 129.489,
     1      130.224, 130.716, 130, 110/
      Data R2(1:25) / 270, 270, 269.998, 269.999126341589, 270, 270.035, 270.092, 270.123925880362, 270.163, 270.092082851331, 270,
     1      269.449, 268.593, 266.543, 265, 266.511, 270, 262.413, 240.138, 195, 199.446, 230, 250.395, 210, 270/
      Data Dfr(1:25) / 0.096, 0.092, 0.081, 0.0708629256967262, 0.063, 0.064, 0.087, 0.101838789463964, 0.12, 0.126961192507398,
     1      0.136, 0.141, 0.138, 0.122, 0.109, 0.1, 0.098, 0.104, 0.105, 0.088, 0.07, 0.061, 0.058, 0.06, 0.1/
      Data Dfv(1:25) / 0.07, 0.03, 0.029, 0.0295631707946263, 0.03, 0.022, 0.014, 0.0144496602867868, 0.015, 0.0280522359513715,
     1      0.045, 0.055, 0.05, 0.049, 0.06, 0.07, 0.02, 0.01, 0.008, 0, 0, 0, 0, 0, 0.07/
      Data l1(1:25) / 0.698, 0.7018, 0.7212, 0.73916514834858, 0.7531, 0.7447, 0.7279, 0.72448258182042, 0.7203, 0.716384329214589,
     1      0.7113, 0.6984, 0.6754, 0.6428, 0.6147, 0.5815, 0.5527, 0.5317, 0.5263, 0.5335, 0.536, 0.5285, 0.5117, 0.5103, 0.6951/
      Data l2(1:25) / 0.4992, 0.5023, 0.5136, 0.524131293859512, 0.5323, 0.5423, 0.5407, 0.538946324881531, 0.5368,
     1      0.537670149063425, 0.5388, 0.5471, 0.5614, 0.5804, 0.599, 0.6218, 0.625, 0.6192, 0.6182, 0.619, 0.6156, 0.6223, 0.6344,
     1      0.6036, 0.4951/
      Data t1(1:25) / 0.4019, 0.4087, 0.4449, 0.477451271929402, 0.5027, 0.4744, 0.4153, 0.387780790448648, 0.3541,
     1      0.34953171741702, 0.3436, 0.35, 0.3634, 0.381, 0.4101, 0.4572, 0.4983, 0.5248, 0.5325, 0.5369, 0.5427, 0.532, 0.511,
     1      0.4869, 0.3982/
      Data t2(1:25) / 0.3446, 0.3464, 0.364, 0.39908554050522, 0.4263, 0.4658, 0.4583, 0.42664391581021, 0.3879, 0.353355082182037,
     1      0.3085, 0.2664, 0.229, 0.2097, 0.2235, 0.2664, 0.2984, 0.3151, 0.3291, 0.3438, 0.3492, 0.3354, 0.2699, 0.2392, 0.348/
 
C     Set constant parameters
      mref = 4.5
      rref = 1.0
      vref = 760.0
      f1 = 0.0
      f3 = 0.1
      V1 = 225.0
      V2 = 300.0
C First check for the PGA case (i.e., specT=0.0) 
      nPer = 25
      if (specT .eq. 0.0) then
         period1 = period(1)
         e0T = e0(1)
         e1T = e1(1)
         e2T = e2(1)
         e3T = e3(1)
         e4T = e4(1)
         e5T = e5(1)
         e6T = e6(1)
         mhT = mh(1)
         c1T = c1(1)
         c2T = c2(1)
         c3T = c3(1)        
         hT = h(1)
         cT = c(1)
         VcT = vc(1)
         phi2T = phi2(1)
         phi3T = phi3(1)
         Dc3GlobalT = DC3global(1)
         Dc3ChinaTrkT = DC3ChinaTrk(1)
         Dc3ItalyJapanT = DC3ItalyJapan(1)
         f4T = f4(1)
         f5T = f5(1)
         f6T = f6(1)
         f7T = f7(1)
         R1T = R1(1)
         R2T = R2(1)
         l1T = l1(1)
         l2T = l2(1)
         t1T = t1(1)
         t2T = t2(1)
         rjbbarT = rjbbar(1)
         DfrT = Dfr(1)
         DfvT = Dfv(1)
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
      write (*,*) 'BSSA (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),e0(count1),e0(count2),
     +                   specT,e0T,iflag)
            call S24_interp (period(count1),period(count2),e1(count1),e1(count2),
     +                   specT,e1T,iflag)
            call S24_interp (period(count1),period(count2),e2(count1),e2(count2),
     +                   specT,e2T,iflag)
            call S24_interp (period(count1),period(count2),e3(count1),e3(count2),
     +                   specT,e3T,iflag)
            call S24_interp (period(count1),period(count2),e4(count1),e4(count2),
     +                   specT,e4T,iflag)
            call S24_interp (period(count1),period(count2),e5(count1),e5(count2),
     +                   specT,e5T,iflag)
            call S24_interp (period(count1),period(count2),e6(count1),e6(count2),
     +                   specT,e6T,iflag)
            call S24_interp (period(count1),period(count2),mh(count1),mh(count2),
     +                   specT,mhT,iflag)
            call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +                   specT,c2T,iflag)
            call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +                   specT,c3T,iflag)
            call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +                   specT,hT,iflag)
            call S24_interp (period(count1),period(count2),c(count1),c(count2),
     +                   specT,cT,iflag)
            call S24_interp (period(count1),period(count2),phi2(count1),phi2(count2),
     +                   specT,phi2T,iflag)
            call S24_interp (period(count1),period(count2),phi3(count1),phi3(count2),
     +                   specT,phi3T,iflag)
            call S24_interp (period(count1),period(count2),DC3Global(count1),DC3Global(count2),
     +                   specT,DC3GlobalT,iflag)
            call S24_interp (period(count1),period(count2),DC3ChinaTrk(count1),DC3ChinaTrk(count2),
     +                   specT,DC3ChinaTrkT,iflag)
            call S24_interp (period(count1),period(count2),DC3ItalyJapan(count1),DC3ItalyJapan(count2),
     +                   specT,DC3ItalyJapanT,iflag)
            call S24_interp (period(count1),period(count2),Vc(count1),Vc(count2),
     +                   specT,VcT,iflag)
            call S24_interp (period(count1),period(count2),f4(count1),f4(count2),
     +                   specT,f4T,iflag)
            call S24_interp (period(count1),period(count2),f5(count1),f5(count2),
     +                   specT,f5T,iflag)
            call S24_interp (period(count1),period(count2),f6(count1),f6(count2),
     +                   specT,f6T,iflag)
            call S24_interp (period(count1),period(count2),f7(count1),f7(count2),
     +                   specT,f7T,iflag)
            call S24_interp (period(count1),period(count2),R1(count1),R1(count2),
     +                   specT,R1T,iflag)
            call S24_interp (period(count1),period(count2),R2(count1),R2(count2),
     +                   specT,R2T,iflag)
            call S24_interp (period(count1),period(count2),l1(count1),l1(count2),
     +                   specT,l1T,iflag)
            call S24_interp (period(count1),period(count2),l2(count1),l2(count2),
     +                   specT,l2T,iflag)
            call S24_interp (period(count1),period(count2),t1(count1),t1(count2),
     +                   specT,t1T,iflag)
            call S24_interp (period(count1),period(count2),t2(count1),t2(count2),
     +                   specT,t2T,iflag)
            call S24_interp (period(count1),period(count2),rjbbar(count1),rjbbar(count2),
     +                   specT,rjbbarT,iflag)
            call S24_interp (period(count1),period(count2),Dfr(count1),Dfr(count2),
     +                   specT,DfrT,iflag)
            call S24_interp (period(count1),period(count2),Dfv(count1),Dfv(count2),
     +                   specT,DfvT,iflag)
 1011 period1 = specT                                                                                                              
C.....Set the mechanism terms based on ftype............
C     Set mechanism term and corresponding Frv and Fnm values.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                    -120 < Rake <  -60
C     -0.5      Normal/Oblique            -150 < Rake < -120
C                                          -60 < Rake <  -30
C       0       Strike-Slip               -180 < Rake < -150
C                                          -30 < Rake <   30
C                                          150 < Rake <  180
C      0.5      Reverse/Oblique             30 < Rake <   60
C                                          120 < Rake <  150
C       1       Reverse                     60 < Rake <  120 
C     Note: Unknown Mechanism is not currently coded.  
      if (ftype .eq. -1.0) then
         mechS = 0.0
         mechN = 1.0
         mechR = 0.0
      elseif (ftype .eq. -0.5) then 
         mechS = 0.0
         mechN = 1.0
         mechR = 0.0
      elseif (ftype .eq. 0.0) then 
         mechS = 1.0
         mechN = 0.0
         mechR = 0.0
      elseif (ftype .eq. 0.5) then
         mechS = 0.0
         mechN = 0.0
         mechR = 1.0
      elseif (ftype .eq. 1.0) then
         mechS = 0.0
         mechN = 0.0
         mechR = 1.0
      endif 
C.....First compute the Reference Rock PGA value...........
C.....This will include the regional dependence for PGA....
C.....MAGNITUDE DEPENDENCE.................................
      if (mag .le. mh(1)) then
         term1 = e1(1) + e2(1)*mechN + e3(1)*mechR +
     1           e4(1)*(mag-mh(1)) + e5(1)*(mag-mh(1))**2.0
      else
         term1 = e1(1) + e2(1)*mechN + e3(1)*mechR +
     1           e6(1)*(mag-mh(1))
      endif
C.....Distance dependence......
      Rp = SQRT( Rbjf*Rbjf+h(1)*h(1) )
C.....Apply Regional term.....
         TERM2 = ( c1(1) + c2(1)*(mag-mref) ) * alog(Rp/rref) + 
     1           c3(1)  * (Rp-rref)
     
      pga4nl = exp(term1+term2)
C.....Now compute the requested ground motion value........
C.....MAGNITUDE DEPENDENCE.................................
      if (mag .le. mhT) then
         term1 = e1T + e2T*mechN + e3T*mechR +
     1           e4T*(mag-mhT) + e5T*(mag-mhT)**2.0
      else
         term1 = e1T + e2T*mechN + e3T*mechR +
     1           e6T*(mag-mhT)
      endif
C.....Distance dependence......
      R = SQRT( Rbjf*Rbjf+hT*hT )
C     Now apply the regional attenuation differnece.
C     Global Case
         TERM2 = ( c1T + c2T*(mag-mref) ) * alog(R/rref) + 
     1        c3T * (R-rref) 
C.....Site Response Term.........
C.....Now compute the site term........
C.....First the linear term......
      if (vs .le. VcT ) then
         flin = cT*alog(Vs/Vref)
      else
         flin = cT*alog(VcT/Vref)
      endif
C.....Next the non-linear term......
      f2 = f4T*(exp(f5T*(min(vs,760.0)-360.0))-exp(f5T*(760.0-360.0))) 
      
C.....Now compute the basin effect term......
C Deviation from ln(Vs30) scaling: bedrock depth (Z1) effect for California.
c      if (basinflag .eq. 1) then 
c     Compute the DeltaZ1 term. Apply the California model for all regions except for Japan. 
          deltaz1 = z10 -
     1           exp(-2.63/4.0 * alog((vs**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))/1000.0        
          fbasin = min(f7T, f6T * deltaZ1)

c      else
c         fbasin = 0.0
c      endif
      TERM3 = flin + f1 + f2*alog((pga4nl+f3)/f3) + fBasin
      
      lnY = term1 + term2 + term3 
      period2 = period1
c     Now compute the sigma value which is a function of magnitude and Vs
C     Tau (Eq. 4.11)
      if (mag .le. 4.5) then
         tau = t1T
      elseif (mag .gt. 4.5 .and. mag .lt. 5.5) then
         tau = t1T + (t2T - t1T)*(mag - 4.5)
      else
         tau = t2T
      endif
      
C     Phi - Magnitude (Eq. 4.12) 
      if (mag .le. 4.5) then
         phi = l1T
      elseif (mag .gt. 4.5 .and. mag .lt. 5.5) then
         phi = l1T + (l2T - l1T)*(mag - 4.5)
      else
         phi = l2T
      endif
C     Phi - Distance (Eq. 4.13)
      if (rbjf .le. R1T) then
          phi = phi
      elseif (rbjf .gt. R1T .and. rbjf .le. R2T) then
          phi = phi + DfrT*( (alog(rbjf/R1T))/(alog(R2T/R1T)) )
      else
          phi = phi + DfrT 
      endif
C     Phi - Vs30 (Eq. 4.14)
      if (vs .ge. V2) then
         phi = phi
      elseif (vs .ge. V1 .and. vs .lt. V2) then
         phi = phi - DfvT*( alog(V2/vs) / alog(V2/V1) )
      else
         phi = phi - DfvT
      endif
      sigma = sqrt (tau**2 + phi**2)
C     Convert ground motion to units of gals.
      lnY = lnY + 6.89
      return
      END
c ---------------------------------------------------------------------            
C ** Chiou and Youngs (NGA West2-2013 Model) Horizontal **
C     Earthquake Spectra Paper:
C        Update of the Chiou and Youngs NGA Model for the
C            Average Horizontal Component of Peak
C            Ground Motion and Response Spectra
C         B. S.J. Chiou and R.R. Youngs
C     Notes:
C        Applicable Range (see Abstract):  
C           3.5 <= M <= 8.5 Strike-slip
C           3.5 <= M <= 8.0 Reverse and Normal
C           Rrup <= 300 km
C           Ztor <= 20 km
C           180m/s <= Vs30m <= 1500m/s
C           Mainshock events only
C        Regional attenuation included based on Regionflag
C             0 = Global
C             1 = Japan/Italy
C             2 = Wenchuan (note only applicable for M7.9 event)
C         Sigma dependent on estimated or measured Vs30m based on 
C             Vs30_Class
C             0 = Estimated Vs30m
C             1 = Measured Vs30m
c ---------------------------------------------------------------------            
      Subroutine S04_CY14_TW_C01 ( m, Rrup, Rbjf, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, Delta, DTor, Ftype, depthvs10, vs30_class,
     3                     hwflag, Rx, regionflag, phi, tau )
C     Last Updated: 8/1/13
      parameter (MAXPER=25)
      REAL Period(MAXPER), C1(MAXPER), C1a(MAXPER), C1b(MAXPER)
      REAL cn(MAXPER), cm(MAXPER), c5(MAXPER), c6(MAXPER), c8(MAXPER)
      REAL c7(MAXPER), c9(MAXPER), gamma1(MAXPER), gamma2(MAXPER)
      REAL phi1(MAXPER), phi2(MAXPER), phi3(MAXPER), phi4(MAXPER)
      REAL phi5(MAXPER), phi6(MAXPER), phi1jp(MAXPER), phi5jp(MAXPER), phi6jp(MAXPER)
      REAL tau1(MAXPER), tau2(MAXPER), sigma1(MAXPER), sigma2(MAXPER)
      REAL sigma3(MAXPER), c9a(MAXPER)
      REAL c3(MAXPER), gm(MAXPER), c9b(MAXPER)
      Real CHM(MAXPER), C1c(MAXPER), C1d(MAXPER), C7b(MAXPER), C8b(MAXPER)
      real C11(MAXPER), C11B(MAXPER)
      real gscaleJapIt(MAXPER), gscaleWen(MAXPER), sigma2Jap(MAXPER)
 
      REAL c1T, c1aT, c1bT, cnT, cmT, c5T, c6T, c7T, c9T, c9aT, c3T, c8T
      REAL gamma1T, gamma2T, phi1T, phi2T, phi3T, phi4T, sigma3T
      REAL phi5T, phi6T, tau1T, tau2T, sigma1T, sigma2T
      real C1cT, C1dT, C7bT, C11T, C11bT, CHMT
      REAL c2, c4, c4a, cRB, pi, d2r, gamma, gmT, c9bT
      real cc, cosdelta, r1, r2, r3, r4, hw, psa_ref, psa, a, b, c
      integer iflag, count1, count2, hwflag, vs30_class, regionflag
      REAL M, RRUP, RBJF, DTOR, Delta, specT, sigma, Ftype, Rx
      REAL period2, lnY, F_RV, F_NM, depthvs10, tau, rkdepth
      real c8a, c8bT, deltaZ1, fd
      real NL0, sigma_NL0, F_Measured, F_Inferred, mz_TOR, deltaZ_TOR, coshM
      real gscaleJapItT, gscaleWenT, sigma2JapT, phi
      real phi1jpT, phi5jpT, phi6jpT
 
 
      Data period(1:25) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2,
     1      3, 4, 5, 7.5, 10/
      Data c1(1:25) / -1.31266420077258, -1.31266420077258, -1.27052159841083, -1.15311329373526, -1.00621910291969,
     1      -0.847745315342867, -0.506698132668687, -0.34011520483206, -0.304795250685725, -0.376848741110481, -0.437416783354322,
     1      -0.536485802809589, -0.705314012467511, -0.884076206174429, -1.25954659234479, -1.57601313090642, -2.13561675748938,
     1      -2.43819675871175, -2.86642751914157, -3.04157367665102, -3.45879762569914, -3.67827599612404, -3.80553683316369,
     1      -4.14696033561859, -4.59281925797601/
      Data c1a(1:25) / 0.18246709014369, 0.18246709014369, 0.178309224403595, 0.158716871514361, 0.136398107463104,
     1      0.11936121227405, 0.107073183963987, 0.159893513444579, 0.222972184210787, 0.297692608763902, 0.311114442939766,
     1      0.310324821840864, 0.369412804286422, 0.434878712929772, 0.489932697682775, 0.52782991756011, 0.459542012425824,
     1      0.363115431495551, 0.234026514862835, 0.173812951157414, 0.064526297612305, -0.0890925475898438, -0.168344291391401,
     1      -0.389030032261896, -0.365470360970178/
      Data c1b(1:25) / -0.167760549207099, -0.167760549207099, -0.151515173695688, -0.163132859539498, -0.162472505617202,
     1      -0.221893645353074, -0.284833283451462, -0.269908988021925, -0.253098159297076, -0.167316655333099, -0.128183304040673,
     1      -0.0815500335830345, -0.0388320763843968, -0.0710297213106874, -0.0831223427861709, -0.117869383232857,
     1      -0.168476132541299, -0.338866311701415, -0.453619201467275, -0.452653695489928, -0.320455666547956,
     1      -0.0163505888400868, -0.0385758573257073, -0.0688076694199481, -0.201948320081607/
      Data c1c(1:25) / 0.108289511741018, 0.108289511741018, 0.117493112308231, 0.166099933928779, 0.236069585545508,
     1      0.287333206939413, 0.290840813361722, 0.140626983543901, -0.00892928682460603, -0.178471416007927, -0.202091399690708,
     1      -0.179864089834803, -0.296685382036618, -0.422875653161692, -0.484879731863843, -0.532238646990341, -0.3777514162866,
     1      -0.181032482609988, -0.0234286549494071, 0.0430130382831375, 0.172098227714921, 0.377617402576038, 0.445213852683411,
     1      0.641855063614217, 0.678886019451774/
      Data c1d(1:25) / 0.305816299049586, 0.305816299049586, 0.287611795011559, 0.323395585491379, 0.380558607840273,
     1      0.499784003610166, 0.562091379441839, 0.407925321289974, 0.393524081511871, 0.305684778605879, 0.219827912926136,
     1      0.13378053808553, 0.0446584872502858, 0.10219487164099, 0.122601530820636, 0.100196434191248, 0.141479092793457,
     1      0.353765525486546, 0.551135413760324, 0.596493950421738, 0.623873999403419, 0.543943402945403, -0.101415672062659,
     1      0.0891998953807712, 0.00297285487244325/
      Data cn(1:25) / 16.0875, 16.0875, 15.7118, 15.8819, 16.4556, 17.6453, 20.1772, 19.9992, 18.7106, 16.6246, 15.3709, 13.7012,
     1      11.2667, 9.1908, 6.5459, 5.2305, 3.7896, 3.3024, 2.8498, 2.5417, 2.1488, 1.8957, 1.7228, 1.5737, 1.5265/
      Data cM(1:25) / 4.9993, 4.9993, 4.9993, 4.9993, 4.9993, 4.9993, 5.0031, 5.0172, 5.0315, 5.0547, 5.0704, 5.0939, 5.1315,
     1      5.167, 5.2317, 5.2893, 5.4109, 5.5106, 5.6705, 5.7981, 5.9983, 6.1552, 6.2856, 6.5428, 6.7415/
      Data c3(1:25) / 1.16592367966216, 1.16592367966216, 1.12830120699638, 1.07978867952166, 1.06523941578527, 1.04218816686773,
     1      1.14153732705161, 1.14747191330801, 1.18716213766583, 1.24109745493194, 1.29977223747474, 1.40201902033596,
     1      1.51446360204558, 1.64340184698047, 1.77923120143032, 1.85464518277339, 2.21522772742916, 2.60354231791946,
     1      2.9855138517209, 3.15987852562695, 3.26159427177443, 3.30665315926803, 3.27930774108803, 3.19724048046314,
     1      3.00085309764975/
      Data c5(1:25) / 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.8305, 7.1333, 7.3621, 7.4365, 7.4972, 7.5416, 7.56,
     1      7.5735, 7.5778, 7.5808, 7.5814, 7.5817, 7.5818, 7.5818, 7.5818, 7.5818, 7.5818, 7.5818/
      Data cHM(1:25) / 3.0956, 3.0956, 3.0963, 3.0974, 3.0988, 3.1011, 3.1094, 3.2381, 3.3407, 3.43, 3.4688, 3.5146, 3.5746,
     1      3.6232, 3.6945, 3.7401, 3.7941, 3.8144, 3.8284, 3.833, 3.8361, 3.8369, 3.8376, 3.838, 3.838/
      Data c6(1:25) / 0.4908, 0.4908, 0.4925, 0.4992, 0.5037, 0.5048, 0.5048, 0.5048, 0.5048, 0.5045, 0.5036, 0.5016, 0.4971,
     1      0.4919, 0.4807, 0.4707, 0.4575, 0.4522, 0.4501, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45/
      Data c7(1:25) / 0.0206913326016277, 0.0206913326016277, 0.0207295305470051, 0.0207758273344614, 0.0213063949326634,
     1      0.0223491127548945, 0.0242401920598548, 0.0280050260478514, 0.0275296793481186, 0.0284827099580447, 0.0288451884510878,
     1      0.0273907953297194, 0.0238321506745443, 0.0213283545610859, 0.0202646590294867, 0.0195915595381851, 0.0166455490057581,
     1      0.0141714984212119, 0.0118838505771454, 0.0104635707293479, 0.00237705091080518, -0.00600975826092364,
     1      -0.0108012350729907, -0.0152387032851066, -0.0180976550452297/
      Data c7b(1:25) / 0.0193395871686318, 0.0193395871686318, 0.0196111041243464, 0.0219350089169555, 0.0238886678030938,
     1      0.0246148494646156, 0.0235295139674962, 0.0138917888720577, 0.0124974976814102, 0.00630765722218034,
     1      0.00482354753413154, 0.00397587341462028, 0.00611253593560759, 0.0087913752708122, 0.0076244392096242,
     1      0.00504899023750188, 0.00672890981765682, 0.00910206491385782, 0.00562559209974137, 0.00605716808870929,
     1      0.00765199592400173, 0.0221861811513163, 0.0312372545353795, 0.0483241134635139, 0.0652336102841609/
      Data c8b(1:25) / 0.4833, 0.4833, 1.2144, 1.6421, 1.9456, 2.181, 2.6087, 2.9122, 3.1045, 3.3399, 3.4719, 3.6434, 3.8787,
     1      4.0711, 4.3745, 4.6099, 5.0376, 5.3411, 5.7688, 6.0723, 6.5, 6.8035, 7.0389, 7.4666, 7.77/
      Data c9(1:25) / 0.9228, 0.9228, 0.9296, 0.9396, 0.9661, 0.9794, 1.026, 1.0177, 1.0008, 0.9801, 0.9652, 0.9459, 0.9196,
     1      0.8829, 0.8302, 0.7884, 0.6754, 0.6196, 0.5101, 0.3917, 0.1244, 0.0086, 0, 0, 0/
      Data c9a(1:25) / 0.1202, 0.1202, 0.1217, 0.1194, 0.1166, 0.1176, 0.1171, 0.1146, 0.1128, 0.1106, 0.115, 0.1208, 0.1208,
     1      0.1175, 0.106, 0.1061, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1/
      Data c9b(1:25) / 6.8607, 6.8607, 6.8697, 6.9113, 7.0271, 7.0959, 7.3298, 7.2588, 7.2372, 7.2109, 7.2491, 7.2988, 7.3691,
     1      6.8789, 6.5334, 6.526, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5/
      Data c11(1:25) / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      Data c11b(1:25) / -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.444,
     1      -0.3539, -0.2688, -0.1793, -0.1428, -0.1138, -0.1062, -0.102, -0.1009, -0.1003, -0.1001, -0.1001, -0.1, -0.1/
      Data gamma1(1:25) / -0.00783493702135279, -0.00783493702135279, -0.00786128761621098, -0.00778079877094438,
     1      -0.00780076858327062, -0.00810910287093641, -0.00912151012877523, -0.0106273581887152, -0.0110305348350741,
     1      -0.0122773873835144, -0.0123798742277027, -0.012084026107517, -0.0113593004922249, -0.0100391491127564,
     1      -0.00770496369836205, -0.00593665980807693, -0.00477597053626469, -0.00578040071664831, -0.00715653723401061,
     1      -0.00843770016421773, -0.00717222257275659, -0.00610799690698496, -0.00476134903941374, -0.00377468396919992,
     1      -0.00158837559327459/
      Data gamma2(1:25) / -0.0153089549093592, -0.0153089549093592, -0.0156946231637086, -0.0180473134509439, -0.0196742447234555,
     1      -0.019990050777154, -0.0178990960171151, -0.0127698636872442, -0.0102565070902433, -0.00547601025969077,
     1      -0.00389061715956447, -0.00242026662929076, -0.000968015013739525, -0.00143340361892399, -0.00209160361537034,
     1      -0.00326463289495902, -0.00426625880121567, -0.00254580133231734, 0.003706273954965, 0.00606790118927974,
     1      0.00770828253872531, 0.00687887371514652, 0.00257091239027114, 0.00221896024249813, -0.00410767355413492/
      Data gm(1:25) / 4.2542, 4.2542, 4.2386, 4.2519, 4.296, 4.3578, 4.5455, 4.7603, 4.8963, 5.0644, 5.1371, 5.188, 5.2164, 5.1954,
     1      5.0899, 4.7854, 4.3304, 4.1667, 4.0029, 3.8949, 3.7928, 3.7443, 3.709, 3.6632, 3.623/
      Data phi1(1:25) / -0.478295601468722, -0.478295601468722, -0.469350764099948, -0.46172967994569, -0.451078399262125,
     1      -0.434858964898445, -0.425178256509881, -0.430145396922332, -0.446790915722155, -0.465667035381067, -0.482521496092297,
     1      -0.493897938945111, -0.492815169820979, -0.532833891414964, -0.589495724586852, -0.62420819116978, -0.787788908491935,
     1      -0.901784713782365, -1.00639980254843, -1.02277834678832, -1.04028777564407, -1.04611389204112, -0.997886189821876,
     1      -0.933873385132495, -0.837232537244328/
      Data phi2(1:25) / -0.1417, -0.1417, -0.1364, -0.1403, -0.1591, -0.1862, -0.2538, -0.2943, -0.3077, -0.3113, -0.3062, -0.2927,
     1      -0.2662, -0.2405, -0.1975, -0.1633, -0.1028, -0.0699, -0.0425, -0.0302, -0.0129, -0.0016, 0, 0, 0/
      Data phi3(1:25) / -4.96041757793564, -4.96041757793564, -4.92276178884602, -4.91251089610306, -4.96513625402504,
     1      -5.04104295651768, -5.16134191152096, -5.18427465050738, -5.16799110487246, -5.14216868405811, -5.1228525972173,
     1      -5.09276768363488, -5.04538203042609, -5.00505091580236, -4.94414555282742, -4.90155669907002, -4.81342512480855,
     1      -4.77429914899606, -4.86562627218628, -5.34080741816899, -6.30453280595018, -6.48707320506911, -6.54311216539423,
     1      -6.59367473267582, -6.59953555531281/
      Data phi4(1:25) / -2.28130316824122, -2.28130316824122, -2.22229626176524, -2.12119730536018, -2.01259817873199,
     1      -1.90429902597891, -1.65759927434112, -1.4668018434511, -1.37369802906703, -1.32250111763273, -1.3277990635296,
     1      -1.36550006882328, -1.4629983153057, -1.5736992160094, -1.79900163048511, -2.01119988561788, -2.46330564055044,
     1      -2.83710591027157, -3.44869787781101, -3.92632479013901, -4.64152301545389, -5.22525279569103, -5.73744267614884,
     1      -6.78200407367658, -7.57134365730054/
      Data phi5(1:25) / 0.0139926199706649, 0.0139926199706649, 0.0137650452723175, 0.0140339515445724, 0.0173393648081361,
     1      0.0221061501659787, 0.0416600809735219, 0.0637396168303911, 0.0508535274577875, 0.041714225815666, 0.0309579228046502,
     1      0.025016098336011, 0.0264586215769839, 0.0225557672812644, 0.0174484052235742, -0.000523452997460243,
     1      -0.0150610705243032, -0.00345810627559283, 0.0631461836018504, 0.0955736653255159, 0.130800857675134,
     1      0.128174100854778, 0.125397135621374, 0.130974459736013, 0.115120841547345/
      Data phi6(1:25) / 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
     1      300, 300, 300, 300/
      Data tau1(1:25) / 0.4, 0.4, 0.40260697107594, 0.406254677901244, 0.409536599100115, 0.412362151356064, 0.417868294570895,
     1      0.421889928589278, 0.424427487574697, 0.427481452707345, 0.429157682348667, 0.431284929241235, 0.434101130455258,
     1      0.436302393283048, 0.439579110622969, 0.441948893180584, 0.445868511544276, 0.448355188391465, 0.451470005001172,
     1      0.453423408395737, 0.455849456889221, 0.457360780855491, 0.458424490145996, 0.460138736069957, 0.461201624142458/
      Data tau2(1:25) / 0.26, 0.26, 0.2637242443942, 0.268935254144634, 0.273623713000165, 0.277660216222949, 0.285526135101279,
     1      0.291271326556111, 0.294896410820996, 0.29925921815335, 0.301653831926667, 0.304692756058908, 0.308715900650368,
     1      0.311860561832925, 0.316541586604241, 0.319926990257977, 0.325526445063252, 0.329078840559236, 0.333528578573103,
     1      0.336319154851053, 0.339784938413173, 0.341943972650702, 0.343463557351422, 0.345912480099939, 0.347430891632083/
      Data sigma1(1:25) / 0.491174072226289, 0.491174072226289, 0.490392102266757, 0.49884629286719, 0.504893375136174,
     1      0.509564336873576, 0.517904200145505, 0.523587992281292, 0.526987000584287, 0.530830818870199, 0.532795897236441,
     1      0.535092985589259, 0.537693444822054, 0.539459005039374, 0.542235827503619, 0.543286483861406, 0.529446500271604,
     1      0.510487062001946, 0.478270612673322, 0.468082689790949, 0.461651268718433, 0.457088102848047, 0.453548636419445,
     1      0.447117215346928, 0.442554049476542/
      Data sigma2(1:25) / 0.376232982006498, 0.376232982006498, 0.376232982006498, 0.384851051972815, 0.390965669923574,
     1      0.395708536077288, 0.404326606043605, 0.410441223994364, 0.414316427806967, 0.419059293960682, 0.421719608502208,
     1      0.42517391191144, 0.429916778065155, 0.433791981877758, 0.439906599828517, 0.444649465982231, 0.453267535948548,
     1      0.459382153899307, 0.468000223865624, 0.468082689790949, 0.461651268718433, 0.457088102848047, 0.453548636419445,
     1      0.447117215346928, 0.442554049476542/
      Data sigma3(1:25) / 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7999, 0.7997, 0.7988, 0.7966, 0.7792,
     1      0.7504, 0.7136, 0.7035, 0.7006, 0.7001, 0.7, 0.7, 0.7/
      Data gscaleJapIt(1:25) / 1.581682, 1.581682, 1.574012, 1.554376, 1.550154, 1.539148, 1.480416, 1.40939, 1.368192, 1.324078,
     1      1.307137, 1.293133, 1.314989, 1.351437, 1.405064, 1.440233, 1.52797, 1.652328, 1.887186, 2.134757, 3.575186, 3.864586,
     1      3.729218, 2.376267, 1.767935/
      Data phi1Jp(1:25) / -0.684621362856666, -0.684621362856666, -0.668104490153807, -0.631364652348356, -0.585522686567435,
     1      -0.545659123699481, -0.468510173737368, -0.498485421565652, -0.560262194908129, -0.645129113225058, -0.698050600056553,
     1      -0.765310551536843, -0.846945325590174, -0.899885267655382, -0.961772613950091, -0.994485940646074, -1.02254530165975,
     1      -1.0001741238606, -0.924505251934862, -0.862638634236657, -0.788170491672007, -0.719494506504419, -0.655991993310966,
     1      -0.520223846723544, -0.406803274628334/
      Data phi5Jp(1:25) / 0.459, 0.459, 0.458, 0.462, 0.453, 0.436, 0.383, 0.375, 0.377, 0.379, 0.38, 0.384, 0.393, 0.408, 0.462,
     1      0.524, 0.658, 0.78, 0.96, 1.11, 1.291, 1.387, 1.433, 1.46, 1.464/
      Data phi6Jp(1:25) / 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800,
     1      800, 800, 800, 800/
      Data gscaleWen(1:25) / 0.75937, 0.75937, 0.760565, 0.76416, 0.767627, 0.773871, 0.795556, 0.79319, 0.776833, 0.743674,
     1      0.721908, 0.692222, 0.657901, 0.636173, 0.604851, 0.550656, 0.358158, 0.200263, 0.035633, 0, 0, 0, 0, 0, 0/
      Data sigma2jap(1:25) / 0.452832136505141, 0.452832136505141, 0.455081474488854, 0.457089607558063, 0.464232248075043,
     1      0.471593304149618, 0.502235058788957, 0.522972224766759, 0.527828908598616, 0.530427671390517, 0.530985974389397,
     1      0.531176947346726, 0.53091866775, 0.53069348684478, 0.530953922960708, 0.531256171294916, 0.530904125925903,
     1      0.53016308922439, 0.527603492956116, 0.516719512970182, 0.491681001472351, 0.468235579988731, 0.451748161403632,
     1      0.416734686252964, 0.375505067353143/
      Data c8(1:25) / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0991, 0.1982, 0.2154, 0.2154, 0.2154, 0.2154, 
     1      0.2154, 0.2154, 0.2154, 0.2154/
  
 
C Find the requested spectral period and corresponding coefficients
      nPer = 25
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T      = c1(1)
         c1aT     = c1a(1)
         c1bT     = c1b(1)
         c1cT     = c1c(1)
         c1dT     = c1d(1)
         c3T      = c3(1)
         cHMT     = CHM(1)
         cnT      = cn(1)
         cmT      = cm(1)
         c5T      = c5(1)
         c6T      = c6(1)
         c7T      = c7(1) 
         c7bT      = c7b(1) 
         c8bT      = c8b(1)
         c9T      = c9(1)
         c9aT      = c9a(1)
         c9bT      = c9b(1)
         c11T      = c11(1)
         c11bT      = c11b(1)
         gamma1T  = gamma1(1)
         gamma2T  = gamma2(1)
         gmT      = gm(1)
         phi1T    = phi1(1)
         phi2T    = phi2(1)
         phi3T    = phi3(1)
         phi4T    = phi4(1)
         phi5T    = phi5(1)
         phi6T    = phi6(1)
         tau1T    = tau1(1)
         tau2T    = tau2(1)
         sigma1T = sigma1(1)
         sigma2T = sigma2(1)
         sigma2JapT = sigma2Jap(1)
         sigma3T = sigma3(1)
         gscaleJapItT = gscaleJapIt(1)
         phi1jpT = phi1jp(1)
         phi5jpT = phi5jp(1)
         phi6jpT = phi6jp(1)
         gscaleWenT = gscaleWen(1)
         c8T = c8(1)
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
      write (*,*) 'Chiou and Youngs (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c1a(count1),c1a(count2),
     +                   specT,c1aT,iflag)
            call S24_interp (period(count1),period(count2),c1b(count1),c1b(count2),
     +                   specT,c1bT,iflag)
            call S24_interp (period(count1),period(count2),c1c(count1),c1c(count2),
     +                   specT,c1cT,iflag)
            call S24_interp (period(count1),period(count2),c1d(count1),c1d(count2),
     +                   specT,c1dT,iflag)
            call S24_interp (period(count1),period(count2),cn(count1),cn(count2),
     +                   specT,cnT,iflag)
            call S24_interp (period(count1),period(count2),cm(count1),cm(count2),
     +                   specT,cmT,iflag)
            call S24_interp (period(count1),period(count2),cHM(count1),cHM(count2),
     +                   specT,cHMT,iflag)
            call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +                   specT,c3T,iflag)
            call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +                   specT,c5T,iflag)
            call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +                   specT,c6T,iflag)
            call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +                   specT,c7T,iflag)
            call S24_interp (period(count1),period(count2),c7b(count1),c7b(count2),
     +                   specT,c7bT,iflag)
            call S24_interp (period(count1),period(count2),c8b(count1),c8b(count2),
     +                   specT,c8bT,iflag)
            call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +                   specT,c9T,iflag)
            call S24_interp (period(count1),period(count2),c9a(count1),c9a(count2),
     +                   specT,c9aT,iflag)
            call S24_interp (period(count1),period(count2),c9b(count1),c9b(count2),
     +                   specT,c9bT,iflag)
            call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +                   specT,c11T,iflag)
            call S24_interp (period(count1),period(count2),c11b(count1),c11b(count2),
     +                   specT,c11bT,iflag)
            call S24_interp (period(count1),period(count2),gamma1(count1),gamma1(count2),
     +                   specT,gamma1T,iflag)
            call S24_interp (period(count1),period(count2),gamma2(count1),gamma2(count2),
     +                   specT,gamma2T,iflag)
            call S24_interp (period(count1),period(count2),gm(count1),gm(count2),
     +                   specT,gmT,iflag)
            call S24_interp (period(count1),period(count2),phi1(count1),phi1(count2),
     +                   specT,phi1T,iflag)
            call S24_interp (period(count1),period(count2),phi2(count1),phi2(count2),
     +                   specT,phi2T,iflag)
            call S24_interp (period(count1),period(count2),phi3(count1),phi3(count2),
     +                   specT,phi3T,iflag)
            call S24_interp (period(count1),period(count2),phi4(count1),phi4(count2),
     +                   specT,phi4T,iflag)
            call S24_interp (period(count1),period(count2),phi5(count1),phi5(count2),
     +                   specT,phi5T,iflag)
            call S24_interp (period(count1),period(count2),phi6(count1),phi6(count2),
     +                   specT,phi6T,iflag)
            call S24_interp (period(count1),period(count2),tau1(count1),tau1(count2),
     +                   specT,tau1T,iflag)
            call S24_interp (period(count1),period(count2),tau2(count1),tau2(count2),
     +                   specT,tau2T,iflag)
            call S24_interp (period(count1),period(count2),sigma1(count1),sigma1(count2),
     +                   specT,sigma1T,iflag)
            call S24_interp (period(count1),period(count2),sigma2(count1),sigma2(count2),
     +                   specT,sigma2T,iflag)
            call S24_interp (period(count1),period(count2),sigma2Jap(count1),sigma2Jap(count2),
     +                   specT,sigma2JapT,iflag)
            call S24_interp (period(count1),period(count2),sigma3(count1),sigma3(count2),
     +                   specT,sigma3T,iflag)
            call S24_interp (period(count1),period(count2),gscaleJapIt(count1),gscaleJapIt(count2),
     +                   specT,gscaleJapItT,iflag)
            call S24_interp (period(count1),period(count2),phi1jp(count1),phi1jp(count2),
     +                   specT,phi1jpT,iflag)
            call S24_interp (period(count1),period(count2),phi5jp(count1),phi5jp(count2),
     +                   specT,phi5jpT,iflag)
            call S24_interp (period(count1),period(count2),phi6jp(count1),phi6jp(count2),
     +                   specT,phi6jpT,iflag)
            call S24_interp (period(count1),period(count2),gscaleWen(count1),gscaleWen(count2),
     +                   specT,gscaleWenT,iflag)
            call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +                   specT,c8T,iflag)
 1011 period1 = specT                                                                                                              
c     Set the fault mechanism term.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
         if (ftype .eq. -1) then
            F_RV = 0.0
            F_NM = 1.0
         elseif (ftype .ge. 0.5) then
            F_RV = 1.0
            F_NM = 0.0
         else
            F_RV = 0.0
            F_NM = 0.0
         endif
C     Constant terms
        c2 = 1.06
        c4 = -2.1
        c4a = -0.5
        cRB = 50.0
C        c8 = 0.2154
        c8a = 0.2695
        pi = atan(1.0)*4.0
        d2r = pi/180.0
        cc = c5T* cosh(c6T * max((M-cHMT),0.0))
        gamma = gamma1T + gamma2T/cosh(max((M-gmT),0.0))
C     Apply Regional scaling factor.
C     Regionflag = 0 Global
C     Regionflag = 1 Japan and Italy
C        Also set sigma2 equal to Japan specific value
C     Regionflag = 2 Wenchuan (note only for M7.9)
c        if (regionflag .eq. 1 ) then
c           gamma = gamma * gscaleJapItT
c           sigma2T = sigma2JapT
c        elseif (regionflag .eq. 2) then
c           gamma = gamma * gscaleWenT        
c        endif
c        cosDELTA = cos(abs(DELTA)*d2r)
c Magnitude scaling
        r1 = c1T + c2 * (M-6.0) +
     1       (c2-c3T)/cnT *
     1             alog(1.0 + exp(-cnT*(M-cMT)))
c Near-field magnitude and distance scaling
        r2 = c4 * alog(Rrup + cc)
c Distance scaling at large distance
        r3 = (c4a-c4)/2.0 *
     1            alog( Rrup*Rrup+cRB*cRB ) +
     1       Rrup * gamma
c Center Z_TOR on the Z_TOR-M relation in Chiou and Youngs (2013)
        if (F_RV.EQ.1) then
          if (M .le. 5.849) then
              mZ_TOR = 2.704*2.704
          else
              mZ_TOR = max(2.704-1.226*(M-5.849), 0.0)
              mZ_TOR = mZ_TOR * mZ_TOR
          endif
        else
          if (M .le. 4.970) then
              mZ_TOR = 2.673*2.673
          else
              mZ_TOR = max(2.673-1.136*(M-4.970), 0.0)
              mZ_TOR = mZ_TOR * mZ_TOR
          endif
        endif
        deltaZ_TOR = Dtor - 2*mZ_TOR
        
c Scaling with other source variables (F_RV, F_NM, deltaZ_TOR, and Dip)
        coshM = cosh(2*max(M-4.5,0.0))
        cosDELTA = cos(DELTA*d2r)
        r4 = (c1aT+c1cT/coshM) * F_RV +
     1       (c1bT+c1dT/coshM) * F_NM +
     1       (c7T +c7bT/coshM) * deltaZ_TOR +
     1       (c11T+c11bT/coshM)* cosDELTA**2
        
c HW effect
        if (HWFlag .eq. 0) then
           hw = 0.0
        else
         hw = c9T * (cosDELTA) * (c9aT+(1-c9aT)
     1        *tanh(abs(Rx)/c9bT)) *
     1          (1.0 - sqrt(Rbjf**2+DTor**2)/(Rrup + 1))
        endif
C     Current version of the code sets cDPP=0 (i.e., no directivity)
c Directivity effect
        cDPP = 0.0
        fd = c8T * exp(-c8aT * (M-c8bT)**2) *
     1       max(0.0, 1.0-max(0.0,Rrup-40.0)/30.0) *
     1       min(max(0.0,M-5.5)/0.8, 1.0) * cDPP
c Predicted median Sa on reference condition (Vs=1130 m/sec)
c        fd = 0.0
        psa_ref = r1+r2+r3+r4+hw+fd
C     Set Phi1, Phi5, and Phi6 term for Japan is region is requested
c        if (regionflag .eq. 1) then
c           phi1T = phi1jpT
c           phi5T = phi5jpT
c           phi6T = phi6jpT
c        endif       
        
c Linear soil amplification
        a = phi1T * min(alog(Vs/1130.0), 0.0)
c Nonlinear soil amplification
        b = phi2T *
     1      (exp(-exp(phi3T)*(min(Vs,1130.0)-360.0)) - exp(-exp(phi3T)*(1130.0-360.0)))
        c = exp(phi4T)
C Deviation from ln(Vs30) scaling: bedrock depth (Z1) effect.
c        if (regionflag .eq. 1) then
c           deltaZ1 = depthvs10*1000.0 -
c     1     exp(-5.23/2.0 * alog((VS**2.0 + 412.0**2.0)/(1360.0**2.0 + 412.0**2.0)))      
c        else
           deltaZ1 = depthvs10*1000.0 -
     1     exp(-2.63 / 4.0 * alog((VS**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))
c        endif
        rkdepth = phi5T * ( 1.0 - exp(-deltaZ1/phi6T ) )
c Sa on soil condition
        psa = psa_ref + (a + b * alog((exp(psa_ref)+c)/c)) + rkdepth
C Compute the sigma term
        NL0 = b * exp(psa_ref)/(exp(psa_ref)+c)
        tau = tau1T +
     1            (tau2T-tau1T)/1.5*(min(max(M,5.0),6.5)-5.0)
        sigma_NL0 = sigma1T +
     1              (sigma2T-sigma1T)/1.5*(min(max(M,5.0),6.5)-5.0)
C     Current code set for Measured Vs30 values (i.e., Vs30class=1)
      if (vs30_class .eq. 0) then
         F_measured = 0.0
         F_Inferred = 1.0
      elseif (vs30_class .eq. 1) then      
         F_measured = 1.0
         F_Inferred = 0.0
      endif
        sigma_NL0 = sigma_NL0 *
     1        sqrt(0.7*F_Measured+F_Inferred*sigma3T+(1+NL0)**2.0)
        sigma = sqrt((tau*(1.0+NL0))**2.0+sigma_NL0**2.0)
      phi = sigma_NL0
      tau = (tau*(1.0+NL0))
      
C     Convert ground motion to units of gals.
      lnY = psa + 6.89
      period2 = period1
      return
      end 
  
C  ***** PEER NGA-West 2 MODELS (2013) **********
c ---------------------------------------------------------------------            
C ** Idriss (NGA-2013) Horizontal **
C     PEER Report 2013/08
C        NGA-West2 Model for Estimating Average Horizontal Values of 
C            Pseudo-Absolute Spectral Accelerations Generated by
C            Crustal Earthquakes
C         I. M. Idriss
C     Notes:
C        Applicable Range (see Abstract):  
C           5 <= M <= 8.5
C           Vs>=450 m/sec
C              for Vs>1200 use Vs=1200
C           Rrup <= 150 km
C        Mechanisms: Strike-slip and Norml (0)
C                    Reverse and Oblique (1)
c ---------------------------------------------------------------------            
      Subroutine S04_I14_TW_C01 ( m, Rrup, ftype, vs30, specT,
     1                     period2, lnY, sigma, iflag )
C     Last Updated: 5/18/13
      implicit none
      integer MAXPER
      parameter (MAXPER=25)
      REAL Period(MAXPER), a1mlt675(MAXPER), a2mlt675(MAXPER), a3mlt675(MAXPER)
      REAL b1mlt675(MAXPER), b2mlt675(MAXPER)
      REAL a1(MAXPER), a2(MAXPER), a3(MAXPER), b1(MAXPER), b2(MAXPER) 
      real gam(MAXPER), phi(MAXPER), xsi(MAXPER), period1
      REAL M, Rrup, Vs30, specT, sigma
      REAL SOF, period2, lnY, PhiT, gamT, XsiT, ftype
      real a1T, a2T, a3T, b1T, b2T, a1mlt675T, a2mlt675T, a3mlt675T, b1mlt675T, b2mlt675T
      integer iflag, count1, count2, nPer, i
 
 
      Data Period(1:25) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2,
     1      3, 4, 5, 7.5, 10/
      Data a1mlt675(1:25) / 2.500233380713, 2.50558343921901, 2.54905389007533, 2.66315218302058, 3.0465390943314,
     1      3.41168122348334, 3.53305607116099, 4.21379608491631, 4.12227689763015, 3.72362656328493, 4.00303950423054,
     1      4.20873194310916, 3.95556883190431, 4.2483082536625, 4.06438168106831, 3.97839067786446, 3.2169039370047,
     1      2.03558750338934, 1.33122659044408, -0.185188948776329, -1.92219595865101, -2.75022102589597, -4.14418500864792,
     1      -5.36713918636102, -7.44125982871721/
      Data a2mlt675(1:25) / 0.501420705570071, 0.50045274486304, 0.49013841997293, 0.476630410483486, 0.364470577366067,
     1      0.260841471360512, 0.2614288686757, 0.220757752258636, 0.258514730282064, 0.334072823230007, 0.329265244826973,
     1      0.325192262976761, 0.350819255594322, 0.329717876309448, 0.379780263480321, 0.419759529827459, 0.603851798024095,
     1      0.856643686843613, 0.960706519770273, 1.14832767985058, 1.33436528125161, 1.29373670045768, 1.37441955857378,
     1      1.31373625595917, 1.41503742549729/
      Data a3mlt675(1:25) / 0.0589, 0.0589, 0.0589, 0.0589, 0.0492134623324272, 0.0417, 0.0527, 0.0442, 0.0391188387593093, 0.0329,
     1      0.0267654491028554, 0.0188, 0.0095, -0.0039, -0.0133, -0.0224, -0.0267, -0.0198, -0.0367, -0.0291, -0.0214, -0.024,
     1      -0.0202, -0.0219, -0.0035/
      Data b1mlt675(1:25) / 2.9935, 2.9935, 2.9935, 2.9935, 2.92192099200299, 2.8664, 2.9406, 3.019, 2.91472377949414, 2.7871,
     1      2.81929551534672, 2.8611, 2.8289, 2.8423, 2.83, 2.856, 2.7544, 2.7339, 2.68, 2.6837, 2.6907, 2.5782, 2.5468, 2.4478,
     1      2.3922/
      Data b2mlt675(1:25) / -0.2287, -0.2287, -0.2287, -0.2287, -0.236077537409605, -0.2418, -0.2513, -0.2516, -0.23900951196997,
     1      -0.2236, -0.223295447827801, -0.2229, -0.22, -0.2284, -0.2318, -0.2337, -0.2392, -0.2398, -0.2417, -0.245, -0.2389,
     1      -0.2514, -0.2541, -0.2593, -0.2586/
      Data xsi(1:25) / -0.418699682752246, -0.418615223959223, -0.410872380718745, -0.403366643966195, -0.390638638085288,
     1      -0.369589338311506, -0.343323354414622, -0.330917990949051, -0.349509622112699, -0.3645285353545, -0.383246394748638,
     1      -0.390867061900157, -0.393160229911773, -0.439607616873627, -0.50735066440335, -0.555860494660592, -0.750091656095736,
     1      -0.872310784811448, -0.967396910488233, -0.979257377149927, -0.987613564785966, -0.997568672611442, -0.953592248112969,
     1      -0.891337622105761, -0.794820984688996/
      Data gam(1:25) / -0.00283013920396674, -0.00282803033499828, -0.00290526852332255, -0.0035024282565208, -0.00571328650969498,
     1      -0.00758878173241265, -0.00846045072079165, -0.00769077769427285, -0.00776930447201177, -0.00738500857491979,
     1      -0.00628375110597102, -0.0046841295810579, -0.00346640983735254, -0.00322976125564237, -0.00218135066244653,
     1      -0.000916435063493218, -0.00195789280161822, -0.00282711858804776, -0.00364343441913228, -0.004785164050821,
     1      -0.00281482371882879, -0.00461977029165467, -0.00480352660804264, -0.00599454233673703, -0.0062195400597961/
      Data phi(1:25) / 0.200739641446067, 0.200813356085181, 0.199656078139181, 0.194723442703274, 0.18902048140083,
     1      0.182402028007246, 0.160818007837235, 0.155060599538377, 0.149501396019003, 0.147001197734822, 0.161553742915627,
     1      0.1945576588434, 0.215937776551892, 0.216160624099592, 0.239974142176991, 0.257520562672612, 0.277554637391253,
     1      0.297363818562307, 0.276583423687761, 0.256719938372776, 0.180567884489731, 0.0750113965684316, 0.0606074683089297,
     1      -0.107929819900768, -0.121318376812063/
      Data a1(1:25) / 9.0138, 9.0408, 9.1338, 8.48609726910026, 7.9837, 7.756, 9.4252, 9.51468239707057, 9.6242, 10.2793352298525,
     1      11.13, 11.3629, 11.7818, 11.6097, 11.4484, 10.9065, 9.8565, 8.3363, 6.8656, 4.1178, 1.8102, 0.0977, -3.0563, -4.4387,
     1      9.0138/
      Data a2(1:25) / -0.0794, -0.0794, -0.0794, -0.0794, -0.142981982713312, -0.1923, -0.1614, -0.1887, -0.133751512954654,
     1      -0.0665, -0.111443199125889, -0.1698, -0.1766, -0.2798, -0.3048, -0.2911, -0.3097, -0.2565, -0.232, -0.1226, 0.1724,
     1      0.3001, 0.4609, 0.6948, 0.8393/
      Data a3(1:25) / 0.0589, 0.0589, 0.0589, 0.0589, 0.0492134623324272, 0.0417, 0.0527, 0.0442, 0.0391188387593093, 0.0329,
     1      0.0267654491028554, 0.0188, 0.0095, -0.0039, -0.0133, -0.0224, -0.0267, -0.0198, -0.0367, -0.0291, -0.0214, -0.024,
     1      -0.0202, -0.0219, -0.0035/
      Data b1(1:25) / 2.9935, 2.9935, 2.9935, 2.9935, 2.88424486584249, 2.7995, 2.8143, 2.8131, 2.63143724413814, 2.4091,
     1      2.44595081283604, 2.4938, 2.3773, 2.3772, 2.3413, 2.3477, 2.2042, 2.1493, 2.0408, 2.0013, 1.9408, 1.7763, 1.703,
     1      1.5212, 1.4195/
      Data b2(1:25) / -0.2287, -0.2287, -0.2287, -0.2287, -0.230502146542804, -0.2319, -0.2326, -0.2211, -0.197043174656907,
     1      -0.1676, -0.167991567078541, -0.1685, -0.1531, -0.1595, -0.1594, -0.1584, -0.1577, -0.1532, -0.147, -0.1439, -0.1278,
     1      -0.1326, -0.1291, -0.122, -0.1145/
 
       
C Find the requested spectral period and corresponding coefficients
      nPer = 25
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         a1T     = a1(1)
         a2T     = a2(1)
         a3T     = a3(1)
         b1T     = b1(1)
         b2T     = b2(1)
         a1mlt675T = a1mlt675(1)
         a2mlt675T = a2mlt675(1)
         a3mlt675T = a3mlt675(1)
         b1mlt675T = b1mlt675(1)
         b2mlt675T = b2mlt675(1)
         xsiT   = xsi(1)
         gamT   = gam(1)
         phiT   = phi(1)
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
      write (*,*) 'Idriss (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),b1(count1),b1(count2),
     +                   specT,b1T,iflag)
            call S24_interp (period(count1),period(count2),b2(count1),b2(count2),
     +                   specT,b2T,iflag)
            call S24_interp (period(count1),period(count2),a1mlt675(count1),a1mlt675(count2),
     +                   specT,a1mlt675T,iflag)
            call S24_interp (period(count1),period(count2),a2mlt675(count1),a2mlt675(count2),
     +                   specT,a2mlt675T,iflag)
            call S24_interp (period(count1),period(count2),a3mlt675(count1),a3mlt675(count2),
     +                   specT,a3mlt675T,iflag)
            call S24_interp (period(count1),period(count2),b1mlt675(count1),b1mlt675(count2),
     +                   specT,b1mlt675T,iflag)
            call S24_interp (period(count1),period(count2),b2mlt675(count1),b2mlt675(count2),
     +                   specT,b2mlt675T,iflag)
            call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +                   specT,phiT,iflag)
            call S24_interp (period(count1),period(count2),gam(count1),gam(count2),
     +                   specT,gamT,iflag)
            call S24_interp (period(count1),period(count2),xsi(count1),xsi(count2),
     +                   specT,xsiT,iflag)
 1011 period1 = specT                                                                                                              
C.....Compute the Ground motion.......
C.....Set the mechanism term.....................
C     Strike-slip and normal events --> SOF = 0
C     Reverse and oblique events     --> SOF = 1
C     Otherwise assume SOF = 0
      if (ftype .gt. 0.0) then
         SOF = 1.0
      else
         SOF = 0.0
      endif
      
      a1T = a1mlt675T + (a2mlt675T - a2T) * 6.75
   
      if (m .le. 6.75) then
         lnY = a1mlt675T + a2mlt675T*m + a3mlt675T*(8.5 - m)**2.0 - (b1mlt675T+b2mlt675T*m) * alog(Rrup+10.0) +
     1         xsiT*alog(Vs30) + gamT*rRup + SOF*phiT
      else
         lnY = a1T + a2T*m + a3T*(8.5 - m)**2.0 - (b1T+b2T*m) * alog(Rrup+10.0) +
     1         xsiT*alog(Vs30) + gamT*rRup + SOF*phiT
      endif
C     Convert ground motion to units of gals.
      lnY = lnY + 6.89
C     Compute Sigma which is Period and magnitude dependent.
C     Note report does not state a limit on sigma for M<5 but 
C     Since model is only applicable for M>=5 a limit is retained
C     for sigma with M<5 equal to M=5 values. 
      if (specT .le. 0.05) then
         if (m .le. 5.0) then
            sigma = 1.18 + 0.035*alog(0.05) - 0.06*5.0
         elseif (m .ge. 7.5) then
            sigma = 1.18 + 0.035*alog(0.05) - 0.06*7.5
         else
            sigma = 1.18 + 0.035*alog(0.05) - 0.06*m
         endif
      elseif (specT .ge. 3.0) then
         if (m .le. 5.0) then
            sigma = 1.18 + 0.035*alog(3.00) - 0.06*5.0
         elseif (m .ge. 7.5) then
            sigma = 1.18 + 0.035*alog(3.00) - 0.06*7.5
         else
            sigma = 1.18 + 0.035*alog(3.00) - 0.06*m
         endif
      else
         if (m .le. 5.0) then
            sigma = 1.18 + 0.035*alog(specT) - 0.06*5.0
         elseif (m .ge. 7.5) then
            sigma = 1.18 + 0.035*alog(specT) - 0.06*7.5
         else
            sigma = 1.18 + 0.035*alog(specT) - 0.06*m
         endif
      endif
      period2 = period1
      return
      end 

c ---------------------------------------------------------------------------            

      subroutine S04_CB14_TW_C01 ( mag, Rrup, Rbjf, Ftype, specT, 
     1                     period2, lnY, sigma, iflag, vs,
     2                     depthtop, D25, Dip, depth, HWflag, Rx, rupwidth, regionflag, phi, tau ) 

C     Last Updated: 5/17/17
C     Coefficients updated from PEER Report version to be consistent with EQ Spectra paper in press. 
C     Minor change to T=5, 7.5, and 10 sec for coefficient C6

      parameter (MAXPER=25)
      REAL Period(MAXPER), C0(MAXPER), C1(MAXPER), C2(MAXPER), C3(MAXPER), C4(MAXPER), C5(MAXPER)
      REAL C6(MAXPER), C7(MAXPER), C8(MAXPER), C9(MAXPER), C10(MAXPER), C11(MAXPER), C12(MAXPER)
      REAL C13(MAXPER), C14(MAXPER), C15(MAXPER), C16(MAXPER), C17(MAXPER), C18(MAXPER), C19(MAXPER)
      REAL A2(MAXPER), h1(MAXPER), h2(MAXPER), h3(MAXPER), h4(MAXPER)
      REAL h5(MAXPER), h6(MAXPER)
      REAL K1(MAXPER), K2(MAXPER), K3(MAXPER)
      REAL C20(MAXPER), DC20CA(MAXPER), DC20JP(MAXPER), DC20CH(MAXPER)
      REAL T1(MAXPER), T2(MAXPER), phi1(MAXPER), Phi2(MAXPER), phic(MAXPER)
      REAL flnAF(MAXPER), rho(MAXPER)

      REAL MAG, RRUP, RBJF, VS, D25, FHWR, FHWM, FHWZ, FHWD, PGAROCK, C, N
      real lnY, ftype, Dip, pgasoil, Rx, R1, R2, f1, f2
      real fhypH, D25_RK
      INTEGER count1, count2, HWFlag, regionflag, iflag

      real c0T, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, c10T, c11T, c12T
      real c13T, c14T, c15T ,c16T, c17T, c18T, c19T, c20T, Dc20CAT, Dc20JPT, Dc20CHT
      real k1T, k2T, k3T, a2T, h1T, h2T, h3T, h4T, h5T, h6T
      real t1T, t2T, phi1T, phi2T, phicT
      real rhoT, flnAFT

      real alpha, tau, depthtop, depth, rupwidth, specT, period2, sigma
      real tau_lnyB, tau_lnPGAB, phi_lnY, phi_lnyB, phi_lnPGAB, phi, sigmatot

C.....MODEL COEFFICIENTS.....................

      Data Period(1:25) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2,
     1      3, 4, 5, 7.5, 10/
      Data c0(1:25) / 0.356117232783958, 0.386758773602206, 0.431142669934497, 0.512003932121085, 0.649473207725169,
     1      0.83030930999764, 1.1884910081139, 1.3219880004909, 1.34639286130098, 1.29509898875361, 1.22285289551775,
     1      1.15727406735158, 0.909093052011905, 0.641421467723059, 0.318531776548646, 0.0272076106161339, -0.843846880903822,
     1      -1.43885245652115, -2.44958793225659, -2.94483330428711, -3.72583189029108, -4.39207737933935, -4.93632451517939,
     1      -5.39658678253602, -6.08834234327286/
      Data c1(1:25) / 0.559096866797525, 0.642869120021438, 0.610036408494032, 0.52214613681899, 0.437792237469586,
     1      0.358623916113537, 0.51578122821281, 0.688524087477354, 0.869930448885324, 1.14455435981339, 1.20493625362397,
     1      1.40329281705379, 1.4904299432736, 1.63350604615105, 1.98913041881546, 2.2979390413256, 2.79162470458416,
     1      2.97138814888335, 3.0269458557659, 2.98088696061972, 2.8427912585638, 2.93520078299344, 3.22558297406824,
     1      3.61157935351029, -4.46691141314152/
      Data c2(1:25) / 0.56879460840871, 0.548250606173165, 0.544503251398314, 0.546550656490723, 0.548607737985833,
     1      0.534084853902167, 0.566967391899305, 0.603886489691274, 0.66348181582985, 0.760818638761807, 0.81016622646059,
     1      0.819292926892119, 0.772660670116403, 0.796485459698103, 0.879449399855572, 0.949861255467957, 1.22026103137552,
     1      1.55566628762785, 1.91095548183813, 2.07554636000273, 2.16930529375473, 2.34782157727094, 2.58402506728671,
     1      2.04361774499736, 1.87859526803749/
      Data c3(1:25) / -0.285486382775716, -0.0545857227354485, -0.0468164681031079, -0.0591023602938874, -0.0881036367623818,
     1      -0.149826837033275, -0.230321700539041, -0.214719346117054, -0.201381168883291, -0.100368099517681,
     1      -0.0722912727125354, 0.00798058403150975, 0.153613771030137, 0.259331343435423, 0.51727423145433, 0.652041670932998,
     1      0.848887051786159, 1.06949114987441, 1.1526912631002, 1.42678748863222, 1.77988512832679, 2.13558950144393,
     1      2.52588186467635, 2.9719228836894, 3.01922347974208/
      Data c4(1:25) / -0.474, -0.474, -0.464, -0.452, -0.446368292053737, -0.442, -0.437, -0.417, -0.401261889962462, -0.382,
     1      -0.384610447190274, -0.388, -0.383, -0.37, -0.301, -0.266, -0.221, -0.123, 0.066, 0.14, 0.313, 0.467, 0.544, 0.559,
     1      0.417/
      Data c5(1:25) / -2.773, -2.773, -2.772, -2.782, -2.78706853715164, -2.791, -2.745, -2.633, -2.55430944981231, -2.458,
     1      -2.44190224232664, -2.421, -2.392, -2.376, -2.303, -2.296, -2.232, -2.158, -2.063, -2.104, -2.051, -1.986, -2.021,
     1      -2.179, -2.244/
      Data c6(1:25) / 0.248, 0.248, 0.247, 0.246, 0.242620975232242, 0.24, 0.227, 0.21, 0.197859172256757, 0.183,
     1      0.182564925468288, 0.182, 0.189, 0.195, 0.185, 0.186, 0.186, 0.169, 0.158, 0.158, 0.148, 0.135, 0.135, 0.165, 0.18/
      Data c7(1:25) / 6.768, 6.753, 6.502, 6.291, 6.30564244066028, 6.317, 6.861, 7.294, 7.62539963136186, 8.031, 8.18501638422618,
     1      8.385, 7.534, 6.99, 7.012, 6.902, 5.522, 5.65, 5.795, 6.632, 6.759, 7.978, 8.538, 8.468, 6.564/
      Data c8(1:25) / 0.0182803849428244, 0.0352717972565958, 0.0336028179581038, 0.016285534576622, -0.00650582531427639,
     1      -0.0345641613176312, -0.0866151455252148, -0.0948362008376924, -0.0766783286040763, -0.0412765896403683,
     1      -0.0216181596004451, 0.0350874237327737, 0.100288157761276, 0.123377479855529, 0.177587576029485, 0.224540830616795,
     1      0.245077396356484, 0.220376520456257, 0.204137207098202, 0.205066753804615, 0.169509715267224, 0.0906713960479872,
     1      0.0509663089420929, -0.0518868818629973, -0.0615381415093379/
      Data c9(1:25) / -0.125687581469247, -0.117608954989698, -0.104761537320645, -0.114299550101683, -0.109461799739078,
     1      -0.14856799063587, -0.194842329812706, -0.224015245432348, -0.197576107268713, -0.120680493944779, -0.09533913761543,
     1      -0.0658759795655963, -0.0210037837890234, -0.0480314832976, -0.0478090629350223, -0.0953914692233888,
     1      -0.171312750492123, -0.314868168006307, -0.361138481474666, -0.332504676330922, -0.193596524311098, 0.128500586959169,
     1      -0.114089568223998, 0.0689814675704926, 0.129699911937625/
      Data c10(1:25) / 0.72, 0.72, 0.73, 0.759, 0.796732443239964, 0.826, 0.815, 0.831, 0.794127856483483, 0.749,
     1      0.755526117975686, 0.764, 0.716, 0.737, 0.738, 0.718, 0.795, 0.556, 0.48, 0.401, 0.206, 0.105, 0, 0, 0/
      Data c11(1:25) / 0.933171905899399, 0.951465896830195, 0.996975452344212, 1.06638085646268, 1.1251338211444,
     1      1.18278442028379, 1.3522774025506, 1.54260504525065, 1.68699538923742, 1.86859564855112, 1.98109812646793,
     1      2.14309564613518, 2.37009022061089, 2.49124376382254, 2.59195688737813, 2.56205311051791, 2.06801479382536,
     1      1.43182748349166, 0.255480399259129, -0.603422596903406, -0.971549580809544, -0.983248384544924, -0.94345132737757,
     1      -0.881625094489307, -0.799650443247734/
      Data c12(1:25) / 2.186, 2.191, 2.189, 2.164, 2.14935755933972, 2.138, 2.446, 2.969, 3.2275546649024, 3.544, 3.61491714866912,
     1      3.707, 3.343, 3.334, 3.544, 3.016, 2.616, 2.47, 2.108, 1.327, 0.601, 0.568, 0.356, 0.075, -0.027/
      Data c13(1:25) / 1.42, 1.416, 1.453, 1.476, 1.51711146800772, 1.549, 1.772, 1.916, 2.02616677026276, 2.161, 2.29326265764056,
     1      2.465, 2.766, 3.011, 3.203, 3.333, 3.054, 2.562, 1.453, 0.657, 0.367, 0.306, 0.268, 0.374, 0.297/
      Data c14(1:25) / -0.0064, -0.007, -0.0167, -0.0422, -0.0557724161504944, -0.0663, -0.0794, -0.0294, 0.0126882028432437,
     1      0.0642, 0.0783834297338236, 0.0968, 0.1441, 0.1597, 0.141, 0.1474, 0.1764, 0.2593, 0.2881, 0.3112, 0.3478, 0.3747,
     1      0.3382, 0.3754, 0.3506/
      Data c15(1:25) / -0.202, -0.207, -0.199, -0.202, -0.279154398863807, -0.339, -0.404, -0.416, -0.411953057418919, -0.407,
     1      -0.365232844955611, -0.311, -0.172, -0.084, 0.085, 0.233, 0.411, 0.479, 0.566, 0.562, 0.534, 0.522, 0.477, 0.321, 0.174/
      Data c16(1:25) / 0.393, 0.39, 0.387, 0.378, 0.331256824046015, 0.295, 0.322, 0.384, 0.398838789463964, 0.417,
     1      0.411344031087739, 0.404, 0.466, 0.528, 0.54, 0.638, 0.776, 0.771, 0.748, 0.763, 0.686, 0.691, 0.67, 0.757, 0.621/
      Data c17(1:25) / 0.0527738108852, 0.0536396004239697, 0.0543820012311646, 0.0556642971677058, 0.057747107463423,
     1      0.0589149176480399, 0.0596691974981739, 0.0590915909800739, 0.0558072637267952, 0.0525013711190904, 0.0528445069680512,
     1      0.0508750262164707, 0.0479946999744034, 0.0474862395788676, 0.0446805684978132, 0.0403541904713576, 0.0338996848064585,
     1      0.033624053570051, 0.0235397196049072, 0.0213354471059939, 0.00744874217515587, 0.000940480289996758,
     1      0.000203896893899771, 0.00788547168296762, 0.0153935990038117/
      Data c18(1:25) / 0.0653258807156129, 0.0464891780106092, 0.0463973807810442, 0.0492094902793388, 0.0537950886076159,
     1      0.0588986857462317, 0.0663145055000436, 0.0671398120179699, 0.0671344775756613, 0.0650210747636337, 0.0627452520792793,
     1      0.05815787225324, 0.0443938406971418, 0.0352962862117813, 0.027447688540921, 0.0210259018767922, 0.0094635438562194,
     1      0.00651597397714056, 0.00392680184009326, -0.0108379604986975, -0.0259582614544221, -0.0478025498454571,
     1      -0.0744211978989785, -0.109831497410489, -0.0897606936410669/
      Data c19(1:25) / 0.00757, 0.00755, 0.00759, 0.0079, 0.00797321220330142, 0.00803, 0.00811, 0.00744, 0.0073140951196997,
     1      0.00716, 0.00703817913112053, 0.00688, 0.00556, 0.00458, 0.00401, 0.00388, 0.0042, 0.00409, 0.00424, 0.00448, 0.00345,
     1      0.00603, 0.00805, 0.0028, 0.00458/
      Data k1(1:25) / 865, 865, 865, 908, 990.222936015443, 1054, 1086, 1032, 962.752315834834, 878, 821.44031087739, 748, 654,
     1      587, 503, 457, 410, 400, 400, 400, 400, 400, 400, 400, 400/
      Data k2(1:25) / -1.186, -1.186, -1.219, -1.273, -1.31411146800772, -1.346, -1.471, -1.624, -1.76204570804354, -1.931,
     1      -2.04281415465008, -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0, 0, 0, 0, 0/
      Data k3(1:25) / 1.839, 1.839, 1.84, 1.841, 1.84212634158925, 1.843, 1.845, 1.847, 1.84924830143393, 1.852, 1.85374029812685,
     1      1.856, 1.861, 1.865, 1.874, 1.883, 1.906, 1.929, 1.974, 2.019, 2.11, 2.2, 2.291, 2.517, 2.744/
      Data a2(1:25) / 0.167, 0.168, 0.166, 0.167, 0.170379024767758, 0.173, 0.198, 0.174, 0.184791846882883, 0.198,
     1      0.200610447190274, 0.204, 0.185, 0.164, 0.16, 0.184, 0.216, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596/
      Data h1(1:25) / 0.241, 0.242, 0.244, 0.246, 0.248815853973132, 0.251, 0.26, 0.259, 0.256751698566066, 0.254,
     1      0.246603732960889, 0.237, 0.206, 0.21, 0.226, 0.217, 0.154, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117/
      Data h2(1:25) / 1.474, 1.471, 1.467, 1.467, 1.45686292569673, 1.449, 1.435, 1.449, 1.45439592344144, 1.461, 1.47100671422938,
     1      1.484, 1.581, 1.586, 1.544, 1.554, 1.626, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616/
      Data h3(1:25) / -0.715, -0.714, -0.711, -0.713, -0.706241950464484, -0.701, -0.695, -0.708, -0.711147622007507, -0.715,
     1      -0.717610447190274, -0.721, -0.787, -0.795, -0.77, -0.77, -0.78, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733,
     1      -0.733, -0.733/
      Data h4(1:25) / 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/
      Data h5(1:25) / -0.337, -0.336, -0.339, -0.338, -0.338, -0.338, -0.347, -0.391, -0.417080296633634, -0.449,
     1      -0.424635826224107, -0.393, -0.339, -0.447, -0.525, -0.407, -0.371, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128,
     1      -0.128, -0.128/
      Data h6(1:25) / -0.27, -0.27, -0.263, -0.259, -0.261252683178505, -0.263, -0.219, -0.201, -0.155134650747747, -0.099,
     1      -0.142072378639526, -0.198, -0.21, -0.121, -0.086, -0.281, -0.285, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756,
     1      -0.756, -0.756/
      Data c20(1:25) / -0.0055, -0.0055, -0.0055, -0.0057, -0.0060379024767758, -0.0063, -0.007, -0.0073, -0.00712013588528528,
     1      -0.0069, -0.00650843292145886, -0.006, -0.0055, -0.0049, -0.0037, -0.0027, -0.0016, -6e-04, 0, 0, 0, 0, 0, 0, 0/
      Data Dc20CA(1:25) / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      Data Dc20JP(1:25) / -0.0035, -0.0035, -0.0034, -0.0035689512383879, -0.0037, -0.0037, -0.0034, -0.00322013588528528, -0.003,
     1      -0.00304350745317124, -0.0031, -0.0033, -0.0035, -0.0034, -0.0034, -0.0032, -0.003, -0.0019, -5e-04, 0, 0, 0, 0, 0,
     1      -0.0035/
      Data Dc20CH(1:25) / 0.0036, 0.0036, 0.0037, 0.0038689512383879, 0.004, 0.0039, 0.0042, 0.0042, 0.0042, 0.00415649254682876,
     1      0.0041, 0.0036, 0.0031, 0.0028, 0.0025, 0.0016, 6e-04, 0, 0, 0, 0, 0, 0, 0, 0.0036/
      Data t1(1:25) / 0.404, 0.417, 0.446, 0.480916589266832, 0.508, 0.504, 0.445, 0.416671401932432, 0.382, 0.363291795136368,
     1      0.339, 0.34, 0.34, 0.356, 0.379, 0.43, 0.47, 0.497, 0.499, 0.5, 0.543, 0.534, 0.523, 0.466, 0.409/
      Data t2(1:25) / 0.325, 0.326, 0.344, 0.362584636222669, 0.377, 0.418, 0.426, 0.408463248815315, 0.387, 0.365681347946093,
     1      0.338, 0.316, 0.3, 0.264, 0.263, 0.326, 0.353, 0.399, 0.4, 0.417, 0.393, 0.421, 0.438, 0.438, 0.322/
      Data phi1(1:25) / 0.734, 0.738, 0.747, 0.76389512383879, 0.777, 0.782, 0.769, 0.769, 0.769, 0.765519403746301, 0.761, 0.744,
     1      0.727, 0.69, 0.663, 0.606, 0.579, 0.541, 0.529, 0.527, 0.521, 0.502, 0.457, 0.441, 0.734/
      Data phi2(1:25) / 0.492, 0.496, 0.503, 0.512573903508648, 0.52, 0.535, 0.543, 0.543, 0.543, 0.546915670785411, 0.552, 0.545,
     1      0.568, 0.593, 0.611, 0.633, 0.628, 0.603, 0.588, 0.578, 0.559, 0.551, 0.546, 0.543, 0.492/
      Data flnaf(1:25) / 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
     1      0.3, 0.3, 0.3, 0.3/
      Data phic(1:25) / 0.166, 0.166, 0.165, 0.163310487616121, 0.162, 0.158, 0.17, 0.174496602867868, 0.18, 0.182610447190274,
     1      0.186, 0.191, 0.198, 0.206, 0.208, 0.221, 0.225, 0.222, 0.226, 0.229, 0.237, 0.237, 0.271, 0.29, 0.166/
      Data rho(1:25) / 1, 0.998, 0.986, 0.958967801857936, 0.938, 0.887, 0.87, 0.872697961720721, 0.876, 0.873389552809726, 0.87,
     1      0.85, 0.819, 0.743, 0.684, 0.562, 0.467, 0.364, 0.298, 0.234, 0.202, 0.184, 0.176, 0.154, 1/
   

      nPer = 25
      c = 1.88
      n = 1.18

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c0T = c0(1)
         c1T = c1(1)
         c2T = c2(1)
         c3T = c3(1)
         c4T = c4(1)
         c5T = c5(1)
         c6T = c6(1)
         c7T = c7(1)
         c8T = c8(1)
         c9T = c9(1)
         c10T = c10(1)
         c11T = c11(1)
         c12T = c12(1)
         c13T = c13(1)
         c14T = c14(1)
         c15T = c15(1)
         c16T = c16(1)
         c17T = c17(1)
         c18T = c18(1)
         c19T = c19(1)

         a2T = a2(1)
         h1T = h1(1)
         h2T = h2(1)
         h3T = h3(1)
         h4T = h4(1)
         h5T = h5(1)
         h6T = h6(1)

         k1T = k1(1)
         k2T = k2(1)
         k3T = k3(1)
         c20T = c20(1)
         Dc20CAT = Dc20CA(1)
         Dc20JPT = Dc20JP(1)
         Dc20CHT = Dc20CH(1)
         
         phi1T = phi1(1)
         phi2T = phi2(1)
         t1T = t1(1)
         t2T = t2(1)
         flnAFT = flnaf(1)
         phicT = phic(1)
         rhoT = rho(1)

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
      write (*,*) 'Campbell&Bozorgnia (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),c0(count1),c0(count2),
     +                   specT,c0T,iflag)
            call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +                   specT,c2T,iflag)
            call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +                   specT,c3T,iflag)
            call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +                   specT,c4T,iflag)
            call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +                   specT,c5T,iflag)
            call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +                   specT,c6T,iflag)
            call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +                   specT,c7T,iflag)
            call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +                   specT,c8T,iflag)
            call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +                   specT,c9T,iflag)
            call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +                   specT,c10T,iflag)
            call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +                   specT,c11T,iflag)
            call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +                   specT,c12T,iflag)
            call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +                   specT,c13T,iflag)
            call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +                   specT,c14T,iflag)
            call S24_interp (period(count1),period(count2),c15(count1),c15(count2),
     +                   specT,c15T,iflag)
            call S24_interp (period(count1),period(count2),c16(count1),c16(count2),
     +                   specT,c16T,iflag)
            call S24_interp (period(count1),period(count2),c17(count1),c17(count2),
     +                   specT,c17T,iflag)
            call S24_interp (period(count1),period(count2),c18(count1),c18(count2),
     +                   specT,c18T,iflag)
            call S24_interp (period(count1),period(count2),c19(count1),c19(count2),
     +                   specT,c19T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),h1(count1),h1(count2),
     +                   specT,h1T,iflag)
            call S24_interp (period(count1),period(count2),h2(count1),h2(count2),
     +                   specT,h2T,iflag)
            call S24_interp (period(count1),period(count2),h3(count1),h3(count2),
     +                   specT,h3T,iflag)
            call S24_interp (period(count1),period(count2),h4(count1),h4(count2),
     +                   specT,h4T,iflag)
            call S24_interp (period(count1),period(count2),h5(count1),h5(count2),
     +                   specT,h5T,iflag)
            call S24_interp (period(count1),period(count2),h6(count1),h6(count2),
     +                   specT,h6T,iflag)

            call S24_interp (period(count1),period(count2),k1(count1),k1(count2),
     +                   specT,k1T,iflag)
            call S24_interp (period(count1),period(count2),k2(count1),k2(count2),
     +                   specT,k2T,iflag)
            call S24_interp (period(count1),period(count2),k3(count1),k3(count2),
     +                   specT,k3T,iflag)

            call S24_interp (period(count1),period(count2),c20(count1),c20(count2),
     +                   specT,c20T,iflag)
            call S24_interp (period(count1),period(count2),Dc20CA(count1),Dc20CA(count2),
     +                   specT,Dc20CAT,iflag)
            call S24_interp (period(count1),period(count2),Dc20JP(count1),Dc20JP(count2),
     +                   specT,Dc20JPT,iflag)
            call S24_interp (period(count1),period(count2),Dc20CH(count1),Dc20CH(count2),
     +                   specT,Dc20CHT,iflag)
     
            call S24_interp (period(count1),period(count2),phi1(count1),phi1(count2),
     +                   specT,phi1T,iflag)
            call S24_interp (period(count1),period(count2),phi2(count1),phi2(count2),
     +                   specT,phi2T,iflag)
            call S24_interp (period(count1),period(count2),t1(count1),t1(count2),
     +                   specT,t1T,iflag)
            call S24_interp (period(count1),period(count2),t2(count1),t2(count2),
     +                   specT,t2T,iflag)
            call S24_interp (period(count1),period(count2),flnAF(count1),flnAF(count2),
     +                   specT,flnAfT,iflag)
            call S24_interp (period(count1),period(count2),phic(count1),phic(count2),
     +                   specT,phicT,iflag)
            call S24_interp (period(count1),period(count2),rho(count1),rho(count2),
     +                   specT,rhoT,iflag)

 1011 period1 = specT                                                                                                              

C.....COMPUTE ROCK PGA VALUE FIRST.........................
C.....MAGNITUDE DEPENDENCE (Eq 3.2)........................
      IF (MAG .LE. 4.5) THEN
         TERM1 = C0(1) + C1(1)*(MAG-4.5)
      elseif (mag .le. 5.5) then
         TERM1 = C0(1) + c2(1)*(MAG-4.5)
      elseif (mag .le. 6.5) then
         TERM1 = C0(1) + c2(1) + c3(1)*(mag-5.5)
      ELSE
         TERM1 = C0(1) + c2(1) + C3(1) + c4(1)*(mag-6.5)
      ENDIF
      
C.....Distance dependence (Eq 3.3).....
      R = SQRT( RRUP*RRUP+C7(1)*C7(1) )
      TERM2 = (C5(1) + C6(1)*MAG)*ALOG(R)

C.....SET UP STYLE OF FAULTING TERMS (Eq 3.4, 3.5, and 3.6).........

C     Set mechanism term and corresponding Frv and Fnm values.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C    -1,-0.5    Normal and NMl/Obl       -150 < Rake < -30.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C       0       Strike-Slip                    Otherwise
      IF (Ftype .EQ. 0.0) THEN
         TERM3 = 0.0
      ELSEIF (Ftype .ge. 0.5) THEN
         if (mag .le. 4.5) then
            TERM3 = 0.0
         elseif (mag .le. 5.5) then
            TERM3 = C8(1)*(mag-4.5)
         else
            TERM3 = c8(1)
         endif
      ELSEIF (Ftype .le. -0.5) THEN
         if (mag .le. 4.5) then
            TERM3 = 0.0
         elseif (mag .le. 5.5) then
            TERM3 = C9(1)*(mag-4.5)
         else
            TERM3 = C9(1)
         endif
      ENDIF

C.....SET UP HANGING WALL TERMS (Eq 3.7)..............
      if (HWflag .eq. 1) then
         R1 = rupwidth*cos(abs(dip)*3.14159/180.0)
         R2 = 62.0*mag - 350.0
         f1 = h1(1) + h2(1)*(Rx/R1) + h3(1)*(Rx/R1)**2.0
         f2 = h4(1) + h5(1)*((Rx-R1)/(R2-R1)) + h6(1)*((Rx-R1)/(R2-R1))**2.0
         if (Rrup .eq. 0.0) then
            fhwrrup = 1.0
         else
            fhwrrup = ((Rrup-Rbjf)/Rrup)         
         endif
         if (Rx .lt. R1) Then
            fhwr = f1*fhwrrup
         else
            fhwr = max(f2,0.0)*fhwrrup
         endif
         if (mag .le. 5.5) then
            fhwm = 0.0         
         elseif (mag .le. 6.5) then
            fhwm = (mag-5.5)*(1.0+a2(1)*(mag-6.5))
         else
            fhwm = 1.0 + a2(1)*(mag-6.5)        
         endif
         if (depthtop .le. 16.66) then
            fhwz = 1.0 - 0.06*depthtop 
         else
            fhwz = 0.0
         endif
         fhwd = (90.0 - dip)/45.0        
         TERM4 = c10(1)*fhwr*fhwm*fhwz*fhwd
      else   
         term4 = 0.0
      endif 

C.....NOW COMPUTE THE SITE CONDITION FACTORS...............
C.....(FOR PGA ROCK, VS=1100, i.e., Vs>k1)
      TERM5_RK = (C11(1) + K2(1)*n)*ALOG( 1100.0/K1(1) )

C.....NOW COMPUTE THE SEDIMENT DEPTH DEPENDENCE (Eq 3.17)............
C     For Rock PGA the D25 value should be set at the recommended value of D25=0.398
      D25_RK = 0.398
      
      TERM6_RK = C14(1)*(D25_RK-1.0)
 
C.....Now compute the hypocentral depth term (Eq 3.21).........
      if (depth .le. 7.0) then
         fhypH = 0.0
      elseif (depth .le. 20.0) then
         fhypH = depth - 7.0      
      else
         fhypH = 13.0
      endif
      if (mag .le. 5.5) then
          term7 = c17(1)*fhypH
      elseif (mag .le. 6.5) then
          term7 = (c17(1) + (c18(1)-c17(1))*(mag-5.5))*fhypH
      else
          term7 = c18(1)*fhypH     
      endif

C.....Compute Rupture Dip term (Eq 3.24)............
      if (mag .le. 4.5) then
          term8 = c19(1)*dip
      elseif (mag .le. 5.5) then
          term8 = c19(1)*(5.5-mag)*dip      
      else 
          term8 =0
      endif

C.....Compute anelastic attenuation term.....
      if (Rrup .le. 80.0) then
         term9 = 0.0
      else      
         term9 = (c20(1)+Dc20CA(1) ) * (Rrup-80.0)

      endif      
      
      PGAROCK = EXP(TERM1+TERM2+TERM3+TERM4+TERM5_RK+TERM6_RK+TERM7+TERM8+TERM9)
C    write(*,*) "fmag  = ", TERM1
C    write(*,*) "fdis  = ", TERM2
C    write(*,*) "fflt  = ", TERM3
C    write(*,*) "fhng  = ", TERM4
C    write(*,*) "fsite = ", TERM5_RK
C    write(*,*) "fsed  = ", TERM6_RK
C    write(*,*) "fhyp  = ", TERM7
C    write(*,*) "fdip  = ", TERM8
C    write(*,*) "fatn  = ", TERM9
            
C.....For PGA Specific Vs30m Value
      if (vs .le. k1(1) ) then
         term5 = c11(1)*alog(vs/k1(1)) + 
     1           k2(1)*(alog(pgarock+c*((vs/k1(1))**n)) - 
     2           alog(pgarock+c))
      else
         term5 = (c11(1) + k2(1)*n)*alog(vs/k1(1))
      endif


C.....NOW COMPUTE THE SEDIMENT DEPTH DEPENDENCE (Eq 3.17)............
C     For Rock PGA the D25 value should be set at the recommended value of D25=0.398
      if (D25 .le. 1.0) then
            TERM6 = C14(1)*(D25-1.0)
      elseif (D25 .GT. 1.0 .AND. D25 .LE. 3.0) then
        TERM6 = 0.0
      elseif (D25 .GT. 3.0) then
         TERM6 = c16(1)*k3(1)*exp(-0.75)*(1.0-exp(-0.25*(D25-3.0)))
      endif

      pgasoil = alog(pgarock) - term5_rk - term6_RK + term5 + term6
      psoil2 = (TERM1+TERM2+TERM3+TERM4+TERM5+TERM6+TERM7+TERM8+TERM9)
c    write(*,*) "PGAROCK = ", PGAROCK
c    write(*,*) "pgasoil = ", pgasoil
c    write(*,*) "psoil2  = ", psoil2


C.....NOW COMPUTE THE GROUND MOTION VALUES.................
C.....MAGNITUDE DEPENDENCE.................................
      IF (MAG .LE. 4.5) THEN
         TERM1 = C0T + C1T*(MAG-4.5)
      elseif (mag .le. 5.5) then
         TERM1 = C0T + c2T*(MAG-4.5)
      elseif (mag .le. 6.5) then
         TERM1 = C0T + c2T + c3T*(mag-5.5)
      ELSE
         TERM1 = C0T + c2T + C3T + c4T*(mag-6.5)
      ENDIF

C.....Distance dependence......
      R = SQRT( RRUP*RRUP+C7T*C7T )
      TERM2 = (C5T + C6T*MAG)*ALOG(R)

C.....SET UP STYLE OF FAULTING TERMS...........

C     Set mechanism term and corresponding Frv and Fnm values.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C    -1,-0.5    Normal and NMl/Obl       -150 < Rake < -30.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C       0       Strike-Slip                    Otherwise
      IF (Ftype .EQ. 0.0) THEN
         TERM3 = 0.0
      ELSEIF (Ftype .ge. 0.5) THEN
         if (mag .le. 4.5) then
            TERM3 = 0.0
         elseif (mag .le. 5.5) then
            TERM3 = C8T*(mag-4.5)
         else
            TERM3 = C8T
         endif
      ELSEIF (Ftype .le. -0.5) THEN
         if (mag .le. 4.5) then
            TERM3 = 0.0
         elseif (mag .le. 5.5) then
            TERM3 = C9T*(mag-4.5)
         else
            TERM3 = C9T
         endif
      ENDIF

C.....SET UP HANGING WALL TERMS................
      if (HWflag .eq. 1) then
         R1 = rupwidth*cos(abs(dip)*3.14159/180.0)
         R2 = 62.0*mag - 350.0
         f1 = h1T + h2T*(Rx/R1) + h3T*(Rx/R1)**2.0
         f2 = h4T + h5T*((Rx-R1)/(R2-R1)) + h6T*((Rx-R1)/(R2-R1))**2.0
         if (Rrup .eq. 0.0) then
            fhwrrup = 1.0
         else
            fhwrrup = ((Rrup-Rbjf)/Rrup)         
         endif
         if (Rx .lt. R1) Then
            fhwr = f1*fhwrrup
         else
            fhwr = max(f2,0.0)*fhwrrup
         endif
         if (mag .le. 5.5) then
            fhwm = 0.0         
         elseif (mag .le. 6.5) then
            fhwm = (mag-5.5)*(1.0+a2T*(mag-6.5))
         else
            fhwm = 1.0 + a2T*(mag-6.5)        
         endif
         if (depthtop .le. 16.66) then
            fhwz = 1.0 - 0.06*depthtop 
         else
            fhwz = 0.0
         endif
         fhwd = (90.0 - dip)/45.0        
         TERM4 = c10T*fhwr*fhwm*fhwz*fhwd
      else   
         term4 = 0.0
      endif 

C.....NOW COMPUTE THE SITE CONDITION FACTORS...............
      IF (VS .LE. K1T ) THEN
         TERM5 = C11T*ALOG( VS/K1T ) +
     1           K2T*( ALOG( PGAROCK+c*( (VS/K1T)**N) ) -
     2           ALOG( PGAROCK+c ) )
      ELSE
         TERM5 = ( C11T+K2T*n )*ALOG( VS/K1T )
      ENDIF

C.....NOW COMPUTE THE SEDIMENT DEPTH DEPENDENCE.............
      IF (D25 .LE. 1.0) THEN
         TERM6 = C14T*(D25-1.0)
      ELSEIF (D25 .GT. 1.0 .AND. D25 .LE. 3.0) THEN
         TERM6 = 0.0
      ELSEIF (D25 .GT. 3.0)  THEN
         TERM6 = c16T*k3T*exp(-0.75)*( 1.0 - exp(-0.25*(D25-3.0)))
      ENDIF
      
C.....Now compute the hypocentral depth term..........
      if (depth .le. 7.0) then
         fhypH = 0.0
      elseif (depth .le. 20.0) then
         fhypH = depth - 7.0      
      else
         fhypH = 13.0
      endif
      if (mag .le. 5.5) then
          term7 = c17T*fhypH
      elseif (mag .le. 6.5) then
          term7 = (c17T + (c18T-c17T)*(mag-5.5))*fhypH
      else
          term7 = c18T*fhypH     
      endif

C.....Compute Rupture Dip term.............
      if (mag .le. 4.5) then
          term8 = c19T*dip
      elseif (mag .le. 5.5) then
          term8 = c19T*(5.5-mag)*dip      
      else 
          term8 = 0
      endif

C.....Compute anelastic attenuation term.....
      if (Rrup .le. 80.0) then
         term9 = 0.0
      else
         term9 = (c20T+Dc20CAT)*(Rrup-80.0)
      endif

      LnY = (TERM1+TERM2+TERM3+TERM4+TERM5+TERM6+TERM7+TERM8+TERM9)

c   write(*,*) "fmag  = ", TERM1
c   write(*,*) "fdis  = ", TERM2
c   write(*,*) "fflt  = ", TERM3
c   write(*,*) "fhng  = ", TERM4
c   write(*,*) "fsite = ", TERM5
c   write(*,*) "fsed  = ", TERM6
c   write(*,*) "fhyp  = ", TERM7
c   write(*,*) "fdip  = ", TERM8
c   write(*,*) "fatn  = ", TERM9
c
c   write(*,*) "LnY = ", LnY
c   write(*,*) "Sa = ", exp(LnY)

C    Check that SA is not less than PGA for T<0.25sec
c     if (specT .lt. 0.25) then
c        if (lnY .lt. pgasoil ) then
c           lnY = pgasoil
c        endif
c     endif

C.....Now compute the sigma value..........
      IF (Vs .LT. k1T) THEN
        alpha = k2T*pgarock*(1/(pgarock  
     &    +c*(Vs/k1T)**n) 
     &    -1/(pgarock + c))
      ELSE
        alpha = 0.0
      ENDIF

      If (Mag.le.4.5) then
      tau_lnyB = t1T
       tau_lnPGAB = t1(1)
      elseif (Mag.lt.5.5) then
      tau_lnyB = t2T + 
     &          (t1T - t2T)*(5.5-mag)
      tau_lnPGAB = t2(1) + 
     &          (t1(1) - t2(1))*(5.5-Mag)
      else
      tau_lnyB = t2T
      tau_lnPGAB = t2(1)
      endif

      tau = SQRT(tau_lnyB**2 +  
     &           (alpha * tau_lnPGAB)**2 +
     &           2.0*alpha*rhoT*tau_lnyB*tau_lnPGAB)

      If (Mag.le.4.5) then
          phi_lny = phi1T
           phi_lnPGAB = phi1(1)
      elseif (Mag.lt.5.5) then
          phi_lny = phi2T + 
     &          (phi1T - phi2T)*(5.5-mag)
          phi_lnPGAB = phi2(1) + 
     &          (phi1(1) - phi2(1))*(5.5-mag)
      else
          phi_lny = phi2T 
          phi_lnPGAB = phi2(1) 
      endif

      phi_lnyB = SQRT(phi_lny**2 - flnAFT**2)

      phi_lnPGAB = SQRT(phi_lnPGAB**2 - flnAF(1)**2)

      phi = SQRT(phi_lny**2 + 
     &           (alpha*phi_lnPGAB)**2 +
     &           2.0*alpha*rhoT*phi_lnyB*phi_lnPGAB)
 
      Sigmatot = SQRT(phi**2 + Tau**2)

      period2 = period1
      
C     Convert ground motion to units of gals.

      lnY = lnY + 6.89
      sigma = sigmaTot

      return
      END

c--------------------Adjust E Version --------------------------------------------------  

c ---------------------------------------------------------------------            
C     *** Akkar, Sandikkaya, and Bommer (2013) *** Adjusted in Taiwan SSHAC Project
c ---------------------------------------------------------------------            
      subroutine S04_ASB14_TW_E02 ( mag, Rbjf, specT, 
     1                     period2, lnY, sigma, iflag, ftype, Vs, phiT, tauT ) 
      implicit none
      integer MAXPER
      parameter (MAXPER=22)
      REAL Period(MAXPER), a1(MAXPER), a3(MAXPER), a4(MAXPER), a8(MAXPER)
      Real a9(MAXPER), b1(MAXPER), b2(MAXPER), phi(MAXPER), tau(MAXPER)
      real specT, a1T, a3T, a4T, a8T, a9T, b1T, b2T, phiT, tauT, period2
      real mag, Rbjf, Ftype, Fn, Fr, Vs, lnY
      real a5, a6, a7, c1, c, n, sigma, period1, pgaref
      INTEGER iFlag, count1, count2, nPer, i
      real a2(MAXPER), a2T
 
 
      Data Period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4/
      Data a1 / 2.03945675856391, 2.0105755214897, 2.09441962806437, 2.19342029563189, 2.31904682103906, 2.45390770646749,
     1      2.74068576984133, 3.02725901648496, 3.14444959319828, 3.27802577633879, 3.21237091996084, 3.11514512674121,
     1      2.95906158340172, 2.7808944616179, 2.37220179517804, 2.15108279768441, 1.3774519534676, 0.898857916754095,
     1      0.2513574194084, -0.248428626839261, -0.937137489057441, -1.4545229162242/
      Data a3 / -0.02807, -0.0274, -0.02715, -0.02403, -0.01797, -0.01248, -0.00532, -0.00925, -0.0149516924364565, -0.02193,
     1      -0.0274510958074301, -0.03462, -0.0467825076620119, -0.05672, -0.07684, -0.0949, -0.12347, -0.14345, -0.17187,
     1      -0.19029, -0.21392, -0.23848/
      Data a4 / -1.23452, -1.23698, -1.25363, -1.27525, -1.30123, -1.32632, -1.35722, -1.38182, -1.37407684986153, -1.3646,
     1      -1.33160829826025, -1.28877, -1.22380239685518, -1.17072, -1.0653, -1.01909, -0.88393, -0.81838, -0.75751, -0.72033,
     1      -0.69085, -0.66482/
      Data a8 / 0.104553265468639, 0.0120117886289134, 0.00735113146215132, -0.000968056992192673, 0.00592806686242958,
     1      0.000506378970806889, 0.0195584168540039, 0.0199008069526318, 0.0711615823806575, 0.0909863018067894,
     1      0.122980430881693, 0.122293459202686, 0.0815956680067303, 0.0826805756231909, 0.0464275875594828, -0.0224107102290375,
     1      -0.0530221693454976, -0.0518991410288072, -0.0838484024943528, -0.129115027502304, -0.149034891946652,
     1      -0.176695513944899/
      Data a9 / 0.248308147882621, 0.226557955393241, 0.227183473276702, 0.224971695893834, 0.231865222854536, 0.233678929149701,
     1      0.231586777371364, 0.24429370127944, 0.250876857540622, 0.239927616870316, 0.243232642934037, 0.252967738138007,
     1      0.244561507444697, 0.23142234226282, 0.196576492581241, 0.162932138512573, 0.163598922854411, 0.16397439996027,
     1      0.136763456368592, 0.10290935123472, 0.020715592797924, -0.0518151815751219/
      Data b1 / -0.50847449458028, -0.538549867425188, -0.530969716151637, -0.51477742792559, -0.4874004233809, -0.454349946166605,
     1      -0.398928162316836, -0.447336306798626, -0.477689733960318, -0.529059553231679, -0.556558326532214, -0.578103741261652,
     1      -0.59884406919878, -0.660851743331629, -0.696131436888087, -0.707234728955517, -0.830868141852288, -0.937681728021304,
     1      -1.00546910992409, -0.977458388978301, -0.98301709636186, -0.97009213939352/
      Data b2 / -0.28846, -0.28685, -0.28241, -0.26842, -0.24759, -0.22385, -0.17525, -0.29293, -0.339056152218589, -0.39551,
     1      -0.417668345900112, -0.44644, -0.452416689285495, -0.4573, -0.43008, -0.37408, -0.28957, -0.28702, -0.24695, -0.17336,
     1      -0.13336, -0.07749/
      Data phi / 0.6201, 0.6215, 0.6266, 0.641, 0.6534, 0.6622, 0.6626, 0.667, 0.672665719613514, 0.6796, 0.673030374571143,
     1      0.6645, 0.661968437319219, 0.6599, 0.6697, 0.6512, 0.6744, 0.6787, 0.7164, 0.7254, 0.6997, 0.6196/
      Data tau / 0.3501, 0.3526, 0.3555, 0.3565, 0.3484, 0.3551, 0.3759, 0.4067, 0.39887591100991, 0.3893, 0.387081119888267,
     1      0.3842, 0.382769116745646, 0.3816, 0.3962, 0.4021, 0.4043, 0.3943, 0.3799, 0.3717, 0.4046, 0.3566/
      Data a2(1:22) / 0.156596, 0.179855, 0.176579, 0.171114, 0.173395, 0.163795, 0.154743, 0.154450, 0.172438, 0.217818, 
     1      0.219640, 0.225918, 0.245985, 0.270563, 0.256632, 0.238360, 0.233643, 0.262352, 0.256889, 0.220679, 0.190910, 
     1      0.080137  /

 
 
C First check for the PGA case (i.e., specT=0.0) 
      nPer = 22
      if (specT .eq. 0.0) then
         period1 = period(1)
         a1T = a1(1)
         a2T = a2(1)
         a3T = a3(1)
         a4T = a4(1)
         a8T = a8(1)
         a9T = a9(1)
         b1T = b1(1)
         b2T = b2(1)
         phiT = phi(1)
         tauT = tau(1)
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
      write (*,*) 'Adjusted Akkar,Sandikkaya&Bommer (2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +                   specT,a9T,iflag)
            call S24_interp (period(count1),period(count2),b1(count1),b1(count2),
     +                   specT,b1T,iflag)
            call S24_interp (period(count1),period(count2),b2(count1),b2(count2),
     +                   specT,b2T,iflag)
            call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +                   specT,phiT,iflag)
            call S24_interp (period(count1),period(count2),tau(count1),tau(count2),
     +                   specT,tauT,iflag)
 1011 period1 = specT                                                                                                              
C.....Set the mechanism terms based on ftype............
C     Set mechanism term and corresponding Frv and Fnm values.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                    -120 < Rake <  -60
C     -0.5      Normal/Oblique            -150 < Rake < -120
C                                          -60 < Rake <  -30
C       0       Strike-Slip               -180 < Rake < -150
C                                          -30 < Rake <   30
C                                          150 < Rake <  180
C      0.5      Reverse/Oblique             30 < Rake <   60
C                                          120 < Rake <  150
C       1       Reverse                     60 < Rake <  120   
      if (ftype .eq. -1.0) then
         Fr = 0.0
         Fn = 1.0
      elseif (ftype .eq. -0.5) then 
         Fr = 0.0
         Fn = 1.0
      elseif (ftype .eq. 0.0) then 
         Fr = 0.0
         Fn = 0.0
      elseif (ftype .eq. 0.5) then
         Fr = 1.0
         Fn = 0.0
      elseif (ftype .eq. 1.0) then
         Fr = 1.0
         Fn = 0.0
      endif 
C     Set frequency independent terms
c      a2 = 0.0029
      a5 = 0.2529
      a6 = 7.5
      a7 = -0.5096
      c1 = 6.75
      c = 2.5
      n = 3.2
C     Compute the PGA for reference Vs=750m/s.
      if (mag .lt. c1 ) then
         pgaref = a1(1) + a2(1)*(mag-c1) + a3(1)*(8.5-mag)**2.0 + 
     1                 (a4(1)+a5*(mag-c1))*alog(sqrt(Rbjf*Rbjf+a6*a6)) + 
     2                  a8(1)*Fn + a9(1)*Fr
      else
         pgaref = a1(1) + a7*(mag-c1) + a3(1)*(8.5-mag)**2.0 + 
     1                 (a4(1)+a5*(mag-c1))*alog(sqrt(Rbjf*Rbjf+a6*a6)) + 
     2                  a8(1)*Fn + a9(1)*Fr      
      endif
      pgaref = exp(pgaref) 
C.....Now compute the ground motion value........
      if (mag .lt. c1 ) then
         lnY = a1T + a2T*(mag-c1) + a3T*(8.5-mag)**2.0 + 
     1                 (a4T+a5*(mag-c1))*alog(sqrt(Rbjf*Rbjf+a6*a6)) + 
     2                  a8T*Fn + a9T*Fr
      else
         lnY = a1T + a7*(mag-c1) + a3T*(8.5-mag)**2.0 + 
     1                 (a4T+a5*(mag-c1))*alog(sqrt(Rbjf*Rbjf+a6*a6)) + 
     2                  a8T*Fn + a9T*Fr      
      endif
C.....Now apply site amplification term......
      if (vs .le. 750.0) then
         lnY = lnY + b1T*alog(Vs/750.0) + 
     1         b2T*alog( (pgaref + c*(Vs/750.0)**n) / ((pgaref+c)*(Vs/750.0)**n) )
      else
         lnY = lnY + b1T*alog( min(Vs,1000.0)/750.0)
      endif
C.....Set Sigma value.........
      sigma = sqrt (phiT*phiT + tauT*tauT)
C     Convert ground motion to units of gals.
      lnY = lnY + 6.89
      period2 = period1
      return
      END


c ---------------------------------------------------------------------            
C     *** Bindi (2014) *** Adjusted in Taiwan SSHAC Project
c ---------------------------------------------------------------------            
              
      Subroutine S04_Bindi14_TW_E02 ( m, jbDist, ftype, specT,
     1                     period2, lnY, sigma, iflag, vs, phiT, tauT )
      implicit none
      integer MAXPER
      parameter (MAXPER=21)
      REAL Period(MAXPER), e1(MAXPER), c1(MAXPER), c2(MAXPER), h(MAXPER)
      REAL c3(MAXPER), b1(MAXPER), b2(MAXPER), b3(MAXPER), gamma(MAXPER)
      REAL sofN(MAXPER), sofR(MAXPER), sofS(MAXPER), phi(MAXPER), tau(MAXPER), sig(MAXPER), sigs2s(MAXPER)
      real e1T, c1T, c2T, hT, c3T, b1T, b2T, b3T, gammaT, sofNT, sofRT, sofST, sigs2sT
      real phiT, tauT, sigT, period1
      real Rref, Mref, Mh, R, Vref, vs
 
      REAL M, jbDist, specT, sigma, termsof
      REAL period2, lnY, ftype
      integer iflag, count1, count2, nPer, i
      real f_D, f_M, f_S
  
 
      Data period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3/
      Data e1 / 0.602962884594764, 0.602962884594764, 0.621561504845865, 0.726381825288257, 0.830096762727687, 0.938766589905189,
     1      1.22302346770631, 1.39715735313782, 1.39063856362439, 1.33235570611944, 1.33153232812115, 1.31633816650828,
     1      1.26708499354147, 1.27951948079044, 1.0912836412466, 0.789633867906657, 0.45939097848205, 0.383689218228279,
     1      -0.238839681034164, -0.775967717654436, -1.62622126714907/
      Data c1 / -1.26358, -1.26358, -1.26358, -1.29088019990866, -1.31025, -1.30237878943392, -1.28882132582112, -1.28178,
     1      -1.23465110534188, -1.17697, -1.14479188763455, -1.10301, -1.08792194443592, -1.10591, -1.09538, -1.05767,
     1      -1.04831378793046, -1.0527, -0.983388, -0.979215, -0.940373/
      Data c2 / 0.220527, 0.220527, 0.220527, 0.234653259429915, 0.244676, 0.239572871690745, 0.229465497782482, 0.219406,
     1      0.202883682422306, 0.182662, 0.161122330083983, 0.133154, 0.118226690912901, 0.108276, 0.101111, 0.112197,
     1      0.122345650356548, 0.103471, 0.109072, 0.163344, 0.227241/
      Data h / 5.20082, 5.20082, 5.20082, 5.0346146046701, 4.91669, 5.09314627212803, 5.50666272693079, 6.12146, 5.95062506384396,
     1      5.74154, 5.55812998041133, 5.31998, 5.16226984536375, 5.12846, 4.95386, 4.43205, 4.18619685740845, 4.41613, 4.56697,
     1      4.58186, 5.74173/
      Data c3 / -0.00125182256184671, -0.00125182256184671, -0.00133789871642764, -0.00145781319421002, -0.00173468008070477,
     1      -0.00246123088762426, -0.00381606405587413, -0.00435108777453373, -0.0047057194777153, -0.00460999914182417,
     1      -0.0043200209052312, -0.00390394298744224, -0.00275425145452962, -0.00158999623017042, -0.000272441818527292,
     1      0.000164300223313353, 0.000502999605010327, 0.000398719751676808, -0.000948565808555795, -0.0024781223824519,
     1      -0.00322395485836136/
      Data b1 / -0.264028732682727, -0.264028732682727, -0.277024473197614, -0.283451417352844, -0.308005433242966,
     1      -0.284467972471899, -0.246398118410714, -0.314858885881003, -0.303807335017873, -0.280546390151647, -0.186148612791909,
     1      -0.0609586208135754, -0.0377024981203794, -0.00172994133464409, 0.00469921518191794, 0.0572336411166358,
     1      0.317560454195245, 0.655123906459053, 0.893461205087385, 0.806835366863026, 0.969592203201732/
      Data b2 / -0.210951794, -0.210951794, -0.210951794, -0.193123478297033, -0.180474087, -0.1608422730442, -0.139406648252206,
     1      -0.173459492, -0.201518836635462, -0.235860699, -0.238616639899499, -0.24219511, -0.285257767788449, -0.319087637,
     1      -0.373226018, -0.377872634, -0.365550259031872, -0.339764851, -0.322055669, -0.341431923, -0.256611596/
      Data b3 / 0, 0, 0, 0, 0, 0, 0, 0, 0.0765190098116516, 0.170170709, 0.246865015286964, 0.346449256, 0.402141459151449,
     1      0.436758946, 0.517683299, 0.356325043, 0.186190244791055, 0.213765783, 0.227350805, 0, 0/
      Data gamma / -0.420010505364829, -0.420010505364829, -0.413347035937856, -0.402883815827122, -0.383986069840125,
     1      -0.361691454954229, -0.32852135916708, -0.33282648363129, -0.345703356657057, -0.372689074563319, -0.391727413988005,
     1      -0.402371664144613, -0.419633565670697, -0.47940685364329, -0.525072668541241, -0.554508489346, -0.715799389410168,
     1      -0.828844831674151, -0.915855863015358, -0.922171987503322, -0.937731571051244/
      Data sofN / -0.0341137359074037, -0.0341137359074037, -0.0392742814171438, -0.0368435752483731, -0.0242782024245241,
     1      -0.0321899120452947, -0.0235686801843898, -0.0462006347636813, -0.00702791404031222, 0.00122126989570676,
     1      0.0312035502188266, 0.0277947946911042, -0.0203490413230642, -0.0236806245379859, -0.0647608232700541,
     1      -0.119553289633279, -0.129207746499647, -0.120818803004753, -0.14582252672715, -0.185290984847367, -0.170081819191526/
      Data sofR / 0.188104514724646, 0.188104514724646, 0.186927382190412, 0.192042043686743, 0.203145831357016, 0.20436206557435,
     1      0.200592230823518, 0.201020548031326, 0.200999016576307, 0.178429474522643, 0.179943044142034, 0.187213636902413,
     1      0.173592950213276, 0.158653833000087, 0.123383881665229, 0.0949571434239225, 0.119199060471125, 0.132643796729745,
     1      0.118405804003951, 0.103287532691844, 0.0373650184045664/
      Data sofS / -0.096283977, -0.096283977, -0.096283977, -0.0966103375235148, -0.096841894, -0.0963783235595202,
     1      -0.0979532926457108, -0.107435167, -0.101695221063626, -0.094670095, -0.0968825621131357, -0.099755355,
     1      -0.0954761145312324, -0.090850567, -0.108074134, -0.098696165, -0.0876704883298241, -0.103411629, -0.096162861,
     1      -0.058467241, -0.005358461/
      Data tau / 0.365351177, 0.365351177, 0.365351177, 0.359897473917039, 0.35602801, 0.372705149349208, 0.396474106671501,
     1      0.390727967, 0.373344954217707, 0.352069866, 0.349217752310341, 0.345514406, 0.347774367936185, 0.363684105,
     1      0.363182142, 0.416353437, 0.435717323627481, 0.44906396, 0.409289105, 0.379829832, 0.406512188/
      Data phi / 0.650148717, 0.650148717, 0.650148717, 0.661984155059029, 0.670381532, 0.670708390598036, 0.675863668060965,
     1      0.695304713, 0.699277470177307, 0.704139732, 0.698434503079109, 0.69102651, 0.69555035047304, 0.684793412, 0.683849352,
     1      0.667650666, 0.665568654301517, 0.685956217, 0.728335296, 0.759659664, 0.723391646/
      Data sigs2s / 0.423581251, 0.423581251, 0.423581251, 0.428228145723766, 0.43152517, 0.443005623288873, 0.463997900223042,
     1      0.479347559, 0.483433172815405, 0.48843356, 0.466739664863982, 0.43857108, 0.431849944902504, 0.400921511, 0.391170063,
     1      0.380133773, 0.387723196000042, 0.462301522, 0.503273721, 0.518887551, 0.477203853/
      Data sig / 0.745772773, 0.745772773, 0.745772773, 0.753543188479092, 0.759056386, 0.767508805717544, 0.783602925138285,
     1      0.797567122, 0.792929650007148, 0.787253843, 0.780874405181787, 0.772590981, 0.777647499576309, 0.775374807,
     1      0.77430871, 0.78683247, 0.795506751916842, 0.819874566, 0.835458462, 0.84932463, 0.829789498/
 
C Find the requested spectral period and corresponding coefficients
      nPer = 21
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         e1T      = e1(1)
         c1T      = c1(1)
         c2T      = c2(1)
         hT       = h(1)
         c3T      = c3(1)
         b1T      = b1(1)
         b2T      = b2(1)
         b3T      = b3(1)
         gammaT   = gamma(1)
         sofNT    = sofN(1)
         sofRT    = sofR(1)
         sofST    = sofS(1)
         sigs2sT  = sigs2s(1)
         phiT     = phi(1)
         tauT     = tau(1)
         sigT     = sig(1)
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
      write (*,*) 'Adjusted Bindi et al. (2013) Horizontal'
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
 1020 call S24_interp (period(count1),period(count2),e1(count1),e1(count2),
     +             specT,e1T,iflag)
      call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +             specT,hT,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),b1(count1),b1(count2),
     +             specT,b1T,iflag)
      call S24_interp (period(count1),period(count2),b2(count1),b2(count2),
     +             specT,b2T,iflag)
      call S24_interp (period(count1),period(count2),b3(count1),b3(count2),
     +             specT,b3T,iflag)
      call S24_interp (period(count1),period(count2),gamma(count1),gamma(count2),
     +             specT,gammaT,iflag)
      call S24_interp (period(count1),period(count2),sofN(count1),sofN(count2),
     +             specT,sofNT,iflag)
      call S24_interp (period(count1),period(count2),sofR(count1),sofR(count2),
     +             specT,sofRT,iflag)
      call S24_interp (period(count1),period(count2),sofS(count1),sofS(count2),
     +             specT,sofST,iflag)
      call S24_interp (period(count1),period(count2),sigs2s(count1),sigs2s(count2),
     +             specT,sigs2sT,iflag)
      call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +             specT,phiT,iflag)
      call S24_interp (period(count1),period(count2),tau(count1),tau(count2),
     +             specT,tauT,iflag)
      call S24_interp (period(count1),period(count2),sig(count1),sig(count2),
     +             specT,sigT,iflag)
   
 1011 period1 = specT                                                                                                              
C     Set Constant Terms
      Mref = 5.5
      Rref = 1.0
      Mh = 6.75
      Vref = 800.0
C     Set the mechanism term. 
      if (ftype .eq. 0 ) then
         termsof = 0
      elseif (ftype .ge. 0.5) then
         termsof = sofRT
      elseif (ftype .le. -0.5) then
         termsof = sofNT
      endif
      R = sqrt (jbdist**2 + hT**2)
      f_D = (c1T+c2T*(M-Mref))*alog(R/Rref) + c3T*(R-Rref)
C     Compute the ground motion for the given spectral period. 
      if (M .le. Mh) then
         f_M = b1T*(M-Mh) + b2T*(M-Mh)**2.0 
      else
         f_M = b3T*(M-Mh)  
      endif
   
      f_S = gammaT*alog(vs/vref)
   
      lnY = e1T + f_D + f_M + f_S + termsof
   
C     Set the sigma value and convert from log10 to Ln units
c      phiT = phiT*alog(10.0)
c      tauT = tauT*alog(10.0)
      sigma = sigT
c      sigs2sT = sigs2sT*alog(10.0)
C     Convert ground motion to units of gals in natural log units.
c      lnY = lnY*alog(10.0)
      lnY = lnY + 6.89
      period2 = period1
      return
      end
 
      
c ------------------------------------------------------------------            
C *** Abrahamson, Silva, and Kamai (NGA-West2 2013) Horizontal ****
C     Earthquake Spectra Paper:
C        Summary of the Abrahamson, Silva, and Kamai NGA-West2 
C            Ground-Motion Relations for Active Crustal REgions
C         N. A. Abrahamson, S. J. Silva, and R. Kamai
C     Notes:
C        Applicable Range (see Abstract):  
C           3 <= M <= 8.5
C           Rrup <= 300 km
C        Regional attenuation included based on Regionflag
C             0 = Global
C             1 = Taiwan
C             2 = China
C             3 = Japan
C         Mainshock and Aftershocks included based on MSASFlag
C             0 = Mainshocks
C             1 = Aftershocks
C         Sigma dependent on estimated or measured Vs30m based on 
C             Vs30_Class
C             0 = Estimated Vs30m
C             1 = Measured Vs30m
c ------------------------------------------------------------------            
      subroutine S04_ASK14_TW_E03 ( mag, dip, fType, fltWidth, rRup, Rjb,  
     1                     vs30, hwflag, lnY, sigma, specT, period2, ztor,
     2                     iflag, vs30_class, z10, Rx, Ry0, regionflag, msasflag,
     1                     phi, tau )
C     Last Updated: 8/1/13
      implicit none
 
      real mag, dip, fType, rRup, rjb, Rx, Ry0, vs30, SA1180,
     1      Z10,  ZTOR, fltWidth, lnSa, sigma, lnY, vs30_rock
      real Fn, Frv, specT, period2, CRjb, phi, tau, z10_rock, SA_rock
      integer hwflag, iflag, vs30_class, regionflag, msasflag
 
c     Vs30 class is to distinguish between the sigma if the Vs30 is measured
c     vs the VS30 being estimated from surface geology.
c         Vs30_class = 0 for estimated
c         Vs30_class = 1 for measured 
 
C     Current version is not programmed for Aftershock cases. 
C       For implementation of Aftershock a new distance metric, CRjb
C       will need to be computed and passed along to this subroutine. 
 
      CRjb = 999.9 
 
C     Set mechanism term and corresponding Frv and Fnm values.     
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
C 
      if ( fType .eq. 1.0 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .eq. 0.5 ) then
        Frv = 1.0
        Fn = 0.0
      elseif ( fType .eq. -1.0 ) then
        Frv = 0.0
        Fn = 1.0
      elseif ( fType .eq. -0.5 ) then
        Frv = 0.0
        Fn = 1.0
      else
        Frv = 0.0
        Fn = 0.0
      endif
 
c     Compute SA1180
      vs30_rock = 1180.
      z10_rock = 0.005
      SA_rock = 0.
      
      call S04_ASK14_TW_E03_model ( mag, dip, fltWidth, ZTOR, Frv, Fn, rRup, rjb, rx, Ry0, 
     1                     vs30_rock, SA_rock, Z10_rock, hwflag, vs30_class,
     2                     specT, lnSa, phi, tau, iflag, regionflag, msasflag, CRjb )
      Sa1180 = exp(lnSa)
 
c     Compute Sa at spectral period for given Vs30
 
      call S04_ASK14_TW_E03_model ( mag, dip, fltWidth, ZTOR, Frv, Fn, rRup, rjb, rx, Ry0, 
     1                     vs30, SA1180, Z10, hwflag, vs30_class,
     2                     specT, lnSa, phi, tau, iflag, regionflag, msasflag, CRjb )
 
c     compute Sa (given the PGA rock value)
      sigma = sqrt( phi**2 + tau**2 )
 
      lnY = lnSa + 6.89
 
      period2 = specT
 
      return
      end
 
c ----------------------------------------------------------------------
      subroutine S04_ASK14_TW_E03_model ( mag, dip, FltWidth, ZTOR, Frv, Fn, rRup, rjb, Rx, Ry0, 
     1                     vs30, Sa1180, Z1, hwflag, vs30_class,
     3                     specT, lnSa, phi, tau, iflag, regionflag, msasflag, CRjb)
 
      implicit none
      
      integer MAXPER     
      parameter (MAXPER=25)
      real Vlin(MAXPER), b(MAXPER), c4(MAXPER), M1(MAXPER), a1(MAXPER)
      real a2(MAXPER), a3(MAXPER), a6(MAXPER), a8(MAXPER), a10(MAXPER)
      real a11(MAXPER), a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER)
      real a17(MAXPER), a43(MAXPER), a44(MAXPER), a45(MAXPER), a46(MAXPER)
      real a25(MAXPER), a28(MAXPER), a29(MAXPER), a31(MAXPER), a36(MAXPER)
      real a37(MAXPER), a38(MAXPER), a39(MAXPER), a40(MAXPER), a41(MAXPER), a42(MAXPER)
      real s1est(MAXPER), s2est(MAXPER), s1msr(MAXPER), s2msr(MAXPER)
      real s3(MAXPER), s4(MAXPER), s5(MAXPER), s6(MAXPER), period(MAXPER)
      real a4(MAXPER), a5, a7
      real VlinT, bT, c4T, M1T, a1T
      real a2T, a3T, a6T, a8T, a10T, a4T
      real a11T, a12T, a13T, a14T, a15T
      real a17T, a43T, a44T, a45T, a46T
      real a25T, a28T, a29T, a31T, a36T
      real a37T, a38T, a39T, a40T, a41T, a42T
      real s1estT, s2estT, s1msrT, s2msrT
      real s3T, s4T, s5T, s6T, c4_mag
      real phiA_est, phiA_msr, period1
      
      real M2
      real lnSa, SA1180, rjb, rRup, Rx, Ry0, dip, mag, vs30
      real HW_taper1, HW_taper2, HW_taper3, HW_taper4, HW_taper5
      real damp_dSA1180, sigAmp, fltWidth
      real f1, f4, f5, f6, f7, f8, f10, f11, fReg, f12, f13
      real Ry1, ZTOR, Frv, Fn, SpecT
      real phiA, phiB, tauA, tauB, phi, tau
      integer vs30_class, hwflag, iflag, nPer, regionflag, msasflag
      real n, c, z1, z1_ref
      real R, V1, Vs30Star, hw_a2, h1, h2, h3, R1, R2, CRjb
      integer count1, count2, i
      real y1, y2, x1, x2, y1z, y2z, x1z, x2z
 
 
      Data Period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4,
     1      5, 7.5, 10/
      Data Vlin / 660, 660, 680, 770, 851.659765220817, 915, 960, 910, 833.557751246246, 740, 674.738820243143, 590, 495, 430, 360,
     1      340, 330, 330, 330, 330, 330, 330, 330, 330, 330/
      Data b / -1.47, -1.47, -1.459, -1.39, -1.2936977941189, -1.219, -1.152, -1.23, -1.39052872238288, -1.587, -1.77190667597776,
     1      -2.012, -2.411, -2.757, -3.278, -3.599, -3.8, -3.5, -2.4, -1, 0, 0, 0, 0, 0/
      Data c4 / 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     1      4.5, 4.5/
      Data M1 / 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75,
     1      6.75, 6.82, 6.92, 7, 7.15, 7.25/
      Data a1 / 0.515752059105666, 0.515752059105666, 0.525965272598101, 0.551232209070966, 0.579804327267391, 0.634993029400247,
     1      0.883852047776952, 1.05080201084861, 1.17945339477789, 1.26593036089864, 1.30919277345635, 1.359730302838,
     1      1.50077522112501, 1.59564384622605, 1.51277593451914, 1.43908031259574, 1.3196502459015, 1.2063286269123,
     1      0.949663235267868, 0.750639640295542, 0.337153706963871, 0.0754357515863762, -0.207939424723418, -1.22376727970931,
     1      -2.08785274897304/
      Data a2 / -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79,
     1      -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.765, -0.634, -0.529/
      Data a3 / 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275,
     1      0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275/
      Data a4(1:25) / 0.0495698, 0.0495698, 0.0390541, 0.0255496, -0.007645, -0.030922, -0.061750, -0.081645, -0.073836, 
     1               -0.062996, -0.051664, -0.044848, 0.0133874, 0.0809235, 0.1084020, 0.1504476, 0.1442501,  
     1               0.2095746, 0.3581460, 0.4368578, 0.3650917, 0.2733890, 0.2247384, 0.0816875, -0.104597/
      Data a6 / 1.23717748598244, 1.23717748598244, 1.21929244416125, 1.18776475236154, 1.1386198050469, 1.06637205570935,
     1      1.02389263942901, 1.07752774917395, 1.1812142280341, 1.27818420722981, 1.33346592714088, 1.41780002102162,
     1      1.55359352494286, 1.66821078559092, 1.82895235653472, 1.92716558892461, 2.23060244069684, 2.47207263637695,
     1      2.65903360278295, 2.76456465596202, 2.81455506509579, 2.86433089087378, 2.78708834917703, 2.27921975558958,
     1      1.41525726709607/
      Data a8 / -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.0181476220075075, -0.022, -0.0254805962536991,
     1      -0.03, -0.038, -0.045, -0.055, -0.065, -0.095, -0.11, -0.124, -0.138, -0.172, -0.197, -0.218, -0.255, -0.285/
      Data a10 / 1.89083492480792, 1.89083492480792, 1.87460135247031, 1.76211350257916, 1.61443183119401, 1.50363187017527,
     1      1.4090869756904, 1.47890076707259, 1.63332017493193, 1.82597122248374, 2.07773336055039, 2.42097852808943,
     1      3.03179342267475, 3.33576148779323, 4.13356936485604, 4.65243364489275, 4.87770413836195, 4.34778741199272,
     1      2.68122643944686, 0.675958400923867, -0.701731284520998, -0.681551288424886, -0.675729625257139, -0.715692117050037,
     1      -0.76096582037326/
      Data a11 / 0.0204898906386533, 0.0204898906386533, 0.0221041041712509, 0.0216820961375558, 0.0348706632203091,
     1      0.0420746347091259, 0.0593187964347741, 0.0817510779375231, 0.114137690473519, 0.128772315620918, 0.137560367749909,
     1      0.139072825872178, 0.137079592197487, 0.151473534772405, 0.0851888969474017, 0.023694753084204, 0.0511085267557186,
     1      0.0590746772575226, 0.0538160326818372, 0.0272210601385645, -0.0290399691465416, -0.084284975099982,
     1      -0.138189866483548, -0.225041479336417, -0.173841111661655/
      Data a12 / -0.0434881906880183, -0.0434881906880183, -0.0406064297893105, -0.0360572370481667, -0.00450128728753472,
     1      0.0175844682739117, 0.0415661398017674, -0.00940639550634408, 0.0377318270533452, 0.0504726922720043,
     1      0.102654347142401, 0.102315899582994, 0.0212337142422029, 0.0226318694825402, 0.00382205344070005, -0.0558206573265485,
     1      -0.102934247138437, -0.117317479408773, -0.144610004682186, -0.194298490614156, -0.209597487336665, -0.219820299126258,
     1      -0.391508524187674, -0.243508338360703, -0.389377983928878/
      Data a13 / 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.58, 0.56, 0.53, 0.5, 0.42, 0.35, 0.2, 0,
     1      0, 0, 0/
      Data a14 / -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.24, -0.19, -0.11, -0.04, 0.07, 0.15,
     1      0.27, 0.35, 0.46, 0.54, 0.61, 0.72, 0.8/
      Data a15 / 0.952886534345006, 0.952886534345006, 0.945485578184931, 0.924040323141569, 0.911895250528834, 0.926379261832576,
     1      0.913774531356106, 0.928600111926654, 0.938643840034865, 0.956069292534134, 0.988813836565905, 0.997918632289772,
     1      0.961361805536286, 0.948242197134087, 0.987444458169033, 0.989386698853271, 0.865419045408074, 0.812275301247927,
     1      0.699815092188431, 0.665629217476703, 0.462283418340137, 0.329492634149113, 0.26821523588858, 0.18175093025376,
     1      0.0425866525919946/
      Data a17 / -0.00797507454784998, -0.00797507454784998, -0.00811713628445259, -0.0088371482560873, -0.00944700779246997,
     1      -0.010056910806905, -0.0108414565925128, -0.0106900501055218, -0.0103707000297666, -0.00916470944805998,
     1      -0.00845544209662438, -0.00738040686016529, -0.00615546315315055, -0.00534167482845196, -0.00386413815130589,
     1      -0.00286182573942047, -0.00312552599477905, -0.00390018562669184, -0.00378727984340934, -0.00485179523035504,
     1      -0.00336928957584865, -0.00237403743233928, -0.00209352593747079, -0.00294860514044089, -0.00344279044006146/
      Data a25 / -0.0015, -0.0015, -0.0015, -0.0016, -0.00182526831785053, -0.002, -0.0027, -0.0033, -0.00338993205735736, -0.0035,
     1      -0.00341298509365752, -0.0033, -0.0029, -0.0027, -0.0023, -0.002, -0.001, -5e-04, -4e-04, -2e-04, 0, 0, 0, 0, 0/
      Data a28 / 0.0025, 0.0025, 0.0024, 0.0023, 0.00252526831785053, 0.0027, 0.0032, 0.0036, 0.00346510191396396, 0.0033,
     1      0.00303895528097257, 0.0027, 0.0024, 0.002, 0.001, 8e-04, 7e-04, 7e-04, 6e-04, 3e-04, 0, 0, 0, 0, 0/
      Data a29 / -0.0034, -0.0034, -0.0033, -0.0034, -0.00334368292053737, -0.0033, -0.0029, -0.0025, -0.0025, -0.0025,
     1      -0.00276104471902743, -0.0031, -0.0036, -0.0039, -0.0048, -0.005, -0.0041, -0.0032, -0.002, -0.0017, -0.002, -0.002,
     1      -0.002, -0.002, -0.002/
      Data a31 / -0.1503, -0.1503, -0.1479, -0.1447, -0.137885633385021, -0.1326, -0.1353, -0.1128, -0.0448563306665158, 0.0383,
     1      0.0553549216431254, 0.0775, 0.0741, 0.2548, 0.2136, 0.1542, 0.0787, 0.0476, -0.0163, -0.1203, -0.2719, -0.2958,
     1      -0.2718, -0.14, -0.0216/
      Data a36 / 0.265, 0.265, 0.255, 0.249, 0.222530972652563, 0.202, 0.126, 0.022, -0.049046325312313, -0.136,
     1      -0.110765677160682, -0.078, 0.037, -0.091, 0.129, 0.31, 0.505, 0.358, 0.131, 0.123, 0.109, 0.135, 0.189, 0.15, 0.092/
      Data a37 / 0.337, 0.337, 0.328, 0.32, 0.302541705366584, 0.289, 0.275, 0.256, 0.213731933042042, 0.162, 0.188974620966168,
     1      0.224, 0.248, 0.203, 0.232, 0.252, 0.208, 0.208, 0.108, 0.068, -0.023, 0.028, 0.031, -0.07, -0.159/
      Data a38 / 0.188, 0.188, 0.184, 0.18, 0.172678779669858, 0.167, 0.173, 0.189, 0.15257751677027, 0.108, 0.111045521721987,
     1      0.115, 0.122, 0.096, 0.123, 0.134, 0.129, 0.152, 0.118, 0.119, 0.093, 0.084, 0.058, 0, -0.05/
      Data a39 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      Data a40 / 0.088, 0.088, 0.088, 0.093, 0.115526831785053, 0.133, 0.186, 0.16, 0.118631253615615, 0.068, 0.0592985093657524,
     1      0.048, 0.055, 0.073, 0.143, 0.16, 0.158, 0.145, 0.131, 0.083, 0.07, 0.101, 0.095, 0.151, 0.124/
      Data a41 / -0.196, -0.196, -0.194, -0.175, -0.127130482456762, -0.09, 0.09, 0.006, -0.0668449664594602, -0.156,
     1      -0.207338794742061, -0.274, -0.248, -0.203, -0.154, -0.159, -0.141, -0.144, -0.126, -0.075, -0.021, 0.072, 0.205,
     1      0.329, 0.301/
      Data a42 / 0.044, 0.044, 0.061, 0.162, 0.324756359647008, 0.451, 0.506, 0.335, 0.146592339836334, -0.084, -0.124897005980964,
     1      -0.178, -0.187, -0.159, -0.023, -0.029, 0.061, 0.062, 0.037, -0.143, -0.028, -0.097, 0.015, 0.299, 0.243/
      Data a43 / -0.149495952993246, -0.149495952993246, -0.151128232484514, -0.168289086295502, -0.152744340369102,
     1      -0.125178594984686, -0.0708221042932599, -0.018687600885894, -0.0488962033737504, -0.102650669997252,
     1      -0.14139394076209, -0.163740294866553, -0.136395538821148, -0.124602153543284, -0.141053741759318, -0.196182748731043,
     1      -0.202510576145646, -0.240547327695437, -0.0826441436634076, -0.0180271258142243, 0.0892524330661483,
     1      0.116738540780062, 0.116530419956562, 0.245400179869857, 0.232382919073963/
      Data a44 / 0.256715289083819, 0.256715289083819, 0.255339396571963, 0.264260006001574, 0.272015226453774, 0.285510349704578,
     1      0.324678146856865, 0.350313882311414, 0.352438749966816, 0.334134637027212, 0.308471519152929, 0.27514124588578,
     1      0.2563436061904, 0.238864930136566, 0.205874034975355, 0.15795474716122, 0.0912683348414232, 0.0841261231937438,
     1      0.138117392538999, 0.207629332605861, 0.253411979944384, 0.23872346346654, 0.238588535068823, 0.202834869900874,
     1      0.134203153980519/
      Data a45 / 0.0457700473517813, 0.0457700473517813, 0.0435318067889647, 0.0501947432901603, 0.0624301316668515,
     1      0.0687884080659366, 0.0647670033588215, 0.0590494940840289, 0.0437403341501233, 0.014042923325033, 0.0156366357333638,
     1      0.0104774389934965, -0.0280510978438352, -0.0303761507159176, 0.0150254526329061, 0.0685318298837249,
     1      0.160653952013054, 0.210407637256458, 0.254117761826276, 0.25349767918346, 0.241497248250486, 0.19997793858101,
     1      0.186738729301699, 0.118521497630021, 0.100176124550952/
      Data a46 / -0.003997315926222, -0.003997315926222, -0.00357648428663343, -0.0109957133438085, -0.0286145009912559,
     1      -0.0514643243277616, -0.0804767575159452, -0.128361501561549, -0.126208108964704, -0.0724493444170797,
     1      -0.0492564870359541, -0.0281250300810575, 0.0719710790597949, 0.112299541489194, 0.189126676107113, 0.217751558860112,
     1      0.272661172941417, 0.282861088845214, 0.240034763108598, 0.269167681227446, 0.232475892491132, 0.223280865347115,
     1      0.219466643291953, 0.17212732187748, 0.143452062783061/
      Data s1est / 0.754, 0.754, 0.76, 0.781, 0.797331953044163, 0.81, 0.81, 0.81, 0.805953057418919, 0.801, 0.795779105619451,
     1      0.789, 0.77, 0.74, 0.699, 0.676, 0.631, 0.609, 0.578, 0.555, 0.548, 0.527, 0.505, 0.457, 0.429/
      Data s2est / 0.52, 0.52, 0.52, 0.52, 0.525631707946263, 0.53, 0.54, 0.55, 0.554496602867868, 0.56, 0.562175372658562, 0.565,
     1      0.57, 0.58, 0.59, 0.6, 0.615, 0.63, 0.64, 0.65, 0.64, 0.63, 0.63, 0.63, 0.63/
      Data s1msr / 0.741, 0.741, 0.747, 0.769, 0.785331953044163, 0.798, 0.798, 0.795, 0.785107473690691, 0.773, 0.764298509365752,
     1      0.753, 0.729, 0.693, 0.644, 0.616, 0.566, 0.541, 0.506, 0.48, 0.472, 0.447, 0.425, 0.378, 0.359/
      Data s2msr / 0.501, 0.501, 0.501, 0.501, 0.50719487874089, 0.512, 0.522, 0.527, 0.523402717705706, 0.519, 0.516824627341438,
     1      0.514, 0.513, 0.519, 0.524, 0.532, 0.548, 0.565, 0.576, 0.587, 0.576, 0.565, 0.568, 0.575, 0.585/
      Data s3 / 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47,
     1      0.47, 0.47, 0.47, 0.47, 0.47, 0.47/
      Data s4 / 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36,
     1      0.36, 0.36, 0.36, 0.36, 0.36, 0.36/
      Data s5 / 0.54, 0.54, 0.54, 0.55, 0.555631707946263, 0.56, 0.57, 0.57, 0.574496602867868, 0.58, 0.584350745317124, 0.59,
     1      0.61, 0.63, 0.66, 0.69, 0.73, 0.77, 0.8, 0.8, 0.8, 0.76, 0.72, 0.67, 0.64/
      Data s6 / 0.63, 0.63, 0.63, 0.63, 0.641263415892527, 0.65, 0.69, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.69, 0.68,
     1      0.66, 0.62, 0.55, 0.52, 0.5, 0.5, 0.5/
 
 
C Find the requested spectral period and corresponding coefficients
      nPer = 25
 
C First check for the PGA, PGV, PGD cases 
      if (specT .eq. 0.0) then
         period1 = period(1)
         a1T = a1(1)
         a2T = a2(1)
         a3T = a3(1)
         a4T = a4(1)
         a6T = a6(1)
         a8T = a8(1)
         M1T = M1(1)
         a10T = a10(1)
         a11T = a11(1)
         a12T = a12(1)
         a13T = a13(1)
         a14T = a14(1)
         a15T = a15(1)
         a17T = a17(1)
         a43T = a43(1)
         a44T = a44(1)
         a45T = a45(1)
         a46T = a46(1)
         a25T = a25(1)
         a28T = a28(1)
         a29T = a29(1)
         a31T = a31(1)
         a36T = a36(1)
         a37T = a37(1)
         a38T = a38(1)
         a39T = a39(1)
         a40T = a40(1)
         a41T = a41(1)
         a42T = a42(1)
         VlinT = Vlin(1)
         bT = b(1)
         c4T = c4(1)
         s1estT = s1est(1)
         s2estT = s2est(1)
         s1msrT = s1msr(1)
         s2msrT = s2msr(1)
         s3T = s3(1)
         s4T = s4(1)
         s5T = s5(1)
         s6T = s6(1)
         goto 1011
      endif
C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Adjusted Abrahamson, Silva, and Kamai (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),M1(count1),M1(count2),
     +                   specT,M1T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a17(count1),a17(count2),
     +                   specT,a17T,iflag)
            call S24_interp (period(count1),period(count2),a43(count1),a43(count2),
     +                   specT,a43T,iflag)
            call S24_interp (period(count1),period(count2),a44(count1),a44(count2),
     +                   specT,a44T,iflag)
            call S24_interp (period(count1),period(count2),a45(count1),a45(count2),
     +                   specT,a45T,iflag)
            call S24_interp (period(count1),period(count2),a46(count1),a46(count2),
     +                   specT,a46T,iflag)
            call S24_interp (period(count1),period(count2),a25(count1),a25(count2),
     +                   specT,a25T,iflag)
            call S24_interp (period(count1),period(count2),a28(count1),a28(count2),
     +                   specT,a28T,iflag)
            call S24_interp (period(count1),period(count2),a29(count1),a29(count2),
     +                   specT,a29T,iflag)
            call S24_interp (period(count1),period(count2),a31(count1),a31(count2),
     +                   specT,a31T,iflag)
            call S24_interp (period(count1),period(count2),a36(count1),a36(count2),
     +                   specT,a36T,iflag)
            call S24_interp (period(count1),period(count2),a37(count1),a37(count2),
     +                   specT,a37T,iflag)
            call S24_interp (period(count1),period(count2),a38(count1),a38(count2),
     +                   specT,a38T,iflag)
            call S24_interp (period(count1),period(count2),a39(count1),a39(count2),
     +                   specT,a39T,iflag)
            call S24_interp (period(count1),period(count2),a40(count1),a40(count2),
     +                   specT,a40T,iflag)
            call S24_interp (period(count1),period(count2),a41(count1),a41(count2),
     +                   specT,a41T,iflag)
            call S24_interp (period(count1),period(count2),a42(count1),a42(count2),
     +                   specT,a42T,iflag)
            call S24_interp (period(count1),period(count2),Vlin(count1),Vlin(count2),
     +                   specT,VlinT,iflag)
            call S24_interp (period(count1),period(count2),b(count1),b(count2),
     +                   specT,bT,iflag)
            call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +                   specT,c4T,iflag)
            call S24_interp (period(count1),period(count2),s1est(count1),s1est(count2),
     +                   specT,s1estT,iflag)
            call S24_interp (period(count1),period(count2),s2est(count1),s2est(count2),
     +                   specT,s2estT,iflag)
            call S24_interp (period(count1),period(count2),s1msr(count1),s1msr(count2),
     +                   specT,s1msrT,iflag)
            call S24_interp (period(count1),period(count2),s2msr(count1),s2msr(count2),
     +                   specT,s2msrT,iflag)
            call S24_interp (period(count1),period(count2),s3(count1),s3(count2),
     +                   specT,s3T,iflag)
            call S24_interp (period(count1),period(count2),s4(count1),s4(count2),
     +                   specT,s4T,iflag)
            call S24_interp (period(count1),period(count2),s5(count1),s5(count2),
     +                   specT,s5T,iflag)
            call S24_interp (period(count1),period(count2),s6(count1),s6(count2),
     +                   specT,s6T,iflag)
 1011 period1 = specT                                                                                                              
C     Constant values
      n = 1.5
      M2 = 5.0
c      a4 = -0.1
      a5 = -0.41
      a7 = 0.0
C     Set C term
      if (period1 .eq. -1.0) then
         c = 2400.0
      else
         c = 2.4
      endif
C     Magnitude dependent taper for C4 (eq. 4.4)
      if (mag .ge. 5.0) then
         c4_mag = c4T
      elseif (mag .ge. 4.0) then
         c4_mag = c4T - (c4T-1.0) * (5.0-mag)
      else
         c4_mag = 1.0
      endif 
     
c     Set distance (eq 4.3)
      R = sqrt(rRup**2 + c4_mag**2)
          
C     Base Model (eq 4.2)
      if ( mag .lt. M2 ) then
        f1 = a1T + a6T*(Mag-M2) + a7*(Mag-M2)**2 + a4T*(M2-M1T) + a8T*(8.5-M2)**2 +
     1                (a2T + a3T*(M2-M1T)) * alog(R) + a17T*Rrup
      elseif ( mag .le. M1T ) then
        f1 = a1T + a4T*(Mag-M1T) + a8T*(8.5-Mag)**2 + (a2T + a3T*(Mag-M1T)) * alog(R) + a17T*Rrup
      else
        f1 = a1T + a5*(Mag-M1T) + a8T*(8.5-Mag)**2 + (a2T + a3T*(Mag-M1T)) * alog(R) + a17T*Rrup
      endif
   
c     style of faulting (eq 4.5 and 4.6) 
      if ( mag .gt. 5. ) then
        f7 = Frv * a11T
        f8 = Fn * a12T
      elseif ( mag .ge. 4. ) then
        f7 = Frv * a11T * (mag-4.)
        f8 = Fn * a12T * (mag-4.)
      else 
        f7 = 0
        f8 = 0
      endif
c     ZTOR (eq 14) 
c     form modified"Extend the upper bound ZTOR to 50km:"
c     2018/10/09  ZTOR modified to 20km:"
      if (ZTOR .lt. 20.) then 
        f6 = a15T * ZTOR/20.0
      else
        f6 = a15T
      endif    
c     Set VS30_star (eq 4.8 and 4.9)
      if ( specT .gt. 3.0 ) then
        V1 = 800.
      elseif ( specT .gt. 0.5 ) then
        V1 = exp( -0.35 * alog(specT/0.5)  + alog(1500.) )
      else
        V1=1500.
      endif      
      if ( vs30 .lt. v1 ) then 
         vs30Star = vs30
      else
      vs30Star = v1
      endif  
c     Compute site amplification (Eq. 4.7)  
      if (vs30 .lt. vLinT) then
        f5 = a10T*alog(vs30Star/vLinT) - bT*alog(c+Sa1180) 
     1              + bT*alog(Sa1180+c*((vs30Star/vLinT)**(n)) )
      else
      f5 = (a10T + bT*n) * alog(vs30Star/vLinT)
      endif
      if (vs30 .eq. 1180.) then
      f5 = (a10T + bT*n) * alog(1180/vLinT)
      endif
   
c     Set Regional z1 reference (eq 4.18)
      if (regionflag .eq. 1) then  
         z1_ref =  exp(-2.629 / 4.0 * alog((vs30**4.0 + 253.299**4.0)/(2491.945**4.0 + 253.299**4.0)))/ 1000.
      elseif (regionflag .eq. 3) then
         z1_ref = exp ( -5.23/2. * alog( (Vs30**2.0 + 412.**2.0)/(1360.**2.0+412.**2.0) ) ) / 1000.
      else
         z1_ref = exp ( -7.67/4. * alog( (Vs30**4.0 + 610.**4.0)/(1360.**4.0+610.**4.0) ) ) / 1000.
      endif 
      
C     Soil Depth Model (eq 4.17)
C     Updated 8/1/13
      if ( vs30 .lt. 150.0 ) then
         y1z = a43T
         y2z = a43T 
         x1z = 50.0
         x2z = 150.0
      elseif ( vs30 .lt. 250.0 ) then
         y1z = a43T
         y2z = a44T 
         x1z = 150.0
         x2z = 250.0
      elseif ( vs30 .lt. 400.0 ) then
         y1z = a44T
         y2z = a45T 
         x1z = 250.0
         x2z = 400.0
      elseif ( vs30 .lt. 700.0 ) then
         y1z = a45T
         y2z = a46T 
         x1z = 400.0
         x2z = 700.0
      else
         y1z = a46T
         y2z = a46T 
         x1z = 700.0
         x2z = 1000.0
      endif
C     Calculation f10 term and set it equal to zero for Vs=1180m/s (i.e., reference condition)
      if (vs30 .eq. 1180.0) then
          f10 = 0.0
      else
         f10 = ( y1z + (Vs30-x1z)*(y2z-y1z)/(x2z-x1z))*alog( (z1 + 0.01) / (z1_ref+0.01) )
      endif
c     Compute HW taper1 (eq 4.11) 
      if ( dip .le. 30. ) then
        HW_taper1 = 60./ 45.
      else
        HW_taper1 = (90.-dip)/45.
      endif
c     Compute HW taper2 (eq. 4.12)
      hw_a2 = 0.2
      if( mag .ge. 6.5 ) then
        HW_taper2 = 1. + hw_a2 * (mag-6.5) 
      elseif ( mag .gt. 5.5 ) then
        HW_taper2 = 1. + HW_a2 * (mag-6.5) - (1.0 - HW_a2)*(mag-6.5)**2
      else
        HW_taper2 = 0.
      endif
c     Compute HW taper 3 (eq. 4.13)
C   April 11, correction by ronnie for HW_taper3 when Rx.gt.R2
      h1 = 0.25
      h2 = 1.5
      h3 = -0.75
      R1 = fltWidth * cos(dip*3.1415926/180.)
      R2 = 3.*R1
      if ( Rx .lt. R1 ) then
        HW_taper3 = h1 + h2*(Rx/R1) + h3*(Rx/R1)**2
      elseif ( Rx . le. R2 ) then
        HW_taper3 = 1. - (Rx-R1)/(R2-R1)
      else
        HW_taper3 = 0.
      endif 
c     Compute HW taper 4 (eq 4.14)
      if ( ZTOR .le. 10. ) then
        HW_taper4 = 1. - (ZTOR**2) / 100.
      else
        HW_taper4 = 0.
      endif
      
c     Compute HW taper 5 (eq. 13)  **** Ry0 version ***
      Ry1 = Rx * tan(20.*3.1415926/180.)
      if ( Ry0 .lt. Ry1 ) then
        HW_taper5 = 1.
      elseif ( Ry0-Ry1 .lt. 5. ) then
        HW_taper5 = 1. - (Ry0-Ry1) / 5.
      else
        HW_taper5 = 0.
      endif
c     Compute HW taper 5 (eq. 4.15b)  **** No Ry0 version ***     
c      if (Rjb .eq. 0. ) then
c        HW_taper5 = 1. 
c      elseif ( Rjb .lt. 30. ) then
c        HW_taper5 = 1 - Rjb/30.
c      else
c        HW_taper5 = 0.
c      endif
c     Hanging wall Model (eq 4.10)
      if ( HWFlag .eq. 1 ) then
        f4 = a13T * HW_taper1 * HW_taper2 * HW_taper3 * HW_taper4 * HW_taper5
      else
        f4 = 0.
      endif
C     Add aftershock factor (eq 4.21)
      if (msasflag .eq. 1) then
         if (CRjb .ge. 15.0) then
             f11 = 0.0
         elseif (CRjb .le. 5.0) then
             f11 = a14T
         else
             f11 = a14T * ( 1.0 - (CRjb - 5.0) /10.0) 
         endif
      elseif (msasflag .eq. 0) then
         f11 = 0.0
      endif 
C     Now apply the regional attenuation differences (eq 4.22)
C     Global No Change
      if (regionflag .eq. 0) then
         freg = 0.0
C     Taiwan
      elseif (regionflag .eq. 1) then
         f12 = a31T * alog(vs30Star/VlinT) 
         freg = f12 + a25T*Rrup 
      elseif (regionflag .eq. 1 .and. vs30 .eq. 1180.0) then
         f12 = a31T * alog(1180/VlinT) 
         freg = f12 + a25T*Rrup 
C     China
      elseif (regionflag .eq. 2) then
         freg = a28T*Rrup 
C     Japan
      elseif (regionflag .eq. 3) then
         freg = a42T + a29T*Rrup 
      endif
C     Set the Sigma Values
      if (regionflag .ne. 3)  then
c     Compute within-event term, phiA, at the surface for linear site response (eq 7.1)
        if (mag .lt. 4.0) then
           phiA_est = s1estT
        elseif (mag .le. 6.0) then
           phiA_est = s1estT + ((s2estT-s1estT)/2.0)*(mag-4.0)
        else
           phiA_est = s2estT
        endif
c     Compute within-event term, phiA, for known Vs30
        if (mag .lt. 4.0) then
           phiA_msr = s1msrT
        elseif (mag .le. 6.0) then
           phiA_msr = s1msrT + ((s2msrT-s1msrT)/2.0)*(mag-4.0)
        else
           phiA_msr = s2msrT
        endif
C     choose phiA by Vs30 class
        if (vs30_class .eq. 0 ) then
      phiA = phiA_est
        elseif (vs30_class .eq. 1) then
      phiA = phiA_msr
        else
      stop 99
        endif
C     Set Sigma values for Japan Region
      else
C calculate phi_A for Japan (eq. 7.3)
        if (Rrup .lt. 30) then
           phiA = s5T        
        elseif (Rrup .le. 80) then
           phiA = s5T + (s6T-s5T)/50*(Rrup-30)
        else
           phiA = s6T
        endif
      endif
   
c     Compute between-event term, tau (eq. 7.2)
      if (mag .lt. 5.0) then
         tauA = s3T
      elseif (mag .le. 7.0) then
         tauA = s3T + ((s4T-s3T)/2.0)*(mag-5.0)
      else
         tauA = s4T
      endif
      tauB = tauA
c     Compute phiB, within-event term with site amp variability removed (eq. 7.7)
c     with fix to model for small mag at long periods - limit sigAmp to be less than phiA
      sigAmp = 0.4
      if (phiA .le. sigAmp) then
        sigAmp = phiA*0.99
      endif
      phiB = sqrt( phiA**2 - sigAmp**2)
      
c     Compute partial derivative of alog(soil amp) w.r.t. alog(Sa1180) (eq. 7.10)
      if ( vs30 .ge. vLinT) then
        dAmp_dSa1180 = 0.
      else
        dAmp_dSa1180 = bT*Sa1180 * ( -1. / (Sa1180+c) 
     1              + 1./ (Sa1180 + c*(vs30/vLinT)**(n)) )
      endif
C     Compute phi, with non-linear effects (eq. 7.8)
      phi = sqrt( phiB**2 * (1. + dAmp_dSa1180)**2 + sigAmp**2 )
C     Compute tau, with non-linear effects (eq. 7.9)
      tau = tauB * (1. + dAmp_dSa1180)
      

c     Compute median ground motion (eq. 1)
      lnSa = f1 + f4 + f5 + f6 + f7 + f8 + f10 + f11 + freg 
   
c       write(*,*) "f1   = " , f1  
c       write(*,*) "f4   = " , f4  
c       write(*,*) "f5   = " , f5  
c       write(*,*) "f6   = " , f6  
c       write(*,*) "f7   = " , f7  
c       write(*,*) "f8   = " , f8  
c       write(*,*) "f10  = " , f10 
c       write(*,*) "f11  = " , f11 
c       write(*,*) "freg = " , freg
c       write(*,*) "lnSa = ", lnSa
c       write(*,*) "Sa = ", exp(lnSa)

      return
      end
 
 

 
c ---------------------------------------------------------------------------            
C     *** Boore, Stewart, Seyhan and Atkinson NGA West 2 (NGA West2-2013) ***
C         Earthquake Spectra Report:
C            NGA-West2 Equations for Predicting PGA, PGV, and 5%-Damped
C                PSA for Shallow Crustal Earthquakes.
C             D. M. Boore, J. P. Stewart, E. Seyhan, and G. M. Atkinson
C     Notes:
C        Applicable Range:
C            3.0 < M < 8.5 (Strike-Slip)
C            3.0 < M < 7.0 (Normal)
C            Distance < 300km
C            150 < Vs < 1500 m/s
C            0.0 < Z1 < 3.0 km
C            Region Flag:
C               0 = Global 
C               1 = China-Turkey
C               2 = Italy-Japan
c ---------------------------------------------------------------------------            
      subroutine S04_BSSA14_TW_E02 ( mag, Rbjf, specT, 
     1        period2, lnY, sigma, iflag, vs, ftype, pga4nl, z10, regionflag, basinflag,
     1        phi, tau ) 
C     Last Updated: 9/16/13
      parameter (MAXPER=24)
      REAL Period(MAXPER), c1(MAXPER), c2(MAXPER), c3(MAXPER)
      real h(MAXPER), DC3ChinaTrk(MAXPER), DC3ItalyJapan(MAXPER), e0(MAXPER)
      real e1(MAXPER), e2(MAXPER), e3(MAXPER), e4(MAXPER)
      real e5(MAXPER), e6(MAXPER), mh(MAXPER), c(MAXPER), Vc(MAXPER)
      real phi2(MAXPER), phi3(MAXPER), f4(MAXPER)
      real l1(MAXPER), l2(MAXPER), t1(MAXPER), t2(MAXPER)
      real f5(MAXPER), rjbbar(MAXPER), Dfr(MAXPER), Dfv(MAXPER)
      real R1(MAXPER), R2(MAXPER), DC3Global(MAXPER), f6(MAXPER), f7(MAXPER)
      real Mref, Rref, Vref, f1, f3, specT
      REAL MAG, RBJF, VS, z10
      real ftype, Rp, R
      INTEGER iFlag, count1, count2, regionflag, basinflag
      real lnY, mechS, mechN, mechR, pga4nl
      real f2, flin, fBasin, phi, tau
      real c1T, c2T, c3T, hT, e0T, e1T, e2T, e3T, e4T
      real e5T, e6T, mhT, cT, VcT, phi2T, phi3T, f4T, l1T, l2T, t1T, t2T
      real deltaz1, f5T, rjbbarT, DfrT, DfvT, R1T, R2T, DC3GlobalT
      real DC3ChinaTrkT, DC3ItalyJapanT, f6T, f7T
 
 
      Data Period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4,
     1      5, 7.5/
      Data e0 / 0.4473, 0.4534, 0.48598, 0.56916, 0.673459231164795, 0.75436, 0.96447, 1.1268, 1.20895293439595, 1.3095,
     1      1.3164611925074, 1.3255, 1.2766, 1.2217, 1.1046, 0.96991, 0.66903, 0.3932, -0.14954, -0.58669, -1.1898, -1.6388,
     1      -1.966, -2.5865/
      Data e1 / -0.0595130663398657, -0.0434648257872447, 0.00723856448506604, 0.0927610149820673, 0.16457754638023,
     1      0.227041716245343, 0.44200643374538, 0.628456255453818, 0.726456018249075, 0.796155204184605, 0.832996288802006,
     1      0.839710890536238, 0.870329523227315, 0.86421353604234, 0.859428598430318, 0.799146696191413, 0.501318810676072,
     1      0.271210383253748, -0.216071094556293, -0.527765325832417, -1.23517562268132, -1.80451795360888, -2.14029683889059,
     1      -2.57205309494597/
      Data e2 / 0.0128644012118545, -0.039657974514646, -0.0432972654368019, -0.0528753047838778, -0.0371997189165065,
     1      -0.0694243260609077, -0.0665458750869888, -0.078126004267113, -0.0381648642482661, -0.015034889910705,
     1      0.0277769958059142, 0.0232791124271304, -0.0132120586353397, -0.0136120469735889, -0.0353660883149616,
     1      -0.089958791754687, -0.108469227753284, -0.0962910050989717, -0.109993355573954, -0.143706248763231,
     1      -0.134469002079297, -0.15369594108218, -0.355158141044777, -0.279997608428904/
      Data e3 / 0.203074215311353, 0.201176182119364, 0.20039108802914, 0.199653683877596, 0.210766879607819, 0.206791720434049,
     1      0.202176201184481, 0.209133744968306, 0.20864732810785, 0.189898726349265, 0.193454117511752, 0.200065319053198,
     1      0.190509627343089, 0.177634411439118, 0.151740289401428, 0.120064967833506, 0.135964753367863, 0.151592453559548,
     1      0.138426843070023, 0.120156010265541, 0.0481440566062522, -0.0190941855692312, -0.112287130179668, -0.215408568323982/
      Data e4 / 0.712266102827589, 0.684842932365268, 0.683879021139178, 0.672977363930926, 0.620065668793878, 0.614099394851088,
     1      0.603475935406328, 0.639500382498947, 0.674037854179504, 0.674489114985437, 0.62843666806635, 0.627594892015755,
     1      0.600149445878001, 0.642710126756053, 0.65270724948819, 0.689679935416504, 1.08672406855482, 1.39585543347496,
     1      1.72422222138858, 1.9015545343119, 2.12964285080081, 2.21837879531008, 2.2105384076795, 1.85456007759323/
      Data e5 / 0.05053, 0.04932, 0.05339, 0.06144, 0.0647739711041879, 0.06736, 0.07355, 0.05523, 0.0114780540956452, -0.04207,
     1      -0.072042284489666, -0.11096, -0.16213, -0.1959, -0.22608, -0.23522, -0.21591, -0.18983, -0.1467, -0.11237, -0.04332,
     1      -0.01464, -0.01486, -0.08161/
      Data e6 / -0.1662, -0.1659, -0.16561, -0.1669, -0.174739337461198, -0.18082, -0.19665, -0.19838, -0.19116744899994, -0.18234,
     1      -0.171976524654611, -0.15852, -0.12784, -0.09286, -0.02319, 0.02912, 0.10829, 0.17895, 0.33896, 0.44788, 0.62694,
     1      0.76303, 0.87314, 1.0121/
      Data Mh / 5.62683852817361, 5.7278287417983, 5.72634382479267, 5.72419805461193, 5.78950009876979, 5.6454778542563,
     1      5.60564354629504, 5.6095882429428, 5.62540554303087, 5.71790190281766, 5.84351085002455, 5.87951468808749,
     1      6.01455506496882, 6.0580925842216, 6.22967481182813, 6.31771461927764, 6.22329198777463, 6.19803864164175,
     1      6.1770931403571, 6.19479909698098, 6.20544074705526, 6.18063946709863, 6.21078404445105, 6.59049042088875/
      Data c1 / -1.134, -1.134, -1.1394, -1.1421, -1.12734492518079, -1.1159, -1.0831, -1.0652, -1.05980407655856, -1.0532,
     1      -1.05646305898784, -1.0607, -1.0773, -1.0948, -1.1243, -1.1459, -1.1777, -1.193, -1.2063, -1.2159, -1.2179, -1.2162,
     1      -1.2189, -1.2543/
      Data c2 / 0.1917, 0.1916, 0.18962, 0.18842, 0.187670982843147, 0.18709, 0.18225, 0.17203, 0.163927121632102, 0.15401,
     1      0.150042120270783, 0.14489, 0.13925, 0.13388, 0.12512, 0.12015, 0.11054, 0.10248, 0.09645, 0.09636, 0.09764, 0.10218,
     1      0.10353, 0.12507/
      Data c3 / -0.00789237345533951, -0.00783583541961706, -0.00776556896476328, -0.00829539839832367, -0.00920134243728047,
     1      -0.00996907638994333, -0.0115155455012644, -0.0118950222605932, -0.0115088807317956, -0.0105015393785244,
     1      -0.00973285233813413, -0.00868358737872155, -0.00699522992186121, -0.00579250321997269, -0.00356117007757743,
     1      -0.00182136049301181, -0.00032734268017595, 1.97247798069236e-05, 0.000714578116996264, 0.000165203729102031,
     1      0.00127217001357917, 0.00138046929786292, 0.00144089026812282, 0.0011090363057596/
      Data h / 4.5, 4.5, 4.5, 4.49, 4.32668046955837, 4.2, 4.04, 4.13, 4.24691167456457, 4.39, 4.48571639697672, 4.61, 4.78, 4.93,
     1      5.16, 5.34, 5.6, 5.74, 6.18, 6.54, 6.93, 7.32, 7.78, 9.48/
      Data Dc3Global / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      Data Dc3ChinaTrk / 0.002858, 0.002816, 0.00278, 0.002765, 0.00287312879256825, 0.002957, 0.002957, 0.002879,
     1      0.00283718159332883, 0.002786, 0.00271029703148205, 0.002612, 0.002444, 0.002196, 0.002107, 0.002348, 0.00269,
     1      0.002921, 0.003039, 0.002923, 0.002616, 0.002605, 0.002604, 0.0026/
      Data Dc3ItalyJapan / -0.00255, -0.002437, -0.00234, -0.002168, -0.00206831876935114, -0.001991, -0.002159, -0.002439,
     1      -0.00255905929657207, -0.002706, -0.00282085967637207, -0.00297, -0.00314, -0.003297, -0.003212, -0.002907, -0.002527,
     1      -0.002089, -0.001518, -0.00117, -0.001188, -0.001083, -0.000571, 0.000385/
      Data c / -0.491005720442681, -0.496858012636404, -0.492693576830718, -0.487762396793923, -0.474108229985153,
     1      -0.454433171424995, -0.422912974621248, -0.432287456269094, -0.447361030268909, -0.477038393985033, -0.495648486633132,
     1      -0.505872247475676, -0.523355295292678, -0.580611412765796, -0.622063967379841, -0.643464365148664, -0.799483990066468,
     1      -0.902512999802672, -0.969173634547672, -0.959488485054134, -0.976840024588964, -0.994362745349477, -0.97477634851446,
     1      -0.890924799294348/
      Data Vc / 1500, 1500.2, 1500.36, 1502.95, 1502.08834868422, 1501.42, 1494, 1479.12, 1462.81082139824, 1442.85,
     1      1420.99185552677, 1392.61, 1356.21, 1308.47, 1252.66, 1203.91, 1147.59, 1109.95, 1072.39, 1009.49, 922.43, 844.48,
     1      793.13, 771.01/
      Data f4 / -0.15, -0.14833, -0.1471, -0.15485, -0.1782103245611, -0.19633, -0.22866, -0.24916, -0.252743792485691, -0.25713,
     1      -0.252539963690434, -0.24658, -0.23574, -0.21912, -0.19582, -0.17041, -0.13866, -0.10521, -0.06794, -0.03614, -0.01358,
     1      -0.00321, -0.00025, -5e-05/
      Data f5 / -0.00701, -0.00701, -0.00728, -0.00735, -0.00685440970072883, -0.00647, -0.00573, -0.0056, -0.0057124150716967,
     1      -0.00585, -0.00597617161419659, -0.00614, -0.00644, -0.0067, -0.00713, -0.00744, -0.00812, -0.00844, -0.00771,
     1      -0.00479, -0.00183, -0.00152, -0.00144, -0.00137/
      Data f6 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.092, 0.367, 0.638, 0.871, 1.135, 1.271, 1.329, 1.329/
      Data f7 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.059, 0.208, 0.309, 0.382, 0.516, 0.629, 0.738, 0.809/
      Data R1 / 110, 111.667, 113.105, 112.133, 104.132595691538, 97.927, 85.989, 79.587, 80.3689592387222, 81.326,
     1      85.4944490883363, 90.907, 97.039, 103.152, 106.018, 105.536, 108.388, 116.388, 125.38, 130.369, 130.365, 129.489,
     1      130.224, 130.716/
      Data R2 / 270, 270, 270, 269.998, 269.999126341589, 270, 270.035, 270.092, 270.123925880362, 270.163, 270.092082851331, 270,
     1      269.449, 268.593, 266.543, 265, 266.511, 270, 262.413, 240.138, 195, 199.446, 230, 250.395/
      Data Dfr / 0.1, 0.096, 0.092, 0.081, 0.0708629256967262, 0.063, 0.064, 0.087, 0.101838789463964, 0.12, 0.126961192507398,
     1      0.136, 0.141, 0.138, 0.122, 0.109, 0.1, 0.098, 0.104, 0.105, 0.088, 0.07, 0.061, 0.058/
      Data Dfv / 0.07, 0.07, 0.03, 0.029, 0.0295631707946263, 0.03, 0.022, 0.014, 0.0144496602867868, 0.015, 0.0280522359513715,
     1      0.045, 0.055, 0.05, 0.049, 0.06, 0.07, 0.02, 0.01, 0.008, 0, 0, 0, 0/
      Data l1 / 0.6951, 0.698, 0.7018, 0.7212, 0.73916514834858, 0.7531, 0.7447, 0.7279, 0.72448258182042, 0.7203,
     1      0.716384329214589, 0.7113, 0.6984, 0.6754, 0.6428, 0.6147, 0.5815, 0.5527, 0.5317, 0.5263, 0.5335, 0.536, 0.5285,
     1      0.5117/
      Data l2 / 0.4951, 0.4992, 0.5023, 0.5136, 0.524131293859512, 0.5323, 0.5423, 0.5407, 0.538946324881531, 0.5368,
     1      0.537670149063425, 0.5388, 0.5471, 0.5614, 0.5804, 0.599, 0.6218, 0.625, 0.6192, 0.6182, 0.619, 0.6156, 0.6223, 0.6344/
      Data t1 / 0.3982, 0.4019, 0.4087, 0.4449, 0.477451271929402, 0.5027, 0.4744, 0.4153, 0.387780790448648, 0.3541,
     1      0.34953171741702, 0.3436, 0.35, 0.3634, 0.381, 0.4101, 0.4572, 0.4983, 0.5248, 0.5325, 0.5369, 0.5427, 0.532, 0.511/
      Data t2 / 0.348, 0.3446, 0.3464, 0.364, 0.39908554050522, 0.4263, 0.4658, 0.4583, 0.42664391581021, 0.3879,
     1      0.353355082182037, 0.3085, 0.2664, 0.229, 0.2097, 0.2235, 0.2664, 0.2984, 0.3151, 0.3291, 0.3438, 0.3492, 0.3354,
     1      0.2699/
      
C     Set constant parameters
      mref = 4.5
      rref = 1.0
      vref = 760.0
      f1 = 0.0
      f3 = 0.1
      V1 = 225.0
      V2 = 300.0
C First check for the PGA case (i.e., specT=0.0) 
      nPer = 24
      if (specT .eq. 0.0) then
         period1 = period(1)
         e0T = e0(1)
         e1T = e1(1)
         e2T = e2(1)
         e3T = e3(1)
         e4T = e4(1)
         e5T = e5(1)
         e6T = e6(1)
         mhT = mh(1)
         c1T = c1(1)
         c2T = c2(1)
         c3T = c3(1)        
         hT = h(1)
         cT = c(1)
         VcT = vc(1)
         phi2T = phi2(1)
         phi3T = phi3(1)
         Dc3GlobalT = DC3global(1)
         Dc3ChinaTrkT = DC3ChinaTrk(1)
         Dc3ItalyJapanT = DC3ItalyJapan(1)
         f4T = f4(1)
         f5T = f5(1)
         f6T = f6(1)
         f7T = f7(1)
         R1T = R1(1)
         R2T = R2(1)
         l1T = l1(1)
         l2T = l2(1)
         t1T = t1(1)
         t2T = t2(1)
         rjbbarT = rjbbar(1)
         DfrT = Dfr(1)
         DfvT = Dfv(1)
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
      write (*,*) 'Adjusted BSSA (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),e0(count1),e0(count2),
     +                   specT,e0T,iflag)
            call S24_interp (period(count1),period(count2),e1(count1),e1(count2),
     +                   specT,e1T,iflag)
            call S24_interp (period(count1),period(count2),e2(count1),e2(count2),
     +                   specT,e2T,iflag)
            call S24_interp (period(count1),period(count2),e3(count1),e3(count2),
     +                   specT,e3T,iflag)
            call S24_interp (period(count1),period(count2),e4(count1),e4(count2),
     +                   specT,e4T,iflag)
            call S24_interp (period(count1),period(count2),e5(count1),e5(count2),
     +                   specT,e5T,iflag)
            call S24_interp (period(count1),period(count2),e6(count1),e6(count2),
     +                   specT,e6T,iflag)
            call S24_interp (period(count1),period(count2),mh(count1),mh(count2),
     +                   specT,mhT,iflag)
            call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +                   specT,c2T,iflag)
            call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +                   specT,c3T,iflag)
            call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +                   specT,hT,iflag)
            call S24_interp (period(count1),period(count2),c(count1),c(count2),
     +                   specT,cT,iflag)
            call S24_interp (period(count1),period(count2),phi2(count1),phi2(count2),
     +                   specT,phi2T,iflag)
            call S24_interp (period(count1),period(count2),phi3(count1),phi3(count2),
     +                   specT,phi3T,iflag)
            call S24_interp (period(count1),period(count2),DC3Global(count1),DC3Global(count2),
     +                   specT,DC3GlobalT,iflag)
            call S24_interp (period(count1),period(count2),DC3ChinaTrk(count1),DC3ChinaTrk(count2),
     +                   specT,DC3ChinaTrkT,iflag)
            call S24_interp (period(count1),period(count2),DC3ItalyJapan(count1),DC3ItalyJapan(count2),
     +                   specT,DC3ItalyJapanT,iflag)
            call S24_interp (period(count1),period(count2),Vc(count1),Vc(count2),
     +                   specT,VcT,iflag)
            call S24_interp (period(count1),period(count2),f4(count1),f4(count2),
     +                   specT,f4T,iflag)
            call S24_interp (period(count1),period(count2),f5(count1),f5(count2),
     +                   specT,f5T,iflag)
            call S24_interp (period(count1),period(count2),f6(count1),f6(count2),
     +                   specT,f6T,iflag)
            call S24_interp (period(count1),period(count2),f7(count1),f7(count2),
     +                   specT,f7T,iflag)
            call S24_interp (period(count1),period(count2),R1(count1),R1(count2),
     +                   specT,R1T,iflag)
            call S24_interp (period(count1),period(count2),R2(count1),R2(count2),
     +                   specT,R2T,iflag)
            call S24_interp (period(count1),period(count2),l1(count1),l1(count2),
     +                   specT,l1T,iflag)
            call S24_interp (period(count1),period(count2),l2(count1),l2(count2),
     +                   specT,l2T,iflag)
            call S24_interp (period(count1),period(count2),t1(count1),t1(count2),
     +                   specT,t1T,iflag)
            call S24_interp (period(count1),period(count2),t2(count1),t2(count2),
     +                   specT,t2T,iflag)
            call S24_interp (period(count1),period(count2),rjbbar(count1),rjbbar(count2),
     +                   specT,rjbbarT,iflag)
            call S24_interp (period(count1),period(count2),Dfr(count1),Dfr(count2),
     +                   specT,DfrT,iflag)
            call S24_interp (period(count1),period(count2),Dfv(count1),Dfv(count2),
     +                   specT,DfvT,iflag)
 1011 period1 = specT                                                                                                              
C.....Set the mechanism terms based on ftype............
C     Set mechanism term and corresponding Frv and Fnm values.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                    -120 < Rake <  -60
C     -0.5      Normal/Oblique            -150 < Rake < -120
C                                          -60 < Rake <  -30
C       0       Strike-Slip               -180 < Rake < -150
C                                          -30 < Rake <   30
C                                          150 < Rake <  180
C      0.5      Reverse/Oblique             30 < Rake <   60
C                                          120 < Rake <  150
C       1       Reverse                     60 < Rake <  120 
C     Note: Unknown Mechanism is not currently coded.  
      if (ftype .eq. -1.0) then
         mechS = 0.0
         mechN = 1.0
         mechR = 0.0
      elseif (ftype .eq. -0.5) then 
         mechS = 0.0
         mechN = 1.0
         mechR = 0.0
      elseif (ftype .eq. 0.0) then 
         mechS = 1.0
         mechN = 0.0
         mechR = 0.0
      elseif (ftype .eq. 0.5) then
         mechS = 0.0
         mechN = 0.0
         mechR = 1.0
      elseif (ftype .eq. 1.0) then
         mechS = 0.0
         mechN = 0.0
         mechR = 1.0
      endif 
C.....First compute the Reference Rock PGA value...........
C.....This will include the regional dependence for PGA....
C.....MAGNITUDE DEPENDENCE.................................
      if (mag .le. mh(1)) then
         term1 = e1(1) + e2(1)*mechN + e3(1)*mechR +
     1           e4(1)*(mag-mh(1)) + e5(1)*(mag-mh(1))**2.0
      else
         term1 = e1(1) + e2(1)*mechN + e3(1)*mechR +
     1           e6(1)*(mag-mh(1))
      endif
C.....Distance dependence......
      Rp = SQRT( Rbjf*Rbjf+h(1)*h(1) )
C.....Apply Regional term.....
         TERM2 = ( c1(1) + c2(1)*(mag-mref) ) * alog(Rp/rref) + 
     1           c3(1)  * (Rp-rref)
     
      pga4nl = exp(term1+term2)
C.....Now compute the requested ground motion value........
C.....MAGNITUDE DEPENDENCE.................................
      if (mag .le. mhT) then
         term1 = e1T + e2T*mechN + e3T*mechR +
     1           e4T*(mag-mhT) + e5T*(mag-mhT)**2.0
      else
         term1 = e1T + e2T*mechN + e3T*mechR +
     1           e6T*(mag-mhT)
      endif
C.....Distance dependence......
      R = SQRT( Rbjf*Rbjf+hT*hT )
C     Now apply the regional attenuation differnece.
C     Global Case
         TERM2 = ( c1T + c2T*(mag-mref) ) * alog(R/rref) + 
     1        c3T * (R-rref) 
C.....Site Response Term.........
C.....Now compute the site term........
C.....First the linear term......
      if (vs .le. VcT ) then
         flin = cT*alog(Vs/Vref)
      else
         flin = cT*alog(VcT/Vref)
      endif
C.....Next the non-linear term......
      f2 = f4T*(exp(f5T*(min(vs,760.0)-360.0))-exp(f5T*(760.0-360.0))) 
      
C.....Now compute the basin effect term......
C Deviation from ln(Vs30) scaling: bedrock depth (Z1) effect for California.
c      if (basinflag .eq. 1) then 
c     Compute the DeltaZ1 term. Apply the California model for all regions except for Japan. 
          deltaz1 = z10 -
     1           exp(-2.63/4.0 * alog((vs**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))/1000.0        
          fbasin = min(f7T, f6T * deltaZ1)

c      else
c         fbasin = 0.0
c      endif
      TERM3 = flin + f1 + f2*alog((pga4nl+f3)/f3) + fBasin
      
      lnY = term1 + term2 + term3 
      period2 = period1
c     Now compute the sigma value which is a function of magnitude and Vs
C     Tau (Eq. 4.11)
      if (mag .le. 4.5) then
         tau = t1T
      elseif (mag .gt. 4.5 .and. mag .lt. 5.5) then
         tau = t1T + (t2T - t1T)*(mag - 4.5)
      else
         tau = t2T
      endif
      
C     Phi - Magnitude (Eq. 4.12) 
      if (mag .le. 4.5) then
         phi = l1T
      elseif (mag .gt. 4.5 .and. mag .lt. 5.5) then
         phi = l1T + (l2T - l1T)*(mag - 4.5)
      else
         phi = l2T
      endif
C     Phi - Distance (Eq. 4.13)
      if (rbjf .le. R1T) then
          phi = phi
      elseif (rbjf .gt. R1T .and. rbjf .le. R2T) then
          phi = phi + DfrT*( (alog(rbjf/R1T))/(alog(R2T/R1T)) )
      else
          phi = phi + DfrT 
      endif
C     Phi - Vs30 (Eq. 4.14)
      if (vs .ge. V2) then
         phi = phi
      elseif (vs .ge. V1 .and. vs .lt. V2) then
         phi = phi - DfvT*( alog(V2/vs) / alog(V2/V1) )
      else
         phi = phi - DfvT
      endif
      sigma = sqrt (tau**2 + phi**2)
C     Convert ground motion to units of gals.
      lnY = lnY + 6.89
      return
      END
c ---------------------------------------------------------------------            
C ** Chiou and Youngs (NGA West2-2013 Model) Horizontal **
C     Earthquake Spectra Paper:
C        Update of the Chiou and Youngs NGA Model for the
C            Average Horizontal Component of Peak
C            Ground Motion and Response Spectra
C         B. S.J. Chiou and R.R. Youngs
C     Notes:
C        Applicable Range (see Abstract):  
C           3.5 <= M <= 8.5 Strike-slip
C           3.5 <= M <= 8.0 Reverse and Normal
C           Rrup <= 300 km
C           Ztor <= 20 km
C           180m/s <= Vs30m <= 1500m/s
C           Mainshock events only
C        Regional attenuation included based on Regionflag
C             0 = Global
C             1 = Japan/Italy
C             2 = Wenchuan (note only applicable for M7.9 event)
C         Sigma dependent on estimated or measured Vs30m based on 
C             Vs30_Class
C             0 = Estimated Vs30m
C             1 = Measured Vs30m
c ---------------------------------------------------------------------            
      Subroutine S04_CY14_TW_E04 ( m, Rrup, Rbjf, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, Delta, DTor, Ftype, depthvs10, vs30_class,
     3                     hwflag, Rx, regionflag, phi, tau )
C     Last Updated: 8/1/13
      parameter (MAXPER=25)
      REAL Period(MAXPER), C1(MAXPER), C1a(MAXPER), C1b(MAXPER)
      REAL cn(MAXPER), cm(MAXPER), c5(MAXPER), c6(MAXPER), c8(MAXPER)
      REAL c7(MAXPER), c9(MAXPER), gamma1(MAXPER), gamma2(MAXPER)
      REAL phi1(MAXPER), phi2(MAXPER), phi3(MAXPER), phi4(MAXPER)
      REAL phi5(MAXPER), phi6(MAXPER), phi1jp(MAXPER), phi5jp(MAXPER), phi6jp(MAXPER)
      REAL tau1(MAXPER), tau2(MAXPER), sigma1(MAXPER), sigma2(MAXPER)
      REAL sigma3(MAXPER), c9a(MAXPER)
      REAL c3(MAXPER), gm(MAXPER), c9b(MAXPER)
      Real CHM(MAXPER), C1c(MAXPER), C1d(MAXPER), C7b(MAXPER), C8b(MAXPER)
      real C11(MAXPER), C11B(MAXPER)
      real gscaleJapIt(MAXPER), gscaleWen(MAXPER), sigma2Jap(MAXPER)
 
      REAL c1T, c1aT, c1bT, cnT, cmT, c5T, c6T, c7T, c9T, c9aT, c3T, c8T
      REAL gamma1T, gamma2T, phi1T, phi2T, phi3T, phi4T, sigma3T
      REAL phi5T, phi6T, tau1T, tau2T, sigma1T, sigma2T
      real C1cT, C1dT, C7bT, C11T, C11bT, CHMT
      REAL c2, c4, c4a, cRB, pi, d2r, gamma, gmT, c9bT
      real cc, cosdelta, r1, r2, r3, r4, hw, psa_ref, psa, a, b, c
      integer iflag, count1, count2, hwflag, vs30_class, regionflag
      REAL M, RRUP, RBJF, DTOR, Delta, specT, sigma, Ftype, Rx
      REAL period2, lnY, F_RV, F_NM, depthvs10, tau, rkdepth
      real c8a, c8bT, deltaZ1, fd
      real NL0, sigma_NL0, F_Measured, F_Inferred, mz_TOR, deltaZ_TOR, coshM
      real gscaleJapItT, gscaleWenT, sigma2JapT, phi
      real phi1jpT, phi5jpT, phi6jpT
 
 
      Data period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4,
     1      5, 7.5, 10/
      Data c1 / -1.29685236824339, -1.29685236824339, -1.25271582686019, -1.11773250980953, -0.956465870903099, -0.799324896281897,
     1      -0.473133817081457, -0.353655570079099, -0.321042638540326, -0.430478115531005, -0.488513320442751, -0.588000246169921,
     1      -0.754765039959403, -0.948897280768392, -1.27976500725669, -1.5161340602198, -2.06982415253695, -2.43130880701782,
     1      -2.91990260094065, -3.14328949115078, -3.59175920286139, -3.73548289693155, -3.8619921765288, -4.27387883486744,
     1      -4.8739084917412/
      Data c1a / 0.132636216298042, 0.132636216298042, 0.121913678741804, 0.0954458190369074, 0.0780343270248127,
     1      0.0797499701125017, 0.126409272490449, 0.237857890187294, 0.298285920791561, 0.358866830018622, 0.36378185586797,
     1      0.356638682110162, 0.412416702936915, 0.491578096184498, 0.470551468547512, 0.397454699266285, 0.407088725126396,
     1      0.370003161641225, 0.296244421779607, 0.249035757405967, 0.131863234230354, -0.0602151515860512, -0.170108793250811,
     1      -0.317888319174529, -0.213657545750405/
      Data c1b / -0.255, -0.255, -0.255, -0.255, -0.255, -0.255, -0.254, -0.253, -0.252, -0.25, -0.248, -0.2449, -0.2382, -0.2313,
     1      -0.2146, -0.1972, -0.162, -0.14, -0.1184, -0.11, -0.104, -0.102, -0.101, -0.101, -0.1/
      Data c1c / 0.374104895476179, 0.374104895476179, 0.394779299081704, 0.450367178792894, 0.502125276282922, 0.508541914329817,
     1      0.428533580148518, 0.244753765147496, 0.117861600910215, -0.0501929400220308, -0.0759243806200689, -0.0673618474199209,
     1      -0.193405747617701, -0.336816868253795, -0.340127470321498, -0.254193921547823, -0.234709477226598, -0.165352574104918,
     1      -0.114381565715004, -0.0730641756720992, 0.0136196696946929, 0.297168386190877, 0.457028397025905, 0.601615377508697,
     1      0.389474218593599/
      Data c1d / 0.579885366931696, 0.579885366931696, 0.571532295941719, 0.542986936430351, 0.539763891586565, 0.52141078447487,
     1      0.522557997517347, 0.491383246687325, 0.497598699249502, 0.538049677027715, 0.569106172547023, 0.582648948364874,
     1      0.543121898449559, 0.526919263871434, 0.566056464373584, 0.508787011694501, 0.453602572723648, 0.386616890673916,
     1      0.318326198639264, 0.302169672767348, 0.311997432161276, 0.377860259670876, 0.0664297305804952, 0.285446228730392,
     1      -0.0136805779652169/
      Data cn / 16.0875, 16.0875, 15.7118, 15.8819, 16.4556, 17.6453, 20.1772, 19.9992, 18.7106, 16.6246, 15.3709, 13.7012,
     1      11.2667, 9.1908, 6.5459, 5.2305, 3.7896, 3.3024, 2.8498, 2.5417, 2.1488, 1.8957, 1.7228, 1.5737, 1.5265/
      Data cM / 4.9993, 4.9993, 4.9993, 4.9993, 4.9993, 4.9993, 5.0031, 5.0172, 5.0315, 5.0547, 5.0704, 5.0939, 5.1315, 5.167,
     1      5.2317, 5.2893, 5.4109, 5.5106, 5.6705, 5.7981, 5.9983, 6.1552, 6.2856, 6.5428, 6.7415/
      Data c3 / 1.23775675850738, 1.23775675850738, 1.2100424798881, 1.17153664416715, 1.1376437803661, 1.08246488917746,
     1      1.10369070050752, 1.13570908666677, 1.24297786137047, 1.30184311690581, 1.36171523129664, 1.44527356984657,
     1      1.594640592675, 1.72776997676668, 1.89169957915458, 1.97790604653059, 2.23537178645213, 2.49270651199425,
     1      2.75122931083948, 2.8970091263776, 3.00181841245375, 3.13027839998167, 3.1358113461517, 3.03809613747909,
     1      2.73950225730737/
      Data c5 / 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.8305, 7.1333, 7.3621, 7.4365, 7.4972, 7.5416, 7.56,
     1      7.5735, 7.5778, 7.5808, 7.5814, 7.5817, 7.5818, 7.5818, 7.5818, 7.5818, 7.5818, 7.5818/
      Data cHM / 3.0956, 3.0956, 3.0963, 3.0974, 3.0988, 3.1011, 3.1094, 3.2381, 3.3407, 3.43, 3.4688, 3.5146, 3.5746, 3.6232,
     1      3.6945, 3.7401, 3.7941, 3.8144, 3.8284, 3.833, 3.8361, 3.8369, 3.8376, 3.838, 3.838/
      Data c6 / 0.4908, 0.4908, 0.4925, 0.4992, 0.5037, 0.5048, 0.5048, 0.5048, 0.5048, 0.5045, 0.5036, 0.5016, 0.4971, 0.4919,
     1      0.4807, 0.4707, 0.4575, 0.4522, 0.4501, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45/
      Data c7 / 0.0156320196168315, 0.0156320196168315, 0.0155659066225557, 0.0146375738629334, 0.0141084359110348,
     1      0.0144675463384911, 0.0152557523823881, 0.017607076212419, 0.018295042859516, 0.0212021767029952, 0.0216044681666273,
     1      0.0216325540991847, 0.0206030525582857, 0.0192274906261871, 0.0203438004744285, 0.0221610151832725, 0.0181240815582568,
     1      0.0150093892866716, 0.0152967789184208, 0.0158948163799267, 0.010176392608306, 0.00462828820896198,
     1      0.00101109863417872, -0.00259603156324936, -0.00775735563319467/
      Data c7b / 0.0373479124617126, 0.0373479124617126, 0.0376647327000708, 0.0410806029766391, 0.0427431654503428,
     1      0.042628675341005, 0.0397954134478023, 0.0344321279195352, 0.0317755358269565, 0.0234082190447402, 0.0220299636347817,
     1      0.0202156172523618, 0.0193661188593883, 0.0198764817033073, 0.0175409091342005, 0.0153015300326434, 0.0202430270775378,
     1      0.023988718198258, 0.0173169162080137, 0.0127766342571745, 0.0135894379296717, 0.0201630919552565, 0.0311749656859344,
     1      0.041123655355196, 0.0482773665664252/
      Data c8b / 0.4833, 0.4833, 1.2144, 1.6421, 1.9456, 2.181, 2.6087, 2.9122, 3.1045, 3.3399, 3.4719, 3.6434, 3.8787, 4.0711,
     1      4.3745, 4.6099, 5.0376, 5.3411, 5.7688, 6.0723, 6.5, 6.8035, 7.0389, 7.4666, 7.77/
      Data c9 / 0.9228, 0.9228, 0.9296, 0.9396, 0.9661, 0.9794, 1.026, 1.0177, 1.0008, 0.9801, 0.9652, 0.9459, 0.9196, 0.8829,
     1      0.8302, 0.7884, 0.6754, 0.6196, 0.5101, 0.3917, 0.1244, 0.0086, 0, 0, 0/
      Data c9a / 0.1202, 0.1202, 0.1217, 0.1194, 0.1166, 0.1176, 0.1171, 0.1146, 0.1128, 0.1106, 0.115, 0.1208, 0.1208, 0.1175,
     1      0.106, 0.1061, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1/
      Data c9b / 6.8607, 6.8607, 6.8697, 6.9113, 7.0271, 7.0959, 7.3298, 7.2588, 7.2372, 7.2109, 7.2491, 7.2988, 7.3691, 6.8789,
     1      6.5334, 6.526, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5/
      Data c11 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      Data c11b / -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.444,
     1      -0.3539, -0.2688, -0.1793, -0.1428, -0.1138, -0.1062, -0.102, -0.1009, -0.1003, -0.1001, -0.1001, -0.1, -0.1/
      Data gamma1 / -0.00640307889116, -0.00640307889116, -0.00638131012042546, -0.00625997127234089, -0.00642538130455206,
     1      -0.00690019669933074, -0.00848542829220626, -0.00990011050213089, -0.0103515868006968, -0.0120997918097204,
     1      -0.0122285761449459, -0.0119797715281222, -0.010902997808741, -0.00964286794696091, -0.00672196790194451,
     1      -0.00368546023588717, -0.00352284967615295, -0.00448724464278282, -0.00541975711124784, -0.00682327611588849,
     1      -0.00573735719392469, -0.00497998238605772, -0.0040599277758327, -0.00299825668266239, -0.0012793982951989/
      Data gamma2 / -0.0201255906298917, -0.0201255906298917, -0.0207779647283925, -0.023413720264442, -0.0245651180616863,
     1      -0.0240875383609615, -0.0201742310538795, -0.0148236450210469, -0.0119579201248049, -0.00572565425101764,
     1      -0.00400723668978553, -0.00255605013467344, -0.00170317615421141, -0.00204042057882375, -0.00358961151902558,
     1      -0.00755476059376392, -0.00812678237342369, -0.00629773345694908, -0.00193473747610703, 0.000540012774488294,
     1      0.00214131525339182, 0.00122238039154492, -0.00215278384414996, -0.00462114847977663, -0.00863524911340661/
      Data gm / 4.2542, 4.2542, 4.2386, 4.2519, 4.296, 4.3578, 4.5455, 4.7603, 4.8963, 5.0644, 5.1371, 5.188, 5.2164, 5.1954,
     1      5.0899, 4.7854, 4.3304, 4.1667, 4.0029, 3.8949, 3.7928, 3.7443, 3.709, 3.6632, 3.623/
      Data phi1 / -0.495117853043077, -0.495117853043077, -0.486937271935439, -0.476934166316566, -0.461711595080219,
     1      -0.444512512093825, -0.423067885770678, -0.440786048766373, -0.454699191928289, -0.483160893946246, -0.498472875843845,
     1      -0.510388096539616, -0.524579693588633, -0.57786783775584, -0.62231921014123, -0.640207000667697, -0.77744017182197,
     1      -0.88016132249611, -0.970289132476145, -0.978205396429725, -0.996982646440291, -0.999903132506138, -0.967067027725684,
     1      -0.880204381174846, -0.802423619903068/
      Data phi2 / -0.1417, -0.1417, -0.1364, -0.1403, -0.1591, -0.1862, -0.2538, -0.2943, -0.3077, -0.3113, -0.3062, -0.2927,
     1      -0.2662, -0.2405, -0.1975, -0.1633, -0.1028, -0.0699, -0.0425, -0.0302, -0.0129, -0.0016, 0, 0, 0/
      Data phi3 / -4.96041757793564, -4.96041757793564, -4.92276178884602, -4.91251089610306, -4.96513625402504, -5.04104295651768,
     1      -5.16134191152096, -5.18427465050738, -5.16799110487246, -5.14216868405811, -5.1228525972173, -5.09276768363488,
     1      -5.04538203042609, -5.00505091580236, -4.94414555282742, -4.90155669907002, -4.81342512480855, -4.77429914899606,
     1      -4.86562627218628, -5.34080741816899, -6.30453280595018, -6.48707320506911, -6.54311216539423, -6.59367473267582,
     1      -6.59953555531281/
      Data phi4 / -2.28130316824122, -2.28130316824122, -2.22229626176524, -2.12119730536018, -2.01259817873199, -1.90429902597891,
     1      -1.65759927434112, -1.4668018434511, -1.37369802906703, -1.32250111763273, -1.3277990635296, -1.36550006882328,
     1      -1.4629983153057, -1.5736992160094, -1.79900163048511, -2.01119988561788, -2.46330564055044, -2.83710591027157,
     1      -3.44869787781101, -3.92632479013901, -4.64152301545389, -5.22525279569103, -5.73744267614884, -6.78200407367658,
     1      -7.57134365730054/
      Data phi5 / 0.123790392326453, 0.123790392326453, 0.122085061000195, 0.124067312794112, 0.131759301635619, 0.142345432287797,
     1      0.166089530281892, 0.189041259587371, 0.179097115599691, 0.158033112380444, 0.138419563263086, 0.117610235180921,
     1      0.113214950210755, 0.108641249845724, 0.109834885229066, 0.0891870961945205, 0.077818193600786, 0.0816052270055997,
     1      0.15609806082559, 0.208065905730029, 0.25476022675178, 0.239629693652289, 0.233671524012721, 0.221987195100417,
     1      0.173323672829449/
      Data phi6 / 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
     1      300, 300, 300/
      Data tau1 / 0.4, 0.4, 0.40260697107594, 0.406254677901244, 0.409536599100115, 0.412362151356064, 0.417868294570895,
     1      0.421889928589278, 0.424427487574697, 0.427481452707345, 0.429157682348667, 0.431284929241235, 0.434101130455258,
     1      0.436302393283048, 0.439579110622969, 0.441948893180584, 0.445868511544276, 0.448355188391465, 0.451470005001172,
     1      0.453423408395737, 0.455849456889221, 0.457360780855491, 0.458424490145996, 0.460138736069957, 0.461201624142458/
      Data tau2 / 0.26, 0.26, 0.2637242443942, 0.268935254144634, 0.273623713000165, 0.277660216222949, 0.285526135101279,
     1      0.291271326556111, 0.294896410820996, 0.29925921815335, 0.301653831926667, 0.304692756058908, 0.308715900650368,
     1      0.311860561832925, 0.316541586604241, 0.319926990257977, 0.325526445063252, 0.329078840559236, 0.333528578573103,
     1      0.336319154851053, 0.339784938413173, 0.341943972650702, 0.343463557351422, 0.345912480099939, 0.347430891632083/
      Data sigma1 / 0.491174072226289, 0.491174072226289, 0.490392102266757, 0.49884629286719, 0.504893375136174,
     1      0.509564336873576, 0.517904200145505, 0.523587992281292, 0.526987000584287, 0.530830818870199, 0.532795897236441,
     1      0.535092985589259, 0.537693444822054, 0.539459005039374, 0.542235827503619, 0.543286483861406, 0.529446500271604,
     1      0.510487062001946, 0.478270612673322, 0.468082689790949, 0.461651268718433, 0.457088102848047, 0.453548636419445,
     1      0.447117215346928, 0.442554049476542/
      Data sigma2 / 0.376232982006498, 0.376232982006498, 0.376232982006498, 0.384851051972815, 0.390965669923574,
     1      0.395708536077288, 0.404326606043605, 0.410441223994364, 0.414316427806967, 0.419059293960682, 0.421719608502208,
     1      0.42517391191144, 0.429916778065155, 0.433791981877758, 0.439906599828517, 0.444649465982231, 0.453267535948548,
     1      0.459382153899307, 0.468000223865624, 0.468082689790949, 0.461651268718433, 0.457088102848047, 0.453548636419445,
     1      0.447117215346928, 0.442554049476542/
      Data sigma3 / 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7999, 0.7997, 0.7988, 0.7966, 0.7792, 0.7504,
     1      0.7136, 0.7035, 0.7006, 0.7001, 0.7, 0.7, 0.7/
      Data gscaleJapIt / 1.581682, 1.581682, 1.574012, 1.554376, 1.550154, 1.539148, 1.480416, 1.40939, 1.368192, 1.324078,
     1      1.307137, 1.293133, 1.314989, 1.351437, 1.405064, 1.440233, 1.52797, 1.652328, 1.887186, 2.134757, 3.575186, 3.864586,
     1      3.729218, 2.376267, 1.767935/
      Data phi1Jp / -0.684621362856666, -0.684621362856666, -0.668104490153807, -0.631364652348356, -0.585522686567435,
     1      -0.545659123699481, -0.468510173737368, -0.498485421565652, -0.560262194908129, -0.645129113225058, -0.698050600056553,
     1      -0.765310551536843, -0.846945325590174, -0.899885267655382, -0.961772613950091, -0.994485940646074, -1.02254530165975,
     1      -1.0001741238606, -0.924505251934862, -0.862638634236657, -0.788170491672007, -0.719494506504419, -0.655991993310966,
     1      -0.520223846723544, -0.406803274628334/
      Data phi5Jp / 0.459, 0.459, 0.458, 0.462, 0.453, 0.436, 0.383, 0.375, 0.377, 0.379, 0.38, 0.384, 0.393, 0.408, 0.462, 0.524,
     1      0.658, 0.78, 0.96, 1.11, 1.291, 1.387, 1.433, 1.46, 1.464/
      Data phi6Jp / 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800,
     1      800, 800, 800/
      Data gscaleWen / 0.75937, 0.75937, 0.760565, 0.76416, 0.767627, 0.773871, 0.795556, 0.79319, 0.776833, 0.743674, 0.721908,
     1      0.692222, 0.657901, 0.636173, 0.604851, 0.550656, 0.358158, 0.200263, 0.035633, 0, 0, 0, 0, 0, 0/
      Data sigma2jap / 0.452832136505141, 0.452832136505141, 0.455081474488854, 0.457089607558063, 0.464232248075043,
     1      0.471593304149618, 0.502235058788957, 0.522972224766759, 0.527828908598616, 0.530427671390517, 0.530985974389397,
     1      0.531176947346726, 0.53091866775, 0.53069348684478, 0.530953922960708, 0.531256171294916, 0.530904125925903,
     1      0.53016308922439, 0.527603492956116, 0.516719512970182, 0.491681001472351, 0.468235579988731, 0.451748161403632,
     1      0.416734686252964, 0.375505067353143/
      Data c8(1:25) / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0991, 0.1982, 0.2154, 0.2154, 0.2154, 0.2154,  
     1                0.2154, 0.2154, 0.2154, 0.2154 /
  
 
C Find the requested spectral period and corresponding coefficients
      nPer = 25
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T      = c1(1)
         c1aT     = c1a(1)
         c1bT     = c1b(1)
         c1cT     = c1c(1)
         c1dT     = c1d(1)
         c3T      = c3(1)
         cHMT     = CHM(1)
         cnT      = cn(1)
         cmT      = cm(1)
         c5T      = c5(1)
         c6T      = c6(1)
         c7T      = c7(1) 
         c7bT      = c7b(1) 
         c8bT      = c8b(1)
         c9T      = c9(1)
         c9aT      = c9a(1)
         c9bT      = c9b(1)
         c11T      = c11(1)
         c11bT      = c11b(1)
         gamma1T  = gamma1(1)
         gamma2T  = gamma2(1)
         gmT      = gm(1)
         phi1T    = phi1(1)
         phi2T    = phi2(1)
         phi3T    = phi3(1)
         phi4T    = phi4(1)
         phi5T    = phi5(1)
         phi6T    = phi6(1)
         tau1T    = tau1(1)
         tau2T    = tau2(1)
         sigma1T = sigma1(1)
         sigma2T = sigma2(1)
         sigma2JapT = sigma2Jap(1)
         sigma3T = sigma3(1)
         gscaleJapItT = gscaleJapIt(1)
         phi1jpT = phi1jp(1)
         phi5jpT = phi5jp(1)
         phi6jpT = phi6jp(1)
         gscaleWenT = gscaleWen(1)
         c8T = c8(1)
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
      write (*,*) 'Adjusted Chiou and Youngs (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c1a(count1),c1a(count2),
     +                   specT,c1aT,iflag)
            call S24_interp (period(count1),period(count2),c1b(count1),c1b(count2),
     +                   specT,c1bT,iflag)
            call S24_interp (period(count1),period(count2),c1c(count1),c1c(count2),
     +                   specT,c1cT,iflag)
            call S24_interp (period(count1),period(count2),c1d(count1),c1d(count2),
     +                   specT,c1dT,iflag)
            call S24_interp (period(count1),period(count2),cn(count1),cn(count2),
     +                   specT,cnT,iflag)
            call S24_interp (period(count1),period(count2),cm(count1),cm(count2),
     +                   specT,cmT,iflag)
            call S24_interp (period(count1),period(count2),cHM(count1),cHM(count2),
     +                   specT,cHMT,iflag)
            call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +                   specT,c3T,iflag)
            call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +                   specT,c5T,iflag)
            call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +                   specT,c6T,iflag)
            call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +                   specT,c7T,iflag)
            call S24_interp (period(count1),period(count2),c7b(count1),c7b(count2),
     +                   specT,c7bT,iflag)
            call S24_interp (period(count1),period(count2),c8b(count1),c8b(count2),
     +                   specT,c8bT,iflag)
            call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +                   specT,c9T,iflag)
            call S24_interp (period(count1),period(count2),c9a(count1),c9a(count2),
     +                   specT,c9aT,iflag)
            call S24_interp (period(count1),period(count2),c9b(count1),c9b(count2),
     +                   specT,c9bT,iflag)
            call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +                   specT,c11T,iflag)
            call S24_interp (period(count1),period(count2),c11b(count1),c11b(count2),
     +                   specT,c11bT,iflag)
            call S24_interp (period(count1),period(count2),gamma1(count1),gamma1(count2),
     +                   specT,gamma1T,iflag)
            call S24_interp (period(count1),period(count2),gamma2(count1),gamma2(count2),
     +                   specT,gamma2T,iflag)
            call S24_interp (period(count1),period(count2),gm(count1),gm(count2),
     +                   specT,gmT,iflag)
            call S24_interp (period(count1),period(count2),phi1(count1),phi1(count2),
     +                   specT,phi1T,iflag)
            call S24_interp (period(count1),period(count2),phi2(count1),phi2(count2),
     +                   specT,phi2T,iflag)
            call S24_interp (period(count1),period(count2),phi3(count1),phi3(count2),
     +                   specT,phi3T,iflag)
            call S24_interp (period(count1),period(count2),phi4(count1),phi4(count2),
     +                   specT,phi4T,iflag)
            call S24_interp (period(count1),period(count2),phi5(count1),phi5(count2),
     +                   specT,phi5T,iflag)
            call S24_interp (period(count1),period(count2),phi6(count1),phi6(count2),
     +                   specT,phi6T,iflag)
            call S24_interp (period(count1),period(count2),tau1(count1),tau1(count2),
     +                   specT,tau1T,iflag)
            call S24_interp (period(count1),period(count2),tau2(count1),tau2(count2),
     +                   specT,tau2T,iflag)
            call S24_interp (period(count1),period(count2),sigma1(count1),sigma1(count2),
     +                   specT,sigma1T,iflag)
            call S24_interp (period(count1),period(count2),sigma2(count1),sigma2(count2),
     +                   specT,sigma2T,iflag)
            call S24_interp (period(count1),period(count2),sigma2Jap(count1),sigma2Jap(count2),
     +                   specT,sigma2JapT,iflag)
            call S24_interp (period(count1),period(count2),sigma3(count1),sigma3(count2),
     +                   specT,sigma3T,iflag)
            call S24_interp (period(count1),period(count2),gscaleJapIt(count1),gscaleJapIt(count2),
     +                   specT,gscaleJapItT,iflag)
            call S24_interp (period(count1),period(count2),phi1jp(count1),phi1jp(count2),
     +                   specT,phi1jpT,iflag)
            call S24_interp (period(count1),period(count2),phi5jp(count1),phi5jp(count2),
     +                   specT,phi5jpT,iflag)
            call S24_interp (period(count1),period(count2),phi6jp(count1),phi6jp(count2),
     +                   specT,phi6jpT,iflag)
            call S24_interp (period(count1),period(count2),gscaleWen(count1),gscaleWen(count2),
     +                   specT,gscaleWenT,iflag)
            call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +                   specT,c8T,iflag)
 1011 period1 = specT                                                                                                              
c     Set the fault mechanism term.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
         if (ftype .eq. -1) then
            F_RV = 0.0
            F_NM = 1.0
         elseif (ftype .ge. 0.5) then
            F_RV = 1.0
            F_NM = 0.0
         else
            F_RV = 0.0
            F_NM = 0.0
         endif
C     Constant terms
        c2 = 1.06
        c4 = -2.1
        c4a = -0.5
        cRB = 50.0
C        c8 = 0.2154
        c8a = 0.2695
        pi = atan(1.0)*4.0
        d2r = pi/180.0
        cc = c5T* cosh(c6T * max((M-cHMT),0.0))
        gamma = gamma1T + gamma2T/cosh(max((M-gmT),0.0))
C     Apply Regional scaling factor.
C     Regionflag = 0 Global
C     Regionflag = 1 Japan and Italy
C        Also set sigma2 equal to Japan specific value
C     Regionflag = 2 Wenchuan (note only for M7.9)
c        if (regionflag .eq. 1 ) then
c           gamma = gamma * gscaleJapItT
c           sigma2T = sigma2JapT
c        elseif (regionflag .eq. 2) then
c           gamma = gamma * gscaleWenT        
c        endif
c        cosDELTA = cos(abs(DELTA)*d2r)
c Magnitude scaling
        r1 = c1T + c2 * (M-6.0) +
     1       (c2-c3T)/cnT *
     1             alog(1.0 + exp(-cnT*(M-cMT)))
c Near-field magnitude and distance scaling
        r2 = c4 * alog(Rrup + cc)
c Distance scaling at large distance
        r3 = (c4a-c4)/2.0 *
     1            alog( Rrup*Rrup+cRB*cRB ) +
     1       Rrup * gamma
c Center Z_TOR on the Z_TOR-M relation in Chiou and Youngs (2013)
        if (F_RV.EQ.1) then
          if (M .le. 5.849) then
              mZ_TOR = 2.704*2.704
          else
              mZ_TOR = max(2.704-1.226*(M-5.849), 0.0)
              mZ_TOR = mZ_TOR * mZ_TOR
          endif
        else
          if (M .le. 4.970) then
              mZ_TOR = 2.673*2.673
          else
              mZ_TOR = max(2.673-1.136*(M-4.970), 0.0)
              mZ_TOR = mZ_TOR * mZ_TOR
          endif
        endif
        deltaZ_TOR = Dtor - 2*mZ_TOR
        
c Scaling with other source variables (F_RV, F_NM, deltaZ_TOR, and Dip)
        coshM = cosh(2*max(M-4.5,0.0))
        cosDELTA = cos(DELTA*d2r)
        r4 = (c1aT+c1cT/coshM) * F_RV +
     1       (c1bT+c1dT/coshM) * F_NM +
     1       (c7T +c7bT/coshM) * deltaZ_TOR +
     1       (c11T+c11bT/coshM)* cosDELTA**2
        
c HW effect
        if (HWFlag .eq. 0) then
           hw = 0.0
        else
         hw = c9T * (cosDELTA) * (c9aT+(1-c9aT)
     1        *tanh(abs(Rx)/c9bT)) *
     1          (1.0 - sqrt(Rbjf**2+DTor**2)/(Rrup + 1))
        endif
C     Current version of the code sets cDPP=0 (i.e., no directivity)
c Directivity effect
        cDPP = 0.0
        fd = c8T * exp(-c8aT * (M-c8bT)**2) *
     1       max(0.0, 1.0-max(0.0,Rrup-40.0)/30.0) *
     1       min(max(0.0,M-5.5)/0.8, 1.0) * cDPP
c Predicted median Sa on reference condition (Vs=1130 m/sec)
c        fd = 0.0
        psa_ref = r1+r2+r3+r4+hw+fd
C     Set Phi1, Phi5, and Phi6 term for Japan is region is requested
c        if (regionflag .eq. 1) then
c           phi1T = phi1jpT
c           phi5T = phi5jpT
c           phi6T = phi6jpT
c        endif       
        
c Linear soil amplification
        a = phi1T * min(alog(Vs/1130.0), 0.0)
c Nonlinear soil amplification
        b = phi2T *
     1      (exp(-exp(phi3T)*(min(Vs,1130.0)-360.0)) - exp(-exp(phi3T)*(1130.0-360.0)))
        c = exp(phi4T)
C Deviation from ln(Vs30) scaling: bedrock depth (Z1) effect.
c        if (regionflag .eq. 1) then
c           deltaZ1 = depthvs10*1000.0 -
c     1     exp(-5.23/2.0 * alog((VS**2.0 + 412.0**2.0)/(1360.0**2.0 + 412.0**2.0)))      
c        else
           deltaZ1 = depthvs10*1000.0 -
     1     exp(-2.63 / 4.0 * alog((VS**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))
c        endif
        rkdepth = phi5T * ( 1.0 - exp(-deltaZ1/phi6T ) )
c Sa on soil condition
        psa = psa_ref + (a + b * alog((exp(psa_ref)+c)/c)) + rkdepth
C Compute the sigma term
        NL0 = b * exp(psa_ref)/(exp(psa_ref)+c)
        tau = tau1T +
     1            (tau2T-tau1T)/1.5*(min(max(M,5.0),6.5)-5.0)
        sigma_NL0 = sigma1T +
     1              (sigma2T-sigma1T)/1.5*(min(max(M,5.0),6.5)-5.0)
C     Current code set for Measured Vs30 values (i.e., Vs30class=1)
      if (vs30_class .eq. 0) then
         F_measured = 0.0
         F_Inferred = 1.0
      elseif (vs30_class .eq. 1) then      
         F_measured = 1.0
         F_Inferred = 0.0
      endif
        sigma_NL0 = sigma_NL0 *
     1        sqrt(0.7*F_Measured+F_Inferred*sigma3T+(1+NL0)**2.0)
        sigma = sqrt((tau*(1.0+NL0))**2.0+sigma_NL0**2.0)
      phi = sigma_NL0
      tau = (tau*(1.0+NL0))
      
C     Convert ground motion to units of gals.
      lnY = psa + 6.89
      period2 = period1
      return
      end 
  
C  ***** PEER NGA-West 2 MODELS (2013) **********
c ---------------------------------------------------------------------            
C ** Idriss (NGA-2013) Horizontal **
C     PEER Report 2013/08
C        NGA-West2 Model for Estimating Average Horizontal Values of 
C            Pseudo-Absolute Spectral Accelerations Generated by
C            Crustal Earthquakes
C         I. M. Idriss
C     Notes:
C        Applicable Range (see Abstract):  
C           5 <= M <= 8.5
C           Vs>=450 m/sec
C              for Vs>1200 use Vs=1200
C           Rrup <= 150 km
C        Mechanisms: Strike-slip and Norml (0)
C                    Reverse and Oblique (1)
c ---------------------------------------------------------------------            
      Subroutine S04_I14_TW_E04 ( m, Rrup, ftype, vs30, specT,
     1                     period2, lnY, sigma, iflag, Ztor )
      implicit none
      integer MAXPER
      parameter (MAXPER=25)
      REAL Period(MAXPER), a1mlt675(MAXPER), a2mlt675(MAXPER), a3mlt675(MAXPER)
      REAL b1mlt675(MAXPER), b2mlt675(MAXPER)
      REAL a1(MAXPER), a2(MAXPER), a3(MAXPER), b1(MAXPER), b2(MAXPER) ,a4(MAXPER)
      real gam(MAXPER), phi(MAXPER), xsi(MAXPER), period1
      REAL M, Rrup, Vs30, specT, sigma
      REAL SOF, period2, lnY, PhiT, gamT, XsiT, ftype, Ztor, a5
      real a1T, a2T, a3T, b1T, b2T, a1mlt675T, a2mlt675T, a3mlt675T, b1mlt675T, b2mlt675T, a4T
      integer iflag, count1, count2, nPer, i
 
 
      Data Period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4,
     1      5, 7.5, 10/
      Data a1mlt675 / 1.71685064584986, 1.71685064584986, 1.77876166629157, 1.92536250277134, 2.40902758641293, 2.859631936143,
     1      2.93099394107436, 3.61835712622066, 3.38293263767327, 3.04156524316998, 3.28039978748316, 3.45038410819511,
     1      3.16616863849086, 3.53647601520311, 3.17395975608739, 3.03648720441387, 2.08200569354534, 0.877460000459805,
     1      0.216213417632509, -1.23963578496273, -2.80372057210227, -3.56608887992888, -4.55737815257571, -5.78609036852066,
     1      -7.22273678696677/
      Data a2mlt675 / 0.528852639669921, 0.528852639669921, 0.516037115360633, 0.492476077090731, 0.357607293500749,
     1      0.236783913545005, 0.228518265363278, 0.191844760129495, 0.249131657210971, 0.330989114828171, 0.335257653142616,
     1      0.346063934766942, 0.404409985365426, 0.388859493480024, 0.461123202179168, 0.50021234991494, 0.695805122055084,
     1      0.943592129322998, 1.03770822530097, 1.20591426565113, 1.3828482598272, 1.34717380241795, 1.38254546549063,
     1      1.29518525654939, 1.31032284619722/
      Data a3mlt675 / 0.0589, 0.0589, 0.0589, 0.0589, 0.0492134623324272, 0.0417, 0.0527, 0.0442, 0.0391188387593093, 0.0329,
     1      0.0267654491028554, 0.0188, 0.0095, -0.0039, -0.0133, -0.0224, -0.0267, -0.0198, -0.0367, -0.0291, -0.0214, -0.024,
     1      -0.0202, -0.0219, -0.0035/
      Data b1mlt675 / 2.9935, 2.9935, 2.9935, 2.9935, 2.92192099200299, 2.8664, 2.9406, 3.019, 2.91472377949414, 2.7871,
     1      2.81929551534672, 2.8611, 2.8289, 2.8423, 2.83, 2.856, 2.7544, 2.7339, 2.68, 2.6837, 2.6907, 2.5782, 2.5468, 2.4478,
     1      2.3922/
      Data b2mlt675 / -0.2287, -0.2287, -0.2287, -0.2287, -0.236077537409605, -0.2418, -0.2513, -0.2516, -0.23900951196997,
     1      -0.2236, -0.223295447827801, -0.2229, -0.22, -0.2284, -0.2318, -0.2337, -0.2392, -0.2398, -0.2417, -0.245, -0.2389,
     1      -0.2514, -0.2541, -0.2593, -0.2586/
      Data xsi / -0.409999751830372, -0.409999751830372, -0.403233611021475, -0.393481052577978, -0.37706088187339,
     1      -0.355814430765091, -0.321300496735889, -0.323366702130062, -0.335833018872878, -0.362629820180182, -0.380383821924299,
     1      -0.38935871479828, -0.407130850107857, -0.467922466005464, -0.514670896613802, -0.543860145153594, -0.70876791754058,
     1      -0.822530157344685, -0.910136480902729, -0.916903880880618, -0.931048836251643, -0.94349528103761, -0.916043672142136,
     1      -0.835331966695324, -0.761357905540046/
      Data gam / -0.00281171552594942, -0.00281171552594942, -0.00289820257356861, -0.00350046984876955, -0.00577140539086469,
     1      -0.00763284577998198, -0.00852343893560511, -0.00773824638190913, -0.00782724759447787, -0.00734771093084663,
     1      -0.00622034033001967, -0.00471131925907849, -0.00360158865115855, -0.0033778158066833, -0.00232228393467311,
     1      -0.000932812535934199, -0.00222731372835322, -0.00290633165497401, -0.00366008106428892, -0.00476789995604717,
     1      -0.00273212517266628, -0.00475317864417355, -0.00509860339188468, -0.00633316100977675, -0.00643937122066843/
      Data phi / 0.173688616576542, 0.173688616576542, 0.174464843218782, 0.176059052545283, 0.182307217395956, 0.186414664900519,
     1      0.18044164297813, 0.18598354274442, 0.175314142592085, 0.151918884140306, 0.146073813578121, 0.157299856583112,
     1      0.161995973984013, 0.151671212694175, 0.1338550576422, 0.120877380466125, 0.145695260398602, 0.158201912612096,
     1      0.153992092718155, 0.148919063082084, 0.0844120327268032, 0.0211462842495985, 0.0199788102248451, -0.116039977911879,
     1      -0.0450237899128813/
      Data a1 / 9.0138, 9.0138, 9.0408, 9.1338, 8.48609726910026, 7.9837, 7.756, 9.4252, 9.51468239707057, 9.6242,
     1      10.2793352298525, 11.13, 11.3629, 11.7818, 11.6097, 11.4484, 10.9065, 9.8565, 8.3363, 6.8656, 4.1178, 1.8102, 0.0977,
     1      -3.0563, -4.4387/
      Data a2 / -0.0794, -0.0794, -0.0794, -0.0794, -0.142981982713312, -0.1923, -0.1614, -0.1887, -0.133751512954654, -0.0665,
     1      -0.111443199125889, -0.1698, -0.1766, -0.2798, -0.3048, -0.2911, -0.3097, -0.2565, -0.232, -0.1226, 0.1724, 0.3001,
     1      0.4609, 0.6948, 0.8393/
      Data a3 / 0.0589, 0.0589, 0.0589, 0.0589, 0.0492134623324272, 0.0417, 0.0527, 0.0442, 0.0391188387593093, 0.0329,
     1      0.0267654491028554, 0.0188, 0.0095, -0.0039, -0.0133, -0.0224, -0.0267, -0.0198, -0.0367, -0.0291, -0.0214, -0.024,
     1      -0.0202, -0.0219, -0.0035/
      Data b1 / 2.9935, 2.9935, 2.9935, 2.9935, 2.88424486584249, 2.7995, 2.8143, 2.8131, 2.63143724413814, 2.4091,
     1      2.44595081283604, 2.4938, 2.3773, 2.3772, 2.3413, 2.3477, 2.2042, 2.1493, 2.0408, 2.0013, 1.9408, 1.7763, 1.703,
     1      1.5212, 1.4195/
      Data b2 / -0.2287, -0.2287, -0.2287, -0.2287, -0.230502146542804, -0.2319, -0.2326, -0.2211, -0.197043174656907, -0.1676,
     1      -0.167991567078541, -0.1685, -0.1531, -0.1595, -0.1594, -0.1584, -0.1577, -0.1532, -0.147, -0.1439, -0.1278, -0.1326,
     1      -0.1291, -0.122, -0.1145/
      data a4 / 0.034345982, 0.034345982, 0.034178643, 0.034303764, 0.033911183, 0.033619895, 0.035247655,   
     1           0.0379925, 0.038355387, 0.037966827, 0.037984881, 0.037240959, 0.035830469, 0.035193632,   
     1           0.033500003, 0.031625775, 0.030589191, 0.030691086, 0.028359981, 0.027695125, 0.021529407,   
     1           0.016492213, 0.014525774, 0.013942019, 0.009132389 / 
       
C Find the requested spectral period and corresponding coefficients
      nPer = 25
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         a1T     = a1(1)
         a2T     = a2(1)
         a3T     = a3(1)
         a4T     = a4(1)
         b1T     = b1(1)
         b2T     = b2(1)
         a1mlt675T = a1mlt675(1)
         a2mlt675T = a2mlt675(1)
         a3mlt675T = a3mlt675(1)
         b1mlt675T = b1mlt675(1)
         b2mlt675T = b2mlt675(1)
         xsiT   = xsi(1)
         gamT   = gam(1)
         phiT   = phi(1)
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
      write (*,*) 'Adjusted Idriss (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),b1(count1),b1(count2),
     +                   specT,b1T,iflag)
            call S24_interp (period(count1),period(count2),b2(count1),b2(count2),
     +                   specT,b2T,iflag)
            call S24_interp (period(count1),period(count2),a1mlt675(count1),a1mlt675(count2),
     +                   specT,a1mlt675T,iflag)
            call S24_interp (period(count1),period(count2),a2mlt675(count1),a2mlt675(count2),
     +                   specT,a2mlt675T,iflag)
            call S24_interp (period(count1),period(count2),a3mlt675(count1),a3mlt675(count2),
     +                   specT,a3mlt675T,iflag)
            call S24_interp (period(count1),period(count2),b1mlt675(count1),b1mlt675(count2),
     +                   specT,b1mlt675T,iflag)
            call S24_interp (period(count1),period(count2),b2mlt675(count1),b2mlt675(count2),
     +                   specT,b2mlt675T,iflag)
            call S24_interp (period(count1),period(count2),phi(count1),phi(count2),
     +                   specT,phiT,iflag)
            call S24_interp (period(count1),period(count2),gam(count1),gam(count2),
     +                   specT,gamT,iflag)
            call S24_interp (period(count1),period(count2),xsi(count1),xsi(count2),
     +                   specT,xsiT,iflag)
 1011 period1 = specT                                                                                                              
C.....Compute the Ground motion.......
C.....Set the mechanism term.....................
C     Strike-slip and normal events --> SOF = 0
C     Reverse and oblique events     --> SOF = 1
C     Otherwise assume SOF = 0
      if (ftype .gt. 0.0) then
         SOF = 1.0
      else
         SOF = 0.0
      endif
      
      a1T = a1mlt675T + (a2mlt675T - a2T) * 6.75
      a5 = 10.0
   
      if (m .le. 6.75) then
         lnY = a1mlt675T + a2mlt675T*m + a3mlt675T*(8.5 - m)**2.0 - (b1mlt675T+b2mlt675T*m) * alog(Rrup+10.0) +
     1         xsiT*alog(Vs30) + gamT*rRup + SOF*phiT + a4T * min(max(Ztor, a5), 50.0)
      else
         lnY = a1T + a2T*m + a3T*(8.5 - m)**2.0 - (b1T+b2T*m) * alog(Rrup+10.0) +
     1         xsiT*alog(Vs30) + gamT*rRup + SOF*phiT + a4T * min(max(Ztor, a5), 50.0)
      endif
C     Convert ground motion to units of gals.
      lnY = lnY + 6.89
C     Compute Sigma which is Period and magnitude dependent.
C     Note report does not state a limit on sigma for M<5 but 
C     Since model is only applicable for M>=5 a limit is retained
C     for sigma with M<5 equal to M=5 values. 
      if (specT .le. 0.05) then
         if (m .le. 5.0) then
            sigma = 1.18 + 0.035*alog(0.05) - 0.06*5.0
         elseif (m .ge. 7.5) then
            sigma = 1.18 + 0.035*alog(0.05) - 0.06*7.5
         else
            sigma = 1.18 + 0.035*alog(0.05) - 0.06*m
         endif
      elseif (specT .ge. 3.0) then
         if (m .le. 5.0) then
            sigma = 1.18 + 0.035*alog(3.00) - 0.06*5.0
         elseif (m .ge. 7.5) then
            sigma = 1.18 + 0.035*alog(3.00) - 0.06*7.5
         else
            sigma = 1.18 + 0.035*alog(3.00) - 0.06*m
         endif
      else
         if (m .le. 5.0) then
            sigma = 1.18 + 0.035*alog(specT) - 0.06*5.0
         elseif (m .ge. 7.5) then
            sigma = 1.18 + 0.035*alog(specT) - 0.06*7.5
         else
            sigma = 1.18 + 0.035*alog(specT) - 0.06*m
         endif
      endif
      period2 = period1
      return
      end 

c ---------------------------------------------------------------------------            

      subroutine S04_CB14_TW_E05 ( mag, Rrup, Rbjf, Ftype, specT, 
     1                     period2, lnY, sigma, iflag, vs,
     2                     depthtop, D25, Dip, depth, HWflag, Rx, rupwidth, regionflag, phi, tau ) 

C     Last Updated: 5/17/17
C     Coefficients updated from PEER Report version to be consistent with EQ Spectra paper in press. 
C     Minor change to T=5, 7.5, and 10 sec for coefficient C6

      parameter (MAXPER=25)
      REAL Period(MAXPER), C0(MAXPER), C1(MAXPER), C2(MAXPER), C3(MAXPER), C4(MAXPER), C5(MAXPER)
      REAL C6(MAXPER), C7(MAXPER), C8(MAXPER), C9(MAXPER), C10(MAXPER), C11(MAXPER), C12(MAXPER)
      REAL C13(MAXPER), C14(MAXPER), C15(MAXPER), C16(MAXPER), C17(MAXPER), C18(MAXPER), C19(MAXPER)
      REAL A2(MAXPER), h1(MAXPER), h2(MAXPER), h3(MAXPER), h4(MAXPER)
      REAL h5(MAXPER), h6(MAXPER)
      REAL K1(MAXPER), K2(MAXPER), K3(MAXPER)
      REAL C20(MAXPER), DC20CA(MAXPER), DC20JP(MAXPER), DC20CH(MAXPER)
      REAL T1(MAXPER), T2(MAXPER), phi1(MAXPER), Phi2(MAXPER), phic(MAXPER)
      REAL flnAF(MAXPER), rho(MAXPER)

      REAL MAG, RRUP, RBJF, VS, D25, FHWR, FHWM, FHWZ, FHWD, PGAROCK, C, N
      real lnY, ftype, Dip, pgasoil, Rx, R1, R2, f1, f2
      real fhypH, D25_RK
      INTEGER count1, count2, HWFlag, regionflag, iflag

      real c0T, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, c10T, c11T, c12T
      real c13T, c14T, c15T ,c16T, c17T, c18T, c19T, c20T, Dc20CAT, Dc20JPT, Dc20CHT
      real k1T, k2T, k3T, a2T, h1T, h2T, h3T, h4T, h5T, h6T
      real t1T, t2T, phi1T, phi2T, phicT
      real rhoT, flnAFT

      real alpha, tau, depthtop, depth, rupwidth, specT, period2, sigma
      real tau_lnyB, tau_lnPGAB, phi_lnY, phi_lnyB, phi_lnPGAB, phi, sigmatot

C.....MODEL COEFFICIENTS.....................

      Data Period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4,
     1      5, 7.5, 10/
      Data c0 / 0.335093047810229, 0.351983138983978, 0.391664799222186, 0.480136356183228, 0.641529297579733, 0.82188120439914,
     1      1.19798152569248, 1.32463106269343, 1.3182111773268, 1.23786207049298, 1.17310035741048, 1.09596095473213,
     1      0.817547341300483, 0.569893713035775, 0.248317116508005, -0.0557086278316989, -0.898682094186014, -1.44551136634166,
     1      -2.44640589145251, -2.92728614327044, -3.74974403104352, -4.35687560387378, -4.90171652535439, -5.41286572183955,
     1      -6.00453776820377/
      Data c1 / 0.723018626766693, 0.76923053586393, 0.746659591778324, 0.705276848354753, 0.667081586168126, 0.589258437601297,
     1      0.585904684521786, 0.734594671263828, 0.87979012976083, 1.03046934140068, 1.07123988953934, 1.15880253917983,
     1      1.26367510891572, 1.39239252903295, 1.62015530388173, 1.69555434775151, 2.14365561022274, 2.43294493581605,
     1      2.48341437559922, 2.51106466415091, 2.59601777463608, 3.29009556992616, 3.05440516654319, 2.42612031645475,
     1      1.2805453253274/
      Data c2 / 0.582382030368516, 0.548276929260374, 0.547314123358948, 0.546705112037665, 0.540056592867299, 0.516455399622689,
     1      0.535858038675047, 0.520460985630314, 0.63206103275597, 0.752652287959249, 0.794904011310226, 0.821406760373195,
     1      0.805438888112464, 0.82157932533524, 0.95711703484122, 1.06774878746493, 1.27846012307265, 1.55669349721454,
     1      1.88847628445144, 2.03025576587569, 2.12029432816302, 2.25744390137846, 2.42438178382811, 1.89305011622601,
     1      1.66640016976192/
      Data c3 / -0.17376145698123, 0.113398391799347, 0.117582734720799, 0.100296134767346, 0.0395073633176415,
     1      -0.00999677104319945, -0.111157212017641, -0.0377898768562306, -0.0567913232402645, 0.0585226979400388,
     1      0.116550362476395, 0.204916066823206, 0.378548290064969, 0.468826328904193, 0.681767427902487, 0.809717324614938,
     1      1.0301365781875, 1.27535046528841, 1.22484767070439, 1.37539927263016, 1.66802552597525, 2.06876921469128,
     1      2.54220014838118, 3.12904314886059, 3.0514290877133/
      Data c4 / -0.474, -0.474, -0.464, -0.452, -0.446368292053737, -0.442, -0.437, -0.417, -0.401261889962462, -0.382,
     1      -0.384610447190274, -0.388, -0.383, -0.37, -0.301, -0.266, -0.221, -0.123, 0.066, 0.14, 0.313, 0.467, 0.544, 0.559,
     1      0.417/
      Data c5 / -2.773, -2.773, -2.772, -2.782, -2.78706853715164, -2.791, -2.745, -2.633, -2.55430944981231, -2.458,
     1      -2.44190224232664, -2.421, -2.392, -2.376, -2.303, -2.296, -2.232, -2.158, -2.063, -2.104, -2.051, -1.986, -2.021,
     1      -2.179, -2.244/
      Data c6 / 0.248, 0.248, 0.247, 0.246, 0.242620975232242, 0.24, 0.227, 0.21, 0.197859172256757, 0.183, 0.182564925468288,
     1      0.182, 0.189, 0.195, 0.185, 0.186, 0.186, 0.169, 0.158, 0.158, 0.148, 0.135, 0.135, 0.165, 0.18/
      Data c7 / 6.768, 6.753, 6.502, 6.291, 6.30564244066028, 6.317, 6.861, 7.294, 7.62539963136186, 8.031, 8.18501638422618,
     1      8.385, 7.534, 6.99, 7.012, 6.902, 5.522, 5.65, 5.795, 6.632, 6.759, 7.978, 8.538, 8.468, 6.564/
      Data c8 / -0.0569703754274623, -0.0100762198849894, -0.0109390716992389, -0.018000792352, -0.024235150820661,
     1      -0.0259397807540723, -0.0345381504830972, -0.0111793533676018, -0.01961170589007, -0.026725795799002,
     1      -0.0271575498062692, 0.006939002639218, 0.0371637790899891, 0.0557042454003613, 0.0588682536305066, 0.063859150342796,
     1      0.136363714907047, 0.13480261954914, 0.130720642420933, 0.133979703679808, 0.0724452288941395, 0.0157343198390952,
     1      0.0042492609773788, -0.0792819595337705, 0.00670187662714426/
      Data c9 / -0.212, -0.214, -0.208, -0.213, -0.230458294633416, -0.244, -0.266, -0.229, -0.220906114837838, -0.211,
     1      -0.190116422477806, -0.163, -0.15, -0.131, -0.159, -0.153, -0.09, -0.105, -0.058, -0.028, 0, 0, 0, 0, 0/
      Data c10 / 0.72, 0.72, 0.73, 0.759, 0.796732443239964, 0.826, 0.815, 0.831, 0.794127856483483, 0.749, 0.755526117975686,
     1      0.764, 0.716, 0.737, 0.738, 0.718, 0.795, 0.556, 0.48, 0.401, 0.206, 0.105, 0, 0, 0/
      Data c11 / 0.9388705546745, 0.948575269196186, 0.993472316704622, 1.06377202850532, 1.12506666052944, 1.18052683161614,
     1      1.35342376723162, 1.5255133424148, 1.67508288178791, 1.84797588155146, 1.96251994563968, 2.1263747186292,
     1      2.34088720935678, 2.44781132569501, 2.57455153248388, 2.56737766230058, 2.09973952793916, 1.46652488636632,
     1      0.290614567997193, -0.566265658925559, -0.927671102759077, -0.935440318255442, -0.905182604795051, -0.824192926747527,
     1      -0.766146783910017/
      Data c12 / 2.186, 2.191, 2.189, 2.164, 2.14935755933972, 2.138, 2.446, 2.969, 3.2275546649024, 3.544, 3.61491714866912,
     1      3.707, 3.343, 3.334, 3.544, 3.016, 2.616, 2.47, 2.108, 1.327, 0.601, 0.568, 0.356, 0.075, -0.027/
      Data c13 / 1.42, 1.416, 1.453, 1.476, 1.51711146800772, 1.549, 1.772, 1.916, 2.02616677026276, 2.161, 2.29326265764056,
     1      2.465, 2.766, 3.011, 3.203, 3.333, 3.054, 2.562, 1.453, 0.657, 0.367, 0.306, 0.268, 0.374, 0.297/
      Data c14 / -0.0064, -0.007, -0.0167, -0.0422, -0.0557724161504944, -0.0663, -0.0794, -0.0294, 0.0126882028432437, 0.0642,
     1      0.0783834297338236, 0.0968, 0.1441, 0.1597, 0.141, 0.1474, 0.1764, 0.2593, 0.2881, 0.3112, 0.3478, 0.3747, 0.3382,
     1      0.3754, 0.3506/
      Data c15 / -0.202, -0.207, -0.199, -0.202, -0.279154398863807, -0.339, -0.404, -0.416, -0.411953057418919, -0.407,
     1      -0.365232844955611, -0.311, -0.172, -0.084, 0.085, 0.233, 0.411, 0.479, 0.566, 0.562, 0.534, 0.522, 0.477, 0.321, 0.174/
      Data c16 / 0.393, 0.39, 0.387, 0.378, 0.331256824046015, 0.295, 0.322, 0.384, 0.398838789463964, 0.417, 0.411344031087739,
     1      0.404, 0.466, 0.528, 0.54, 0.638, 0.776, 0.771, 0.748, 0.763, 0.686, 0.691, 0.67, 0.757, 0.621/
      Data c17 / 0.0645383794856222, 0.0657962877798967, 0.0664619533475307, 0.0671935057558888, 0.0675290599349501,
     1      0.0682878362543006, 0.0667416957763312, 0.0674023297085519, 0.0659884697784512, 0.0651695636645009, 0.0657845186239398,
     1      0.0642396102150384, 0.0631366384677199, 0.0604817336964761, 0.0569417635049725, 0.0532085062942838, 0.0457272916435991,
     1      0.0439923447205022, 0.0344872404738479, 0.0319850189076888, 0.0224746623950968, 0.012911468121748, 0.0161784929461332,
     1      0.0263137791585881, 0.0285809717700681/
      Data c18 / 0.0587992463139723, 0.031551652901079, 0.0310646621980837, 0.0317510822112648, 0.0344353254061649,
     1      0.0359754970453913, 0.0396525794316812, 0.0380584543140422, 0.0421865523728087, 0.0422123038395977, 0.040654194566602,
     1      0.0381399423031852, 0.0302700235769742, 0.0247168737031733, 0.0236645391552053, 0.019780083453734, 0.0146412584489922,
     1      0.0135373369523852, 0.0322768710375305, 0.0314997066099382, 0.0193706012367015, -0.00969879107064587,
     1      -0.0432513570690934, -0.0941747977724573, -0.0823344155286154/
      Data c19 / 0.00757, 0.00755, 0.00759, 0.0079, 0.00797321220330142, 0.00803, 0.00811, 0.00744, 0.0073140951196997, 0.00716,
     1      0.00703817913112053, 0.00688, 0.00556, 0.00458, 0.00401, 0.00388, 0.0042, 0.00409, 0.00424, 0.00448, 0.00345, 0.00603,
     1      0.00805, 0.0028, 0.00458/
      Data c20 / -0.00634813007182598, -0.00634726499913435, -0.00625681860523333, -0.00650954432760824, -0.00685797038214189,
     1      -0.00732353722115189, -0.00859294805144714, -0.00946052167518335, -0.00948362932067619, -0.00878015959373452,
     1      -0.00826846734104271, -0.00751759793134993, -0.00705019712409618, -0.0066510144444269, -0.00500810118556177,
     1      -0.00344975891987256, -0.00344009816618098, -0.00347600014315439, -0.00436750799911943, -0.00536369094637663,
     1      -0.00277649854773078, -0.00175430069690324, -0.000972176857147457, -0.00122656065080374, -0.00177728087421456/
      Data k1 / 865, 865, 865, 908, 990.222936015443, 1054, 1086, 1032, 962.752315834834, 878, 821.44031087739, 748, 654, 587, 503,
     1      457, 410, 400, 400, 400, 400, 400, 400, 400, 400/
      Data k2 / -1.186, -1.186, -1.219, -1.273, -1.31411146800772, -1.346, -1.471, -1.624, -1.76204570804354, -1.931,
     1      -2.04281415465008, -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0, 0, 0, 0, 0/
      Data k3 / 1.839, 1.839, 1.84, 1.841, 1.84212634158925, 1.843, 1.845, 1.847, 1.84924830143393, 1.852, 1.85374029812685, 1.856,
     1      1.861, 1.865, 1.874, 1.883, 1.906, 1.929, 1.974, 2.019, 2.11, 2.2, 2.291, 2.517, 2.744/
      Data a2 / 0.167, 0.168, 0.166, 0.167, 0.170379024767758, 0.173, 0.198, 0.174, 0.184791846882883, 0.198, 0.200610447190274,
     1      0.204, 0.185, 0.164, 0.16, 0.184, 0.216, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596/
      Data h1 / 0.241, 0.242, 0.244, 0.246, 0.248815853973132, 0.251, 0.26, 0.259, 0.256751698566066, 0.254, 0.246603732960889,
     1      0.237, 0.206, 0.21, 0.226, 0.217, 0.154, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117/
      Data h2 / 1.474, 1.471, 1.467, 1.467, 1.45686292569673, 1.449, 1.435, 1.449, 1.45439592344144, 1.461, 1.47100671422938,
     1      1.484, 1.581, 1.586, 1.544, 1.554, 1.626, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616/
      Data h3 / -0.715, -0.714, -0.711, -0.713, -0.706241950464484, -0.701, -0.695, -0.708, -0.711147622007507, -0.715,
     1      -0.717610447190274, -0.721, -0.787, -0.795, -0.77, -0.77, -0.78, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733,
     1      -0.733, -0.733/
      Data h4 / 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/
      Data h5 / -0.337, -0.336, -0.339, -0.338, -0.338, -0.338, -0.347, -0.391, -0.417080296633634, -0.449, -0.424635826224107,
     1      -0.393, -0.339, -0.447, -0.525, -0.407, -0.371, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128/
      Data h6 / -0.27, -0.27, -0.263, -0.259, -0.261252683178505, -0.263, -0.219, -0.201, -0.155134650747747, -0.099,
     1      -0.142072378639526, -0.198, -0.21, -0.121, -0.086, -0.281, -0.285, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756,
     1      -0.756, -0.756/
      Data Dc20CA / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      Data Dc20JP / -0.0035, -0.0035, -0.0035, -0.0034, -0.0035689512383879, -0.0037, -0.0037, -0.0034, -0.00322013588528528,
     1      -0.003, -0.00304350745317124, -0.0031, -0.0033, -0.0035, -0.0034, -0.0034, -0.0032, -0.003, -0.0019, -5e-04, 0, 0, 0,
     1      0, 0/
      Data Dc20CH / 0.0036, 0.0036, 0.0036, 0.0037, 0.0038689512383879, 0.004, 0.0039, 0.0042, 0.0042, 0.0042, 0.00415649254682876,
     1      0.0041, 0.0036, 0.0031, 0.0028, 0.0025, 0.0016, 6e-04, 0, 0, 0, 0, 0, 0, 0/
      Data t1 / 0.409, 0.404, 0.417, 0.446, 0.480916589266832, 0.508, 0.504, 0.445, 0.416671401932432, 0.382, 0.363291795136368,
     1      0.339, 0.34, 0.34, 0.356, 0.379, 0.43, 0.47, 0.497, 0.499, 0.5, 0.543, 0.534, 0.523, 0.466/
      Data t2 / 0.322, 0.325, 0.326, 0.344, 0.362584636222669, 0.377, 0.418, 0.426, 0.408463248815315, 0.387, 0.365681347946093,
     1      0.338, 0.316, 0.3, 0.264, 0.263, 0.326, 0.353, 0.399, 0.4, 0.417, 0.393, 0.421, 0.438, 0.438/
      Data phi1 / 0.734, 0.734, 0.738, 0.747, 0.76389512383879, 0.777, 0.782, 0.769, 0.769, 0.769, 0.765519403746301, 0.761, 0.744,
     1      0.727, 0.69, 0.663, 0.606, 0.579, 0.541, 0.529, 0.527, 0.521, 0.502, 0.457, 0.441/
      Data phi2 / 0.492, 0.492, 0.496, 0.503, 0.512573903508648, 0.52, 0.535, 0.543, 0.543, 0.543, 0.546915670785411, 0.552, 0.545,
     1      0.568, 0.593, 0.611, 0.633, 0.628, 0.603, 0.588, 0.578, 0.559, 0.551, 0.546, 0.543/
      Data flnaf / 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
     1      0.3, 0.3, 0.3/
      Data phic / 0.166, 0.166, 0.166, 0.165, 0.163310487616121, 0.162, 0.158, 0.17, 0.174496602867868, 0.18, 0.182610447190274,
     1      0.186, 0.191, 0.198, 0.206, 0.208, 0.221, 0.225, 0.222, 0.226, 0.229, 0.237, 0.237, 0.271, 0.29/
      Data rho / 1, 1, 0.998, 0.986, 0.958967801857936, 0.938, 0.887, 0.87, 0.872697961720721, 0.876, 0.873389552809726, 0.87,
     1      0.85, 0.819, 0.743, 0.684, 0.562, 0.467, 0.364, 0.298, 0.234, 0.202, 0.184, 0.176, 0.154/
   

      nPer = 25
      c = 1.88
      n = 1.18

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c0T = c0(1)
         c1T = c1(1)
         c2T = c2(1)
         c3T = c3(1)
         c4T = c4(1)
         c5T = c5(1)
         c6T = c6(1)
         c7T = c7(1)
         c8T = c8(1)
         c9T = c9(1)
         c10T = c10(1)
         c11T = c11(1)
         c12T = c12(1)
         c13T = c13(1)
         c14T = c14(1)
         c15T = c15(1)
         c16T = c16(1)
         c17T = c17(1)
         c18T = c18(1)
         c19T = c19(1)

         a2T = a2(1)
         h1T = h1(1)
         h2T = h2(1)
         h3T = h3(1)
         h4T = h4(1)
         h5T = h5(1)
         h6T = h6(1)

         k1T = k1(1)
         k2T = k2(1)
         k3T = k3(1)
         c20T = c20(1)
         Dc20CAT = Dc20CA(1)
         Dc20JPT = Dc20JP(1)
         Dc20CHT = Dc20CH(1)
         
         phi1T = phi1(1)
         phi2T = phi2(1)
         t1T = t1(1)
         t2T = t2(1)
         flnAFT = flnaf(1)
         phicT = phic(1)
         rhoT = rho(1)

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
      write (*,*) 'Campbell&Bozorgnia (NGA West2-2013) Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),c0(count1),c0(count2),
     +                   specT,c0T,iflag)
            call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +                   specT,c2T,iflag)
            call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +                   specT,c3T,iflag)
            call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +                   specT,c4T,iflag)
            call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +                   specT,c5T,iflag)
            call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +                   specT,c6T,iflag)
            call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +                   specT,c7T,iflag)
            call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +                   specT,c8T,iflag)
            call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +                   specT,c9T,iflag)
            call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +                   specT,c10T,iflag)
            call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +                   specT,c11T,iflag)
            call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +                   specT,c12T,iflag)
            call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +                   specT,c13T,iflag)
            call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +                   specT,c14T,iflag)
            call S24_interp (period(count1),period(count2),c15(count1),c15(count2),
     +                   specT,c15T,iflag)
            call S24_interp (period(count1),period(count2),c16(count1),c16(count2),
     +                   specT,c16T,iflag)
            call S24_interp (period(count1),period(count2),c17(count1),c17(count2),
     +                   specT,c17T,iflag)
            call S24_interp (period(count1),period(count2),c18(count1),c18(count2),
     +                   specT,c18T,iflag)
            call S24_interp (period(count1),period(count2),c19(count1),c19(count2),
     +                   specT,c19T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),h1(count1),h1(count2),
     +                   specT,h1T,iflag)
            call S24_interp (period(count1),period(count2),h2(count1),h2(count2),
     +                   specT,h2T,iflag)
            call S24_interp (period(count1),period(count2),h3(count1),h3(count2),
     +                   specT,h3T,iflag)
            call S24_interp (period(count1),period(count2),h4(count1),h4(count2),
     +                   specT,h4T,iflag)
            call S24_interp (period(count1),period(count2),h5(count1),h5(count2),
     +                   specT,h5T,iflag)
            call S24_interp (period(count1),period(count2),h6(count1),h6(count2),
     +                   specT,h6T,iflag)

            call S24_interp (period(count1),period(count2),k1(count1),k1(count2),
     +                   specT,k1T,iflag)
            call S24_interp (period(count1),period(count2),k2(count1),k2(count2),
     +                   specT,k2T,iflag)
            call S24_interp (period(count1),period(count2),k3(count1),k3(count2),
     +                   specT,k3T,iflag)

            call S24_interp (period(count1),period(count2),c20(count1),c20(count2),
     +                   specT,c20T,iflag)
            call S24_interp (period(count1),period(count2),Dc20CA(count1),Dc20CA(count2),
     +                   specT,Dc20CAT,iflag)
            call S24_interp (period(count1),period(count2),Dc20JP(count1),Dc20JP(count2),
     +                   specT,Dc20JPT,iflag)
            call S24_interp (period(count1),period(count2),Dc20CH(count1),Dc20CH(count2),
     +                   specT,Dc20CHT,iflag)
     
            call S24_interp (period(count1),period(count2),phi1(count1),phi1(count2),
     +                   specT,phi1T,iflag)
            call S24_interp (period(count1),period(count2),phi2(count1),phi2(count2),
     +                   specT,phi2T,iflag)
            call S24_interp (period(count1),period(count2),t1(count1),t1(count2),
     +                   specT,t1T,iflag)
            call S24_interp (period(count1),period(count2),t2(count1),t2(count2),
     +                   specT,t2T,iflag)
            call S24_interp (period(count1),period(count2),flnAF(count1),flnAF(count2),
     +                   specT,flnAfT,iflag)
            call S24_interp (period(count1),period(count2),phic(count1),phic(count2),
     +                   specT,phicT,iflag)
            call S24_interp (period(count1),period(count2),rho(count1),rho(count2),
     +                   specT,rhoT,iflag)

 1011 period1 = specT                                                                                                              

C.....COMPUTE ROCK PGA VALUE FIRST.........................
C.....MAGNITUDE DEPENDENCE (Eq 3.2)........................
      IF (MAG .LE. 4.5) THEN
         TERM1 = C0(1) + C1(1)*(MAG-4.5)
      elseif (mag .le. 5.5) then
         TERM1 = C0(1) + c2(1)*(MAG-4.5)
      elseif (mag .le. 6.5) then
         TERM1 = C0(1) + c2(1) + c3(1)*(mag-5.5)
      ELSE
         TERM1 = C0(1) + c2(1) + C3(1) + c4(1)*(mag-6.5)
      ENDIF
      
C.....Distance dependence (Eq 3.3).....
      R = SQRT( RRUP*RRUP+C7(1)*C7(1) )
      TERM2 = (C5(1) + C6(1)*MAG)*ALOG(R)

C.....SET UP STYLE OF FAULTING TERMS (Eq 3.4, 3.5, and 3.6).........

C     Set mechanism term and corresponding Frv and Fnm values.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C    -1,-0.5    Normal and NMl/Obl       -150 < Rake < -30.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C       0       Strike-Slip                    Otherwise
      IF (Ftype .EQ. 0.0) THEN
         TERM3 = 0.0
      ELSEIF (Ftype .ge. 0.5) THEN
         if (mag .le. 4.5) then
            TERM3 = 0.0
         elseif (mag .le. 5.5) then
            TERM3 = C8(1)*(mag-4.5)
         else
            TERM3 = c8(1)
         endif
      ELSEIF (Ftype .le. -0.5) THEN
         if (mag .le. 4.5) then
            TERM3 = 0.0
         elseif (mag .le. 5.5) then
            TERM3 = C9(1)*(mag-4.5)
         else
            TERM3 = C9(1)
         endif
      ENDIF

C.....SET UP HANGING WALL TERMS (Eq 3.7)..............
      if (HWflag .eq. 1) then
         R1 = rupwidth*cos(abs(dip)*3.14159/180.0)
         R2 = 62.0*mag - 350.0
         f1 = h1(1) + h2(1)*(Rx/R1) + h3(1)*(Rx/R1)**2.0
         f2 = h4(1) + h5(1)*((Rx-R1)/(R2-R1)) + h6(1)*((Rx-R1)/(R2-R1))**2.0
         if (Rrup .eq. 0.0) then
            fhwrrup = 1.0
         else
            fhwrrup = ((Rrup-Rbjf)/Rrup)         
         endif
         if (Rx .lt. R1) Then
            fhwr = f1*fhwrrup
         else
            fhwr = max(f2,0.0)*fhwrrup
         endif
         if (mag .le. 5.5) then
            fhwm = 0.0         
         elseif (mag .le. 6.5) then
            fhwm = (mag-5.5)*(1.0+a2(1)*(mag-6.5))
         else
            fhwm = 1.0 + a2(1)*(mag-6.5)        
         endif
         if (depthtop .le. 16.66) then
            fhwz = 1.0 - 0.06*depthtop 
         else
            fhwz = 0.0
         endif
         fhwd = (90.0 - dip)/45.0        
         TERM4 = c10(1)*fhwr*fhwm*fhwz*fhwd
      else   
         term4 = 0.0
      endif 

C.....NOW COMPUTE THE SITE CONDITION FACTORS...............
C.....(FOR PGA ROCK, VS=1100, i.e., Vs>k1)
      TERM5_RK = (C11(1) + K2(1)*n)*ALOG( 1100.0/K1(1) )

C.....NOW COMPUTE THE SEDIMENT DEPTH DEPENDENCE (Eq 3.17)............
C     For Rock PGA the D25 value should be set at the recommended value of D25=0.398
      D25_RK = 0.398
      
      TERM6_RK = C14(1)*(D25_RK-1.0)
 
C.....Now compute the hypocentral depth term (Eq 3.21).........
      if (depth .le. 7.0) then
         fhypH = 0.0
      elseif (depth .le. 20.0) then
         fhypH = depth - 7.0      
      else
         fhypH = 13.0
      endif
      if (mag .le. 5.5) then
          term7 = c17(1)*fhypH
      elseif (mag .le. 6.5) then
          term7 = (c17(1) + (c18(1)-c17(1))*(mag-5.5))*fhypH
      else
          term7 = c18(1)*fhypH     
      endif

C.....Compute Rupture Dip term (Eq 3.24)............
      if (mag .le. 4.5) then
          term8 = c19(1)*dip
      elseif (mag .le. 5.5) then
          term8 = c19(1)*(5.5-mag)*dip      
      else 
          term8 =0
      endif

C.....Compute anelastic attenuation term.....
      if (Rrup .le. 80.0) then
         term9 = 0.0
      else      
         term9 = (c20(1)+Dc20CA(1) ) * (Rrup-80.0)

      endif      
      
      PGAROCK = EXP(TERM1+TERM2+TERM3+TERM4+TERM5_RK+TERM6_RK+TERM7+TERM8+TERM9)
C    write(*,*) "fmag  = ", TERM1
C    write(*,*) "fdis  = ", TERM2
C    write(*,*) "fflt  = ", TERM3
C    write(*,*) "fhng  = ", TERM4
C    write(*,*) "fsite = ", TERM5_RK
C    write(*,*) "fsed  = ", TERM6_RK
C    write(*,*) "fhyp  = ", TERM7
C    write(*,*) "fdip  = ", TERM8
C    write(*,*) "fatn  = ", TERM9
            
C.....For PGA Specific Vs30m Value
      if (vs .le. k1(1) ) then
         term5 = c11(1)*alog(vs/k1(1)) + 
     1           k2(1)*(alog(pgarock+c*((vs/k1(1))**n)) - 
     2           alog(pgarock+c))
      else
         term5 = (c11(1) + k2(1)*n)*alog(vs/k1(1))
      endif


C.....NOW COMPUTE THE SEDIMENT DEPTH DEPENDENCE (Eq 3.17)............
C     For Rock PGA the D25 value should be set at the recommended value of D25=0.398
      if (D25 .le. 1.0) then
            TERM6 = C14(1)*(D25-1.0)
      elseif (D25 .GT. 1.0 .AND. D25 .LE. 3.0) then
        TERM6 = 0.0
      elseif (D25 .GT. 3.0) then
         TERM6 = c16(1)*k3(1)*exp(-0.75)*(1.0-exp(-0.25*(D25-3.0)))
      endif

      pgasoil = alog(pgarock) - term5_rk - term6_RK + term5 + term6
      psoil2 = (TERM1+TERM2+TERM3+TERM4+TERM5+TERM6+TERM7+TERM8+TERM9)
c    write(*,*) "PGAROCK = ", PGAROCK
c    write(*,*) "pgasoil = ", pgasoil
c    write(*,*) "psoil2  = ", psoil2


C.....NOW COMPUTE THE GROUND MOTION VALUES.................
C.....MAGNITUDE DEPENDENCE.................................
      IF (MAG .LE. 4.5) THEN
         TERM1 = C0T + C1T*(MAG-4.5)
      elseif (mag .le. 5.5) then
         TERM1 = C0T + c2T*(MAG-4.5)
      elseif (mag .le. 6.5) then
         TERM1 = C0T + c2T + c3T*(mag-5.5)
      ELSE
         TERM1 = C0T + c2T + C3T + c4T*(mag-6.5)
      ENDIF

C.....Distance dependence......
      R = SQRT( RRUP*RRUP+C7T*C7T )
      TERM2 = (C5T + C6T*MAG)*ALOG(R)

C.....SET UP STYLE OF FAULTING TERMS...........

C     Set mechanism term and corresponding Frv and Fnm values.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C    -1,-0.5    Normal and NMl/Obl       -150 < Rake < -30.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C       0       Strike-Slip                    Otherwise
      IF (Ftype .EQ. 0.0) THEN
         TERM3 = 0.0
      ELSEIF (Ftype .ge. 0.5) THEN
         if (mag .le. 4.5) then
            TERM3 = 0.0
         elseif (mag .le. 5.5) then
            TERM3 = C8T*(mag-4.5)
         else
            TERM3 = C8T
         endif
      ELSEIF (Ftype .le. -0.5) THEN
         if (mag .le. 4.5) then
            TERM3 = 0.0
         elseif (mag .le. 5.5) then
            TERM3 = C9T*(mag-4.5)
         else
            TERM3 = C9T
         endif
      ENDIF

C.....SET UP HANGING WALL TERMS................
      if (HWflag .eq. 1) then
         R1 = rupwidth*cos(abs(dip)*3.14159/180.0)
         R2 = 62.0*mag - 350.0
         f1 = h1T + h2T*(Rx/R1) + h3T*(Rx/R1)**2.0
         f2 = h4T + h5T*((Rx-R1)/(R2-R1)) + h6T*((Rx-R1)/(R2-R1))**2.0
         if (Rrup .eq. 0.0) then
            fhwrrup = 1.0
         else
            fhwrrup = ((Rrup-Rbjf)/Rrup)         
         endif
         if (Rx .lt. R1) Then
            fhwr = f1*fhwrrup
         else
            fhwr = max(f2,0.0)*fhwrrup
         endif
         if (mag .le. 5.5) then
            fhwm = 0.0         
         elseif (mag .le. 6.5) then
            fhwm = (mag-5.5)*(1.0+a2T*(mag-6.5))
         else
            fhwm = 1.0 + a2T*(mag-6.5)        
         endif
         if (depthtop .le. 16.66) then
            fhwz = 1.0 - 0.06*depthtop 
         else
            fhwz = 0.0
         endif
         fhwd = (90.0 - dip)/45.0        
         TERM4 = c10T*fhwr*fhwm*fhwz*fhwd
      else   
         term4 = 0.0
      endif 

C.....NOW COMPUTE THE SITE CONDITION FACTORS...............
      IF (VS .LE. K1T ) THEN
         TERM5 = C11T*ALOG( VS/K1T ) +
     1           K2T*( ALOG( PGAROCK+c*( (VS/K1T)**N) ) -
     2           ALOG( PGAROCK+c ) )
      ELSE
         TERM5 = ( C11T+K2T*n )*ALOG( VS/K1T )
      ENDIF

C.....NOW COMPUTE THE SEDIMENT DEPTH DEPENDENCE.............
      IF (D25 .LE. 1.0) THEN
         TERM6 = C14T*(D25-1.0)
      ELSEIF (D25 .GT. 1.0 .AND. D25 .LE. 3.0) THEN
         TERM6 = 0.0
      ELSEIF (D25 .GT. 3.0)  THEN
         TERM6 = c16T*k3T*exp(-0.75)*( 1.0 - exp(-0.25*(D25-3.0)))
      ENDIF
      
C.....Now compute the hypocentral depth term..........
      if (depth .le. 7.0) then
         fhypH = 0.0
      elseif (depth .le. 20.0) then
         fhypH = depth - 7.0      
      else
         fhypH = 13.0
      endif
      if (mag .le. 5.5) then
          term7 = c17T*fhypH
      elseif (mag .le. 6.5) then
          term7 = (c17T + (c18T-c17T)*(mag-5.5))*fhypH
      else
          term7 = c18T*fhypH     
      endif

C.....Compute Rupture Dip term.............
      if (mag .le. 4.5) then
          term8 = c19T*dip
      elseif (mag .le. 5.5) then
          term8 = c19T*(5.5-mag)*dip      
      else 
          term8 = 0
      endif

C.....Compute anelastic attenuation term.....
      if (Rrup .le. 80.0) then
         term9 = 0.0
      else
         term9 = (c20T+Dc20CAT)*(Rrup-80.0)
      endif

      LnY = (TERM1+TERM2+TERM3+TERM4+TERM5+TERM6+TERM7+TERM8+TERM9)

C   write(*,*) "fmag  = ", TERM1
C   write(*,*) "fdis  = ", TERM2
C   write(*,*) "fflt  = ", TERM3
C   write(*,*) "fhng  = ", TERM4
C   write(*,*) "fsite = ", TERM5
C   write(*,*) "fsed  = ", TERM6
C   write(*,*) "fhyp  = ", TERM7
C   write(*,*) "fdip  = ", TERM8
C   write(*,*) "fatn  = ", TERM9

C   write(*,*) "LnY = ", LnY
C   write(*,*) "Sa = ", exp(LnY)

C    Check that SA is not less than PGA for T<0.25sec
c     if (specT .lt. 0.25) then
c        if (lnY .lt. pgasoil ) then
c           lnY = pgasoil
c        endif
c     endif

C.....Now compute the sigma value..........
      IF (Vs .LT. k1T) THEN
        alpha = k2T*pgarock*(1/(pgarock  
     &    +c*(Vs/k1T)**n) 
     &    -1/(pgarock + c))
      ELSE
        alpha = 0.0
      ENDIF

      If (Mag.le.4.5) then
      tau_lnyB = t1T
       tau_lnPGAB = t1(1)
        elseif (Mag.lt.5.5) then
         tau_lnyB = t2T + 
     &          (t1T - t2T)*(5.5-mag)
         tau_lnPGAB = t2(1) + 
     &          (t1(1) - t2(1))*(5.5-Mag)
        else
         tau_lnyB = t2T
         tau_lnPGAB = t2(1)
        endif

      tau = SQRT(tau_lnyB**2 +  
     &           (alpha * tau_lnPGAB)**2 +
     &           2.0*alpha*rhoT*tau_lnyB*tau_lnPGAB)

      If (Mag.le.4.5) then
             phi_lny = phi1T
           phi_lnPGAB = phi1(1)
      elseif (Mag.lt.5.5) then
          phi_lny = phi2T + 
     &          (phi1T - phi2T)*(5.5-mag)
             phi_lnPGAB = phi2(1) + 
     &          (phi1(1) - phi2(1))*(5.5-mag)
      else
             phi_lny = phi2T 
          phi_lnPGAB = phi2(1) 
      endif

      phi_lnyB = SQRT(phi_lny**2 - flnAFT**2)

      phi_lnPGAB = SQRT(phi_lnPGAB**2 - flnAF(1)**2)

      phi = SQRT(phi_lny**2 + 
     &           (alpha*phi_lnPGAB)**2 +
     &           2.0*alpha*rhoT*phi_lnyB*phi_lnPGAB)
    
      Sigmatot = SQRT(phi**2 + Tau**2)

      period2 = period1
      
C     Convert ground motion to units of gals.

      lnY = lnY + 6.89
      sigma = sigmaTot

      return
      END

c ------------------------------------------------------------------            
C *** BCHydro Subduction (06/2010 - adjustefd Model Version E ) Horizontal ***********
c ------------------------------------------------------------------            
      subroutine S04_AGA16_TW_F10 ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, Ztor, depth, disthypo )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT, Ztor
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C
c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)
      
C     Compute Rock PGA
      call S04_AGA16_TW_F10_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, faba, Ztor, depth, disthypo )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S04_AGA16_TW_F10_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, Ztor, depth, disthypo )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
c ----------------------------------------------------------------------
      subroutine S04_AGA16_TW_F10_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, Ztor, depth, disthypo )

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=25)
      real a1(MAXPER), a2(MAXPER),dC1_itf(MAXPER) ,dC1_itb(MAXPER),
     1     a6(MAXPER), a10(MAXPER), a11(MAXPER), a11a(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a4a(MAXPER), a7(MAXPER), a8(MAXPER),
     1     a15(MAXPER), a16(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo,
     1     mag, a3, a4, a5, a9 ,a7T , a8T , a15T, a16T
      real a1T, a2T, a6T, a11aT, Ztor, a4aT
      real a10T, a11T, a12T, a13T, a14T, sigsT, sigtT, dC1_itfT, dC1_itbT
      real vLinT, b_soilT, sumgm, Ftype, tau, period1
      integer count1, count2, iflag
      real n, c, c4, c1, deltac1, faba, R, testmag, VsStar, depth, specT
      real base, fmag, fdepth, fsite, fbac

      Data Period(1:25) / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 
     1            0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10/
      Data a1(1:25) / 4.35844293849408, 4.40448385002841, 4.42235073045333, 4.51174431025914, 4.6184756364235, 4.73673816849304,
     1      5.1478261284161, 5.32838098916143, 5.43888858346181, 5.54992742949254, 5.53025724823576, 5.50209585586371,
     1      5.42171283115248, 5.1925226061065, 4.79362720821157, 4.36876337362239, 3.6146095430292, 3.14131767795372,
     1      2.41835445173874, 1.83392549350599, 1.09559787881207, 0.493492034594101, -0.0677174706304278, -1.28769323712076,
     1      -2.3529378451964/
      Data a2(1:25) / -1.35, -1.35, -1.35, -1.37212535246749, -1.3878235398683, -1.4, -1.45, -1.45, -1.45, -1.45,
     1      -1.42824627341438, -1.4, -1.35, -1.28, -1.18, -1.08, -0.91, -0.85, -0.77, -0.71, -0.64, -0.58, -0.54, -0.46, -0.4/
      Data a4a(1:25) / 0.490233295747377, 0.501501976617984, 0.487171144951632, 0.47215353028238, 0.427780842429563,
     1      0.384611647919739, 0.344528924478493, 0.328259792218493, 0.330061667882538, 0.383801593508627, 0.412846285048564,
     1      0.438633622467472, 0.509559046177354, 0.611238357609338, 0.821165221077697, 0.962541509212795, 1.23227360622235,
     1      1.35657920015304, 1.46310710338033, 1.39583793060953, 1.14163314305871, 0.828359653814432, 0.593068337437621,
     1      0.635183263258333, 0.334324407572139/
      Data a6(1:25) / 0, -1.2742310369559e-14, -1.04835518780313e-12, -1.38313281100888e-11, -8.62519073094022e-11,
     1      -3.56734213274011e-10, -4.70651720770408e-09, -2.93495449747076e-08, -9.3620855471152e-08, -3.87178855055009e-07,
     1      -8.58410410760826e-07, -2.41301821846127e-06, -9.95823573287351e-06, -3.15658328071745e-05, -0.000187819964423078,
     1      -0.000663253434661024, -0.00261135688215098, -0.00327342900048902, -0.00342678315937337, -0.0034378736922942,
     1      -0.00343983874290077, -0.00343997413984613, -0.00343999374744544, -0.00343999952608246, -0.00343999992400273/
      Data a10(1:25) / 2.50238267428753, 2.45354985208117, 2.45712506053231, 2.48314913393394, 2.51571572736229, 2.54835614651406,
     1      2.63796917560156, 2.66920174985232, 2.65980457395154, 2.65947307783084, 2.50370194429171, 2.28268053526727,
     1      1.97567149316003, 1.74543072521505, 1.46284974570948, 1.18569521344312, 0.759064884501228, 0.472238540201132,
     1      0.105031115220643, 0.0204724942547022, -0.0563468339410934, -0.0796952941502796, -0.0661674056310574,
     1      -0.0727359608786985, 0.132949225428625/
      Data a11(1:25) / 0.0209645459852199, 0.0211185450539438, 0.020907239199861, 0.0206861754405617, 0.0201992952920322,
     1      0.0200863107159024, 0.0213038771564475, 0.0210506222491249, 0.0212386849173846, 0.021974823626841, 0.0218096375422394,
     1      0.0214497986870152, 0.0209887711470591, 0.0202089910414603, 0.0216695638774524, 0.0211315710615843, 0.0200840312848651,
     1      0.0183419906458299, 0.0169140417758707, 0.0152229673546218, 0.0101721058379849, 0.00751807891402711,
     1      0.00662090116843197, 0.00548671441898383, 0.00608785305167983/
      Data a11a(1:25) / 0.0447925753474416, 0.0454773789670112, 0.0453951448647483, 0.0453555996274611, 0.0460678354504749,
     1      0.0463680385857737, 0.0485682360784842, 0.0516565897099376, 0.05242805478271, 0.0517743795770695, 0.0525758958487427,
     1      0.0513712553416057, 0.0525281422211252, 0.0514609710721228, 0.0465114906308642, 0.0433918508645454, 0.0326061440066766,
     1      0.0253674408634384, 0.0210352047730289, 0.0168769997670818, 0.00818548076208565, 0.0016268028317338,
     1      0.00586285954723722, 0.00358594148454458, -0.01151742409672/
      Data a12(1:25) / 0.943715435182963, 0.947291261826117, 0.952369004446095, 1.04329022070683, 1.12605577830093,
     1      1.19877799169209, 1.34969490713516, 1.50771872006435, 1.64463099510739, 1.80153619895198, 1.89875583891682,
     1      2.05409914544544, 2.2556751074265, 2.39309551778298, 2.51270985429006, 2.54116064397387, 2.12606260531648,
     1      1.58081589874553, 0.445249721776344, -0.380429301270945, -0.666755078187764, -0.622309935348892, -0.611246872553534,
     1      -0.504394400041745, -0.487936169204039/
      Data a13(1:25) / -0.0135, -0.0135, -0.0135, -0.0136327521148049, -0.0137269412392098, -0.0138, -0.0142, -0.0145,
     1      -0.0148597282294294, -0.0153, -0.0156915670785411, -0.0162, -0.0172, -0.0183, -0.0206, -0.0231, -0.0296, -0.0363,
     1      -0.0493, -0.061, -0.0798, -0.0935, -0.098, -0.098, -0.098/
      Data a14(1:25) / -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.378246273414381, -0.35, -0.31, -0.28, -0.23,
     1      -0.19, -0.12, -0.07, 0, 0, 0, 0, 0, 0, 0/
      Data Vlin(1:25) / 865.1, 865.1, 865.1, 948.468328097495, 1007.61909822376, 1053.5, 1085.7, 1032.5, 962.847621576726, 877.6,
     1      821.301355596418, 748.2, 654.3, 587.1, 503, 456.6, 410.5, 400, 400, 400, 400, 400, 400, 400, 400/
      Data b_soil(1:25) / -1.186, -1.186, -1.186, -1.25680112789596, -1.30703532757857, -1.346, -1.471, -1.624, -1.76204570804354,
     1      -1.931, -2.04281415465008, -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0, 0, 0, 0, 0/
      Data a7(1:25) / 1.0988, 1.0988, 1.0988, 1.16730009123934, 1.21590167943226, 1.2536, 1.4175, 1.3997, 1.38103909809835, 1.3582,
     1      1.27405658556683, 1.1648, 0.994, 0.8821, 0.7046, 0.5799, 0.3687, 0.1746, -0.082, -0.2821, -0.4466, -0.4344, -0.4368,
     1      -0.4433, -0.4828/
      Data a8(1:25) / -1.42, -1.42, -1.42, -1.52177662135044, -1.59398828339419, -1.65, -1.8, -1.8, -1.75053736845345, -1.69,
     1      -1.60298509365752, -1.49, -1.3, -1.18, -0.98, -0.82, -0.54, -0.34, -0.05, 0.12, 0.3, 0.3, 0.3, 0.3, 0.3/
      Data a15(1:25) / 0.9969, 0.9969, 0.9969, 1.04384999793601, 1.07716155160054, 1.103, 1.2732, 1.3042, 1.28432501532402, 1.26,
     1      1.24390224232664, 1.223, 1.16, 1.05, 0.8, 0.662, 0.48, 0.33, 0.31, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3/
      Data a16(1:25) / -1, -1, -1, -1.07965126888296, -1.13616474352589, -1.18, -1.36, -1.36, -1.33302038279279, -1.3,
     1      -1.27824627341438, -1.25, -1.17, -1.06, -0.78, -0.62, -0.34, -0.14, 0, 0, 0, 0, 0, 0, 0/
      data dC1_itf / 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.14, 0.1, 0.04, 0,
     1      -0.06, -0.1, -0.2, -0.2, -0.2, -0.2, -0.2/
      data dC1_itb / -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3,
     1      -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3/
      data sigs / 0.542501761866534, 0.542540910308792, 0.542470428435724, 0.54743357078777, 0.560177645077094,
     1      0.580620973270876, 0.612843122153575, 0.620992129840455, 0.614264074899902, 0.598009282996332, 0.597938396739046,
     1      0.580250693487794, 0.573346409650175, 0.568315874508338, 0.560071611406591, 0.557525852693131, 0.562785431289318,
     1      0.551292654342151, 0.553950674471151, 0.56095249410327, 0.544349805621334, 0.509976706337188, 0.480872687826229,
     1      0.436932291079628, 0.4169542202548/
      data sigt/ 0.281527305592527, 0.277210251610433, 0.276724083277789, 0.278797924083198, 0.282199550169249,
     1      0.288521475585358, 0.301771906960262, 0.303128517793669, 0.310997206167293, 0.311309354242882, 0.306684702250237,
     1      0.306967620360309, 0.318131596803571, 0.312252665780728, 0.288987070059188, 0.282316238059581, 0.303740388310018,
     1      0.322857833694789, 0.313381866051055, 0.33245875554816, 0.3640745587897, 0.37531046943687, 0.380147961119054,
     1      0.307726020206149, 0.311840210476944/

C Constant parameters            
      n = 1.18
      c = 1.88
      a3 = 0.1
      a4 = 0.9
      a5 = 0.0
      a9 = 0.4
      c4 = 10.0
      c1 = 7.8
C      a7 = 0
C      a8 = 0
C      a15 = 0
C      a16 = 0
      
C Find the requested spectral period and corresponding coefficients
      nPer = 25

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a6T = a6(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a11aT = a11a(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         dC1_itfT = dC1_itf(i1)
         dC1_itbT = dC1_itb(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         a4aT = a4a(i1)
         a7T = a7(i1)
         a8T = a8(i1)         
         a15T = a15(i1)
         a16T = a16(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
 
C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'AGA16_TW_F10 Subduction Horizontal'
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
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a11a(count1),a11a(count2),
     +                   specT,a11aT,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),dC1_itf(count1),dC1_itf(count2),
     +                   specT,dC1_itfT,iflag)
            call S24_interp (period(count1),period(count2),dC1_itb(count1),dC1_itb(count2),
     +                   specT,dC1_itbT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)
            call S24_interp (period(count1),period(count2),a4a(count1),a4a(count2),
     +                   specT,a4aT,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,a16T,iflag)

 1011 period1 = specT                                                                                                              

C 2018/09/04 coef a6 revised
C       if( a6T > 0.0 ) then
C          a6T = 0.0
C       endif    
C 2018/09/11 cancelled
C     Compute the R term and base model based on either Rupture Distance 
c         (Interface events) of Hypocentral distance (Intraslab events). 
      if (ftype .eq. 0.0) then
         deltaC1 = dC1_itfT
         R = rRup + c4*exp( (mag-6.0)*a9 ) 
         base = a1T + a4*deltaC1 + (a2T + a14T*ftype + a3*(mag - 7.8))*alog(R) + a6T*rRup + a10T*ftype
      elseif (ftype .eq. 1.0) then
         deltaC1 = dC1_itbT
         R = disthypo + c4*exp( (mag-6.0)*a9 ) 
         base = a1T + a4*deltaC1 + (a2T + a14T*ftype + a3*(mag - 7.8))*alog(R) + a6T*disthypo + a10T*ftype
      else
         write (*,*) 'AGA16_TW_F10 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif
      
C     Base model for Magnitude scaling.      
      testmag = (7.8 + deltaC1)
C      if (mag .le. testmag ) then
C         fmag = a4*(mag-testmag) + a13T*(10.0-mag)**2.0
C      else
C         fmag = a5*(mag-testmag) + a13T*(10.0-mag)**2.0
C      endif      

C     The fmag term has an additional branch for linear Mag scaling
C     new  coefficient a4a      
      if (mag .le. 5.5 ) then
         fmag = a4aT*(mag-5.5) + a4*(5.5-testmag) + a13T*(10.0-mag)**2.0
      elseif (mag .le. testmag) then
         fmag = a4*(mag-testmag)+ a13T*(10.0-mag)**2.0
      else
        fmag = a5*(mag-testmag)+ a13T*(10.0-mag)**2.0
      endif           
      
C     Depth Scaling
      if (ftype .eq. 0.0) then
        fdepth = a11aT*(Ztor - 20.0)
      elseif (ftype .eq. 1.0) then
C        fdepth = a11T*(depth - 60.0 )
        fdepth =  a11T*(min(depth, 80.0) -60.0 )
      else
         write (*,*) 'AGA16_TW_F10 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif

C     Forearc/Backarc scaling      
      if (ftype .eq. 1) then
         fbac =  (a7T + a8T*alog(max(disthypo,85.0)/40.0))*faba
      elseif (ftype .eq. 0) then   
         fbac =  (a15T + a16T*alog(max(rRup,100.0)/40.0))*faba
      endif 

C     Site Response 
      if (vs30 .ge. 1000.0) then
          VsStar = 1000.0
      else
          VsStar = vs30
      endif
       
      if (vs30 .ge. VlinT) then
         fsite = a12T*alog(VsStar/vLinT) + b_soilT*n*alog(VsStar/vLinT)
      else
         fsite = a12T*alog(VsStar/vLinT) - b_soilT*alog(pgarock + c) +
     1          b_soilT*alog(pgarock + c*(VsStar/vlinT)**n)     
      endif

      sumgm = base + fmag + fdepth + fsite
c      write(*,*) "deltaC1 = ", deltaC1
c      write(*,*) "testmag = ", testmag
c      write(*,*) "base = ", base
c      write(*,*) "fmag = ", fmag
c      write(*,*) "fdepth = ", fdepth
c      write(*,*) "fsite = ", fsite
c      write(*,*) "lnYSa = ", sumgm
c      write(*,*) "Sa = ", exp(sumgm) 
   
C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = sumgm

      return
      end

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

c ------------------------------------------------------------------            
C *** Chao2018 (Crustal and Subduction - Model) Horizontal ***********
c ------------------------------------------------------------------            
      subroutine S04_Chao2018 ( mag, dist, ftype, lnY, sigma, specT, vs, Ztor, Z10,           
     1            vs30_class, attenName, period2, iflag, sourcetype, phi, tau, msasflag )         

      implicit none
C    2018/08/31 revised
      real mag, dip, fType, dist, vs, SA1180,
     1      Z10,  ZTOR, fltWidth, lnSa, sigma, lnY, vs30_rock, sourcetype
      real Fn, Frv, specT, period2, CRjb, phi, tau, z10_rock, SA_rock
      integer hwflag, iflag, vs30_class, regionflag, msasflag
      character*80 attenName                                                    

C     Set the reference spectrum.                
c     sourcetype = 0 for crustal
c                  1 for Subduction 
c     Vs30_class = 0 for estimated
c     Vs30_class = 1 for measured 
C     Mainshock and Aftershocks included based on MSASFlag
C         0 = Mainshocks
C         1 = Aftershocks

c     Compute SA1180
      vs30_rock = 1180.
      z10_rock = 0.004541444
      SA_rock = 0.
      
         call S04_Chaoetal2018 ( mag, dist, ftype, sigma, specT, vs30_rock, Ztor, z10_rock,
     1             SA_rock, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag)
      Sa1180 = exp(lnSa)

c     Compute Sa at spectral period for given Vs30

         call S04_Chaoetal2018 ( mag, dist, ftype, sigma, specT, vs, Ztor, Z10,
     1             sa1180, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag )

C     Convert ground motion to units of gals.
      lnY = lnSa + 6.89

      period2 = specT

      return
      end
c -------------------------------------------------------------------           
C **** Chao et al. 2018 (SSHAC model) *************
c -------------------------------------------------------------------           

      subroutine S04_Chaoetal2018 ( mag, dist, ftype, sigma, specT, vs, Ztor, Z10,           
     1            sa1180, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag )                                   

      implicit none
      
      integer MAXPER                                                                            
      parameter (MAXPER=21)                                                     
      real ftype, dist, mag, lnSa, sigma, specT, lnYref, vs, Ztor, Z10, period1
      real period(MAXPER), c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER), c11(MAXPER)
      real c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER), c16(MAXPER), c17(MAXPER)
      real c18(MAXPER), c19(MAXPER), c20(MAXPER), c21(MAXPER), c22(MAXPER), c23(MAXPER)
      real c24(MAXPER), c25(MAXPER), c26(MAXPER), c27(MAXPER), taucr1(MAXPER), taucr2(MAXPER)
      real tausb1(MAXPER), tausb2(MAXPER), phisscr1(MAXPER), phisscr2(MAXPER), phisssb1(MAXPER)
      real phisssb2(MAXPER), arfacr(MAXPER), arfasb(MAXPER), phis2s(MAXPER)
      character*80 attenName                                                    
      integer nper, count1, count2, C11flag, C23flag, C29flag, iflag, C10flag, C13flag
      integer vs30_class, n, i, msasflag
      integer Fcr, Fsb, Fcrss, Fcrno, Fcrro, Fsbintra, Fsbinter, Fas, Fkuo17, Fks17, Frf, Fmanila 
      real Mc, Mref, Mmax, Rrupref, Vs30ref, Zref, sourcetype
      real c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, c10T, c11T, c12T, c13T, c14T, c15T
      real c16T, c17T, c18T, c19T, c20T, c21T, c22T, c23T, c24T, c25T, c26T, c27T
      real taucr1T, taucr2T, tausb1T, tausb2T, phisscr1T, phisscr2T, phisssb1T, phisssb2T
      real arfacrT, arfasbT, phis2sT, phi, tau, fm, SA1180, Z10ref
      real Ssource, Spath, Ssite, Ssitelin, Ssitenon, Sztor, Smag, Sgeom, Sanel
      real taucr, tausb, phisscr, phisssb, phiss, sigmass
      real c28(MAXPER), c29(MAXPER), c28T, c29T, c30(MAXPER), c30T, h
                                                                                
      data Period / 0, -2, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75,
     1              1, 1.5, 2, 3, 4, 5 / 
      data c1 / -1.042666, 1.895411, -1.042360, -1.009040, -0.962601, -0.870563, -0.785856, -0.727776, 
     1          -0.674517, -0.668684, -0.682285, -0.707153, -0.782246, -0.867588, -1.100953, -1.335849,  
     1          -1.769889, -2.158271, -2.821190, -3.354776, -3.867019 / 
      data c2 / -1.137852, 1.735667, -1.137811, -1.101142, -1.049407, -0.939899, -0.832404, -0.758714,  
     1          -0.705514, -0.716938, -0.752775, -0.800309, -0.908875, -1.016989, -1.278817, -1.515900,  
     1          -1.926198, -2.280583, -2.876528, -3.355129, -3.810868 / 
      data c3 / -1.170990, 1.587449, -1.171510, -1.130461, -1.066827, -0.931421, -0.802566, -0.723218,  
     1          -0.690970, -0.741863, -0.819872, -0.908821, -1.075545, -1.217121, -1.508285, -1.746629,  
     1          -2.151367, -2.498264, -3.089388, -3.563977, -4.047973 / 
      data c4 / -1.314590, 1.882708, -1.298813, -1.262369, -1.214610, -1.092116, -0.964326, -0.894881,  
     1          -0.836942, -0.792193, -0.763280, -0.744597, -0.723252, -0.734591, -0.883649, -1.088087,  
     1          -1.526738, -1.947331, -2.666122, -3.258068, -3.809328 / 
      data c5 / -0.457036, 2.359015, -0.436214, -0.384335, -0.311696, -0.131992, 0.074716, 0.200368, 0.286613,  
     1          0.254600, 0.173035, 0.087029, -0.073676, -0.220276, -0.567236, -0.881354, -1.445217, -1.943705,  
     1          -2.722921, -3.349385, -3.792075 / 
      data c6 / -0.124087, -0.144577, -0.123275, -0.124986, -0.127909, -0.127668, -0.125831, -0.122221,  
     1          -0.115139, -0.105684, -0.095876, -0.088255, -0.081799, -0.084061, -0.105026, -0.125539,  
     1          -0.137664, -0.125372, -0.085938, -0.049451, 0.012311 / 
      data c7 / 0.198357, -0.038842, 0.191344, 0.195147, 0.198967, 0.227624, 0.260255, 0.289816, 0.348032,  
     1          0.370540, 0.361912, 0.339973, 0.258518, 0.169337, 0.025880, -0.053324, -0.113254, -0.138977,  
     1          -0.193895, -0.258713, -0.326743 / 
      data c8 / 0.676721, 1.251174, 0.684339, 0.657072, 0.610194, 0.544252, 0.548393, 0.602698, 0.739429,  
     1          0.875915, 0.986450, 1.078964, 1.233393, 1.354611, 1.564898, 1.701387, 1.863505, 1.953011,  
     1          2.043382, 2.088718, 2.105341 / 
      data c9 / 0.633834, 0.766604, 0.640018, 0.606922, 0.558754, 0.535818, 0.603777, 0.687135, 0.829456,  
     1          0.904226, 0.958892, 0.984563, 1.007287, 1.025986, 1.013332, 1.002811, 0.962813, 0.979659,  
     1          1.052087, 1.157880, 1.147008 / 
      data c10 / -0.135342, -0.210085, -0.136868, -0.131414, -0.122038, -0.108827, -0.109545, -0.120401,  
     1          -0.147755, -0.175139, -0.197287, -0.215793, -0.246679, -0.270922, -0.312979, -0.339111,  
     1          -0.357148, -0.352155, -0.323695, -0.294607, -0.250719 / 
      data c11 / -0.001320, -0.000002, -0.000002, -0.006396, -0.008323, -0.005923, -0.001799, -0.010118,  
     1          -0.080308, -0.163717, -0.238633, -0.297686, -0.346393, -0.334133, -0.229548, -0.131664,  
     1          -0.032403, -0.006790, -0.000014, -0.000002, -0.000001 / 
      data c12 / -0.000034, -0.000014, -0.000001, -0.000004, -0.000005, -0.000016, -0.000069, -0.009043,  
     1          -0.105380, -0.189907, -0.238782, -0.268913, -0.261917, -0.203064, -0.102740, -0.049805,  
     1          -0.011622, -0.002702, -0.000015, -0.000001, -0.000001 / 
      data c13 / -0.000038, -0.291529, -0.000001, -0.000004, -0.000005, -0.000018, -0.000080, -0.000082,  
     1          -0.011888, -0.070362, -0.139438, -0.217492, -0.359123, -0.475375, -0.670648, -0.775270,  
     1          -0.856757, -0.798850, -0.639826, -0.461969, -0.465749 / 
      data c14 / 0.032555, 0.015510, 0.032494, 0.033028, 0.034285, 0.037352, 0.039960, 0.040658, 0.038553,  
     1          0.034446, 0.030226, 0.026384, 0.020530, 0.016445, 0.010777, 0.008159, 0.005621, 0.004032,  
     1          0.000622, -0.002970, -0.008948 / 
      data c15 / 0.017948, 0.016982, 0.018059, 0.018141, 0.018067, 0.017682, 0.017907, 0.018645, 0.020842,  
     1          0.022069, 0.022125, 0.021736, 0.020462, 0.018980, 0.016322, 0.014953, 0.012425, 0.009349,  
     1          0.004309, -0.000454, 0.000502 / 
      data c16 / 0.006634, 0.003036, 0.006626, 0.007008, 0.007575, 0.008618, 0.009183, 0.009159, 0.008404,  
     1          0.007375, 0.006330, 0.005344, 0.003788, 0.002623, 0.001242, 0.000937, 0.000745, 0.000647,  
     1          0.000241, 0.000051, -0.000294 / 
      data c17 / -1.687472, -1.547786, -1.692461, -1.730349, -1.772040, -1.821470, -1.810850, -1.759544,  
     1          -1.643299, -1.551105, -1.488492, -1.446763, -1.398596, -1.375008, -1.347383, -1.333650,  
     1          -1.311795, -1.298714, -1.284844, -1.281901, -1.241532 / 
      data c18 / -1.525058, -1.549392, -1.539999, -1.565127, -1.589292, -1.642560, -1.682785, -1.676278,  
     1          -1.637058, -1.599726, -1.564986, -1.538416, -1.509663, -1.488699, -1.436447, -1.401868,  
     1          -1.333535, -1.279867, -1.210371, -1.152781, -1.113843 / 
      data c19 / 0.393879, 0.303285, 0.388200, 0.395410, 0.410623, 0.421673, 0.405707, 0.381442, 0.338659,  
     1          0.299124, 0.268677, 0.244569, 0.210459, 0.192077, 0.179792, 0.184209, 0.197091, 0.210361,  
     1          0.233088, 0.252812, 0.264010 / 
      data c20 / 0.192940, 0.223730, 0.188931, 0.201543, 0.217245, 0.203093, 0.152740, 0.113630, 0.072366,  
     1          0.055866, 0.044190, 0.040593, 0.041938, 0.044384, 0.071700, 0.099141, 0.153764, 0.199871,  
     1          0.253750, 0.282276, 0.298510 / 
      data c21 / -0.003430, -0.000817, -0.003338, -0.003155, -0.003190, -0.003740, -0.004703, -0.005494,  
     1          -0.006073, -0.005732, -0.005091, -0.004380, -0.003179, -0.002346, -0.001304, -0.000879,  
     1          -0.000564, -0.000475, -0.000459, -0.000460, -0.000879 / 
      data c22 / -0.003964, -0.000992, -0.003844, -0.003962, -0.004180, -0.004364, -0.004364, -0.004447,  
     1          -0.004295, -0.003755, -0.003171, -0.002613, -0.001680, -0.001075, -0.000534, -0.000357,  
     1          -0.000460, -0.000784, -0.001435, -0.002141, -0.002623 / 
      data c23 / -2.405229, -6.913488, -2.378641, -2.418630, -2.313209, -2.067232, -1.854300, -1.693198,  
     1          -1.479479, -1.358600, -1.314082, -1.320243, -1.391252, -1.483500, -1.558715, -1.404591,  
     1          -0.845899, -0.436331, -0.023390, 0.000000, 0.000000 / 
      data c24 / -0.478715, -0.672424, -0.477706, -0.470430, -0.451950, -0.416240, -0.402645, -0.410470,  
     1          -0.446615, -0.481032, -0.513790, -0.542807, -0.596688, -0.648563, -0.742180, -0.798604,  
     1          -0.842403, -0.849752, -0.839143, -0.822511, -0.797141 / 
      data c25 / 0.063345, 0.095748, 0.063485, 0.064655, 0.068531, 0.079314, 0.084622, 0.082002, 0.070590,  
     1          0.062519, 0.061465, 0.064345, 0.074007, 0.083563, 0.103379, 0.118159, 0.139920, 0.153222, 
     1           0.160180, 0.155343, 0.142666 / 
      data c26 / -0.604023, 0.384203, -0.598271, -0.547902, -0.480919, -0.338469, -0.228587, -0.197662,  
     1          -0.245076, -0.344219, -0.449905, -0.550299, -0.725882, -0.874283, -1.154718, -1.353035,  
     1          -1.618376, -1.771092, -1.893009, -1.912582, -1.838861 / 
      data c27 / -0.679356, 0.383524, -0.673352, -0.626423, -0.567634, -0.437480, -0.333207, -0.303548, 
     1           -0.352291, -0.441969, -0.532257, -0.615964, -0.760168, -0.886623, -1.140546, -1.331239,  
     1          -1.584174, -1.727634, -1.829259, -1.833069, -1.746378 / 
      data c28 / -0.650493, 0.305582, -0.644569, -0.596370, -0.529756, -0.377565, -0.258657, -0.227982,  
     1          -0.292081, -0.413965, -0.533268, -0.639714, -0.817367, -0.960564, -1.232582, -1.435782,  
     1          -1.715330, -1.875696, -1.998774, -2.017655, -1.941254 / 
      data c29  / -0.434706511, -0.34661917, -0.447102839, -0.415446676, -0.366999977, -0.366148593, -0.469248166,   
     1          -0.568918055, -0.719037119, -0.790626985, -0.840621112, -0.860236343, -0.873968925, -0.885430175,   
     1          -0.827526363, -0.742505362, -0.57420172, -0.49097237, -0.423101316, -0.426861488, -0.338104658  /
      data c30  / -0.444373594, -0.232141657, -0.460428073, -0.448496325, -0.436361774, -0.492546506, -0.54896049,   
     1          -0.60832263, -0.675866324, -0.689506592, -0.695043713, -0.678658902, -0.636360976, -0.607212602,   
     1          -0.503495503, -0.43304903, -0.313708518, -0.28272403, -0.308813684, -0.394492845, -0.370973146 /
      data taucr1 / 0.367027, 0.439871, 0.366713, 0.367381, 0.365814, 0.360738, 0.364067, 0.376999,  
     1          0.420380, 0.470706, 0.514142, 0.545754, 0.580926, 0.595221, 0.586341, 0.567362, 0.539173,  
     1          0.521031, 0.510998, 0.516268, 0.543015 / 
      data taucr2 / 0.315219, 0.374796, 0.315430, 0.319633, 0.327049, 0.345572, 0.360231, 0.362292, 0.339502,  
     1          0.307997, 0.285602, 0.273645, 0.271198, 0.286852, 0.345657, 0.391395, 0.441542, 0.462329,  
     1          0.475059, 0.480575, 0.470049 / 
      data tausb1 / 0.272750, 0.330733, 0.271058, 0.271673, 0.274006, 0.280698, 0.291948, 0.305046, 0.334109,  
     1          0.369176, 0.401260, 0.428927, 0.466550, 0.482919, 0.477346, 0.450142, 0.405823, 0.366654,  
     1          0.325573, 0.294383, 0.340797 / 
      data tausb2 / 0.532703, 0.566695, 0.536421, 0.557253, 0.577443, 0.601923, 0.599473, 0.579989, 0.528112,  
     1          0.495023, 0.475578, 0.464138, 0.466708, 0.482046, 0.524715, 0.574218, 0.638071, 0.666590,  
     1          0.648573, 0.596551, 0.507714 / 
      data phisscr1 / 0.530015, 0.558388, 0.530680, 0.521767, 0.513496, 0.502312, 0.503270, 0.517247, 0.553158,  
     1          0.581637, 0.599292, 0.608899, 0.610791, 0.600179, 0.561837, 0.526088, 0.477128, 0.448854,  
     1          0.421072, 0.408714, 0.406221 / 
      data phisscr2 / 0.433766, 0.435606, 0.435202, 0.444543, 0.455987, 0.470816, 0.466840, 0.451769,  
     1          0.427517, 0.417334, 0.416486, 0.422015, 0.436118, 0.448613, 0.471981, 0.484933, 0.492704,  
     1          0.488214, 0.469356, 0.447792, 0.429747 / 
      data phisssb1 / 0.430634, 0.478613, 0.431917, 0.426561, 0.421787, 0.407957, 0.406584, 0.419057,  
     1          0.454812, 0.477921, 0.489639, 0.495038, 0.494487, 0.491319, 0.486749, 0.481815, 0.481313,  
     1          0.479424, 0.475432, 0.466465, 0.454073 / 
      data phisssb2 / 0.495377, 0.485676, 0.494439, 0.498452, 0.505789, 0.519604, 0.525237, 0.521053, 
     1           0.511584, 0.503919, 0.496864, 0.491227, 0.487558, 0.486774, 0.488927, 0.493538, 0.488128, 
     1           0.477724, 0.446953, 0.414410, 0.368507 / 
      data phis2s / 0.342398, 0.273837, 0.342667, 0.348735, 0.363743, 0.407580, 0.443690, 0.455723,  
     1          0.439179, 0.411759, 0.389287, 0.372044, 0.352699, 0.344485, 0.341557, 0.347037, 0.357107, 
     1           0.363360, 0.369654, 0.372377, 0.375576 / 

c Set attenuation name                                                            
c     Sourcetype = 0 Crustal
c     Sourcetype = 1 Subduction 
                                                                       
C Find the requested spectral period and corresponding coefficients
      nper = 21

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
        c1T = c1(1)
        c2T = c2(1)
        c3T = c3(1)
        c4T = c4(1)
        c5T = c5(1)
        c6T = c6(1)
        c7T = c7(1)
        c8T = c8(1)
        c9T = c9(1)
        c10T = c10(1)
        c11T = c11(1)
        c12T = c12(1)
        c13T = c13(1)
        c14T = c14(1)
        c15T = c15(1)
        c16T = c16(1)
        c17T = c17(1)
        c18T = c18(1)
        c19T = c19(1)
        c20T = c20(1)
        c21T = c21(1)
        c22T = c22(1)
        c23T = c23(1)
        c24T = c24(1)
        c25T = c25(1)
        c26T = c26(1)
        c27T = c27(1)
        c28T = c28(1)
        c29T = c29(1)        
        c30T = c30(1)        
        taucr1T = taucr1(1)
        taucr2T = taucr2(1)
        tausb1T = tausb1(1)
        tausb2T = tausb2(1)
        phisscr1T = phisscr1(1)
        phisscr2T = phisscr2(1)
        phisssb1T = phisssb1(1)
        phisssb2T = phisssb2(1)
        phis2sT = phis2s(1)
       goto 1011
C   Function Form for PGV Regression     
       elseif (specT .eq. -2.0 .or. specT .eq. -1.0) then
         period1 = period(2)
         c1T = c1(2)
         c2T = c2(2)
         c3T = c3(2)
         c4T = c4(2)
         c5T = c5(2)
         c6T = c6(2)
         c7T = c7(2)
         c8T = c8(2)
         c9T = c9(2)
         c10T = c10(2)
         c11T = c11(2)
         c12T = c12(2)
         c13T = c13(2)
         c14T = c14(2)
         c15T = c15(2)
         c16T = c16(2)
         c17T = c17(2)
         c18T = c18(2)
         c19T = c19(2)
         c20T = c20(2)
         c21T = c21(2)
         c22T = c22(2)
         c23T = c23(2)
         c24T = c24(2)
         c25T = c25(2)
         c26T = c26(2)
         c27T = c27(2)
         c28T = c28(2)
         c29T = c29(2)        
         c30T = c30(2)        
         taucr1T = taucr1(2)
         taucr2T = taucr2(2)
         tausb1T = tausb1(2)
         tausb2T = tausb2(2)
         phisscr1T = phisscr1(2)
         phisscr2T = phisscr2(2)
         phisssb1T = phisssb1(2)
         phisssb2T = phisssb2(2)
         phis2sT = phis2s(2)
         goto 1011      
       endif
C Now loop over the spectral period range of the attenuation relationship.
         do i=3,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1010 
            endif
         enddo
        
      write (*,*) 
      write (*,*) 'Chao et al. (2018) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1010    call S24_interp (period(count1),period(count2),c1(count1),c1(count2), 
     +                specT,c1T,iflag)
         call S24_interp (period(count1),period(count2),c2(count1),c2(count2), 
     +                specT,c2T,iflag)
         call S24_interp (period(count1),period(count2),c3(count1),c3(count2), 
     +                specT,c3T,iflag)
         call S24_interp (period(count1),period(count2),c4(count1),c4(count2), 
     +                specT,c4T,iflag)
         call S24_interp (period(count1),period(count2),c5(count1),c5(count2), 
     +                specT,c5T,iflag)
         call S24_interp (period(count1),period(count2),c6(count1),c6(count2), 
     +                specT,c6T,iflag)
         call S24_interp (period(count1),period(count2),c7(count1),c7(count2), 
     +                 specT,c7T,iflag)
         call S24_interp (period(count1),period(count2),c8(count1),c8(count2), 
     +                 specT,c8T,iflag)
         call S24_interp (period(count1),period(count2),c9(count1),c9(count2), 
     +                 specT,c9T,iflag)
         call S24_interp (period(count1),period(count2),c10(count1),c10(count2), 
     +                 specT,c10T,iflag)
         call S24_interp (period(count1),period(count2),c11(count1),c11(count2), 
     +                 specT,c11T,iflag)
         call S24_interp (period(count1),period(count2),c12(count1),c12(count2), 
     +                 specT,c12T,iflag)
         call S24_interp (period(count1),period(count2),c13(count1),c13(count2), 
     +                 specT,c13T,iflag)
         call S24_interp (period(count1),period(count2),c14(count1),c14(count2), 
     +                 specT,c14T,iflag)
         call S24_interp (period(count1),period(count2),c15(count1),c15(count2), 
     +                 specT,c15T,iflag)
         call S24_interp (period(count1),period(count2),c16(count1),c16(count2), 
     +                 specT,c16T,iflag)
         call S24_interp (period(count1),period(count2),c17(count1),c17(count2), 
     +                 specT,c17T,iflag)
         call S24_interp (period(count1),period(count2),c18(count1),c18(count2), 
     +                 specT,c18T,iflag)
         call S24_interp (period(count1),period(count2),c19(count1),c19(count2), 
     +                 specT,c19T,iflag)
         call S24_interp (period(count1),period(count2),c20(count1),c20(count2), 
     +                 specT,c20T,iflag)
         call S24_interp (period(count1),period(count2),c21(count1),c21(count2), 
     +                 specT,c21T,iflag)
         call S24_interp (period(count1),period(count2),c22(count1),c22(count2), 
     +                 specT,c22T,iflag)
         call S24_interp (period(count1),period(count2),c23(count1),c23(count2), 
     +                 specT,c23T,iflag)
         call S24_interp (period(count1),period(count2),c24(count1),c24(count2), 
     +                 specT,c24T,iflag)
         call S24_interp (period(count1),period(count2),c25(count1),c25(count2), 
     +                 specT,c25T,iflag)
         call S24_interp (period(count1),period(count2),c26(count1),c26(count2), 
     +                 specT,c26T,iflag)
         call S24_interp (period(count1),period(count2),c27(count1),c27(count2), 
     +                 specT,c27T,iflag)
         call S24_interp (period(count1),period(count2),taucr1(count1),taucr1(count2), 
     +                 specT,taucr1T,iflag)
         call S24_interp (period(count1),period(count2),taucr2(count1),taucr2(count2), 
     +                 specT,taucr2T,iflag)
         call S24_interp (period(count1),period(count2),tausb1(count1),tausb1(count2), 
     +                 specT,tausb1T,iflag)
         call S24_interp (period(count1),period(count2),tausb2(count1),tausb2(count2), 
     +                 specT,tausb2T,iflag)
         call S24_interp (period(count1),period(count2),phisscr1(count1),phisscr1(count2), 
     +                 specT,phisscr1T,iflag)
         call S24_interp (period(count1),period(count2),phisscr2(count1),phisscr2(count2), 
     +                 specT,phisscr2T,iflag)
         call S24_interp (period(count1),period(count2),phisssb1(count1),phisssb1(count2), 
     +                 specT,phisssb1T,iflag)
         call S24_interp (period(count1),period(count2),phisssb2(count1),phisssb2(count2), 
     +                 specT,phisssb2T,iflag)
         call S24_interp (period(count1),period(count2),phis2s(count1),phis2s(count2), 
     +                specT,phis2sT,iflag)
         call S24_interp (period(count1),period(count2),c28(count1),c28(count2), 
     +                 specT,c28T,iflag)
         call S24_interp (period(count1),period(count2),c29(count1),c29(count2), 
     +                 specT,c29T,iflag)
         call S24_interp (period(count1),period(count2),c30(count1),c30(count2), 
     +                 specT,c30T,iflag)

  
 1011 period1 = specT

C      h = 10.0
      n = 2.0
      Mc = 7.1
      Mref = 5.5
      Mmax = 8
      Rrupref = 0.0
      Vs30ref = 760.0
    
C     Set the reference spectrum.                
c     sourcetype = 0 for crustal
c                  1 for Subduction 
c     Vs30_class = 0 for estimated
c     Vs30_class = 1 for measured 

      Fcr=0
      Fsb=0
      Fcrss = 0
      Fcrno = 0
      Fcrro = 0
      Fsbintra = 0
      Fsbinter = 0
      Fas = 0
      Fkuo17 = 0
      Fks17 = 0
      Frf = 0
      Fmanila = 0
      C11flag = 0
      C23flag = 0
      C29flag = 0
      C13flag = 0 
      C10flag = 0
   
      if (sourcetype .eq. 0.0 ) then
       Fcr = 1
       Zref = 0
         if(ftype .gt. 0) then
              Fcrro = 1
           elseif(ftype .lt. 0) then
              Fcrno = 1
           else
              Fcrss = 1
         endif
      elseif (sourcetype .eq. 1.0 ) then
        Fsb = 1
         if(ftype .eq. 0) then
              Fsbinter = 1
              Zref = 0
         elseif(ftype .eq. 1) then
              Fsbintra = 1
              Zref = 35
         endif
      endif
   
C     Add aftershock factor 
      if (msasflag .eq. 1) then
           Fas = 1
      endif 

C     choose Site ref by Vs30 class
        if (vs30_class .eq. 0 ) then
         Fks17 = 1
        elseif (vs30_class .eq. 1) then
         Fkuo17 = 1
        endif

      lnYref = c1T*Fcrro + c2T*Fcrss + c3T*Fcrno + c4T*Fsbinter + c5T*Fsbintra +
     &         c6T*Fas + c7T*Fmanila + c26T*Fkuo17 + c27T*Fks17 + c28T*Frf

C     Set Source scaling term 
     
      if(mag .LE. 5 ) then  
       C11flag=1
      endif
      if(mag .GE. Mc ) then  
       C29flag=1
      endif
      if(mag .GE. 7.6 ) then  
       C10flag=1
      endif   
      if(mag .LE. 6 ) then  
       C13flag=1
      endif   
      if (sourcetype .eq. 0.0 ) then
        Smag = c8T*(mag - Mref) + c10T*(mag - Mref)**2 
     1         - c10T*(mag-7.6)**2*C10flag + c11T*(5.0-mag)*C11flag  
      elseif (sourcetype .eq. 1.0 ) then
        Smag = c9T*(mag - Mref) + c29T*Fsbinter*(Mag-Mc)*c29flag + c30T*Fsbintra*(Mag-Mc)*c29flag 
     1        + c12T*(5.0-mag)*C11flag + c13T*(6.0-mag)*C13flag
      endif

      Sztor = c14T * Fcr *(Ztor-Zref) + c15T * Fsbinter * (Ztor-Zref) + c16T * Fsbintra * (Ztor-Zref)    
      Ssource = Smag + Sztor

C     Set Path scaling term

      h = 10.0*Fcr +10.0*Fsbinter*exp(0.3*(mag-7.1)*C29flag) + 10.0*Fsbintra*exp(0.2*(mag-7.1)*C29flag)
   
      if (sourcetype .eq. 0.0 ) then
          Sgeom = (c17T + c19T*(min(mag,Mmax)- Mref )) * alog(SQRT(dist**2 + h**2)/SQRT(Rrupref**2 + h**2))
      elseif (sourcetype .eq. 1.0 ) then
          Sgeom = (c18T + c20T*(min(mag, Mc )- Mref )) * alog(SQRT(dist**2 + h**2)/SQRT(Rrupref**2 + h**2))
      endif

      Sanel = c21T*Fcr*(dist-Rrupref) + c22T*Fsb*(dist-Rrupref)
      Spath = Sgeom + Sanel 
    
C     Set Site scaling term 
    
      Z10ref = exp((-4.08/2.0)*alog((vs**2.0+355.4**2.0)/(1750**2.0+355.4**2.0)))
      if (Z10 .gt. 0.) then
        Ssitelin = c24T * alog(vs/vs30ref) + c25T*alog(Z10*1000/Z10ref)
      else
        Ssitelin = c24T * alog(vs/vs30ref)
      endif

      if(vs .LT. vs30ref ) then  
           C23flag=1
      endif
     
      Ssitenon = c23T * C23flag * (-1.5*alog(vs/vs30ref)-alog(SA1180+2.4)+alog(SA1180+2.4*(vs/vs30ref)**1.5))  
      Ssite = Ssitenon + Ssitelin

      lnSa =  lnYref + Ssource + Spath + Ssite                                        
   
C      write(*,*) "lnYref = ", lnYref
C      write(*,*) "Ssource = ", Ssource
C      write(*,*) "--Smag = ", Smag
C      write(*,*) "--Sztor = ", Sztor
C      write(*,*) "Spath = ", Spath
C      write(*,*) "--Sgeom = ", Sgeom
C      write(*,*) "--Sanel = ", Sanel
C      write(*,*) "Ssite = ", Ssite
C      write(*,*) "--Ssitelin = ", Ssitelin
C      write(*,*) "--Ssitenon = ", Ssitenon
C      write(*,*) "lnSa = ", lnSa
C      write(*,*) "Sa = ", exp(lnSa)

   
C     Set the event-specific residual term
 
      fm = 0.5*(min(6.5, max(4.5, mag))-4.5)
      
      taucr = taucr1T + (taucr2T - taucr1T)*fm
      tausb = tausb1T + (tausb2T - tausb1T)*fm 
      
      tau = taucr*Fcr + tausb*Fsb   
   
C     Set Site-specific residual term

      

C     Set Recoed-specific residual term

      phisscr = phisscr1T + (phisscr2T -phisscr1T)*fm
      phisssb = phisssb1T + (phisssb2T -phisssb1T)*fm

      phiss = phisscr*Fcr + phisssb*Fsb
      
      phi=(phis2sT**2+phiss**2)**0.5
      sigma=(tau**2+phi**2)**0.5
      sigmass=(tau**2+phiss**2)**0.5



c       write(*,*) "Y(gal) = ", exp(lnSa)

      return                                                                    
      end       
          
c ------------------------------------------------------------------            
C *** Phung2018 Crust and Subduction Model- Horizontal ***********
c ------------------------------------------------------------------            

      Subroutine S04_PhungCrust2018 ( m, Rrup, Rbjf, specT, period2, lnY, sigma, iflag, 
     1                     vs, Delta, DTor, Ftype, depthvs10, vs30_class,
     2                       regionflag, phi, tau, HWflag, Rx )

      implicit none
      
      integer MAXPER, i, nPer
      parameter (MAXPER=25)
      REAL Period(MAXPER), C1(MAXPER), C1a(MAXPER), C1b(MAXPER), C1c(MAXPER), C1d(MAXPER)
      REAL cn(MAXPER), cm(MAXPER), c3(MAXPER), c5(MAXPER), c6(MAXPER)
      REAL c7(MAXPER), C7b(MAXPER), C11(MAXPER), C11b(MAXPER), CHM(MAXPER)
      REAL phi1(MAXPER), phi2(MAXPER), phi3(MAXPER), phi4(MAXPER), phi5(MAXPER)
      REAL sigma1inf(MAXPER), sigma2inf(MAXPER)
      REAL tau1(MAXPER), tau2(MAXPER), sigma1(MAXPER), sigma2(MAXPER),tau0(MAXPER)
      REAL sigma3(MAXPER), c8(MAXPER), c8b(MAXPER)
      REAL cg1CA(MAXPER), cg1JP(MAXPER), cg10(MAXPER), dp(MAXPER)
      Real tauT1(MAXPER), phiT1(MAXPER), c9(MAXPER), c9a(MAXPER), c9b(MAXPER)
      real phiss(MAXPER), phis2s(MAXPER) 
      real phiss1M(MAXPER), phiss2M(MAXPER)
      real phi1CA(MAXPER), phi1JP(MAXPER), phi10(MAXPER), cg1(MAXPER), cg2(MAXPER), cg3(MAXPER)
      real vs, phi6
      real Finferred, Fmeasured, dDPP
      REAL c1T, c1aT, c1bT, c1cT, c1dT,cnT, cmT, c5T, c6T, c3T, c9T, c9aT, c9bT
      REAL phi1T, phi2T, phi3T, phi4T, sigma3T, sigma1T, sigma2T
      REAL phi5T, tau1T, tau2T, tauT, tau0T, phi6T
      real c7T, c7bT, c11T, c11bT, cHMT
      real cg1CAT, cg1JPT, cg10T, sigma1infT, sigma2infT
      real phi1CAT, phi1JPT, phi10T
      REAL c2, c4, c4a, cRB, pi, d2r, term14, term15, term16, NL0
      REAL term1, term2, term3, term5, term4, term6, term8, term9, term10, term12, term11
      REAL phissT, phis2sT, sigma1meaT, sigma2meaT, phiss1MT, phiss2MT
      real CNS, cosdelta, psa_ref, psa, cg1T, cg2T, cg3T, dpT
      integer iflag, count1, count2, vs30_class, regionflag, msasflag, HWflag
      REAL M, RRUP, DTOR, Delta, specT, sigma, Ftype, Rbjf, Rx
      REAL period2, lnY, F_RV, F_NM, tau, phi, rkdepth
      real c8T, c8a, c8bT, fd, lnpsa_ref, lnpsa, sa
      real sigmaNL0, F_Measured, F_Inferred, mz_TOR, deltaZ_TOR, coshM
      real period1 ,Ez1, term7, deltaZ1, depthvs10, delc5

C     Mainshock and Aftershocks included based on MSASFlag
C         0 = Mainshocks
C         1 = Aftershocks
C
C     regionflag  
C           = 1 for Taiwan
C           = 0 for global
C
C     vs30_class     Note
C     -------------------------
C      0         estimated
C      1         measured
C

      data period / 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 
     1              0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10 /
      data c1 / -1.654642, -1.574253, -1.520116, -1.424583, -1.301324, -1.164997, -0.847910, -0.672079, 
     1          -0.593532, -0.578521, -0.608837, -0.672703, -0.799766, -0.975575, -1.287988, -1.547417,  
     1          -2.106085, -2.510412, -3.140080, -3.403277, -3.730282, -3.972268, -4.045129, -4.301415, -4.423677 /
      data c3 / 1.543196, 1.513556, 1.476841, 1.434544, 1.386104, 1.321224, 1.249835, 1.314168, 1.396081,  
     1          1.483981, 1.526734, 1.607787, 1.761750, 1.914659, 2.125763, 2.275589, 2.598148, 2.759899,  
     1          2.885389, 2.982531, 3.039484, 3.075429, 3.098058, 3.104144, 3.087479 /
      data cn / 12.148668, 12.148668, 12.248034, 12.533784, 12.991897, 13.650754, 15.714475, 16.772622,  
     1          16.775630, 16.186798, 15.843144, 15.014671, 12.696431, 10.449811, 6.802217, 4.410694,  
     1          3.406400, 3.161200, 2.807800, 2.463100, 2.211100, 1.966800, 1.667100, 1.573700, 1.526500 /
      data cm / 4.883828, 5.386356, 5.455551, 5.509632, 5.512664, 5.533465, 5.566513, 5.513233, 5.483484,  
     1          5.487482, 5.474835, 5.435991, 5.412552, 5.379773, 5.345804, 5.320560, 5.369402, 5.445281,  
     1          5.624212, 5.751105, 6.021403, 6.156573, 6.311369, 6.588973, 6.923637 /
      data c5 / 6.455100, 6.455100, 6.455100, 6.455100, 6.455100, 6.455100, 6.455100, 6.830500, 7.133300,  
     1          7.362100, 7.436500, 7.497200, 7.541600, 7.560000, 7.573500, 7.577800, 7.580800, 7.581400,  
     1          7.581700, 7.581800, 7.581800, 7.581800, 7.581800, 7.581800, 7.581800 /
      data c6 / 0.490800, 0.490800, 0.492500, 0.499200, 0.503700, 0.504800, 0.504800, 0.504800, 0.504800,  
     1          0.504500, 0.503600, 0.501600, 0.497100, 0.491900, 0.480700, 0.470700, 0.457500, 0.452200,  
     1          0.450100, 0.450000, 0.450000, 0.450000, 0.450000, 0.450000, 0.450000 /
      data cHM  / 3.095600, 3.095600, 3.096300, 3.097400, 3.098800, 3.101100, 3.109400, 3.238100, 3.340700,  
     1          3.430000, 3.468800, 3.514600, 3.574600, 3.623200, 3.694500, 3.740100, 3.794100, 3.814400,  
     1          3.828400, 3.833000, 3.836100, 3.836900, 3.837600, 3.838000, 3.838000 /
      data c7 / 0.008581, 0.008035, 0.007593, 0.007250, 0.007006, 0.006860, 0.007008, 0.007247, 0.007456,  
     1          0.007703, 0.007799, 0.007823, 0.008071, 0.008396, 0.009275, 0.010166, 0.012793, 0.013762,  
     1          0.013900, 0.012559, 0.009184, 0.004797, 0.001068, -0.004234, -0.006203 /
      data c7b / 0.020232, 0.021034, 0.021639, 0.022052, 0.022284, 0.022341, 0.021712, 0.020031, 0.018585,  
     1          0.016544, 0.015413, 0.014411, 0.013238, 0.011958, 0.009469, 0.005800, -0.003683, -0.008131,  
     1          -0.010287, -0.008563, -0.003059, 0.003920, 0.013064, 0.027920, 0.041953 /
      data c1a / 0.139286, 0.134926, 0.130854, 0.127072, 0.123583, 0.120391, 0.113706, 0.108868, 0.106295,  
     1          0.104509, 0.104665, 0.107605, 0.118083, 0.130738, 0.153830, 0.179466, 0.158850, 0.145609,  
     1          0.137311, 0.117810, 0.076184, 0.038194, 0.038194, 0.038194, 0.038194 /
      data c1c / 0.126967, 0.136204, 0.142178, 0.146345, 0.156255, 0.171596, 0.177475, 0.167603, 0.155744,  
     1          0.126802, 0.104491, 0.075057, 0.043282, 0.026574, -0.014114, -0.064611, 0.035510, 0.051493,  
     1          0.076539, 0.093041, 0.115673, 0.147164, 0.147164, 0.147164, 0.147164 /
      data c1b / 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  
     1          0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,  
     1          0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000 /
      data c1d / -0.043679, -0.043261, -0.043154, -0.043339, -0.043803, -0.044544, -0.047735, -0.053268,  
     1          -0.059612, -0.069063, -0.074707, -0.093710, -0.130636, -0.167306, -0.238958, -0.292486,  
     1          -0.391164, -0.459741, -0.531023, -0.531023, -0.531023, -0.531023, -0.531023, -0.531023, -0.531023 /
      data c11 / 0.016453, -0.108037, -0.102072, -0.104638, -0.105159, -0.096947, -0.079174, -0.120807,  
     1          -0.127655, -0.123958, -0.120235, -0.128555, -0.104990, -0.125335, -0.131459, -0.102607,  
     1          -0.072843, -0.072287, -0.143270, -0.171096, -0.269172, -0.321537, -0.344322, -0.379467, -0.478011 /
      data c11b  / -0.168281, 0.195952, 0.181778, 0.163170, 0.142063, 0.098054, 0.046297, 0.173997, 0.209295,  
     1          0.217340, 0.218819, 0.262936, 0.231024, 0.270344, 0.306056, 0.272617, 0.265158, 0.303895,  
     1          0.443286, 0.520454, 0.817527, 1.015932, 0.892205, 0.864364, 1.443598 /
      data cg1 / -0.007288, -0.007609, -0.007926, -0.008239, -0.008547, -0.008850, -0.009598, -0.009985,  
     1          -0.010161, -0.010176, -0.009969, -0.009506, -0.008808, -0.008145, -0.006990, -0.006140,  
     1          -0.004769, -0.004051, -0.003342, -0.002973, -0.002487, -0.002123, -0.001764, -0.001079, -0.000742 /
      data cg2 / -0.006969, -0.007127, -0.007249, -0.007328, -0.007362, -0.007361, -0.007052, -0.005719,  
     1          -0.004365, -0.002650, -0.002000, -0.001255, -0.000750, -0.000447, -0.000247, -0.000417,  
     1          -0.001131, -0.001741, -0.002428, -0.002706, -0.004107, -0.005776, -0.007748, -0.009142, -0.012633 /
      data cg3 / 4.222069, 4.225635, 4.230342, 4.236182, 4.250189, 4.303123, 4.446127, 4.610835, 4.723497,  
     1          4.878141, 4.981707, 5.066411, 5.219865, 5.328220, 5.201762, 5.187932, 4.877209, 4.639751,  
     1          4.571204, 4.425117, 3.621904, 3.486264, 3.277906, 3.074948, 3.074948 /
      data dp / -6.785205, -6.750648, -6.716179, -6.681799, -6.647507, -6.613303, -6.528180, -6.443608,  
     1          -6.376345, -6.276109, -6.209723, -6.110798, -5.947666, -5.786703, -5.471253, -5.164376,  
     1          -4.434204, -3.755966, -2.550267, -1.536848, -0.052838, 0.000000, 0.000000, 0.000000, 0.000000 /
      data phi1 / -0.516116, -0.516935, -0.508693, -0.495823, -0.476815, -0.454879, -0.430231, -0.459892,  
     1          -0.479365, -0.509853, -0.531335, -0.551465, -0.580586, -0.626209, -0.677679, -0.702586,  
     1          -0.819168, -0.897580, -1.090857, -1.042582, -0.991035, -0.955669, -0.915233, -0.802877, -0.699613 /
      data phi2 / -0.141700, -0.141700, -0.136400, -0.140300, -0.159100, -0.186200, -0.253800, -0.294300,  
     1          -0.307700, -0.311300, -0.306200, -0.292700, -0.266200, -0.240500, -0.197500, -0.163300,  
     1          -0.102800, -0.069900, -0.042500, -0.030200, -0.012900, -0.001600, 0.000000, 0.000000, 0.000000 /
      data phi3 / -0.007010, -0.007010, -0.007279, -0.007354, -0.006977, -0.006467, -0.005734, -0.005604,  
     1          -0.005696, -0.005845, -0.005959, -0.006141, -0.006439, -0.006704, -0.007125, -0.007435,  
     1          -0.008120, -0.008444, -0.007707, -0.004792, -0.001828, -0.001523, -0.001440, -0.001369, -0.001361 /
      data phi4 / 0.102151, 0.102151, 0.108360, 0.119888, 0.133641, 0.148927, 0.190596, 0.230662, 0.253169,  
     1          0.266468, 0.265060, 0.255253, 0.231541, 0.207277, 0.165464, 0.133828, 0.085153, 0.058595,  
     1          0.031787, 0.019716, 0.009643, 0.005379, 0.003223, 0.001134, 0.000515 /
      data phi5 / 0.131065, 0.131615, 0.131659, 0.138990, 0.149910, 0.167293, 0.199578, 0.230268, 0.221342,  
     1          0.193848, 0.168562, 0.145600, 0.127889, 0.116581, 0.087546, 0.055917, 0.016064, 0.004084,  
     1          0.083716, 0.151694, 0.229009, 0.223222, 0.217675, 0.206050, 0.175716 /
      data c9 / 0.922800, 0.922800, 0.929600, 0.939600, 0.966100, 0.979400, 1.026000, 1.017700, 1.000800,  
     1          0.980100, 0.965200, 0.945900, 0.919600, 0.882900, 0.830200, 0.788400, 0.675400, 0.619600,  
     1          0.510100, 0.391700, 0.124400, 0.008600, 0.000000, 0.000000, 0.000000 /
      data c9a / 0.120200, 0.120200, 0.121700, 0.119400, 0.116600, 0.117600, 0.117100, 0.114600, 0.112800,  
     1          0.110600, 0.115000, 0.120800, 0.120800, 0.117500, 0.106000, 0.106100, 0.100000, 0.100000,  
     1          0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000 /
      data c9b / 6.860700, 6.860700, 6.869700, 6.911300, 7.027100, 7.095900, 7.329800, 7.258800, 7.237200,  
     1          7.210900, 7.249100, 7.298800, 7.369100, 6.878900, 6.533400, 6.526000, 6.500000, 6.500000,  
     1          6.500000, 6.500000, 6.500000, 6.500000, 6.500000, 6.500000, 6.500000 /
      data c8 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0991, 0.1982,  
     1          0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154 /
      data c8b / 0.483300, 0.483300, 1.214400, 1.642100, 1.945600, 2.181000, 2.608700, 2.912200, 3.104500,  
     1          3.339900, 3.471900, 3.643400, 3.878700, 4.071100, 4.374500, 4.609900, 5.037600, 5.341100,  
     1          5.768800, 6.072300, 6.500000, 6.803500, 7.038900, 7.466600, 7.770000 /
      data cg1CA / -0.006490, -0.006470, -0.006536, -0.007179, -0.007676, -0.008257, -0.009532, -0.010166,  
     1          -0.010054, -0.009805, -0.009030, -0.008466, -0.007534, -0.006392, -0.005351, -0.004422,  
     1          -0.002395, -0.002620, -0.002223, -0.000956, -0.000992, -0.001407, -0.000876, -0.000876, -0.000876 /
      data phi1CA / -0.577526, -0.576647, -0.567620, -0.540828, -0.525764, -0.504658, -0.472272, -0.515688,  
     1          -0.536177, -0.547802, -0.604585, -0.639294, -0.705116, -0.748564, -0.799633, -0.879300,  
     1          -0.970139, -0.979947, -1.184797, -1.154081, -0.959930, -0.907932, -0.846230, -0.711399, -0.555313 /
      data cg1JP / -0.011222, -0.011180, -0.011131, -0.011552, -0.012203, -0.012473, -0.012366, -0.012286, 
     1           -0.012403, -0.012407, -0.012416, -0.012184, -0.011760, -0.011313, -0.010318, -0.009411,  
     1          -0.007265, -0.006083, -0.006874, -0.007552, -0.007297, -0.005609, -0.004063, -0.001333, -0.000396 /
      data phi1JP / -0.443000, -0.443292, -0.429683, -0.387328, -0.310552, -0.242082, -0.117036, -0.214463,  
     1          -0.310582, -0.439369, -0.523364, -0.632758, -0.769013, -0.855314, -0.912410, -0.965488,  
     1          -1.037790, -1.030000, -1.118098, -0.960000, -0.762240, -0.659832, -0.587533, -0.407331, -0.364009 /
      data cg10 / -0.005243, -0.005818, -0.006354, -0.006849, -0.007301, -0.007709, -0.008550, -0.009063,  
     1          -0.009177, -0.008970, -0.008694, -0.008187, -0.007169, -0.006289, -0.005024, -0.004042,  
     1          -0.002378, -0.001848, -0.001019, -0.001042, -0.001134, -0.001182, -0.001085, -0.001085, -0.001085 /
      data phi10 / -0.537767, -0.537380, -0.527902, -0.503857, -0.494849, -0.477957, -0.453645, -0.485480,  
     1          -0.502936, -0.515412, -0.561167, -0.601685, -0.677544, -0.732538, -0.771490, -0.814426,  
     1          -0.888440, -0.876327, -1.080612, -1.060215, -0.904488, -0.873702, -0.850298, -0.749727, -0.598953 /
      data tau1 / 0.400000, 0.400000, 0.402600, 0.406300, 0.409500, 0.412400, 0.417900, 0.421900, 0.424400,  
     1          0.427500, 0.429200, 0.431300, 0.434100, 0.436300, 0.439600, 0.441900, 0.445900, 0.448400,  
     1          0.451500, 0.453400, 0.455800, 0.457400, 0.458400, 0.460100, 0.461200 /
      data tau2 / 0.260000, 0.260000, 0.263700, 0.268900, 0.273600, 0.277700, 0.285500, 0.291300, 0.294900,  
     1          0.299300, 0.301700, 0.304700, 0.308700, 0.311900, 0.316500, 0.319900, 0.325500, 0.329100,  
     1          0.333500, 0.336300, 0.339800, 0.341900, 0.343500, 0.345900, 0.347400 /
      data sigma1 / 0.491200, 0.491200, 0.490400, 0.498800, 0.504900, 0.509600, 0.517900, 0.523600, 0.527000,  
     1          0.530800, 0.532800, 0.535100, 0.537700, 0.539500, 0.542200, 0.543300, 0.529400, 0.510500,  
     1          0.478300, 0.468100, 0.461700, 0.457100, 0.453500, 0.447100, 0.442600 /
      data sigma2 / 0.376200, 0.376200, 0.376200, 0.384900, 0.391000, 0.395700, 0.404300, 0.410400, 0.414300,  
     1          0.419100, 0.421700, 0.425200, 0.429900, 0.433800, 0.439900, 0.444600, 0.453300, 0.459400,  
     1          0.468000, 0.468100, 0.461700, 0.457100, 0.453500, 0.447100, 0.442600 /
      data sigma3 / 0.800000, 0.800000, 0.800000, 0.800000, 0.800000, 0.800000, 0.800000, 0.800000, 0.800000,  
     1          0.800000, 0.800000, 0.800000, 0.799900, 0.799700, 0.798800, 0.796600, 0.779200, 0.750400,  
     1          0.713600, 0.703500, 0.700600, 0.700100, 0.700000, 0.700000, 0.700000 /
      data tau0  / 0.372401, 0.376514, 0.376953, 0.388775, 0.403902, 0.421147, 0.443051, 0.439374, 0.424113,  
     1          0.404915, 0.393373, 0.375874, 0.370636, 0.374764, 0.405444, 0.432212, 0.442190, 0.473377,  
     1          0.474720, 0.468637, 0.454880, 0.483916, 0.507364, 0.575778, 0.562976 /
      data phiss / 0.465403, 0.463342, 0.463153, 0.470111, 0.477682, 0.482355, 0.483524, 0.478191, 0.470999,  
     1          0.465492, 0.465018, 0.465715, 0.468827, 0.475041, 0.483249, 0.494461, 0.493867, 0.487026,  
     1          0.479004, 0.471925, 0.458172, 0.456369, 0.453753, 0.445327, 0.443534 /
      data phis2s / 0.317008, 0.317230, 0.317076, 0.324896, 0.337890, 0.355358, 0.389659, 0.400733, 0.396010,  
     1          0.376335, 0.367016, 0.353613, 0.337756, 0.335367, 0.333318, 0.334606, 0.340700, 0.348831,  
     1          0.379082, 0.390649, 0.389996, 0.387436, 0.381651, 0.374743, 0.346654 /


C Find the requested spectral period and corresponding coefficients
      nPer = 25
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T  =   c1(1)
         c3T  =   c3(1)
         cnT  =   cn(1)
         cmT  =   cm(1)
         c5T  =   c5(1)
         c6T  =   c6(1)
         cHMT  =   cHM(1)
         c7T  =   c7(1)
         c7bT  =   c7b(1)
         c1aT  =   c1a(1)
         c1cT  =   c1c(1)
         c1bT  =   c1b(1)
         c1dT  =   c1d(1)
         c11T  =   c11(1)
         c11bT  =   c11b(1)
         cg1T  =   cg1(1)
         cg2T  =   cg2(1)
         cg3T  =   cg3(1)
         dpT   =   dp(1)
         cg1CAT  =   cg1CA(1)
         phi1CAT  =   phi1CA(1)
         cg1JPT  =   cg1JP(1)
         phi1JPT  =   phi1JP(1)
         cg10T  =   cg10(1)
         phi10T  =   phi10(1)

         phi1T  =   phi1(1)
         phi2T  =   phi2(1)
         phi3T  =   phi3(1)
         phi4T  =   phi4(1)
         phi5T  =   phi5(1)
         c9T   =   c9(1)
         c9aT  =   c9a(1)
         c9bT  =   c9b(1)
         c8T   =   c8(1)
         c8bT  =   c8b(1)

         tau1T  =   tau1(1)
         tau2T  =   tau2(1)
         sigma1T  =   sigma1(1)
         sigma2T  =   sigma2(1)
         sigma3T  =   sigma3(1)
         phissT   =   phiss(1)
         phis2sT  =   phis2s(1)
         tau0T  =   tau0(1)

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
      write (*,*) 'Phung et al. 2018 horizontal'
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
 1020       call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +                   specT,c3T,iflag)
            call S24_interp (period(count1),period(count2),cn(count1),cn(count2),
     +                   specT,cnT,iflag)
            call S24_interp (period(count1),period(count2),cm(count1),cm(count2),
     +                   specT,cmT,iflag)
            call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +                   specT,c5T,iflag)
            call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +                   specT,c6T,iflag)
            call S24_interp (period(count1),period(count2),cHM(count1),cHM(count2),
     +                   specT,cHMT,iflag)
            call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +                   specT,c7T,iflag)
            call S24_interp (period(count1),period(count2),c7b(count1),c7b(count2),
     +                   specT,c7bT,iflag)

            call S24_interp (period(count1),period(count2),c1a(count1),c1a(count2),
     +                   specT,c1aT,iflag)
            call S24_interp (period(count1),period(count2),c1b(count1),c1b(count2),
     +                   specT,c1bT,iflag)
            call S24_interp (period(count1),period(count2),c1c(count1),c1c(count2),
     +                   specT,c1cT,iflag)
            call S24_interp (period(count1),period(count2),c1d(count1),c1d(count2),
     +                   specT,c1dT,iflag)
            call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +                   specT,c11T,iflag)
            call S24_interp (period(count1),period(count2),c11b(count1),c11b(count2),
     +                   specT,c11bT,iflag)

            call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +                   specT,c8T,iflag)
            call S24_interp (period(count1),period(count2),c8b(count1),c8b(count2),
     +                   specT,c8bT,iflag)
            call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +                   specT,c9T,iflag)
            call S24_interp (period(count1),period(count2),c9a(count1),c9a(count2),
     +                   specT,c9aT,iflag)
            call S24_interp (period(count1),period(count2),c9b(count1),c9b(count2),
     +                   specT,c9bT,iflag)

  
            call S24_interp (period(count1),period(count2),cg1(count1),cg1(count2),
     +                   specT,cg1T,iflag)
            call S24_interp (period(count1),period(count2),cg2(count1),cg2(count2),
     +                   specT,cg2T,iflag)
            call S24_interp (period(count1),period(count2),cg3(count1),cg3(count2),
     +                   specT,cg3T,iflag)

            call S24_interp (period(count1),period(count2),dp(count1),dp(count2),
     +                   specT,dpT,iflag)

            call S24_interp (period(count1),period(count2),cg1CA(count1),cg1CA(count2),
     +                   specT,cg1CAT,iflag)
             call S24_interp (period(count1),period(count2),cg1JP(count1),cg1JP(count2),
     +                   specT,cg1JPT,iflag)
             call S24_interp (period(count1),period(count2),cg10(count1),cg10(count2),
     +                   specT,cg10T,iflag)
 
     
            call S24_interp (period(count1),period(count2),phi1(count1),phi1(count2),
     +                   specT,phi1T,iflag)
            call S24_interp (period(count1),period(count2),phi2(count1),phi2(count2),
     +                   specT,phi2T,iflag)
            call S24_interp (period(count1),period(count2),phi3(count1),phi3(count2),
     +                   specT,phi3T,iflag)
            call S24_interp (period(count1),period(count2),phi4(count1),phi4(count2),
     +                   specT,phi4T,iflag)
            call S24_interp (period(count1),period(count2),phi5(count1),phi5(count2),
     +                   specT,phi5T,iflag)

            call S24_interp (period(count1),period(count2),phi1CA(count1),phi1CA(count2),
     +                   specT,phi1CAT,iflag)
            call S24_interp (period(count1),period(count2),phi1JP(count1),phi1JP(count2),
     +                   specT,phi1JPT,iflag)

  
            call S24_interp (period(count1),period(count2),phiss(count1),phiss(count2),
     +                   specT,phissT,iflag)
            call S24_interp (period(count1),period(count2),phis2s(count1),phis2s(count2),
     +                   specT,phis2sT,iflag)
            call S24_interp (period(count1),period(count2),tau1(count1),tau1(count2),
     +                   specT,tau1T,iflag)
            call S24_interp (period(count1),period(count2),tau2(count1),tau2(count2),
     +                   specT,tau2T,iflag)
            call S24_interp (period(count1),period(count2),tau0(count1),tau0(count2),
     +                   specT,tau0T,iflag)

            call S24_interp (period(count1),period(count2),sigma1(count1),sigma1(count2),
     +                   specT,sigma1T,iflag)
            call S24_interp (period(count1),period(count2),sigma2(count1),sigma2(count2),
     +                   specT,sigma2T,iflag)
            call S24_interp (period(count1),period(count2),sigma2(count1),sigma2(count2),
     +                   specT,sigma2T,iflag)


 1011 period1 = specT                                                                                                              

c     Set the fault mechanism term.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
         if (ftype .eq. -1) then
            F_RV = 0.0
            F_NM = 1.0
         elseif (ftype .ge. 0.5) then
            F_RV = 1.0
            F_NM = 0.0
         else
            F_RV = 0.0
            F_NM = 0.0
         endif

C     Constant terms
        c2 = 1.06
        c4 = -2.1
        c4a = -0.5
        cRB = 50.0
        phi6 = 300
        c8a = 0.2695
  
      if(regionflag .eq. 1) then
C for Taiwan
        cg1T = cg1T
        phi1T = phi1T

      elseif(regionflag .eq. 2) then
C for California
        cg1T = cg1CAT
        phi1T = phi1CAT

      elseif(regionflag .eq. 3) then
C for others
        cg1T = cg10T
        phi1T = phi10T

        elseif(regionflag .eq. 4) then
C for Japan
        cg1T = cg1JPT
        phi1T = phi1JPT
  
      endif

C     Current code set for Measured Vs30 values (i.e., Vs30class=1)
      if (vs30_class .eq. 0) then
         Fmeasured = 0.0
         FInferred = 1.0

      elseif (vs30_class .eq. 1) then      
         Fmeasured = 1.0
         FInferred = 0.0

      endif       
  
c Center Z_TOR on the Z_TOR-M relation
        if (F_RV.EQ.1) then

            mZ_TOR = max(3.5384-2.60 * max(M-5.8530,0.0),0.0)
            mZ_TOR = mZ_TOR * mZ_TOR

        else
            mZ_TOR = max(2.7482-1.7639*max(M-5.5210,0.0),0.0)
            mZ_TOR = mZ_TOR * mZ_TOR
        endif
  
c        if (Z_TOR .EQ. -999) Z_TOR = mZ_TOR
        deltaZ_TOR = Dtor - mZ_TOR
   
c Reference motion  
  
        pi = atan(1.0)*4.0
        d2r = pi/180.0
        term1 = c1T
  
c Magnitude scaling
        term6 = c2 * (M-6.0) 
        term7 = (c2-c3T)/cnT * alog(1.0 + exp(cnT*(cMT-M)))  

c Near-field magnitude and distance scaling
       if (Dtor > 20 .and. M < 7) then
      delc5 = dpT * max(DTor/50.0-20.0/50.0,0.0)
      else
      delc5 =0.0
       endif
        
        CNS = (c5T + delc5) * cosh(c6T * max((M-cHMT),0.0))  
          
        term8 = c4 * alog(Rrup + CNS)                 

c Distance scaling at large distance
        term9 = (c4a-c4) * alog( sqrt(Rrup*Rrup+cRB*cRB) )  
        term10 = (cg1T + cg2T/cosh(max((M-cg3T),0.0)))*Rrup  


c Scaling with other source variables (F_RV, F_NM, deltaZ_TOR, and Dip)
        coshM = cosh(2*max(M-4.5,0.0))
        cosDELTA = cos(DELTA*d2r)
        term2 = (c1aT+c1cT/coshM) * F_RV 
        term3 = (c1bT+c1dT/coshM) * F_NM 
        term4 = (c7T +c7bT/coshM) * deltaZ_TOR 
        term5 = (c11T+c11bT/coshM)* cosDELTA**2   

c HW effect 
        if (HWFlag .eq. 0) then
           term12 = 0.0
        else
         term12 = c9T * HWFlag *(cosDELTA) * (c9aT+(1-c9aT)
     1        *tanh(abs(Rx)/c9bT)) *
     1          (1.0 - sqrt(Rbjf**2+DTor**2)/(Rrup + 1))
        endif

C     Current version of the code sets dDPP=0 (i.e., no directivity)
c Directivity effect
        dDPP = 0.0
        term11 = c8T * exp(-c8a * (M-c8bT)**2) *
     1       max(0.0, 1.0-max(0.0,Rrup-40.0)/30.0) *
     1       min(max(0.0,M-5.5)/0.8, 1.0) * dDPP
       
c Predicted median Sa on reference condition (Vs=1130 m/sec)
        lnpsa_ref = term1+term2+term3+term5+term4+term6+term7+term8+term9+term10+term11+term12
        psa_ref = exp(lnpsa_ref)
  
c Linear soil amplification
        term14 = phi1T * min(alog(Vs/1130.0), 0.0)   

c Nonlinear soil amplification
        term15 = phi2T *
     1      (exp(phi3T*(min(Vs,1130.0)-360.0)) - exp(phi3T*(1130.0-360.0)))*
     1      alog((psa_ref+phi4T)/phi4T)

C Deviation from ln(Vs30) scaling: bedrock depth (Z1) effect.
        Ez1 = exp(-2.63/4.0 * alog((VS**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))
        deltaZ1 = depthvs10*1000.0 - Ez1
C     1    exp(-2.63/4.0 * alog((VS**4.0 + 253.0**4.0)/(2492.0**4.0 + 253.0**4.0)))
 
        if (regionflag .eq. 0) then
            term16 = 0.0
         elseif (regionflag .eq. 1) then
            term16 = phi5T*( 1.0 -exp(-deltaZ1/phi6))
        endif
  
c Sa on soil condition
        lnpsa = lnpsa_ref + term14 + term15 + term16
        sa = exp(lnpsa_ref + term14 + term15 + term16)
        psa = psa_ref * exp(term14 + term15 + term16)

c        write(*,*) "term1 = " , term1
c        write(*,*) "term2 = " , term2
c        write(*,*) "term3 = " , term3
c        write(*,*) "term4 = " , term4
c        write(*,*) "term5 = " , term5
c        write(*,*) "term6 = " , term6
c        write(*,*) "term7 = " , term7
c        write(*,*) "term8 = " , term8
c        write(*,*) "term9 = " , term9
c        write(*,*) "term10 = ", term10
c        write(*,*) "term14 = ", term14
c        write(*,*) "term15 = ", term15
c        write(*,*) "term16 = ", term16
c        write(*,*) "Ez1 = " , Ez1
c        write(*,*) "deltaZ1 = " , deltaZ1
c        write(*,*) "lnpsa_ref = ", lnpsa_ref
c        write(*,*) "lnpsa = ", lnpsa
c        write(*,*) "psa = ", psa
  
C Compute the sigma term
C Variance Model-1 from CY14 with phi1 from Taiwan


       NL0=phi2T*(exp(phi3T*(min(Vs,1130.0)-360.0))-exp(phi3T*(1130.0-360.0)))
     1    *(psa_ref/(psa_ref+phi4T))
  
       sigmaNL0 = (sigma1T+(sigma2T - sigma1T)/1.5*(min(max(M,5.0),6.5)-5.0))*
     1           sqrt((sigma3T*Finferred + 0.7* Fmeasured) + (1.0+NL0)**2.0)

       tau = tau1T +(tau2T-tau1T)/1.5*(min(max(M,5.0),6.5)-5.0)

       sigma = sqrt((1+NL0)**2.0*(tau)**2.0+sigmaNL0**2.0)

C Variance Model-2 for Taiwan

      sigma = sqrt(tau0T**2 + phissT**2 + phis2sT**2)
      phi = sqrt(phissT**2 + phis2sT**2)

C     Convert ground motion to units of gals.
      lnY = lnpsa + 6.89
      period2 = period1

      return
      end 

c ------------------------------------------------------------------            
C *** Adjusted BCHydro model by Phung and Loh ***********
c ------------------------------------------------------------------            

      subroutine S04_PhungSub2018 ( mag, rRup, vs30, Z10, ZTor, lnY, sigma,  
     2                     specT, period2, iflag, regionflag, ftype )

      implicit none
     
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=21)
      real period(MAXPER), a5(MAXPER), a13(MAXPER), Mref(MAXPER), a2(MAXPER), a14(MAXPER), 
     1     dela1(MAXPER), dela4(MAXPER), a6jp(MAXPER), a12jp(MAXPER), a8jp(MAXPER)
      real phisstj(MAXPER),  phis2stj(MAXPER), tau0(MAXPER)

      real a1tw(MAXPER), a4tw(MAXPER),  a7(MAXPER), a6tw(MAXPER), a12tw(MAXPER), a8(MAXPER),
     1     a11(MAXPER), a10(MAXPER),
     1     phisstw(MAXPER), phis2stw(MAXPER), tautw(MAXPER), phitw(MAXPER)

      real sigma, lnSa, pgaRock, vs30, rRup, disthypo, mag 

      real periodT, a5T, a13T, MrefT, a2T, a14T, dela1T, dela4T, a6jpT, a12jpT, a8jpT
      real phisstjT,  phis2stjT, tau0T, a1twT, a4twT,  a7T, a6twT, a12twT, a8T
      real a11T, a10T, phisstwT, phis2stwT, tautwT, phitwT
   
      real Ez1, fz10, fmag, frup, fsite, fztor, fevt
      real period1, a3, Z10, ZTor, a9, d, b12, lnY, Fs, a11si, a11ss, phiss, phis2s, a1, a4,a6,a12
      integer count1, count2, iflag, regionflag
      real n, c, c4, c1, faba, R, depth, specT, tau, phi, ftype, period2


      data period  /0, 0.01, 0.02, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 1, 1.5, 2, 
     1              2.5, 3, 4, 5 /
      data a5 /  0.03849929, 0.04033665, 0.04190178, 0.04509359, 0.04623140, 0.04819708, 0.04325090, 0.03692059, 
     1            0.06597319, 0.06197944, 0.06979644, 0.08783791, 0.09612877, 0.10612877, 0.22744484, 0.16136621,  
     1            0.22767232, 0.27153377, 0.28822087, 0.32589322, 0.30383949 /  
      data a13 /  -0.0256568, -0.0259617, -0.0262528, -0.0270426, -0.0276048, -0.0280794, -0.0287650, -0.0291017,  
     1            -0.0290970, -0.0287552, -0.0269993, -0.0235859, -0.0180673, -0.0150673, -0.0031849, -0.0031849,  
     1            -0.0031849, -0.0031849, -0.0031849, -0.0031849, -0.0031849 /  
      data Mref /  7.68, 7.68, 7.68, 7.71, 7.77, 7.77, 7.78, 7.72, 7.62, 7.54, 7.42, 7.38, 7.36, 7.32, 7.25, 7.25,  
     1            7.25, 7.25, 7.25, 7.25, 7.25 /     
      data a2  /  -1.552846733, -1.554174269, -1.555152194, -1.556049687, -1.554562252, -1.551165488, -1.539140832,  
     1            -1.520764226, -1.489051706, -1.464118878, -1.414761429, -1.383170353, -1.360022278, -1.313716982,  
     1            -1.236841977, -1.100570482, -0.990254902, -0.896093506, -0.818199517, -0.730697376, -0.734817372 /  
      data a14 /  -0.011876681, -0.012409284, -0.016872732, -0.08510905, -0.118005772, -0.171218187, -0.124720279,  
     1            -0.120958201, -0.116255248, -0.077408811, -0.054966213, -0.034173086, -0.06069315, -0.039053473,  
     1            0.017806808, -0.005705423, 0.053155037, 0.068765677, 0.071577687, 0.042405486, 0.054712361 /  
      data dela1 /  1.141899742, 1.152006702, 1.154339185, 1.515303155, 1.904431142, 1.945456526, 1.787100626,  
     1            1.562515125, 1.356740101, 1.206013896, 0.760110718, 0.431629072, 0.214072689, -0.01782956, -0.204991951, 
     1             -0.382378342, -0.352611734, -0.228719047, -0.16534756, 0.010400184, 0.135306871 /  
      data dela4 /  0.328613796, 0.352192886, 0.367677006, 0.452541112, 0.513739193, 0.499522237, 0.45427803,  
     1            0.363869484, 0.314270529, 0.282854636, 0.176621233, 0.077343234, 0.044354701, -0.012346815, -0.083175191, 
     1             -0.226959428, -0.206168741, -0.163340854, -0.154799226, -0.112793712, -0.003105582 /  
      data a6jp /  -0.006794362, -0.006817094, -0.006816137, -0.007285287, -0.007702243, -0.007674043, -0.00782682,  
     1            -0.007547403, -0.007323965, -0.006976346, -0.006143907, -0.005504091, -0.004739568, -0.004285202,  
     1            -0.003957479, -0.003379651, -0.003469497, -0.003409644, -0.003614492, -0.003749858, -0.003243671 /  
      data a12jp /  -0.7516020, -0.7500948, -0.7307185, -0.4831132, -0.3413025, -0.4948081, -0.8669192, -1.0634892,  
     1            -1.1789740, -1.2253631, -1.2073943, -1.1299835, -1.0859780, -1.0233555, -0.9766258, -0.9437327, -0.8880212,  
     1            -0.8546554, -0.7803988, -0.6937169, -0.6499105 /  
      data a8jp /  0.002650526, 0.002810000, 0.002376243, 0.000548389, 0.003562477, 0.006709533, 0.007165476, 0.004341874,  
     1            0.004771685, 0.005747192, 0.005294191, 0.002624524, 0.00028237, -0.001824156, -0.003649632, -0.005143864,  
     1            -0.005623872, -0.004503918, -0.005225851, -0.006465579, -0.004203115 /  
      data tau0 /  0.426469333, 0.424670673, 0.429099403, 0.477262493, 0.516365961, 0.512863781, 0.461620132, 0.441284014,  
     1            0.434871493, 0.417105574, 0.412231943, 0.396422531, 0.403512982, 0.409058414, 0.428087961, 0.440164103, 
     1            0.451402135, 0.461009777, 0.457681279, 0.470193371, 0.461834624 /  
      data phisstj /  0.420489356, 0.420220542, 0.418935068, 0.419356601, 0.41012804, 0.420066336, 0.433030041, 0.446059203,  
     1            0.456165064, 0.459407794, 0.452177558, 0.443446716, 0.4378372, 0.446095737, 0.44128784, 0.42093392,  
     1            0.425797447, 0.419978946, 0.413207757, 0.36936182, 0.349936999 /  
      data phis2stj /  0.364038777, 0.364065617, 0.364403496, 0.417921293, 0.469674634, 0.469076147, 0.437340846, 0.397161802, 
     1             0.384511586, 0.373773022, 0.372643333, 0.369406327, 0.397812841, 0.419424271, 0.408446403, 0.420315802,  
     1            0.416797894, 0.404296599, 0.362727332, 0.354026172, 0.300771328 /  
      data a1tw /  4.481424147, 4.500413862, 4.524812489, 4.684812623, 4.81707322, 4.943259845, 5.09134277, 5.090645572,  
     1            4.96517978, 4.84672774, 4.616816523, 4.435687208, 4.265016233, 3.916493523, 3.163705631, 2.254253444,  
     1            1.320828841, 0.510914193, -0.121327899, -1.029553279, -1.391905014 /   
      data a4tw /  0.441987425, 0.442328411, 0.436082809, 0.363261281, 0.319469336, 0.325896968, 0.350561168, 0.401101385,  
     1            0.440779304, 0.486141867, 0.593885499, 0.719248494, 0.848115514, 0.96522535, 1.174894012, 1.360979471,  
     1            1.38307024, 1.382803719, 1.391730855, 1.36799257, 1.379913313 /  
      data a7 /  0.681875399, 0.679916748, 0.697566565, 1.036117503, 1.229932174, 1.534436374, 1.263659339, 1.177341019,  
     1            1.046187707, 0.783076472, 0.56945524, 0.437054838, 0.497762038, 0.262897537, -0.126754198, -0.121313834, 
     1            -0.496676074, -0.513466165, -0.600511124, -0.425097306, -0.528599974 /  
      data a6tw /  -0.000639314, -0.000607826, -0.000577165, -0.000490096, -0.00042309, -0.000361045, -0.000251542,  
     1            -0.000161014, -8.90E-05, -3.54E-05, 1.45E-05, -4.21E-05, -7.26E-05, -0.000119483, -0.000199116,  
     1            -0.000362196, -0.000611716, -0.000869846, -0.00106641, -0.001185004, -0.0009885 /  
      data a12tw /  -0.4528715, -0.4516550, -0.4403449, -0.2766783, -0.2833841, -0.3205012, -0.4471684, -0.5552021,  
     1            -0.6466667, -0.7124316, -0.7599690, -0.7702118, -0.8037457, -0.8730668, -0.9821700, -1.0045641, -0.9337591,  
     1            -0.9174852, -0.9334706, -0.8808471, -0.9343411 /  
      data a8 / -0.074484253, -0.074417449, -0.075013122, -0.105054464, -0.116912997, -0.116149188, -0.107454957,  
     1            -0.091894335, -0.071460235, -0.056184877, -0.006874342, 0.040262618, 0.071577385, 0.090332895, 0.125873184, 
     1             0.157605114, 0.159496266, 0.149672181, 0.129557502, 0.105217306, 0.103793858 /  
      data a10 / 0.016025291, 0.017193978, 0.01828222, 0.020842842, 0.022162011, 0.022757257, 0.021797461, 0.020180594,  
     1            0.018556649, 0.016978648, 0.014555899, 0.012627818, 0.011191399, 0.009211979, 0.006851124, 0.003814084,  
     1            0.001733925, 0, 0, 0, 0 /  
      data a11 /  0.014951807, 0.014930723, 0.01491224, 0.01487185, 0.014854753, 0.014852358, 0.014893477, 0.015004213, 
     1             0.015194298, 0.01540766, 0.015952307, 0.016437613, 0.01652538, 0.016212382, 0.015784785, 0.01399451, 
     1             0.011927777, 0.009749305, 0.007785629, 0.00494863, 0.003408571 /  
      data tautw /  0.352252822, 0.349216437, 0.344782755, 0.355375576, 0.380828528, 0.388526115, 0.368100637, 0.368643954,  
     1            0.375291842, 0.365828808, 0.383635154, 0.378521294, 0.369787955, 0.375567859, 0.375762655, 0.39959416,  
     1            0.411391314, 0.428877733, 0.432913094, 0.435941037, 0.415321618 /  
      data phisstw /  0.406623005, 0.406228152, 0.404481748, 0.39813076, 0.387738619, 0.39878984, 0.420231185, 0.435630455, 
     1            0.443766425, 0.447739335, 0.440001589, 0.434476747, 0.427216336, 0.437722887, 0.430194359, 0.403473557,  
     1            0.417620732, 0.414874982, 0.40474639, 0.366802753, 0.331482122 /  
      data phis2stw /  0.342245626, 0.342125322, 0.34182419, 0.392290454, 0.441948427, 0.450860943, 0.409095063, 0.374262742, 
     1             0.352612548, 0.342556195, 0.34961642, 0.351721641, 0.380666872, 0.396790672, 0.384546494, 0.38643462,  
     1            0.383847635, 0.37473218, 0.351739376, 0.343596247, 0.348349141 /  
  
C Constant parameters            

      c4 = 10
      a3 = 0.1
      a9 = 0.25

C     regionflag     Note
C     -------------------------
C      0         for Japan+Taiwan
C      1         for Taiwan
C

C Find the requested spectral period and corresponding coefficients
      nPer = 21

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a5T =        a5(i1)        
         a13T =       a13(i1)       
         MrefT =      Mref(i1)      
         a2T =        a2(i1)        
         a14T =       a14(i1)       
         dela1T =     dela1(i1)     
         dela4T =     dela4(i1)     
         a6jpT =      a6jp(i1)      
         a12jpT =     a12jp(i1)     
         a8jpT   =     a8jp(i1)          
         phisstjT =   phisstj(i1)   
         phis2stjT =  phis2stj(i1)  
         tau0T =      tau0(i1)      
         a1twT =        a1tw(i1)        
         a4twT =        a4tw(i1)        
         a7T =        a7(i1)        
         a6twT =        a6tw(i1)        
         a12twT =       a12tw(i1)       
         a8T   =       a8(i1)            
         a11T =       a11(i1)       
         a10T =       a10(i1)       
         phisstwT =   phisstw(i1)   
         phis2stwT =  phis2stw(i1)  
         tautwT =     tautw(i1)     
         phitwT =     phitw(i1)     

         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Phung et al. Subduction (2018 Model) Horizontal'
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
 1020 call S24_interp (period(count1),period(count2),a5(count1), a5(count2),
     +                 specT, a5T, iflag)
      call S24_interp (period(count1),period(count2),a13(count1), a13(count2),
     +                 specT, a13T, iflag)
      call S24_interp (period(count1),period(count2),Mref(count1), Mref(count2),
     +                 specT, MrefT, iflag)
      call S24_interp (period(count1),period(count2),a2(count1), a2(count2),
     +                 specT, a2T, iflag)
      call S24_interp (period(count1),period(count2),a14(count1), a14(count2),
     +                 specT, a14T, iflag)

      call S24_interp (period(count1),period(count2),dela1(count1), dela1(count2),
     +                 specT, dela1T, iflag)
      call S24_interp (period(count1),period(count2),dela4(count1), dela4(count2),
     +                 specT, dela4T, iflag)

      call S24_interp (period(count1),period(count2),a6jp(count1), a6jp(count2),
     +                 specT, a6jpT, iflag)
      call S24_interp (period(count1),period(count2),a12jp(count1), a12jp(count2),
     +                 specT, a12jpT, iflag)
      call S24_interp (period(count1),period(count2),a8jp(count1), a8jp(count2),
     +                 specT, a8jpT, iflag)
  
      call S24_interp (period(count1),period(count2),phisstj(count1), phisstj(count2),
     +                 specT, phisstjT, iflag)
      call S24_interp (period(count1),period(count2),phis2stj(count1), phis2stj(count2),
     +                 specT, phis2stjT, iflag)       
      call S24_interp (period(count1),period(count2),tau0(count1), tau0(count2),
     +                 specT, tau0T, iflag)

      call S24_interp (period(count1),period(count2),a1tw(count1), a1tw(count2),
     +                 specT, a1twT, iflag)
      call S24_interp (period(count1),period(count2),a4tw(count1), a4tw(count2),
     +                 specT, a4twT, iflag)
      call S24_interp (period(count1),period(count2),a7(count1), a7(count2),
     +                 specT, a7T, iflag)
  
      call S24_interp (period(count1),period(count2),a6tw(count1), a6tw(count2),
     +                 specT, a6twT, iflag)
      call S24_interp (period(count1),period(count2),a12tw(count1), a12tw(count2),
     +                 specT, a12twT, iflag)
      call S24_interp (period(count1),period(count2),a8(count1), a8(count2),
     +                 specT, a8T, iflag)

      call S24_interp (period(count1),period(count2),a10(count1), a10(count2),
     +                 specT, a10T, iflag)
      call S24_interp (period(count1),period(count2),a11(count1), a11(count2),
     +                 specT, a11T, iflag)

      call S24_interp (period(count1),period(count2),tautw(count1), tautw(count2),
     +                 specT, tautwT, iflag)

      call S24_interp (period(count1),period(count2),phisstw(count1), phisstw(count2),
     +                 specT, phisstwT, iflag)
      call S24_interp (period(count1),period(count2),phis2stw(count1), phis2stw(count2),
     +                 specT, phis2stwT, iflag)       


 1011 period1 = specT                                                                                                              

C     Regional term
      if(ftype .eq. 0.0) then 
        fevt = 0.0
      elseif(ftype .eq. 1.0) then 
       fevt = 1.0
      endif

C  Regional term and  Basin Depth term
      if(regionflag .eq. 0) then
       
        a1 = a1twT + dela1T
        a4 = a4twT + dela4T
        a6 = a6jpT
        a12 = a12jpT

        tau = tau0T
        phiss = phisstjT
        phis2s = phis2stjT

      elseif(regionflag .eq. 1) then
       
        a1 = a1twT
        a4 = a4twT
        a6 = a6twT
        a12 = a12twT

        tau = tautwT
        phiss = phisstwT
        phis2s = phis2stwT
  
      endif

C     Magnitude Scaling
      if (mag .le. MrefT ) then
        fmag = a4*(mag-MrefT) + a13T*(10.0-mag)**2.0
      else
        fmag = a5T*(mag-MrefT) + a13T*(10.0-mag)**2.0
      endif 
   
C     Ztor Scaling        
      if  (ftype .eq. 0.0 ) then
         fztor = a10T *(min(Ztor,40.0)-20)
      elseif (ftype .eq. 1.0 ) then
         fztor = a11T *(min(Ztor,80.0)-40)
      endif
      
C     Path Scaling
       R = rRup + c4*exp( (mag-6.0)*a9 ) 
       frup = a1 + a7T*fevt +(a2T + a14T*fevt + a3*(mag - 7.8))*alog(R) + a6*rRup 
     
C     Site Effect
       fsite = a12*min(alog(vs30/760.0),0.0)

C   Basin Depth term
      if(regionflag .eq. 1) then
       
        Ez1 = exp(-4.06/2.0 * alog((vs30**2.0 + 352.7**2.0)/(1750.0**2.0 + 352.7**2.0)))
        fz10 = a8T*(min(alog(Z10*1000.0/Ez1),0.0))     
  
      else
  
        Ez1 = exp(-5.23/2.0 * alog((vs30**2.0 + 412.39**2.0)/(1360.0**2.0 + 412.39**2.0)))
        fz10 = a8jpT*(min(alog(Z10*1000.0/Ez1),0.0))     
  
      endif       


       lnSa = fmag + frup + fztor + fsite + fz10 
   
C     Set sigma values to return
C       tau = tau1T
       sigma = sqrt(tau**2+phiSS**2+phiS2S**2)
    
c     write(*,*) "fz10 = ", fz10
c      write(*,*) "fmag = ", fmag
c     write(*,*) "X = ", frup
c     write(*,*) "fsite = ", fsite
c     write(*,*) "fztor = ", fztor
c     write(*,*) "lnSa = ", lnSa
c     write(*,*) "Sa = ", exp(lnSa)
 
C     Convert ground motion to units of gals.
      lnY = lnSa + 6.89
      period2 = period1
      return
      END
   
c ------------------------------------------------------------------            
C *** Chao2020 (Crustal and Subduction - Model) Horizontal ***********
c ------------------------------------------------------------------            
      subroutine S04_Chao2020 ( mag, dist, ftype, lnY, sigma, specT, vs, Ztor, Z10,           
     1            vs30_class, attenName, period2, iflag, sourcetype, phi, tau, msasflag )         

      implicit none
C    2018/08/31 revised
      real mag, dip, fType, dist, vs, SA1180,
     1      Z10,  ZTOR, fltWidth, lnSa, sigma, lnY, vs30_rock, sourcetype
      real Fn, Frv, specT, period2, CRjb, phi, tau, z10_rock, SA_rock
      integer hwflag, iflag, vs30_class, regionflag, msasflag
      character*80 attenName                                                    

C     Set the reference spectrum.                
c     sourcetype = 0 for crustal
c                  1 for Subduction 
c     Vs30_class = 0 for estimated
c     Vs30_class = 1 for measured 
C     Mainshock and Aftershocks included based on MSASFlag
C         0 = Mainshocks
C         1 = Aftershocks

c     Compute SA1180
      vs30_rock = 1180.
      z10_rock = 0.004541444
      SA_rock = 0.
      
         call S04_Chaoetal2018 ( mag, dist, ftype, sigma, specT, vs30_rock, Ztor, z10_rock,
     1             SA_rock, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag)
      Sa1180 = exp(lnSa)

c     Compute Sa at spectral period for given Vs30

         call S04_Chaoetal2018 ( mag, dist, ftype, sigma, specT, vs, Ztor, Z10,
     1             sa1180, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag )

C     Convert ground motion to units of gals.
      lnY = lnSa + 6.89

      period2 = specT

      return
      end
c -------------------------------------------------------------------           
C **** Chao et al. 2020 (SSHAC model) *************
c -------------------------------------------------------------------           

      subroutine S04_Chaoetal2020 ( mag, dist, ftype, sigma, specT, vs, Ztor, Z10,           
     1            sa1180, vs30_class, attenName, iflag, sourcetype, phi, tau, lnSa, msasflag )                                   

      implicit none
      
      integer MAXPER                                                                            
      parameter (MAXPER=21)                                                     
      real ftype, dist, mag, lnSa, sigma, specT, lnYref, vs, Ztor, Z10, period1
      real period(MAXPER), c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER), c11(MAXPER)
      real c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER), c16(MAXPER), c17(MAXPER)
      real c18(MAXPER), c19(MAXPER), c20(MAXPER), c21(MAXPER), c22(MAXPER), c23(MAXPER)
      real c24(MAXPER), c25(MAXPER), c26(MAXPER), c27(MAXPER), taucr1(MAXPER), taucr2(MAXPER)
      real tausb1(MAXPER), tausb2(MAXPER), phisscr1(MAXPER), phisscr2(MAXPER), phisssb1(MAXPER)
      real phisssb2(MAXPER), arfacr(MAXPER), arfasb(MAXPER), phis2s(MAXPER)
      character*80 attenName                                                    
      integer nper, count1, count2, C11flag, C23flag, C29flag, iflag, C10flag, C13flag
      integer vs30_class, n, i, msasflag
      integer Fcr, Fsb, Fcrss, Fcrno, Fcrro, Fsbintra, Fsbinter, Fas, Fkuo17, Fks17, Frf, Fmanila 
      real Mc, Mref, Mmax, Rrupref, Vs30ref, Zref, sourcetype
      real c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, c10T, c11T, c12T, c13T, c14T, c15T
      real c16T, c17T, c18T, c19T, c20T, c21T, c22T, c23T, c24T, c25T, c26T, c27T
      real taucr1T, taucr2T, tausb1T, tausb2T, phisscr1T, phisscr2T, phisssb1T, phisssb2T
      real arfacrT, arfasbT, phis2sT, phi, tau, fm, SA1180, Z10ref
      real Ssource, Spath, Ssite, Ssitelin, Ssitenon, Sztor, Smag, Sgeom, Sanel
      real taucr, tausb, phisscr, phisssb, phiss, sigmass
      real c28(MAXPER), c29(MAXPER), c28T, c29T, c30(MAXPER), c30T, h
                                                                                
      data Period / 0, -2, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75,
     1              1, 1.5, 2, 3, 4, 5 / 
      data c1 / -1.042666, 1.895411, -1.042360, -1.009040, -0.962601, -0.870563, -0.785856, -0.727776, 
     1          -0.674517, -0.668684, -0.682285, -0.707153, -0.782246, -0.867588, -1.100953, -1.335849,  
     1          -1.769889, -2.158271, -2.821190, -3.354776, -3.867019 / 
      data c2 / -1.137852, 1.735667, -1.137811, -1.101142, -1.049407, -0.939899, -0.832404, -0.758714,  
     1          -0.705514, -0.716938, -0.752775, -0.800309, -0.908875, -1.016989, -1.278817, -1.515900,  
     1          -1.926198, -2.280583, -2.876528, -3.355129, -3.810868 / 
      data c3 / -1.170990, 1.587449, -1.171510, -1.130461, -1.066827, -0.931421, -0.802566, -0.723218,  
     1          -0.690970, -0.741863, -0.819872, -0.908821, -1.075545, -1.217121, -1.508285, -1.746629,  
     1          -2.151367, -2.498264, -3.089388, -3.563977, -4.047973 / 
      data c4 / -1.314590, 1.882708, -1.298813, -1.262369, -1.214610, -1.092116, -0.964326, -0.894881,  
     1          -0.836942, -0.792193, -0.763280, -0.744597, -0.723252, -0.734591, -0.883649, -1.088087,  
     1          -1.526738, -1.947331, -2.666122, -3.258068, -3.809328 / 
      data c5 / -0.457036, 2.359015, -0.436214, -0.384335, -0.311696, -0.131992, 0.074716, 0.200368, 0.286613,  
     1          0.254600, 0.173035, 0.087029, -0.073676, -0.220276, -0.567236, -0.881354, -1.445217, -1.943705,  
     1          -2.722921, -3.349385, -3.792075 / 
      data c6 / -0.124087, -0.144577, -0.123275, -0.124986, -0.127909, -0.127668, -0.125831, -0.122221,  
     1          -0.115139, -0.105684, -0.095876, -0.088255, -0.081799, -0.084061, -0.105026, -0.125539,  
     1          -0.137664, -0.125372, -0.085938, -0.049451, 0.012311 / 
      data c7 / 0.198357, -0.038842, 0.191344, 0.195147, 0.198967, 0.227624, 0.260255, 0.289816, 0.348032,  
     1          0.370540, 0.361912, 0.339973, 0.258518, 0.169337, 0.025880, -0.053324, -0.113254, -0.138977,  
     1          -0.193895, -0.258713, -0.326743 / 
      data c8 / 0.676721, 1.251174, 0.684339, 0.657072, 0.610194, 0.544252, 0.548393, 0.602698, 0.739429,  
     1          0.875915, 0.986450, 1.078964, 1.233393, 1.354611, 1.564898, 1.701387, 1.863505, 1.953011,  
     1          2.043382, 2.088718, 2.105341 / 
      data c9 / 0.633834, 0.766604, 0.640018, 0.606922, 0.558754, 0.535818, 0.603777, 0.687135, 0.829456,  
     1          0.904226, 0.958892, 0.984563, 1.007287, 1.025986, 1.013332, 1.002811, 0.962813, 0.979659,  
     1          1.052087, 1.157880, 1.147008 / 
      data c10 / -0.135342, -0.210085, -0.136868, -0.131414, -0.122038, -0.108827, -0.109545, -0.120401,  
     1          -0.147755, -0.175139, -0.197287, -0.215793, -0.246679, -0.270922, -0.312979, -0.339111,  
     1          -0.357148, -0.352155, -0.323695, -0.294607, -0.250719 / 
      data c11 / -0.001320, -0.000002, -0.000002, -0.006396, -0.008323, -0.005923, -0.001799, -0.010118,  
     1          -0.080308, -0.163717, -0.238633, -0.297686, -0.346393, -0.334133, -0.229548, -0.131664,  
     1          -0.032403, -0.006790, -0.000014, -0.000002, -0.000001 / 
      data c12 / -0.000034, -0.000014, -0.000001, -0.000004, -0.000005, -0.000016, -0.000069, -0.009043,  
     1          -0.105380, -0.189907, -0.238782, -0.268913, -0.261917, -0.203064, -0.102740, -0.049805,  
     1          -0.011622, -0.002702, -0.000015, -0.000001, -0.000001 / 
      data c13 / -0.000038, -0.291529, -0.000001, -0.000004, -0.000005, -0.000018, -0.000080, -0.000082,  
     1          -0.011888, -0.070362, -0.139438, -0.217492, -0.359123, -0.475375, -0.670648, -0.775270,  
     1          -0.856757, -0.798850, -0.639826, -0.461969, -0.465749 / 
      data c14 / 0.032555, 0.015510, 0.032494, 0.033028, 0.034285, 0.037352, 0.039960, 0.040658, 0.038553,  
     1          0.034446, 0.030226, 0.026384, 0.020530, 0.016445, 0.010777, 0.008159, 0.005621, 0.004032,  
     1          0.000622, -0.002970, -0.008948 / 
      data c15 / 0.017948, 0.016982, 0.018059, 0.018141, 0.018067, 0.017682, 0.017907, 0.018645, 0.020842,  
     1          0.022069, 0.022125, 0.021736, 0.020462, 0.018980, 0.016322, 0.014953, 0.012425, 0.009349,  
     1          0.004309, -0.000454, 0.000502 / 
      data c16 / 0.006634, 0.003036, 0.006626, 0.007008, 0.007575, 0.008618, 0.009183, 0.009159, 0.008404,  
     1          0.007375, 0.006330, 0.005344, 0.003788, 0.002623, 0.001242, 0.000937, 0.000745, 0.000647,  
     1          0.000241, 0.000051, -0.000294 / 
      data c17 / -1.687472, -1.547786, -1.692461, -1.730349, -1.772040, -1.821470, -1.810850, -1.759544,  
     1          -1.643299, -1.551105, -1.488492, -1.446763, -1.398596, -1.375008, -1.347383, -1.333650,  
     1          -1.311795, -1.298714, -1.284844, -1.281901, -1.241532 / 
      data c18 / -1.525058, -1.549392, -1.539999, -1.565127, -1.589292, -1.642560, -1.682785, -1.676278,  
     1          -1.637058, -1.599726, -1.564986, -1.538416, -1.509663, -1.488699, -1.436447, -1.401868,  
     1          -1.333535, -1.279867, -1.210371, -1.152781, -1.113843 / 
      data c19 / 0.393879, 0.303285, 0.388200, 0.395410, 0.410623, 0.421673, 0.405707, 0.381442, 0.338659,  
     1          0.299124, 0.268677, 0.244569, 0.210459, 0.192077, 0.179792, 0.184209, 0.197091, 0.210361,  
     1          0.233088, 0.252812, 0.264010 / 
      data c20 / 0.192940, 0.223730, 0.188931, 0.201543, 0.217245, 0.203093, 0.152740, 0.113630, 0.072366,  
     1          0.055866, 0.044190, 0.040593, 0.041938, 0.044384, 0.071700, 0.099141, 0.153764, 0.199871,  
     1          0.253750, 0.282276, 0.298510 / 
      data c21 / -0.003430, -0.000817, -0.003338, -0.003155, -0.003190, -0.003740, -0.004703, -0.005494,  
     1          -0.006073, -0.005732, -0.005091, -0.004380, -0.003179, -0.002346, -0.001304, -0.000879,  
     1          -0.000564, -0.000475, -0.000459, -0.000460, -0.000879 / 
      data c22 / -0.003964, -0.000992, -0.003844, -0.003962, -0.004180, -0.004364, -0.004364, -0.004447,  
     1          -0.004295, -0.003755, -0.003171, -0.002613, -0.001680, -0.001075, -0.000534, -0.000357,  
     1          -0.000460, -0.000784, -0.001435, -0.002141, -0.002623 / 
      data c23 / -2.405229, -6.913488, -2.378641, -2.418630, -2.313209, -2.067232, -1.854300, -1.693198,  
     1          -1.479479, -1.358600, -1.314082, -1.320243, -1.391252, -1.483500, -1.558715, -1.404591,  
     1          -0.845899, -0.436331, -0.023390, 0.000000, 0.000000 / 
      data c24 / -0.478715, -0.672424, -0.477706, -0.470430, -0.451950, -0.416240, -0.402645, -0.410470,  
     1          -0.446615, -0.481032, -0.513790, -0.542807, -0.596688, -0.648563, -0.742180, -0.798604,  
     1          -0.842403, -0.849752, -0.839143, -0.822511, -0.797141 / 
      data c25 / 0.063345, 0.095748, 0.063485, 0.064655, 0.068531, 0.079314, 0.084622, 0.082002, 0.070590,  
     1          0.062519, 0.061465, 0.064345, 0.074007, 0.083563, 0.103379, 0.118159, 0.139920, 0.153222, 
     1           0.160180, 0.155343, 0.142666 / 
      data c26 / -0.604023, 0.384203, -0.598271, -0.547902, -0.480919, -0.338469, -0.228587, -0.197662,  
     1          -0.245076, -0.344219, -0.449905, -0.550299, -0.725882, -0.874283, -1.154718, -1.353035,  
     1          -1.618376, -1.771092, -1.893009, -1.912582, -1.838861 / 
      data c27 / -0.679356, 0.383524, -0.673352, -0.626423, -0.567634, -0.437480, -0.333207, -0.303548, 
     1           -0.352291, -0.441969, -0.532257, -0.615964, -0.760168, -0.886623, -1.140546, -1.331239,  
     1          -1.584174, -1.727634, -1.829259, -1.833069, -1.746378 / 
      data c28 / -0.650493, 0.305582, -0.644569, -0.596370, -0.529756, -0.377565, -0.258657, -0.227982,  
     1          -0.292081, -0.413965, -0.533268, -0.639714, -0.817367, -0.960564, -1.232582, -1.435782,  
     1          -1.715330, -1.875696, -1.998774, -2.017655, -1.941254 / 
      data c29  / -0.434706511, -0.34661917, -0.447102839, -0.415446676, -0.366999977, -0.366148593, -0.469248166,   
     1          -0.568918055, -0.719037119, -0.790626985, -0.840621112, -0.860236343, -0.873968925, -0.885430175,   
     1          -0.827526363, -0.742505362, -0.57420172, -0.49097237, -0.423101316, -0.426861488, -0.338104658  /
      data c30  / -0.444373594, -0.232141657, -0.460428073, -0.448496325, -0.436361774, -0.492546506, -0.54896049,   
     1          -0.60832263, -0.675866324, -0.689506592, -0.695043713, -0.678658902, -0.636360976, -0.607212602,   
     1          -0.503495503, -0.43304903, -0.313708518, -0.28272403, -0.308813684, -0.394492845, -0.370973146 /
      data taucr1 / 0.367027, 0.439871, 0.366713, 0.367381, 0.365814, 0.360738, 0.364067, 0.376999,  
     1          0.420380, 0.470706, 0.514142, 0.545754, 0.580926, 0.595221, 0.586341, 0.567362, 0.539173,  
     1          0.521031, 0.510998, 0.516268, 0.543015 / 
      data taucr2 / 0.315219, 0.374796, 0.315430, 0.319633, 0.327049, 0.345572, 0.360231, 0.362292, 0.339502,  
     1          0.307997, 0.285602, 0.273645, 0.271198, 0.286852, 0.345657, 0.391395, 0.441542, 0.462329,  
     1          0.475059, 0.480575, 0.470049 / 
      data tausb1 / 0.272750, 0.330733, 0.271058, 0.271673, 0.274006, 0.280698, 0.291948, 0.305046, 0.334109,  
     1          0.369176, 0.401260, 0.428927, 0.466550, 0.482919, 0.477346, 0.450142, 0.405823, 0.366654,  
     1          0.325573, 0.294383, 0.340797 / 
      data tausb2 / 0.532703, 0.566695, 0.536421, 0.557253, 0.577443, 0.601923, 0.599473, 0.579989, 0.528112,  
     1          0.495023, 0.475578, 0.464138, 0.466708, 0.482046, 0.524715, 0.574218, 0.638071, 0.666590,  
     1          0.648573, 0.596551, 0.507714 / 
      data phisscr1 / 0.530015, 0.558388, 0.530680, 0.521767, 0.513496, 0.502312, 0.503270, 0.517247, 0.553158,  
     1          0.581637, 0.599292, 0.608899, 0.610791, 0.600179, 0.561837, 0.526088, 0.477128, 0.448854,  
     1          0.421072, 0.408714, 0.406221 / 
      data phisscr2 / 0.433766, 0.435606, 0.435202, 0.444543, 0.455987, 0.470816, 0.466840, 0.451769,  
     1          0.427517, 0.417334, 0.416486, 0.422015, 0.436118, 0.448613, 0.471981, 0.484933, 0.492704,  
     1          0.488214, 0.469356, 0.447792, 0.429747 / 
      data phisssb1 / 0.430634, 0.478613, 0.431917, 0.426561, 0.421787, 0.407957, 0.406584, 0.419057,  
     1          0.454812, 0.477921, 0.489639, 0.495038, 0.494487, 0.491319, 0.486749, 0.481815, 0.481313,  
     1          0.479424, 0.475432, 0.466465, 0.454073 / 
      data phisssb2 / 0.495377, 0.485676, 0.494439, 0.498452, 0.505789, 0.519604, 0.525237, 0.521053, 
     1           0.511584, 0.503919, 0.496864, 0.491227, 0.487558, 0.486774, 0.488927, 0.493538, 0.488128, 
     1           0.477724, 0.446953, 0.414410, 0.368507 / 
      data phis2s / 0.342398, 0.273837, 0.342667, 0.348735, 0.363743, 0.407580, 0.443690, 0.455723,  
     1          0.439179, 0.411759, 0.389287, 0.372044, 0.352699, 0.344485, 0.341557, 0.347037, 0.357107, 
     1           0.363360, 0.369654, 0.372377, 0.375576 / 

c Set attenuation name                                                            
c     Sourcetype = 0 Crustal
c     Sourcetype = 1 Subduction 
                                                                       
C Find the requested spectral period and corresponding coefficients
      nper = 21

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
        c1T = c1(1)
        c2T = c2(1)
        c3T = c3(1)
        c4T = c4(1)
        c5T = c5(1)
        c6T = c6(1)
        c7T = c7(1)
        c8T = c8(1)
        c9T = c9(1)
        c10T = c10(1)
        c11T = c11(1)
        c12T = c12(1)
        c13T = c13(1)
        c14T = c14(1)
        c15T = c15(1)
        c16T = c16(1)
        c17T = c17(1)
        c18T = c18(1)
        c19T = c19(1)
        c20T = c20(1)
        c21T = c21(1)
        c22T = c22(1)
        c23T = c23(1)
        c24T = c24(1)
        c25T = c25(1)
        c26T = c26(1)
        c27T = c27(1)
        c28T = c28(1)
        c29T = c29(1)        
        c30T = c30(1)        
        taucr1T = taucr1(1)
        taucr2T = taucr2(1)
        tausb1T = tausb1(1)
        tausb2T = tausb2(1)
        phisscr1T = phisscr1(1)
        phisscr2T = phisscr2(1)
        phisssb1T = phisssb1(1)
        phisssb2T = phisssb2(1)
        phis2sT = phis2s(1)
       goto 1011
C   Function Form for PGV Regression     
       elseif (specT .eq. -2.0 .or. specT .eq. -1.0) then
         period1 = period(2)
         c1T = c1(2)
         c2T = c2(2)
         c3T = c3(2)
         c4T = c4(2)
         c5T = c5(2)
         c6T = c6(2)
         c7T = c7(2)
         c8T = c8(2)
         c9T = c9(2)
         c10T = c10(2)
         c11T = c11(2)
         c12T = c12(2)
         c13T = c13(2)
         c14T = c14(2)
         c15T = c15(2)
         c16T = c16(2)
         c17T = c17(2)
         c18T = c18(2)
         c19T = c19(2)
         c20T = c20(2)
         c21T = c21(2)
         c22T = c22(2)
         c23T = c23(2)
         c24T = c24(2)
         c25T = c25(2)
         c26T = c26(2)
         c27T = c27(2)
         c28T = c28(2)
         c29T = c29(2)        
         c30T = c30(2)        
         taucr1T = taucr1(2)
         taucr2T = taucr2(2)
         tausb1T = tausb1(2)
         tausb2T = tausb2(2)
         phisscr1T = phisscr1(2)
         phisscr2T = phisscr2(2)
         phisssb1T = phisssb1(2)
         phisssb2T = phisssb2(2)
         phis2sT = phis2s(2)
         goto 1011      
       endif
C Now loop over the spectral period range of the attenuation relationship.
         do i=3,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1010 
            endif
         enddo
        
      write (*,*) 
      write (*,*) 'Chao et al. (2018) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1010    call S24_interp (period(count1),period(count2),c1(count1),c1(count2), 
     +                specT,c1T,iflag)
         call S24_interp (period(count1),period(count2),c2(count1),c2(count2), 
     +                specT,c2T,iflag)
         call S24_interp (period(count1),period(count2),c3(count1),c3(count2), 
     +                specT,c3T,iflag)
         call S24_interp (period(count1),period(count2),c4(count1),c4(count2), 
     +                specT,c4T,iflag)
         call S24_interp (period(count1),period(count2),c5(count1),c5(count2), 
     +                specT,c5T,iflag)
         call S24_interp (period(count1),period(count2),c6(count1),c6(count2), 
     +                specT,c6T,iflag)
         call S24_interp (period(count1),period(count2),c7(count1),c7(count2), 
     +                 specT,c7T,iflag)
         call S24_interp (period(count1),period(count2),c8(count1),c8(count2), 
     +                 specT,c8T,iflag)
         call S24_interp (period(count1),period(count2),c9(count1),c9(count2), 
     +                 specT,c9T,iflag)
         call S24_interp (period(count1),period(count2),c10(count1),c10(count2), 
     +                 specT,c10T,iflag)
         call S24_interp (period(count1),period(count2),c11(count1),c11(count2), 
     +                 specT,c11T,iflag)
         call S24_interp (period(count1),period(count2),c12(count1),c12(count2), 
     +                 specT,c12T,iflag)
         call S24_interp (period(count1),period(count2),c13(count1),c13(count2), 
     +                 specT,c13T,iflag)
         call S24_interp (period(count1),period(count2),c14(count1),c14(count2), 
     +                 specT,c14T,iflag)
         call S24_interp (period(count1),period(count2),c15(count1),c15(count2), 
     +                 specT,c15T,iflag)
         call S24_interp (period(count1),period(count2),c16(count1),c16(count2), 
     +                 specT,c16T,iflag)
         call S24_interp (period(count1),period(count2),c17(count1),c17(count2), 
     +                 specT,c17T,iflag)
         call S24_interp (period(count1),period(count2),c18(count1),c18(count2), 
     +                 specT,c18T,iflag)
         call S24_interp (period(count1),period(count2),c19(count1),c19(count2), 
     +                 specT,c19T,iflag)
         call S24_interp (period(count1),period(count2),c20(count1),c20(count2), 
     +                 specT,c20T,iflag)
         call S24_interp (period(count1),period(count2),c21(count1),c21(count2), 
     +                 specT,c21T,iflag)
         call S24_interp (period(count1),period(count2),c22(count1),c22(count2), 
     +                 specT,c22T,iflag)
         call S24_interp (period(count1),period(count2),c23(count1),c23(count2), 
     +                 specT,c23T,iflag)
         call S24_interp (period(count1),period(count2),c24(count1),c24(count2), 
     +                 specT,c24T,iflag)
         call S24_interp (period(count1),period(count2),c25(count1),c25(count2), 
     +                 specT,c25T,iflag)
         call S24_interp (period(count1),period(count2),c26(count1),c26(count2), 
     +                 specT,c26T,iflag)
         call S24_interp (period(count1),period(count2),c27(count1),c27(count2), 
     +                 specT,c27T,iflag)
         call S24_interp (period(count1),period(count2),taucr1(count1),taucr1(count2), 
     +                 specT,taucr1T,iflag)
         call S24_interp (period(count1),period(count2),taucr2(count1),taucr2(count2), 
     +                 specT,taucr2T,iflag)
         call S24_interp (period(count1),period(count2),tausb1(count1),tausb1(count2), 
     +                 specT,tausb1T,iflag)
         call S24_interp (period(count1),period(count2),tausb2(count1),tausb2(count2), 
     +                 specT,tausb2T,iflag)
         call S24_interp (period(count1),period(count2),phisscr1(count1),phisscr1(count2), 
     +                 specT,phisscr1T,iflag)
         call S24_interp (period(count1),period(count2),phisscr2(count1),phisscr2(count2), 
     +                 specT,phisscr2T,iflag)
         call S24_interp (period(count1),period(count2),phisssb1(count1),phisssb1(count2), 
     +                 specT,phisssb1T,iflag)
         call S24_interp (period(count1),period(count2),phisssb2(count1),phisssb2(count2), 
     +                 specT,phisssb2T,iflag)
         call S24_interp (period(count1),period(count2),phis2s(count1),phis2s(count2), 
     +                specT,phis2sT,iflag)
         call S24_interp (period(count1),period(count2),c28(count1),c28(count2), 
     +                 specT,c28T,iflag)
         call S24_interp (period(count1),period(count2),c29(count1),c29(count2), 
     +                 specT,c29T,iflag)
         call S24_interp (period(count1),period(count2),c30(count1),c30(count2), 
     +                 specT,c30T,iflag)

  
 1011 period1 = specT

C      h = 10.0
      n = 2.0
      Mc = 7.1
      Mref = 5.5
      Mmax = 8
      Rrupref = 0.0
      Vs30ref = 760.0
    
C     Set the reference spectrum.                
c     sourcetype = 0 for crustal
c                  1 for Subduction 
c     Vs30_class = 0 for estimated
c     Vs30_class = 1 for measured 

      Fcr=0
      Fsb=0
      Fcrss = 0
      Fcrno = 0
      Fcrro = 0
      Fsbintra = 0
      Fsbinter = 0
      Fas = 0
      Fkuo17 = 0
      Fks17 = 0
      Frf = 0
      Fmanila = 0
      C11flag = 0
      C23flag = 0
      C29flag = 0
      C13flag = 0 
      C10flag = 0
   
      if (sourcetype .eq. 0.0 ) then
       Fcr = 1
       Zref = 0
         if(ftype .gt. 0) then
              Fcrro = 1
           elseif(ftype .lt. 0) then
              Fcrno = 1
           else
              Fcrss = 1
         endif
      elseif (sourcetype .eq. 1.0 ) then
        Fsb = 1
         if(ftype .eq. 0) then
              Fsbinter = 1
              Zref = 0
         elseif(ftype .eq. 1) then
              Fsbintra = 1
              Zref = 35
         endif
      endif
   
C     Add aftershock factor 
      if (msasflag .eq. 1) then
           Fas = 1
      endif 

C     choose Site ref by Vs30 class
        if (vs30_class .eq. 0 ) then
         Fks17 = 1
        elseif (vs30_class .eq. 1) then
         Fkuo17 = 1
        endif

      lnYref = c1T*Fcrro + c2T*Fcrss + c3T*Fcrno + c4T*Fsbinter + c5T*Fsbintra +
     &         c6T*Fas + c7T*Fmanila + c26T*Fkuo17 + c27T*Fks17 + c28T*Frf

C     Set Source scaling term 
     
      if(mag .LE. 5 ) then  
       C11flag=1
      endif
      if(mag .GE. Mc ) then  
       C29flag=1
      endif
      if(mag .GE. 7.6 ) then  
       C10flag=1
      endif   
      if(mag .LE. 6 ) then  
       C13flag=1
      endif   
      if (sourcetype .eq. 0.0 ) then
        Smag = c8T*(mag - Mref) + c10T*(mag - Mref)**2 
     1         - c10T*(mag-7.6)**2*C10flag + c11T*(5.0-mag)*C11flag  
      elseif (sourcetype .eq. 1.0 ) then
        Smag = c9T*(mag - Mref) + c29T*Fsbinter*(Mag-Mc)*c29flag + c30T*Fsbintra*(Mag-Mc)*c29flag 
     1        + c12T*(5.0-mag)*C11flag + c13T*(6.0-mag)*C13flag
      endif

      Sztor = c14T * Fcr *(Ztor-Zref) + c15T * Fsbinter * (Ztor-Zref) + c16T * Fsbintra * (Ztor-Zref)    
      Ssource = Smag + Sztor

C     Set Path scaling term

      h = 10.0*Fcr +10.0*Fsbinter*exp(0.3*(mag-7.1)*C29flag) + 10.0*Fsbintra*exp(0.2*(mag-7.1)*C29flag)
   
      if (sourcetype .eq. 0.0 ) then
          Sgeom = (c17T + c19T*(min(mag,Mmax)- Mref )) * alog(SQRT(dist**2 + h**2)/SQRT(Rrupref**2 + h**2))
      elseif (sourcetype .eq. 1.0 ) then
          Sgeom = (c18T + c20T*(min(mag, Mc )- Mref )) * alog(SQRT(dist**2 + h**2)/SQRT(Rrupref**2 + h**2))
      endif

      Sanel = c21T*Fcr*(dist-Rrupref) + c22T*Fsb*(dist-Rrupref)
      Spath = Sgeom + Sanel 
    
C     Set Site scaling term 
    
      Z10ref = exp((-4.08/2.0)*alog((vs**2.0+355.4**2.0)/(1750**2.0+355.4**2.0)))
      if (Z10 .gt. 0.) then
        Ssitelin = c24T * alog(vs/vs30ref) + c25T*alog(Z10*1000/Z10ref)
      else
        Ssitelin = c24T * alog(vs/vs30ref)
      endif

      if(vs .LT. vs30ref ) then  
           C23flag=1
      endif
     
      Ssitenon = c23T * C23flag * (-1.5*alog(vs/vs30ref)-alog(SA1180+2.4)+alog(SA1180+2.4*(vs/vs30ref)**1.5))  
      Ssite = Ssitenon + Ssitelin

      lnSa =  lnYref + Ssource + Spath + Ssite                                        
   
C      write(*,*) "lnYref = ", lnYref
C      write(*,*) "Ssource = ", Ssource
C      write(*,*) "--Smag = ", Smag
C      write(*,*) "--Sztor = ", Sztor
C      write(*,*) "Spath = ", Spath
C      write(*,*) "--Sgeom = ", Sgeom
C      write(*,*) "--Sanel = ", Sanel
C      write(*,*) "Ssite = ", Ssite
C      write(*,*) "--Ssitelin = ", Ssitelin
C      write(*,*) "--Ssitenon = ", Ssitenon
C      write(*,*) "lnSa = ", lnSa
C      write(*,*) "Sa = ", exp(lnSa)

   
C     Set the event-specific residual term
 
      fm = 0.5*(min(6.5, max(4.5, mag))-4.5)
      
      taucr = taucr1T + (taucr2T - taucr1T)*fm
      tausb = tausb1T + (tausb2T - tausb1T)*fm 
      
      tau = taucr*Fcr + tausb*Fsb   
   
C     Set Site-specific residual term

      

C     Set Recoed-specific residual term

      phisscr = phisscr1T + (phisscr2T -phisscr1T)*fm
      phisssb = phisssb1T + (phisssb2T -phisssb1T)*fm

      phiss = phisscr*Fcr + phisssb*Fsb
      
      phi=(phis2sT**2+phiss**2)**0.5
      sigma=(tau**2+phi**2)**0.5
      sigmass=(tau**2+phiss**2)**0.5



c       write(*,*) "Y(gal) = ", exp(lnSa)

      return                                                                    
      end       