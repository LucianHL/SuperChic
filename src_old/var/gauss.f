      subroutine gaussinit
      implicit none
      integer nkt,nphi,nb,nphib,nphi4,nkt4

      include 'nsurv.f'
      include 'gaussvars.f'

      nphi4=s2int*4
      nkt4=s2int*4
      nphi=s2int*2
      nkt=s2int*2

      nphi=s2int*4
      nkt=s2int*4

      nb=s2int*2
      nphib=s2int

c      print*,nphi,nkt,s2int

      if(nphi.eq.96)then
         call gauss96(xiphi,wiphi)
      elseif(nphi.eq.64)then
         call gauss64(xiphi,wiphi)
      elseif(nphi.eq.32)then
         call gauss32(xiphi,wiphi)
      elseif(nphi.eq.24)then
         call gauss24(xiphi,wiphi)
      elseif(nphi.eq.16)then
         call gauss16(xiphi,wiphi)
      elseif(nphi.eq.12)then
         call gauss12(xiphi,wiphi)
      elseif(nphi.eq.8)then
         call gauss8(xiphi,wiphi)
      elseif(nphi.eq.6)then
         call gauss6(xiphi,wiphi)
      elseif(nphi.eq.4)then
         call gauss4(xiphi,wiphi)
      else
         print*,'Gaussian points not avaiable for choice of s2int: STOP'
         STOP 1
      endif

      if(nkt.eq.96)then
         call gauss96(xikt,wikt)
      elseif(nkt.eq.64)then
         call gauss64(xikt,wikt)
      elseif(nkt.eq.32)then
         call gauss32(xikt,wikt)
      elseif(nkt.eq.24)then
         call gauss24(xikt,wikt)
      elseif(nkt.eq.16)then
         call gauss16(xikt,wikt)
      elseif(nkt.eq.12)then
         call gauss12(xikt,wikt)
      elseif(nkt.eq.8)then
         call gauss8(xikt,wikt)
      elseif(nkt.eq.6)then
         call gauss6(xikt,wikt)
      elseif(nkt.eq.4)then
         call gauss4(xikt,wikt)
      else
         print*,'Gaussian points not avaiable for choice of s2int: STOP'
         STOP 1
      endif

      if(nphi4.eq.96)then
         call gauss96(xiphi4,wiphi4)
      elseif(nphi4.eq.64)then
         call gauss64(xiphi4,wiphi4)
      elseif(nphi4.eq.32)then
         call gauss32(xiphi4,wiphi4)
      elseif(nphi4.eq.24)then
         call gauss24(xiphi4,wiphi4)
      elseif(nphi4.eq.16)then
         call gauss16(xiphi4,wiphi4)
      elseif(nphi4.eq.12)then
         call gauss12(xiphi4,wiphi4)
      elseif(nphi4.eq.8)then
         call gauss8(xiphi4,wiphi4)
      elseif(nphi4.eq.6)then
         call gauss6(xiphi4,wiphi4)
      elseif(nphi4.eq.4)then
         call gauss4(xiphi4,wiphi4)
      else
         print*,'Gaussian points not avaiable for choice of s2int: STOP'
         STOP 1
      endif

      if(nkt4.eq.96)then
         call gauss96(xikt4,wikt4)
      elseif(nkt4.eq.64)then
         call gauss64(xikt4,wikt4)
      elseif(nkt4.eq.32)then
         call gauss32(xikt4,wikt4)
      elseif(nkt4.eq.24)then
         call gauss24(xikt4,wikt4)
      elseif(nkt4.eq.16)then
         call gauss16(xikt4,wikt4)
      elseif(nkt4.eq.12)then
         call gauss12(xikt4,wikt4)
      elseif(nkt4.eq.8)then
         call gauss8(xikt4,wikt4)
      elseif(nkt4.eq.6)then
         call gauss6(xikt4,wikt4)
      elseif(nkt4.eq.4)then
         call gauss4(xikt4,wikt4)
      else
         print*,'Gaussian points not avaiable for choice of s2int: STOP'
         STOP 1
      endif

      if(nb.eq.96)then
         call gauss96(xib,wib)
      elseif(nb.eq.64)then
         call gauss64(xib,wib)
      elseif(nb.eq.32)then
         call gauss32(xib,wib)
      elseif(nb.eq.24)then
         call gauss24(xib,wib)
      elseif(nb.eq.16)then
         call gauss16(xib,wib)
      elseif(nb.eq.12)then
         call gauss12(xib,wib)
      elseif(nb.eq.8)then
         call gauss8(xib,wib)
      elseif(nb.eq.6)then
         call gauss6(xib,wib)
      elseif(nb.eq.4)then
         call gauss4(xib,wib)
      else
         print*,'Gaussian points not avaiable for choice of s2int: STOP'
         STOP 1
      endif

      if(nphib.eq.96)then
         call gauss96(xiphib,wiphib)
      elseif(nphib.eq.64)then
         call gauss64(xiphib,wiphib)
      elseif(nphib.eq.32)then
         call gauss32(xiphib,wiphib)
      elseif(nphib.eq.24)then
         call gauss24(xiphib,wiphib)
      elseif(nphib.eq.16)then
         call gauss16(xiphib,wiphib)
      elseif(nphib.eq.12)then
         call gauss12(xiphib,wiphib)
      elseif(nphib.eq.8)then
         call gauss8(xiphib,wiphib)
      elseif(nphib.eq.6)then
         call gauss6(xiphib,wiphib)
      elseif(nphib.eq.4)then
         call gauss4(xiphib,wiphib)
      else
         print*,'Gaussian points not avaiable for choice of s2int: STOP'
         STOP 1
      endif

      return
      end

      SUBROUTINE gauss96(xi,wi)
      implicit none
      double precision X(48),W(48),xi(100),wi(100)
      integer i

      X( 1)=   0.01627674484960183561
      X( 2)=   0.04881298513604856015
      X( 3)=   0.08129749546442434360
      X( 4)=   0.11369585011066471632
      X( 5)=   0.14597371465489567682
      X( 6)=   0.17809688236761733390
      X( 7)=   0.21003131046056591064
      X( 8)=   0.24174315616383866556
      X( 9)=   0.27319881259104774468
      X(10)=   0.30436494435449495954
      X(11)=   0.33520852289262397655
      X(12)=   0.36569686147231213885
      X(13)=   0.39579764982890709712
      X(14)=   0.42547898840729897474
      X(15)=   0.45470942216774136446
      X(16)=   0.48345797392059470382
      X(17)=   0.51169417715466604391
      X(18)=   0.53938810832435567233
      X(19)=   0.56651041856139533470
      X(20)=   0.59303236477757022282
      X(21)=   0.61892584012546672523
      X(22)=   0.64416340378496526886
      X(23)=   0.66871831004391424358
      X(24)=   0.69256453664216964528
      X(25)=   0.71567681234896561582
      X(26)=   0.73803064374439816819
      X(27)=   0.75960234117664555964
      X(28)=   0.78036904386743123629
      X(29)=   0.80030874413913884180
      X(30)=   0.81940031073792957139
      X(31)=   0.83762351122818502758
      X(32)=   0.85495903343459936363
      X(33)=   0.87138850590929436968
      X(34)=   0.88689451740241818933
      X(35)=   0.90146063531585023110
      X(36)=   0.91507142312089592706
      X(37)=   0.92771245672230655266
      X(38)=   0.93937033975275308073
      X(39)=   0.95003271778443564022
      X(40)=   0.95968829144874048809
      X(41)=   0.96832682846326217918
      X(42)=   0.97593917458513455843
      X(43)=   0.98251726356301274934
      X(44)=   0.98805412632962202890
      X(45)=   0.99254390032376081654
      X(46)=   0.99598184298720747465
      X(47)=   0.99836437586317963722
      X(48)=   0.99968950388322870559
      W( 1)=   0.03255061449236316962
      W( 2)=   0.03251611871386883307
      W( 3)=   0.03244716371406427668
      W( 4)=   0.03234382256857594104
      W( 5)=   0.03220620479403026124
      W( 6)=   0.03203445623199267876
      W( 7)=   0.03182875889441101874
      W( 8)=   0.03158933077072719007
      W( 9)=   0.03131642559686137819
      W(10)=   0.03101033258631386231
      W(11)=   0.03067137612366917839
      W(12)=   0.03029991542082762553
      W(13)=   0.02989634413632842385
      W(14)=   0.02946108995816795100
      W(15)=   0.02899461415055528410
      W(16)=   0.02849741106508543861
      W(17)=   0.02797000761684838950
      W(18)=   0.02741296272602931385
      W(19)=   0.02682686672559184485
      W(20)=   0.02621234073567250055
      W(21)=   0.02557003600534944960
      W(22)=   0.02490063322248370695
      W(23)=   0.02420484179236479915
      W(24)=   0.02348339908592633665
      W(25)=   0.02273706965832950717
      W(26)=   0.02196664443874448477
      W(27)=   0.02117293989219144572
      W(28)=   0.02035679715433347898
      W(29)=   0.01951908114014518992
      W(30)=   0.01866067962741165898
      W(31)=   0.01778250231604547316
      W(32)=   0.01688547986424539715
      W(33)=   0.01597056290256253144
      W(34)=   0.01503872102699521608
      W(35)=   0.01409094177231515264
      W(36)=   0.01312822956696188190
      W(37)=   0.01215160467108866759
      W(38)=   0.01116210209983888144
      W(39)=   0.01016077053500880978
      W(40)=   0.00914867123078384552
      W(41)=   0.00812687692569928101
      W(42)=   0.00709647079115442616
      W(43)=   0.00605854550423662775
      W(44)=   0.00501420274292825661
      W(45)=   0.00396455433844564804
      W(46)=   0.00291073181793626202
      W(47)=   0.00185396078894924657
      W(48)=   0.00079679206555731759


      DO I=1,48
      XI(I)=-X(49-I)
      WI(I)=W(49-I)
      XI(I+48)=X(I)
      WI(I+48)=W(I)
      enddo
      RETURN
      END

      SUBROUTINE gauss64(xi,wi)
      implicit none
      double precision X(32),W(32),xi(100),wi(100)
      integer i

      W( 1)=  0.0486909570091397
      W( 2)=  0.0485754674415034
      W( 3)=  0.0483447622348030
      W( 4)=  0.0479993885964583
      W( 5)=  0.0475401657148303
      W( 6)=  0.0469681828162100
      W( 7)=  0.0462847965813144
      W( 8)=  0.0454916279274181
      W( 9)=  0.0445905581637566
      W(10)=  0.0435837245293235
      W(11)=  0.0424735151236536
      W(12)=  0.0412625632426235
      W(13)=  0.0399537411327203
      W(14)=  0.0385501531786156
      W(15)=  0.0370551285402400
      W(16)=  0.0354722132568824
      W(17)=  0.0338051618371416
      W(18)=  0.0320579283548516
      W(19)=  0.0302346570724025
      W(20)=  0.0283396726142595
      W(21)=  0.0263774697150547
      W(22)=  0.0243527025687109
      W(23)=  0.0222701738083833
      W(24)=  0.0201348231535302
      W(25)=  0.0179517157756973
      W(26)=  0.0157260304760247
      W(27)=  0.0134630478967186
      W(28)=  0.0111681394601311
      W(29)=  0.0088467598263639
      W(30)=  0.0065044579689784
      W(31)=  0.0041470332605625
      W(32)=  0.0017832807216964

      X( 1)=0.0243502926634244
      X( 2)=0.0729931217877990
      X( 3)=0.1214628192961206
      X( 4)=0.1696444204239928
      X( 5)=0.2174236437400071
      X( 6)=0.2646871622087674
      X( 7)=0.3113228719902110
      X( 8)=0.3572201583376681
      X( 9)=0.4022701579639916
      X(10)=0.4463660172534641
      X(11)=0.4894031457070530
      X(12)=0.5312794640198946
      X(13)=0.5718956462026340
      X(14)=0.6111553551723933
      X(15)=0.6489654712546573
      X(16)=0.6852363130542333
      X(17)=0.7198818501716109
      X(18)=0.7528199072605319
      X(19)=0.7839723589433414
      X(20)=0.8132653151227975
      X(21)=0.8406292962525803
      X(22)=0.8659993981540928
      X(23)=0.8893154459951141
      X(24)=0.9105221370785028
      X(25)=0.9295691721319396
      X(26)=0.9464113748584028
      X(27)=0.9610087996520538
      X(28)=0.9733268277899110
      X(29)=0.9833362538846260
      X(30)=0.9910133714767443
      X(31)=0.9963401167719553
      X(32)=0.9993050417357722

      DO I=1,32
      XI(I)=-X(33-I)
      WI(I)=W(33-I)
      XI(I+32)=X(I)
      WI(I+32)=W(I)
      enddo

      return
      end

      SUBROUTINE gauss32(xi,wi)
      implicit none
      double precision X(16),W(16),xi(100),wi(100)
      integer i

      X(1)=0.048307665687738316235
      X(2)=0.144471961582796493485
      X(3)=0.239287362252137074545
      X(4)=0.331868602282127649780
      X(5)=0.421351276130635345364
      X(6)=0.506899908932229390024
      X(7)=0.587715757240762329041
      X(8)=0.663044266930215200975
      X(9)=0.732182118740289680387
      X(10)=0.794483795967942406963
      X(11)=0.849367613732569970134
      X(12)=0.896321155766052123965
      X(13)=0.934906075937739689171
      X(14)=0.964762255587506430774
      X(15)=0.985611511545268335400
      X(16)=0.997263861849481563545
      W(1)=0.096540088514727800567
      W(2)=0.095638720079274859419
      W(3)=0.093844399080804565639
      W(4)=0.091173878695763884713
      W(5)=0.087652093004403811143
      W(6)=0.083311924226946755222
      W(7)=0.078193895787070306472
      W(8)=0.072345794108848506225
      W(9)=0.065822222776361846838
      W(10)=0.058684093478535547145
      W(11)=0.050998059262376176196
      W(12)=0.042835898022226680657
      W(13)=0.034273862913021433103
      W(14)=0.025392065309262059456
      W(15)=0.016274394730905670605
      W(16)=0.007018610009470096600
      DO I=1,16
      XI(I)=-X(17-I)
      WI(I)=W(17-I)
      XI(I+16)=X(I)
      WI(I+16)=W(I)
      enddo
      RETURN
      END

      SUBROUTINE gauss24(xi,wi)
      implicit none
      double precision X(12),W(12),xi(100),wi(100)
      integer i

      W(1)=0.1279381953467522
      W(2)=0.1258374563468283
      W(3)=0.1216704729278034
      W(4)=0.1155056680537256
      W(5)=0.1074442701159656
      W(6)=0.0976186521041139
      W(7)=0.0861901615319533
      W(8)=0.0733464814110803
      W(9)=0.0592985849154368
      W(10)=0.0442774388174198
      W(11)=0.0285313886289337
      W(12)=0.0123412297999872

      X(1)=0.0640568928626056
      X(2)=0.1911188674736163
      X(3)=0.3150426796961634
      X(4)=0.4337935076260451
      X(5)=0.5454214713888396
      X(6)=0.6480936519369755
      X(7)=0.7401241915785544
      X(8)=0.8200019859739029
      X(9)=0.8864155270044011
      X(10)=0.9382745520027328
      X(11)=0.9747285559713095
      X(12)=0.9951872199970213

      DO I=1,12
      XI(I)=-X(13-I)
      WI(I)=W(13-I)
      XI(I+12)=X(I)
      WI(I+12)=W(I)
      enddo
      RETURN
      END

      SUBROUTINE gauss16(xi,wi)
      implicit none
      double precision X(16),W(16),xi(100),wi(100)
      integer i

      X(1)=0.0950125098376374
      X(2)=0.2816035507792589
      X(3)=0.4580167776572274
      X(4)=0.6178762444026438
      X(5)=0.7554044083550030
      X(6)=0.8656312023878318
      X(7)=0.9445750230732326
      X(8)=0.9894009349916499
      W(1)=0.1894506104550685
      W(2)=0.1826034150449236
      W(3)=0.1691565193950025
      W(4)=0.1495959888165767
      W(5)=0.1246289712555339
      W(6)=0.0951585116824928
      W(7)=0.0622535239386479
      W(8)=0.0271524594117541

      DO I=1,8
      XI(I)=-X(9-I)
      WI(I)=W(9-I)
      XI(I+8)=X(I)
      WI(I+8)=W(I)
      enddo

      RETURN
      END

      SUBROUTINE gauss12(xi,wi)
      implicit none
      double precision X(12),W(12),xi(100),wi(100)
      integer i

       W(1)= 0.2491470458134028
       W(2)= 0.2334925365383548
       W(3)= 0.2031674267230659
       W(4)= 0.1600783285433462
      W(5)= 0.1069393259953184
      W(6)= 0.0471753363865118

      X(1)=0.1252334085114689
      X(2)=0.3678314989981802
      X(3)=0.5873179542866175
      X(4)=0.7699026741943047
      X(5)=0.9041172563704749
      X(6)=0.9815606342467192

      DO I=1,6
      XI(I)=-X(7-I)
      WI(I)=W(7-I)
      XI(I+6)=X(I)
      WI(I+6)=W(I)
      enddo

      RETURN
      END

      SUBROUTINE gauss8(xi,wi)
      implicit none
      double precision X(8),W(8),xi(100),wi(100)
      integer i

      X(1)=0.1834346424956498
      X(2)=0.5255324099163290
      X(3)=0.7966664774136267
      X(4)=0.9602898564975363

      W(1)=0.3626837833783620
      W(2)=0.3137066458778873
      W(3)=0.2223810344533745
      W(4)=0.1012285362903763

      DO I=1,4
      XI(I)=-X(5-I)
      WI(I)=W(5-I)
      XI(I+4)=X(I)
      WI(I+4)=W(I)
      enddo

      RETURN
      END

      SUBROUTINE gauss6(xi,wi)
      implicit none
      double precision X(6),W(6),xi(100),wi(100)
      integer i

      W(1)=0.3607615730481386
      W(2)=0.4679139345726910
      W(3)=0.1713244923791704

      X(1)=0.6612093864662645
      X(2)=0.2386191860831969
      X(3)=0.9324695142031521

      DO I=1,3
      XI(I)=-X(4-I)
      WI(I)=W(4-I)
      XI(I+3)=X(I)
      WI(I+3)=W(I)
      enddo

      RETURN
      END

      SUBROUTINE gauss4(xi,wi)
      implicit none
      double precision X(6),W(6),xi(100),wi(100)
      integer i

      W(1)=0.6521451548625461
      W(2)=0.3478548451374538
      X(1)=0.3399810435848563
      X(2)=0.8611363115940526

      DO I=1,2
      XI(I)=-X(3-I)
      WI(I)=W(3-I)
      XI(I+2)=X(I)
      WI(I+2)=W(I)
      enddo

      RETURN
      END
