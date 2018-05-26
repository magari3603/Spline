!  Spline.f
!
!  関数:
!  Spline - コンソール・アプリケーションのエントリーポイント。
!

!****************************************************************************
!
!  プログラム: Spline
!
!  目的:  コンソール・アプリケーションのエントリーポイント。
!
!****************************************************************************

      program Spline
c  ******************* Parameters ************************
      implicit real*8 (a-h,o-z)
      parameter(NonZeroMax=5000)
      dimension XY(0:NonZeroMax,2),x(0:NonZeroMax),V(0:NonZeroMax)
      dimension abcd(0:NonZeroMax,4),h(0:NonZeroMax)
      dimension dt(0:NonZeroMax)
      character*50 outdt,inptfl
      
      real*8, allocatable :: Ab(:,:),OutData(:,:)
      allocate (Ab(NonZeroMax,NonZeroMax+1),OutData(NonZeroMax*100,2))
      
c *********************************************************      
      
      inptfl='50A1300.txt'
      outdt='New50A1300.txt'
      iDiveNum=2
      
c ********************************************************* 
      
      do 5 i=1,NonZeroMax*100
        OutData(i,1)=0.D0
        OutData(i,2)=0.D0
5     continue
      
      do 200 i=0,1000
        XY(i,1)=0
        XY(i,2)=0
        V(i)=0
        x(i)=0
        h(i)=0
        do 300 j=1,1000
          if(i.ne.0)then
            Ab(i,j)=0.D0
          endif
300     continue
200   continue
      open (2,file=inptfl,status='unknown')
        read(2,*)EddyCoef,HysCoef
        read(2,*)WeightPer,nb
        do 10 i=0,nb
          read(2,*)XY(i,1),XY(i,2)
10      continue
     
      iMatrixNo=nb
c ------------------ Making matrix ------------------
      
      do 100 i=0,iMatrixNo-1
        h(i)=XY(i+1,1)-XY(i,1)
100   continue
      do 400 i=1,iMatrixNo-1
        V(i)=6.D0*(XY(i+1,2)-XY(i,2))/h(i)
        V(i)=V(i)-6.D0*(XY(i,2)-XY(i-1,2))/h(i-1)
400   continue
      
      do 500 i=1,iMatrixNo-1
        do 600 j=1,iMatrixNo-1
          if(i.eq.j)then
            Ab(i,j)=2*(h(i-1)+h(i))
          elseif((i+1).eq.j)then
            Ab(i,j)=h(i)
          elseif((i-1).eq.j)then
             Ab(i,j)=h(i-1)
          endif
600     continue
500   continue
      
      do 900 i=1,iMatrixNo-1
         Ab(i,iMatrixNo)=V(i)
900   continue
C ----------- Gauss Jordan method -----------
      do 1100 k=1,iMatrixNo-1
        pivot=1.D0/Ab(k,k)
        
        do 1200 j=k+1,iMatrixNo
          Ab(k,j)=Ab(k,j)*pivot
1200    continue
        do 1300 i=1,iMatrixNo-1
          if(i.ne.k)then
            do 1400 j=k+1,iMatrixNo
              Ab(i,j)=Ab(i,j)-Ab(i,k)*Ab(k,j)
1400        continue
          endif
1300    continue
1100  continue
          
c ---------------------- x ------------------
      do 1500 i=1,iMatrixNo-1
         x(i)=Ab(i,iMatrixNo)
1500  continue
      
c  -------------- aj,bj,cj,dj ---------------
      do 1600 i=0,iMatrixNo-1
        abcd(i,1)=((1.D0/6.D0)*(x(i+1)-x(i)))/(XY(i+1,1)-XY(i,1))
        abcd(i,2)=(1.D0/2.D0)*x(i)
        abcd(i,3)=(XY(i+1,2)-XY(i,2))/(XY(i+1,1)-XY(i,1))-
     +    (1.D0/6.D0)*(2*x(i)+x(i+1))*(XY(i+1,1)-XY(i,1))
        abcd(i,4)=XY(i,2)
1600  continue
c ------------------ make list ---------------
      do 1700 i=0,iMatrixNo-1
        dt(i)=abs(XY(i+1,1)-XY(i,1))
        dt(i)=dt(i)/real(iDiveNum)
1700  continue
      iCou=0
      do 1900 i=0,iMatrixNo-1
        do 1800 j=1,iDiveNum
          iCou=iCou+1
          OutData(iCou+1,1)=OutData(iCou,1)+dt(i)
1800    continue
1900  continue     
      
c ---------------------------------------------
      iCou=1
      do 2000 i=0,iMatrixNo-1
        do 2100 j=1,iDiveNum
          OutData(iCou,2)=abcd(i,1)*(OutData(iCou,1)-XY(i,1))**3
     +    +abcd(i,2)*(OutData(iCou,1)-XY(i,1))**2
     +    +abcd(i,3)*(OutData(iCou,1)-XY(i,1))+abcd(i,4)
          iCou=iCou+1
2100    continue
2000  continue
      
      open (10,file=outdt,status='unknown')
      write(10,10000)EddyCoef,HysCoef
      write(10,20000)WeightPer,nb*iDiveNum
      do 2200 i=1,iDiveNum*iMatrixNo
        write(10,10000)OutData(i,1),OutData(i,2)
2200  continue
        write(10,10000)OutData(i,1),XY(iMatrixNo,2)
      close(10)
c 
      
10000 format(2ES16.8)
20000 format(ES16.8,i6)
1000  format(10ES12.3)
      
      
      continue
      end

