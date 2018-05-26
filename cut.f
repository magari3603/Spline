c     ++++++++++++++++ cutting program for macec ++++++++++++++++
      subroutine cut(nx,ny,ncur,nn,ne,x,y,ijk,mat,matno,
     +               per,sigma,
     +               xrgn1,xrgn2,yrgn1,yrgn2,
     +               perx1,perx2,pery1,pery2,
     +               curx1,curx2,cury1,cury2,ncrpm,ejj,
     +               sss,cx,cy)
c
      dimension x(nn),y(nn),ijk(ne,3),mat(ne),per(matno),sigma(matno)
      dimension curx1(ncur),curx2(ncur),cury1(ncur),cury2(ncur)
      dimension ncrpm(ncur),ejj(ncur)
c
c     local valuable
      dimension cx(nx+1),cy(ny+1)
c  **************   fem cutting program ************
      pai=.3141592653589794          e+01       
c     resion
      xrgn1=-1.0
      xrgn2=1.0
      yrgn1=-1.0
      yrgn2=1.0
c
c     position of current
      ncrpm(1)=1
      curx1(1)=0.3
      curx2(1)=0.4
      cury1(1)=-0.1
      cury2(1)=0.1
      ncrpm(2)=-1
      curx1(2)=-0.4
      curx2(2)=-0.3
      cury1(2)=-0.1
      cury2(2)=0.1
c
c     position of iron
      perx1=-0.2
      perx2=0.2
      pery1=-0.4
      pery2=0.4
c
c     --------- cut of sigma region --------------
      cx(1)=xrgn1
      cx(2)=-0.4
      cx(3)=-0.3
      cx(4)=-0.2
      cx(5)=-0.1
      cx(6)=-0.0
      cx(7)=0.1
      cx(8)=0.2
      cx(9)=0.3
      cx(10)=0.4
      cx(11)=xrgn2
c
      cy(1)=yrgn1
      cy(2)=-0.4
      cy(3)=-0.3
      cy(4)=-0.2
      cy(5)=-0.1
      cy(6)=-0.0
      cy(7)=0.1
      cy(8)=0.2
      cy(9)=0.3
      cy(10)=0.4
      cy(11)=yrgn2
c  *************** model data  *****************
c     current	eii[A.T]
      eii=1.0
c     current density ejj[A/m*m]
      do 1500 i=1,ncur 
        ss=(curx2(i)-curx1(i))*(cury2(i)-cury1(i))
        ejj(i)=eii/ss	 
1500  continue
      pai=3.141592654
c
c     mat(i)=1: air
c     mat(i)=2: iron
c     mat(i)>10: current
c
      per(1)=4*pai*.0000001
      per(2)=1000*4*pai*.0000001
      per(3)=4*pai*.0000001
      per(4)=4*pai*.0000001
      do 1510 i=11,20
        per(i)=4*pai*.0000001
1510  continue
c
      sigma(1)=0.
      sigma(2)=0.
      sigma(3)=0.
      sigma(4)=0.
      do 1520 i=11,20
        sigma(i)=0.
1520  continue

c  ************* cutting finite elements ***********
      node=0
      do 10420 i=1,nx+1
        do 10410 j=1,ny+1
          node=node+1
          x(node)=cx(i)
          y(node)=cy(j)
10410 continue
10420 continue
      iele=0
      do 10980 i=1,nx
        do 10970 j=1,ny
        if ((nx/2-i+SSS)*(ny/2-j+SSS).ge.0) then
            iele=iele+1
            ijk(iele,1)=j+(ny+1)*(i-1)
            ijk(iele,2)=j+(ny+1)*i
            ijk(iele,3)=j+1+(ny+1)*i
            iele=iele+1
            ijk(iele,1)=j+1+(ny+1)*(i-1)
            ijk(iele,2)=j+(ny+1)*(i-1)
            ijk(iele,3)=j+1+(ny+1)*i
           else
            iele=iele+1
            ijk(iele,1)=j+(ny+1)*(i-1)
            ijk(iele,2)=j+(ny+1)*i
            ijk(iele,3)=j+1+(ny+1)*(i-1)
            iele=iele+1
            ijk(iele,1)=j+(ny+1)*i
            ijk(iele,2)=j+1+(ny+1)*i
            ijk(iele,3)=j+1+(ny+1)*(i-1)
           end if
10970   continue
10980 continue
c
c  ******inmat(ele) ******************************
      do 11193 i=1,ne
        mat(i)=1
11193 continue
        do 12183 j=1,ne
          do 12161 l=1,3
            x1=x(ijk(j,l))
            y1=y(ijk(j,l))
            if (x1.lt.perx1-sss) goto 12183
            if (y1.lt.pery1-sss) goto 12183
            if (x1.gt.perx2+sss) goto 12183
            if (y1.gt.pery2+sss) goto 12183
12161     continue
          mat(j)=2
12183   continue
c
      do 13182 i=1,ncur
        do 13183 j=1,ne
          do 13161 l=1,3
            x1=x(ijk(j,l))
            y1=y(ijk(j,l))
            if (x1.lt.curx1(i)-sss) goto 13183
            if (y1.lt.cury1(i)-sss) goto 13183
            if (x1.gt.curx2(i)+sss) goto 13183
            if (y1.gt.cury2(i)+sss) goto 13183
13161     continue
          mat(j)=10+i
13183   continue
13182 continue
      return
      end

