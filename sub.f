c     *******************************************************
c     *               boundary condition                    *
c     *******************************************************
      subroutine boundary(nn,x,y,a,u,
     +                    xrgn1,xrgn2,yrgn1,yrgn2,
     +                    sss)
      dimension x(nn),y(nn)
      dimension a(nn,nn),u(nn)
c
      bc=0.
c
      do 1000 i=1,nn
        if ((x(i).lt.xrgn1+sss).or.(x(i).gt.xrgn2-sss)
     x  .or.(y(i).lt.yrgn1+sss).or.(y(i).gt.yrgn2-sss)) then
           do 2000 j=1,nn
             a(i,j)=0.
             a(j,i)=0.
2000       continue
           a(i,i)=1.0
           u(i)=0.
        end if
1000  continue
      return
      end
c  ****************************************************
c  *                   making matrix                  *
c  ****************************************************
      subroutine matrix(ncur,nn,ne,x,y,ijk,mat,matno,
     +                  per,sigma,a,u,
     +                  ncrpm,ejj)
      dimension x(nn),y(nn),ijk(ne,3),mat(ne),per(matno),sigma(matno)
      dimension ncrpm(ncur),ejj(ncur)
      dimension a(nn,nn),u(nn)
c
c     local valuable
      real bb(3),cc(3)
c
c  ---------- clear matix -------------
      do 12304 i=1,nn
        do 12304 j=1,nn
          a(i,j)=0.
12304 continue
      do 12307 i=1,nn
        u(i)=0.
12307 continue
c  ----------  making matrix ---------------
      do 12970 i=1,ne
        bb(1)=y(ijk(i,2))-y(ijk(i,3))
        bb(2)=y(ijk(i,3))-y(ijk(i,1))
        bb(3)=y(ijk(i,1))-y(ijk(i,2))
        cc(1)=x(ijk(i,3))-x(ijk(i,2))
        cc(2)=x(ijk(i,1))-x(ijk(i,3))
        cc(3)=x(ijk(i,2))-x(ijk(i,1))
        aa=bb(2)*cc(3)-bb(3)*cc(2)
        dd=aa/2
        do 12960 l=1,3
          do 12950 m=1,3
            te=(bb(l)*bb(m)+cc(l)*cc(m))/(per(mat(i))*4*dd)
            a(ijk(i,l),ijk(i,m))=a(ijk(i,l),ijk(i,m))+te
12950     continue
12960   continue
12970 continue
c  ------------- power matrix ------------------
      do 13210 i=1,ne
	    if (mat(i).le.10) goto 13210
c   current number
        n=mat(i)-10
c
        bb(1)=y(ijk(i,2))-y(ijk(i,3))
        bb(2)=y(ijk(i,3))-y(ijk(i,1))
        bb(3)=y(ijk(i,1))-y(ijk(i,2))
        cc(1)=x(ijk(i,3))-x(ijk(i,2))
        cc(2)=x(ijk(i,1))-x(ijk(i,3))
        cc(3)=x(ijk(i,2))-x(ijk(i,1))
c
        aa=bb(2)*cc(3)-bb(3)*cc(2)
        dd=aa/2
        do 13220 j=1,3
          k=ijk(i,j)
          u(k)=u(k)+dd*ejj(n)*ncrpm(n)/3.
13220   continue
13210 continue
      return
      end

c     *******************************************************
c     *                    gauss method                     *
c     *******************************************************
      subroutine gauss(nn,a,u)
      dimension a(nn,nn),u(nn)
c
      do 1000 m=1,nn
        write(*,100)m
100     format(i6)
        p=a(m,m)
        p1=1/p
        do 2000 j=m,nn
          a(m,j)=a(m,j)*p1
2000    continue
        u(m)=u(m)*p1
        do 1000 i=1,nn
          if (i.eq.m) goto 1000
          q=a(i,m)
          if (q.eq.0) goto 1000
          do 3000 j=m,nn
            a(i,j)=a(i,j)-q*a(m,j)
3000      continue
          u(i)=u(i)-q*u(m)
1000  continue
      return
      end
c     *******************************************************
c     *                    output data                      *
c     *******************************************************
      subroutine out(ncur,nn,ne,x,y,ijk,mat,matno,per,sigma,u,ncrpm,
     +               xrgn1,xrgn2,yrgn1,yrgn2,
     +               perx1,perx2,pery1,pery2,
     +               curx1,curx2,cury1,cury2)
      dimension x(nn),y(nn),ijk(ne,3),mat(ne),per(matno),sigma(matno)
      dimension u(nn)
      dimension ncrpm(ncur)
      dimension curx1(ncur),curx2(ncur),cury1(ncur),cury2(ncur)
c
c     jFileType=-1
c     jFileType=-1: Linear Static
c     jFileType=0: j-omega
c
c     jFileType>0: linear problem
c     jFileType<0: nonlinear problem
c
      ndam=0
      ndam2=1
      jFileType=-1
      write(1,300)nn,ne,ncur,ndam,ndam,matno,jFileType,ndam,ndam2
      do 20 i=1,nn
        write(1,150)ndam,x(i),y(i)
20    continue
      do 30 i=1,ne
        write(1,300)mat(i),ijk(i,1),ijk(i,2),ijk(i,3)
30    continue
      do 45 i=1,matno
        write(1,200)per(i),sigma(i)
45    continue
      do 48 i=1,ncur
        write(1,300)ncrpm(i)
48    continue
      do 50 i=1,nn
        write(1,100)u(i)
50    continue
c
      write(1,400)xrgn1,xrgn2,yrgn1,yrgn2
      write(1,400)perx1,perx2,pery1,pery2
      write(1,300)ncur
      do 60 i=1,ncur
        write(1,400)curx1(i),curx2(i),cury1(i),cury2(i)
60    continue
100   format(e16.7)
150   format(i6,2e16.7)
200   format(2e16.7)
300   format(9i8)
400   format(4e16.7)
      return
      end
