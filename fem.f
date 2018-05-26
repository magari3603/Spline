      program femgss
c  **************************************************************
c
c  *                Analisys for 2D Static fields               *
c  
c  *              using Gauss Elimination  method               *
c
c  *                      named "femgss"                        *
c
c  *            written by K.Yamazaki  march,1992               *
c
c  **************************************************************
c
c  ******************************************************
c  *                control routine                     *
c  ******************************************************
      parameter(nx=10,ny=10,ncur=2)
      parameter(nn=(nx+1)*(ny+1),ne=2*nx*ny)
      parameter(sss=1.0e-7)
      parameter(matno=20)
c  ---Array of mesh
      dimension x(nn),y(nn),ijk(ne,3),mat(ne),per(matno),sigma(matno)
      dimension curx1(ncur),curx2(ncur),cury1(ncur),cury2(ncur)
      dimension ncrpm(ncur),ejj(ncur)
      dimension cx(nx+1),cy(ny+1)
c  ---Array of matrix
      dimension a(nn,nn),u(nn)
c
      character*50 outdt
c
      write(*,3)
3     format( ' electromagnetic output file name is  ')
      read(5,*)outdt
      open (2,file='itdata.txt',status='unknown')
      open (20,file='itval.txt',status='unknown')
      open (1,file=outdt,status='unknown')
c
            call cut(nx,ny,ncur,nn,ne,x,y,ijk,mat,matno,
     +               per,sigma,
     +               xrgn1,xrgn2,yrgn1,yrgn2,
     +               perx1,perx2,pery1,pery2,
     +               curx1,curx2,cury1,cury2,ncrpm,ejj,
     +               sss,cx,cy)
c
            call matrix(ncur,nn,ne,x,y,ijk,mat,matno,
     +                  per,sigma,a,u,
     +                  ncrpm,ejj)
c
            call boundary(nn,x,y,a,u,
     +                    xrgn1,xrgn2,yrgn1,yrgn2,
     +                    sss)
c
            call gauss(nn,a,u)
c
            call out(ncur,nn,ne,x,y,ijk,mat,matno,per,sigma,u,ncrpm,
     +               xrgn1,xrgn2,yrgn1,yrgn2,
     +               perx1,perx2,pery1,pery2,
     +               curx1,curx2,cury1,cury2)
c
10000   continue
      close(1)
      stop
      end
