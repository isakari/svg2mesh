!> svgファイルに描いた絵を境界要素に変換する
program main
  use module_mesher
  implicit none

  !> bin に置いた svg の名前
  !> svg のサイズは 150px x 150px 程度にしておく
  character(len=16) :: fn_svg = "test.svg"

  !> svg のサイズ (px x px)
  integer,parameter :: n(2) = [200, 200]

  ! local variables
  integer :: i, j, k
  character(len=2) :: ctmp(8)
  real(8),allocatable :: chi0(:,:) !< 特性関数
  real(8),allocatable :: chi(:,:) !< 特性関数にガウシアンフィルタをかけたもの
  real(8),allocatable :: lsf(:,:) !< レベルセット関数
  
  ! svg を png に変換
  call system(trim("convert " // trim(fn_svg) // " tmp.png"))

  ! png を txt に変換
  call system("convert tmp.png -type Grayscale tmp.txt")

  ! txt を整形
  call system('sed s/","/", "/g tmp.txt > tmp2.txt')
  call system('sed s/"("/"( "/g tmp2.txt > tmp.txt')

  ! find indcolour
  open(1,file="tmp.txt")
  read(1,*) !空読み
  
  ! load the characteristic function
  allocate(chi0(n(1),n(2)))
  chi0(:,:)=0
  do j=1,n(2)
     do i=1,n(1)
        read(1,*) (ctmp(k),k=1,7)
        if(ctmp(7).ne."#F")then
           chi0(i,n(2)-j+1)=1
        end if
     end do
  end do
  close(1)

  ! 特性にガウシアンフィルタを 10回 かける
  allocate(chi(n(1),n(2)))
  do k=1,10
     chi=chi0
     do j=2,n(2)-1
        do i=2,n(1)-1
           chi(i,j)=(4.d0*chi0(i,j)&
                +2.d0*chi0(i  ,j+1)&
                +2.d0*chi0(i  ,j-1)&
                +2.d0*chi0(i-1,  j)&
                +2.d0*chi0(i+1,  j)&
                +chi0(i-1,j-1)&
                +chi0(i-1,j+1)&
                +chi0(i+1,j-1)&
                +chi0(i+1,j+1))/16.d0
        end do
     end do
     chi0=chi
  end do

  
  ! 特性関数をレベルセット関数に
  allocate(lsf(0:n(1),0:n(2)))
  lsf(0,0)=chi(1,1)
  lsf(n(1),0)=chi(n(1),1)
  lsf(0,n(2))=chi(1,n(2))
  lsf(n(1),n(2))=chi(n(1),n(2))
  do i=1,n(1)-1
     lsf(i,0)=(chi(i,1)+chi(i+1,1))*0.5d0
     lsf(i,n(2))=(chi(i,n(2))+chi(i+1,n(2)))*0.5d0
  end do
  do i=1,n(2)-1
     lsf(0,i)=(chi(1,i)+chi(1,i+1))*0.5d0
     lsf(n(1),i)=(chi(n(1),i)+chi(n(1),i+1))*0.5d0
  end do
  do j=1,n(2)-1
     do i=1,n(1)-1
        lsf(i,j)=(chi(i,j)+chi(i,j+1)+chi(i+1,j)+chi(i+1,j+1))*0.25d0
     end do
  end do
  lsf(:,:)=lsf(:,:)*2.d0-1.d0

  lsf=lsf-epsilon(1.d0)
  
  open(1,file="lsf.txt")
  do j=0,n(2)
     do i=0,n(1)
        write(1,*) dble(i),dble(j),lsf(i,j)
     end do
  end do
  close(1)

  call lsf2elm(n(1),n(2),[-1.d0,-1.d0],[1.d0,1.d0],1,trim("lsf.txt"),trim("mesh.txt"),trim("mesh_gp.txt"))
  
  ! ! 中間ファイルを削除
  call system("rm -f *tmp*")
  call system("rm -f lsf.txt orig.gp *mod a.out")

end program main
