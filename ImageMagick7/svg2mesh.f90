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
  real(8) :: png
  character(len=8) :: ctmp(9)
  real(8),allocatable :: chi(:,:) !< 特性関数
  real(8),allocatable :: lsf(:,:) !< レベルセット関数
  
  ! svg を png に変換
  call system(trim("magick " // trim(fn_svg) // " test.png"))

  ! png を txt に変換
  call system("magick test.png -type Grayscale tmp1.txt")

  ! txt を整形
  call system("awk -f fix_gray.awk tmp1.txt > tmp2.txt")
  call system('sed s/","/" "/g tmp2.txt > tmp3.txt') ! カンマを削除
  call system('sed s/"("/"( "/g tmp3.txt > tmp4.txt') ! 始カッコの直後にスペースを挿入
  call system('sed s/")"/" )"/g tmp4.txt > tmp5.txt') ! 終カッコの直前にスペースを挿入

  ! find indcolour
  open(1,file="tmp5.txt")
  read(1,*) !空読み
  
  ! load the characteristic function
  allocate(chi(n(1),n(2)))
  chi(:,:)=0
  do j=1,n(2)
     do i=1,n(1)
        read(1,*) (ctmp(k),k=1,8), png
        chi(i,n(2)-j+1)=(100.d0-png)/100.d0
     end do
  end do
  close(1)

  open(1,file="chi.txt")
  open(2,file="dom.txt")
  do j=1,n(2)
     do i=1,n(1)
        write(1,*) i,j,chi(i,j)
        if(chi(i,j)>0)then
           write(2,*) (i-100)/100.d0, (j-100)/100.d0
        end if
     end do
     write(1,*) 
  end do
  close(1)
  close(2)
        
  ! ! 特性関数をレベルセット関数に
  allocate(lsf(0:n(1)-2,0:n(2)-2))
  do j=1,n(2)-1
     do i=1,n(1)-1
        lsf(i-1,j-1)=(chi(i,j)+chi(i,j+1)+chi(i+1,j)+chi(i+1,j+1))*0.25d0
     end do
  end do
  lsf(:,:)=2*lsf(:,:)-1.d0
  
  open(1,file="phi.txt")
  do j=0,n(2)-2
     do i=0,n(1)-2
        write(1,*) dble(i),dble(j),lsf(i,j)
     end do
  end do
  close(1)

  call lsf2elm(n(1)-2,n(2)-2,[-1.d0,-1.d0],[1.d0,1.d0],1,trim("phi.txt"),trim("mesh.txt"),trim("mesh_gp.txt"))
  
  ! ! 中間ファイルを削除
  call system("rm -f *tmp*")
  call system("rm -f lsf.txt orig.gp *mod a.out")

end program main
