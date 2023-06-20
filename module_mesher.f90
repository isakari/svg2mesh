module module_mesher
   implicit none
   private
   public lsf2elm

   real(8), allocatable :: an(:,:), at(:,:)

contains

   subroutine improve_mesh(nnode,nelement,nodepos,element,idiv,dmin,dmax,xsize,ysize,eps1,file_elm,file_elm_gp)
      implicit none
      real(8) :: dmin(2), dmax(2)
      integer :: xsize, ysize

      real(8), parameter :: preset_shortestedge=0.2d0
      real(8), parameter :: dFerguson=0.5d0

      integer :: idiv

      integer nnode,nelement
      integer element(2,nelement)
      real(8) nodepos(2,nnode)
      real(8) eps1
      character(len=8) file_elm
      character(len=11) :: file_elm_gp

      real(8), allocatable::elml(:) !elml(i): 要素iの長さ
      integer, allocatable::elm(:,:) !elm(i,j): j番節点のi番要素
      real(8), allocatable::atnd(:,:) !atnd(i,j): j番節点における接線のx_i成分

      integer i,j,k,l,k1,k2

      real(8) xcnt(2)

      real(8) au
      real(8) f1,f2,f3,f4
      real(8) anorm
      real(8) xi,xj,yi,yj

      integer, allocatable::ndiv(:)
      real(8) shortestedge,longestedge

      real(8), allocatable::nodeposnew(:,:)
      integer, allocatable::elementnew(:,:)

      real(8) :: scalex, scaley

      !短い辺を削除
103   continue

      allocate(elml(nelement))
      do i=1,nelement !要素長
         elml(i)=sqrt(dot_product(&
            nodepos(:,element(2,i))-nodepos(:,element(1,i)),&
            nodepos(:,element(2,i))-nodepos(:,element(1,i))))
      end do

      allocate(elm(2,nnode)) !elm(i,j): j番節点のi番要素を返す配列を作る
      do i=1,nnode
         do j=1,nelement
            if(i.eq.element(1,j))then
               elm(1,i)=j
               goto 100
            else if(i.eq.element(2,j))then
               elm(2,i)=j
               goto 100
            end if
100         continue
         end do
      end do

      allocate(atnd(2,nnode)) !atnd(i,j): j番節点における接線のx_i成分
      do i=1,nnode
         atnd(:,i)=at(:,elm(1,i))/elml(elm(1,i))+at(:,elm(2,i))/elml(elm(2,i))
         !atnd(:,i)=at(:,elm(1,i))*elml(elm(1,i))+at(:,elm(2,i))*elml(elm(2,i))
         anorm=1.d0/sqrt(dot_product(atnd(:,i),atnd(:,i)))
         atnd(:,i)=atnd(:,i)*anorm
      end do

!!$  do i=1,nnode !接線チェック
!!$! gnuplot でplot "fort.102" w lp, "fort.88"w lpでチェックできる
!!$     write(102,*) nodepos(:,i)
!!$     write(102,*) nodepos(:,i)+atnd(:,i)
!!$     write(102,*)
!!$     write(102,*)
!!$  end do

      au=0.5d0
      f1=2.d0*au**3-3.d0*au**2+1.d0
      f2=-2.d0*au**3+3.d0*au**2
      f3=au**3-2.d0*au**2+au
      f4=au**3-au**2
      do i=1,nelement
         if(elml(i).le.preset_shortestedge)then
            !Ferguson曲線を作成
            xcnt(:)=nodepos(:,element(1,i))*f1&
               +nodepos(:,element(2,i))*f2&
               +atnd(:,element(1,i))*f3&
               +atnd(:,element(2,i))*f4
            nodepos(:,element(1,i))=xcnt(:) !要素iの両端の点をferguson曲線の中点に移動させる。
            nodepos(:,element(2,i))=xcnt(:) !要素iの両端の点をferguson曲線の中点に移動させる。
         end if
      end do


      ! 重複する節点の統合 ↑で2点を1点にmergeしたので、重複が現れる。これを取り除く。
      do i=1,nnode-1
         xi=nodepos(1,i)
         yi=nodepos(2,i)
         do j=i+1,nnode
            xj=nodepos(1,j)
            yj=nodepos(2,j)
            if(abs(xi-xj)<eps1.and.abs(yi-yj)<eps1) then ! iとjが一致したら
               do k=j,nnode-1
                  nodepos(:,k)=nodepos(:,k+1) ! j+1 以降の点も一つずつ前にずらす
               end do
               do k=1,nelement
                  do l=1,2
                     if(element(l,k).eq.j) then
                        element(l,k)=i ! 要素の節点jをiに置換し
                     else if(element(l,k).gt.j) then
                        element(l,k)=element(l,k)-1 ! j+1以降の点を一つずつ前にずらす
                     end if
                  end do
               end do
               nnode=nnode-1 ! 節点が一つ減った
            end if
         end do
      end do

      ! 同じ点からなる要素を削除
101   continue
      do i=1,nelement
         if(element(1,i).eq.element(2,i)) then
            do j=i,nelement-1
               element(:,j)=element(:,j+1)
            end do
            nelement=nelement-1
            goto 101
         end if
      end do


      ! ! 直線を削除
112   continue
      do i=1,nelement
         do j=i+1,nelement
            if((element(1,i).eq.element(2,j)).and.(element(2,i).eq.element(1,j))) then

               k1=min(element(1,j),element(2,j))
               k2=max(element(1,j),element(2,j))
               do k=k1+1,nnode
                  if((k.gt.k1).and.(k.lt.k2))then
                     nodepos(:,k-1)=nodepos(:,k)
                  else if(k.gt.k2) then
                     nodepos(:,k-2)=nodepos(:,k)
                  end if
               end do
               do k=1,nelement
                  do l=1,2
                     if((element(l,k).gt.k1).and.(element(l,k).lt.k2)) then
                        element(l,k)=element(l,k)-1
                     else if(element(l,k).gt.k2) then
                        element(l,k)=element(l,k)-2
                     end if
                  end do
               end do
               nnode=nnode-2

               do k=i+1,nelement
                  if((k.gt.i).and.(k.lt.j)) then
                     element(:,k-1)=element(:,k)
                  else if(k.gt.j) then
                     element(:,k-2)=element(:,k)
                  end if
               end do
               nelement=nelement-2
               goto 112
            end if
         end do
      end do

      deallocate(at)
      allocate(at(2,nelement))
      do i=1,nelement
         anorm=(nodepos(1,element(2,i))-nodepos(1,element(1,i)))**2+(nodepos(2,element(2,i))-nodepos(2,element(1,i)))**2
         anorm=1d0/sqrt(anorm)
         at(1,i)=(nodepos(1,element(2,i))-nodepos(1,element(1,i)))*anorm
         at(2,i)=(nodepos(2,element(2,i))-nodepos(2,element(1,i)))*anorm
      enddo
      deallocate(elm,atnd)

      deallocate(elml)
      allocate(elml(nelement))
      do i=1,nelement !要素長
         elml(i)=sqrt(dot_product(&
            nodepos(:,element(2,i))-nodepos(:,element(1,i)),&
            nodepos(:,element(2,i))-nodepos(:,element(1,i))))
      end do

      shortestedge=minval(elml)
      longestedge=maxval(elml)
!!$  write(*,*) "# Mesh with short edge is eliminated"
!!$  write(*,*) '## nnode=',nnode !統合後の節点数
!!$  write(*,*) '## nelement=',nelement !統合後の要素数
!!$  write(*,*) "## (Longest edge)",longestedge
!!$  write(*,*) "## (Shortest edge)=",shortestedge
!!$  write(*,*) "## (Longest edge)/(Shortest edge)=",longestedge/shortestedge

      if(shortestedge.le.preset_shortestedge) then
         deallocate(elml)
         goto 103
      end if

!!$  do i=1,nelement !要素改善後の境界要素の絵を描く
!!$     do j=1,2
!!$        write(87,*) nodepos(:,element(j,i))
!!$     end do
!!$     write(87,*)
!!$  end do

      allocate(ndiv(nelement))
      do i=1,nelement
         ndiv(i)=int(dble(idiv)*elml(i)/shortestedge) !各要素を何等分するかを決定
         if(ndiv(i).le.1) ndiv(i)=2
      end do
      allocate(elementnew(2,sum(ndiv)))
      allocate(nodeposnew(2,sum(ndiv)))
      nodeposnew(:,1:nnode)=nodepos(:,1:nnode)

      allocate(atnd(2,nnode)) !atnd(i,j): j番節点の接線ベクトルのx_i成分
      allocate(elm(2,nnode)) !elm(i,j): j番節点のi番要素を返す配列を作る
      do i=1,nnode
         do j=1,nelement
            if(i.eq.element(1,j))then
               elm(1,i)=j
               goto 102
            elseif(i.eq.element(2,j))then
               elm(2,i)=j
               goto 102
            end if
102         continue
         end do
      end do
      do i=1,nnode
         !atnd(:,i)=at(:,elm(1,i))*elml(elm(1,i))+at(:,elm(2,i))*elml(elm(2,i))
         atnd(:,i)=at(:,elm(1,i))/elml(elm(1,i))+at(:,elm(2,i))/elml(elm(2,i))
         anorm=1.d0/sqrt(dot_product(atnd(:,i),atnd(:,i)))
         atnd(:,i)=atnd(:,i)*anorm
      end do
      ! do i=1,nnode !接線チ�����ック
      !    !! gnuplot でplot "fort.104" w lp, "fort.88"w lpでチェックできる
      !    write(104,*) nodepos(:,i)
      !    write(104,*) nodepos(:,i)+atnd(:,i)
      !    write(104,*)
      !    write(104,*)
      ! end do
      ! close(104)
      do i=1,nelement
         !     ferguson曲線を作成
         do j=1,ndiv(i)-1
            au=dble(j)/dble(ndiv(i))
            !        if(ndiv(i).le.2*idiv)then
            if(elml(i).le.dFerguson&
               .or.dot_product(atnd(:,element(1,i)),at(:,i)).lt.0.8d0&
               .or.dot_product(atnd(:,element(2,i)),at(:,i)).lt.0.8d0&
               )then
               f1=1.d0-au
               f2=au
               f3=0.d0
               f4=0.d0
            else
               f1=2.d0*au**3-3.d0*au**2+1.d0
               f2=-2.d0*au**3+3.d0*au**2
               f3=au**3-2.d0*au**2+au
               f4=au**3-au**2
            end if
            xcnt(:)=nodepos(:,element(1,i))*f1&
               +nodepos(:,element(2,i))*f2&
               +atnd(:,element(1,i))*f3&
               +atnd(:,element(2,i))*f4
            !nodeを追加
            nnode=nnode+1
            nodeposnew(:,nnode)=xcnt(:)
            elementnew(1,sum(ndiv(1:i-1))+j)=nnode-1
            elementnew(2,sum(ndiv(1:i-1))+j)=nnode
         end do
         !要素を追加
         elementnew(1,sum(ndiv(1:i-1))+1)=element(1,i)
         elementnew(1,sum(ndiv(1:i)))=nnode
         elementnew(2,sum(ndiv(1:i)))=element(2,i)
      end do

      deallocate(elml)
      allocate(elml(sum(ndiv)))
      do i=1,sum(ndiv) !要素長を計算
         elml(i)=sqrt(dot_product(&
            nodeposnew(:,elementnew(2,i))-nodeposnew(:,elementnew(1,i)),&
            nodeposnew(:,elementnew(2,i))-nodeposnew(:,elementnew(1,i))))
      end do
      shortestedge=minval(elml)
      longestedge=maxval(elml)
!!$  write(*,*) '## nnode=',sum(ndiv)
!!$  write(*,*) '## nelement=',sum(ndiv)
!!$  write(*,*) "## (Longest edge)",longestedge
!!$  write(*,*) "## (Shortest edge)=",shortestedge
!!$  write(*,*) "## (Longest edge)/(Shortest edge)=",longestedge/shortestedge


      scalex=(dmax(1)-dmin(1))/dble(xsize)
      scaley=(dmax(2)-dmin(2))/dble(ysize)
      nodeposnew(1,:)=nodeposnew(1,:)*scalex+dmin(1)
      nodeposnew(2,:)=nodeposnew(2,:)*scaley+dmin(2)

      open(107,file=file_elm_gp)
      do i=1,sum(ndiv) !要素改善後の境界要素の絵を描く
         do j=1,2
            write(107,*) nodeposnew(:,elementnew(j,i)),elementnew(j,i)
         end do
         write(107,*)
      end do
      close(107)

      ! BEMに渡すため、メッシュデータをfileに書き込み
      open(1,file=file_elm)
      write(1,*) sum(ndiv)
      do i=1,sum(ndiv)
         write(1,*) i,nodeposnew(:,i)
      end do
      write(1,*) sum(ndiv)
      do i=1,sum(ndiv)
         write(1,*) i,elementnew(1,i),elementnew(2,i)
      end do
      close(1)

      deallocate(elml,ndiv)
   end subroutine improve_mesh

   !=================================================================================================
   subroutine lsf2elm(xsize,ysize,dmin,dmax,idiv,file_lsf,file_elm,file_elm_gp)
      implicit none
      integer, intent(in) :: xsize, ysize, idiv
      real(8), intent(in) :: dmin(2), dmax(2)
      character(len=7), intent(in) :: file_lsf
      character(len=8), intent(in) :: file_elm
      character(len=11),intent(in) :: file_elm_gp

      real(8), parameter :: eps1=1.0d-8

      integer itmp,ix,iy,i,j,k,l
      real(8) :: tmp
      real(8), allocatable::lsfunc(:,:)
      integer, allocatable::element(:,:)
      real(8), allocatable::nodepos(:,:)
      real(8), allocatable::refvec(:,:)
      !refvec:「voxelの中心から外側に向かう」or「レベルセット関数の値が最大の節点からvoxelの中心に向かう」ベクトル(cnt_nn=2のとき)
      !refvec:voxelの中心から外側に向かうベクトル (cnt_nn=4のとき)
      !refvec:固定設計領域から外側に向かうベクトル (cnt_z>0のとき)

      real(8) lsftmp(0:3)

      integer cnt_z,cnt_p,cnt_m
      real(8) au,av,f1,f2,ww

      integer nnode,nelement,cnt_max,cnt_nn

      real(8) lsfs
      real(8) xi,yi,xj,yj

      real(8) anorm

      integer imaxlsf,iminlsf
      real(8) dmaxlsf,dminlsf

      real(8), allocatable::elml(:) !elml(i): 要素iの長さ

      ! レベルセット関数の読み込み
      open(10,file=file_lsf)
      allocate(lsfunc(0:xsize,0:ysize))
      do iy=0,ysize
         do ix=0,xsize
            read(10,*) tmp,tmp,lsfunc(ix,iy)
            if(abs(lsfunc(ix,iy)).le.5.d-3) lsfunc(ix,iy)=sign(1.d0,lsfunc(ix,iy))*5.d-3
         end do
      end do
      close(10)

      allocate(nodepos(2,xsize*ysize*4),element(2,xsize*ysize*2),refvec(2,xsize*ysize*2))
      nodepos(:,:)=0.d0
      element(:,:)=0

      nnode=0
      nelement=0
      cnt_max=0
      do iy=0,ysize-1    !voxelの左下の節点のloop
         do ix=0,xsize-1 !voxelの左下の節点のloop
            cnt_nn=0
            lsftmp(0)=lsfunc(ix,  iy  ) !voxelの左下の節点のlevel set関数
            lsftmp(1)=lsfunc(ix+1,iy  ) !voxelの右下の節点のlevel set関数
            lsftmp(2)=lsfunc(ix+1,iy+1) !voxelの右上の節点のlevel set関数
            lsftmp(3)=lsfunc(ix,  iy+1) !voxelの左上の節点のlevel set関数

            cnt_z=0 !voxerl内でレベルセット関数がゼロの節点の数
            cnt_p=0 !voxerl内でレベルセット関数が正の節点の数
            cnt_m=0 !voxerl内でレベルセット関数が負の節点の数
            do i=0,3
               if(dabs(lsftmp(i)).le.0.d0) then
                  cnt_z=cnt_z+1
               elseif(lsftmp(i).ge.0.d0) then
                  cnt_p=cnt_p+1
               else
                  cnt_m=cnt_m+1
               end if
            end do

            if((cnt_z.eq.1.and.cnt_m.eq.3))cycle !当該voxelの中に要素は無い
            if(cnt_m.eq.4) cycle !これも上と同じような...?

            do i=0,3 !voxel内の節点loop
               if(i==0.or.i==3) au=dble(ix)
               if(i==1.or.i==2) au=dble(ix+1)
               if(i==0.or.i==1) av=dble(iy)
               if(i==2.or.i==3) av=dble(iy+1) !(au,av)が考えている節点の座標となる

               j=mod(i+1,4) !一つ先の(local)節点番号

               f1=lsftmp(i)
               f2=lsftmp(j)

               if (((f1>=0.d0).and.(f2<0.d0)).or.((f1<0.d0).and.(f2>=0.d0))) then ! 辺の両端の符号が違ったら
                  ww=abs(lsftmp(i)/(lsftmp(i)-lsftmp(j))) ! 節点iからの距離 lsftmp(i):x=(lsftmp(i)-lsftmp(j)):1
                  if( dabs( lsftmp(i) )<eps1 ) ww=0.0d0 !必要か？
                  if( ww>0.95d0 ) ww=0.95d0
                  if( ww<0.05d0 ) ww=0.05d0
                  ! 節点座標計算
                  if( i==0 ) au=au+ww
                  if( i==1 ) av=av+ww
                  if( i==2 ) au=au-ww
                  if( i==3 ) av=av-ww
                  ! 節点カウント
                  cnt_nn=cnt_nn+1 !voxel内の境界要素の節点の数
                  nnode=nnode+1
                  nodepos(1,nnode)=au
                  nodepos(2,nnode)=av
               end if
            end do

            ! 要素追加
            if(cnt_nn.eq.2)then
               nelement=nelement+1
               element(1,nelement)=nnode-1
               element(2,nelement)=nnode
               imaxlsf=0; dmaxlsf=lsftmp(0)
               iminlsf=0; dminlsf=lsftmp(0)
               do i=1,3
                  if(lsftmp(i).gt.dmaxlsf) then
                     imaxlsf=i
                     dmaxlsf=lsftmp(i)
                  end if
                  if(lsftmp(i).lt.dminlsf) then
                     iminlsf=i
                     dminlsf=lsftmp(i)
                  end if
               end do
               if(cnt_m.eq.1)then
                  refvec(1,nelement)=dble(min(mod(iminlsf,3),1))-0.5d0
                  refvec(2,nelement)=dble(iminlsf/2)-0.5d0
               elseif(cnt_p.eq.1)then
                  refvec(1,nelement)=0.5d0-dble(min(mod(imaxlsf,3),1))
                  refvec(2,nelement)=0.5d0-dble(imaxlsf/2)
               elseif(cnt_p.eq.2.and.cnt_m.eq.2)then
                  refvec(1,nelement)=dble(min(mod(iminlsf,3),1))-0.5d0
                  refvec(2,nelement)=dble(iminlsf/2)-0.5d0
               else
                  write(*,*) cnt_p,cnt_m,cnt_z,ix,iy,lsftmp(:)
                  stop
               end if
            else if( cnt_nn.eq.4 )then
               lsfs=lsftmp(0)+lsftmp(1)+lsftmp(2)+lsftmp(3)
               if(lsfs.ge.0.0d0.and.lsftmp(0).ge.0.0d0)then
                  nelement=nelement+1
                  element(1,nelement)=nnode-3
                  element(2,nelement)=nnode-2
                  refvec(1,nelement)=1.d0
                  refvec(2,nelement)=-1.d0
                  nelement=nelement+1
                  element(1,nelement)=nnode-1
                  element(2,nelement)=nnode
                  refvec(1,nelement)=-1.d0
                  refvec(2,nelement)=1.d0
               elseif(lsfs.le.0.d0.and.lsftmp(0).le.0.d0)then
                  nelement=nelement+1
                  element(1,nelement)=nnode-3
                  element(2,nelement)=nnode-2
                  refvec(1,nelement)=-1.d0
                  refvec(2,nelement)=1.d0
                  nelement=nelement+1
                  element(1,nelement)=nnode-1
                  element(2,nelement)=nnode
                  refvec(1,nelement)=1.d0
                  refvec(2,nelement)=-1.d0
               elseif(lsfs.ge.0.d0.and.lsftmp(0).le.0.d0)then
                  nelement=nelement+1
                  element(1,nelement)=nnode-3
                  element(2,nelement)=nnode
                  refvec(1,nelement)=-1.d0
                  refvec(2,nelement)=-1.d0
                  nelement=nelement+1
                  element(1,nelement)=nnode-2
                  element(2,nelement)=nnode-1
                  refvec(1,nelement)=1.d0
                  refvec(2,nelement)=1.d0
               elseif(lsfs.le.0.d0.and.lsftmp(0).ge.0.d0)then
                  nelement=nelement+1
                  element(1,nelement)=nnode-3
                  element(2,nelement)=nnode
                  refvec(1,nelement)=1.d0
                  refvec(2,nelement)=1.d0
                  nelement=nelement+1
                  element(1,nelement)=nnode-2
                  element(2,nelement)=nnode-1
                  refvec(1,nelement)=-1.d0
                  refvec(2,nelement)=-1.d0
               end if
            end if
         end do
      end do

!!$  write(*,*) "# nnode (before integration)=",nnode !統合前の節点数
!!$  write(*,*) "# nelement (before integration)",nelement !統合前の要素数


      ! 重複する節点の統合
      do i=1,nnode-1
         xi=nodepos(1,i)
         yi=nodepos(2,i)
         do j=i+1,nnode
            xj=nodepos(1,j)
            yj=nodepos(2,j)
            if(abs(xi-xj)<eps1.and.abs(yi-yj)<eps1) then ! iとjが一致したら
               do k=j,nnode-1
                  nodepos(:,k)=nodepos(:,k+1) ! j+1 以降の点も一つずつ前にずらす
               end do
               do k=1,nelement
                  do l=1,2
                     if(element(l,k).eq.j) then
                        element(l,k)=i ! 要素の節点jをiに置換し
                     else if(element(l,k).gt.j) then
                        element(l,k)=element(l,k)-1 ! j+1以降の点を一つずつ前にずらす
                     end if
                  end do
               end do
               nnode=nnode-1 ! 節点が一つ減った
            end if
         end do
      end do

      allocate(an(2,nelement),at(2,nelement))
      ! 法線を作る
      do i=1,nelement
         anorm=(nodepos(1,element(2,i))-nodepos(1,element(1,i)))**2+(nodepos(2,element(2,i))-nodepos(2,element(1,i)))**2
         anorm=1d0/sqrt(anorm)
         at(1,i)=(nodepos(1,element(2,i))-nodepos(1,element(1,i)))*anorm
         at(2,i)=(nodepos(2,element(2,i))-nodepos(2,element(1,i)))*anorm
         an(1,i)= at(2,i)
         an(2,i)=-at(1,i)
         if(an(1,i)*refvec(1,i)+an(2,i)*refvec(2,i).ge.0.d0)then ! 向きが揃っていなければ
            itmp=element(1,i) !要素番号を入れ替え
            element(1,i)=element(2,i)
            element(2,i)=itmp
            an(1,i)=-an(1,i)  !法線をひっくり返す
            an(2,i)=-an(2,i)
            at(1,i)=-at(1,i)  !接線もひっくり返す
            at(2,i)=-at(2,i)
         end if
      enddo

      ! 要素チェック
      allocate(elml(nelement))
      do i=1,nelement !要素長
         elml(i)=sqrt(dot_product(&
            nodepos(:,element(2,i))-nodepos(:,element(1,i)),&
            nodepos(:,element(2,i))-nodepos(:,element(1,i))))
      end do
!!$  write(*,*) "# Original mesh"
!!$  write(*,*) '## nnode=',nnode !統合後の節点数
!!$  write(*,*) '## nelement=',nelement !統合後の要素数
!!$  write(*,*) "## (Longest edge)=",maxval(elml)
!!$  write(*,*) "## (Shortest edge)=",minval(elml)
!!$  write(*,*) "## (Longest edge)/(Shortest edge)=",maxval(elml)/minval(elml)
      deallocate(elml)

      open(88,file="orig.gp")
      do i=1,nelement
         do j=1,2
            write(88,*) nodepos(:,element(j,i))
         end do
         write(88,*)
         write(88,*)
      end do
      close(88)

      ! 要素改善
      call improve_mesh(nnode,nelement,nodepos,element,idiv,dmin,dmax,xsize,ysize,eps1,file_elm,file_elm_gp)

      deallocate(lsfunc)
      deallocate(element,nodepos,refvec)
      deallocate(an,at)
   end subroutine lsf2elm

end module module_mesher
