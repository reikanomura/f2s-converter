! ----------------------------------------------------------
module fem_types

  implicit none

! -------------------------------
! 表面メッシュの構造体
! -------------------------------
  type :: surface_mesh
     integer :: node       ! 節点数
     integer :: nelem       ! 三角形要素数
     double precision, pointer :: xy(:,:)    ! 節点座標 (3, mnode)
     integer, pointer :: ncall(:,:) ! グローバル要素節点番号 (3, nelem)
     integer, pointer :: nc(:,:) ! ローカル要素節点番号 (3, nelem)
  end type surface_mesh
end module
! ----------------------------------------------------------
program fluid2solid_force
  use fem_types
  implicit none

  type(surface_mesh) :: fluid   ! 流体側表面
  type(surface_mesh) :: solid   ! 固体側表面

  ! 流体側節点トラクション t_f [N/m^2], force***.txtの中身
  double precision, allocatable :: tf(:,:)

  ! 固体側節点トラクション t_s [N/m^2] （補間される結果）
  double precision, allocatable :: ts(:,:)

  ! 固体側の点PのIDを入力すると
  ! 点Pを取り囲む流体三角形要素の経常関数が出力される配列
  double precision, allocatable :: Nabc(:,:)

  ! 固体側節点の等価節点力 [N] (3 成分, node)
  double precision, allocatable :: force_s(:,:)     

  ! 固体側節点の合力 [N] (3 成分, node)
  double precision, allocatable :: F_global(:,:)     
  double precision :: xP(3), xi(3), xj(3), xk(3)
  double precision :: lambda(3)

  integer, allocatable :: melem(:), mnode(:)

  integer :: i, j, k, n, m
  integer :: na, nb, nc, ni, nj, nk
  integer :: ielem, inode
  integer :: istep

  character(50) :: astep, filename


! 流体側の連成面情報の読み込みe_f,  x_f,  (a,b,c)
open(10, file = 'mangrove-fld-surface.txt', status='old')
  read(10,*) fluid%nelem
  write(6,*) fluid%nelem
  allocate(fluid%ncall(3,fluid%nelem))
  read(10,*) ((fluid%ncall(j,i), j = 1, 3), i = 1, fluid%nelem)
  read(10,*) fluid%node
  write(6,*) fluid%node
  allocate(fluid%xy(3,fluid%node), fluid%nc(3,fluid%nelem))
  read(10,*) ((fluid%xy(j,i), j = 1, 3), i=1, fluid%node)
  read(10,*) ((fluid%nc(j,i), j = 1, 3), i=1, fluid%nelem)
close(10)

! 固体側の連成面情報の読み込みe_s,  x_s,  (i,j,k)
open(20, file = 'mangrove-sld-f-surface.txt', status='old')
  read(20,*) solid%nelem
  write(6,*) solid%nelem
  allocate(solid%nc(3,solid%nelem), solid%ncall(3,solid%nelem))
  read(20,*) ((solid%ncall(j,i), j = 1, 3), i = 1, solid%nelem)
  read(20,*) solid%node
  write(6,*) solid%node
  allocate(solid%xy(3,solid%node))
  read(20,*) ((solid%xy(j,i), j = 1, 3), i=1, solid%node)
  read(20,*) ((solid%nc(j,i), j = 1, 3), i=1, solid%nelem)
close(20)

allocate(tf(3,fluid%node))
allocate(ts(3,solid%node))
allocate(force_s(3,solid%node),F_global(3,solid%node))
allocate(Nabc(3,solid%node))
allocate(melem(solid%node), mnode(solid%node), source=0)

! 固体側の節点が流体側のどの三角形要素内にあるのか判定
do n = 1, solid%node
    xP(:) = solid%xy(:,n)

    if (xP(3) .le. 1.0e-5)then

      !xp(3) = z が0.0に近いということは，底面に接している節点なので
      !等価節点力（外力条件）ではなく，強制変位0.0を３方向に与える節点
      !よってielem = 0を入れて特別に分類する
      melem(n) = 0

    else

      ! xPがどの三角形要素の中にあるか確かめるサブルーチン
      ! ielemにはxPを含んでいる流体要素のローカルIDが入る
      call find_fluid_triangle_for_point(fluid%node, fluid%nelem, fluid%nc, fluid%xy, &
          			       xP, ielem, inode, lambda)

      if (ielem == 0) then
         ! 見つからなかった場合の処理
         ! mnode(n)に，コピー元になる流体側節点ID inodeを格納するyy
         mnode(n) = inode
         ! melem(n)には，フラッグで-1（無効値,NaN等の意味）を入れる．
         melem(n) = -1
         write(6,'(a,i7,a,i7)') 'Solid node:', n, ' <-->  Fluid node:',inode!, xP(1), xP(2), xP(3)
      else

         Nabc(1,n) = lambda(1)  ! N_a = u 
         Nabc(2,n) = lambda(2)  ! N_b = v
         Nabc(3,n) = lambda(3)  ! N_c = 1 - u -v
         melem(n) = ielem

      end if
    end if
end do

!write(6,*) minval(melem), maxval(melem)

call output_surf (solid%node, solid%nelem, solid%xy, solid%nc, melem)

do istep = 1, 300
 

  !==============================================
  ! cal_force_p.f90で出力したforce****.txtを指定します
  !write(filename,'(a,i4.4,a4)') './H0100-3.43/force', istep, '.txt'
  !open(30, file=filename, status='old', action='read')
  ! ***  
  !     tfを読み込むコードを書く 
  !                           ***
  !close(30)
   !==============================================

  do m = 1, solid%nelem
   ni = solid%nc(1,m)
   nj = solid%nc(2,m)
   nk = solid%nc(3,m)
   xi(:) = solid%xy(:,ni)
   xj(:) = solid%xy(:,nj)
   xk(:) = solid%xy(:,nk)
   !==============================================
   ! 武田さん 書き換え部分
    !*** 固体三角形メッシュe_sの面積A_sを求める***
   !==============================================

     do i = 1, 3
       !==============================================
        ! 武田さん 書き換え部分
        !# melem(solid%nc(i,m)) > 0のとき
        !    固体節点位置solid%xy(:,n)での流体力ts(n)を
        !    N_a, N_b, N_cで補間しトラクションベクトルを求める
        !# melem(solid%nc(i,m)) =-1のとき
        !    固体節点nには，流体節点inode(solid%nc(i,m))の値を
        !　　コピーしてトラクションベクトルを求める
       !==============================================

       !==============================================
        ! 武田さん 書き換え部分
        ! tsが求められたら等価節点力force_sを求める 
       !==============================================

       !==============================================
       ! 武田さん 書き換え部分
       ! 等価節点力からグローバル合力F_globalを求める
       !==============================================
     end do

   
     !==============================================
     ! 松原先生に協力してもらう部分
     ! tree-****-bc.cml形式(****はistepなど)で
     ! 流体からの作用力を各ステップでの境界条件として出力
     ! 等価節点力からグローバル合力F_glonalを求める
     ! melem(n)  ==  0 : 固体節点nには強制変位(0, 0, 0)を与える
     ! melem(n) .ne. 0 : 固体節点nには，等価節点力F_global=(fx, fy, fz)を
     !                 与える
     !==============================================
   

  end do
  !write(6,*) istep 
end do



end program
!-----------------------------------------------------------
  subroutine find_fluid_triangle_for_point(node, nelem, elem, xy, xP, ielem, inode, lambda)
!-----------------------------------------------------------
! 与えられた点 P が、どの流体側三角形の内部にあるか
! 探索するサブルーチン
! 入力:
!   mesh_f  : 流体側表面メッシュ
!   xP(3)   : 固体側節点座標
!
! 出力:
!   ielem : 見つかった流体側三角形要素番号（0なら見つからず）
!   lambda  : バリセンター (3)  = (lambda_a, lambda_b, lambda_c)
    !type(surface_mesh), intent(in) :: mesh_f
    integer, intent(in) :: nelem, node, elem(3,nelem)
    double precision, intent(in)  :: xP(3), xy(3,node)
    integer,  intent(out) :: ielem, inode
    double precision, intent(out) :: lambda(3)

    integer :: m, na, nb, nc
    double precision :: xA(3), xB(3), xC(3)
    double precision :: u, v
    double precision :: reps, dx, dy, dz
    logical  :: inside

    ielem = 0
    lambda(:) = 0.0d0


    do m = 1, nelem
       na = elem(1,m)
       nb = elem(2,m)
       nc = elem(3,m)

       xA(:) = xy(:,na)
       xB(:) = xy(:,nb)
       xC(:) = xy(:,nc)

       call barycentric_3d(xA, xB, xC, xP, u, v, inside)

       if (inside) then
          lambda(1) = 1.0d0 - u - v   ! lambda_a
          lambda(2) = u               ! lambda_b
          lambda(3) = v               ! lambda_c
          ielem   = m
          return
       end if
    end do

    ! ここに来たら、どの三角形にも属さなかったことになるので
    ! 最も距離の近い流体節点inodeを突き止めて，その値を
    ! 固体節点にコピーすることにする
    inode = 0
    reps = 1.0d0
    do n = 1, node

       dx = xP(1) - xy(1,n)
       dy = xP(2) - xy(2,n)
       dz = xP(3) - xy(3,n)
       rr = dsqrt(dx*dx + dy*dy + dz*dz)
       if (rr < reps) then
          reps = rr
          inode = n
       end if

    end do
    !write(6,*) reps
  end subroutine find_fluid_triangle_for_point

!-----------------------------------------------------------
  subroutine barycentric_3d(xA, xB, xC, xP, u, v, inside)
!-----------------------------------------------------------
  ! 3D 三角形 ABC と点 P に対して、 u,v を求める
  ! （lambda_b = u, lambda_c = v, lambda_a = 1-u-v）
  ! 戻り値 inside = .true. なら三角形内部または辺上
    double precision, intent(in)  :: xA(3), xB(3), xC(3), xP(3)
    double precision, intent(out) :: u, v
    logical, intent(out)  :: inside

    double precision :: v0(3), v1(3), v2(3)
    double precision :: dot00, dot01, dot02, dot11, dot12
    double precision :: denom
    double precision :: lambda_a, lambda_b, lambda_c
    double precision, parameter :: eps = 1.0e-14

    v0(:) = xB(:) - xA(:)
    v1(:) = xC(:) - xA(:)
    v2(:) = xP(:) - xA(:)

    dot00 = sum(v0*v0)
    dot01 = sum(v0*v1)
    dot02 = sum(v0*v2)
    dot11 = sum(v1*v1)
    dot12 = sum(v1*v2)

    denom = dot00*dot11 - dot01*dot01

    if (abs(denom) < eps) then
       ! 三角形がつぶれている
       u = 0.0d0
       v = 0.0d0
       inside = .false.
       return
    end if

    u = ( dot11*dot02 - dot01*dot12 ) / denom
    v = ( dot00*dot12 - dot01*dot02 ) / denom

    lambda_a = 1.0d0 - u - v
    lambda_b = u
    lambda_c = v

    if (abs(lambda_a) < eps) lambda_a = 0.0d0
    if (abs(lambda_b) < eps) lambda_b = 0.0d0
    if (abs(lambda_c) < eps) lambda_c = 0.0d0


    if ( lambda_a >= -eps .and. lambda_b >= -eps .and. lambda_c >= -eps ) then
       inside = .true.
    else
       inside = .false.
    end if
  end subroutine barycentric_3d


 ! ---------------------------------------------------------------------
subroutine output_vtk (respath, istep, node, nelem, xy, nc, force, vec)
! ---------------------------------------------------------------------
implicit none
character(50), intent(in) :: respath
integer, intent(in) :: istep, node, nelem
integer, intent(in) :: nc(3,nelem)
double precision, intent(in) :: xy(3,node)
double precision, intent(in) :: force(3,node), vec(3,node)

character(80) :: astep

integer :: i

write(astep, '(a,i4.4,a)') trim(adjustl(respath)), istep,'.vtk'
open(50, file = astep)
  write(50,'(a)')'# vtk DataFile Version 2.0'
  write(50,'(a)')'Unstructured Grid Example'
  write(50,'(a)')'ASCII'
  write(50,'(a)')'     '

  write(50,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(50,'(a,i10,a)')'POINTS',node,' float'
  do i = 1, node
    write(50,'(3e15.6)')xy(1,i), xy(2,i), xy(3,i)
  end do

  write(50,'(a,2i10)')'CELLS   ',nelem, (nelem)*4
  do i = 1, nelem
    write(50,'(5i7)')3,nc(1,i)-1, nc(2,i)-1, nc(3,i)-1
  end do

  write(50,'(a,i10)')'CELL_TYPES    ',nelem
  do i = 1, nelem
    write(50,'(i7)')5
  end do

  write(50,'(a)')'     '
  write(50,'(a,i10)')'POINT_DATA   ',node
  write(50,'(a)')'VECTORS force float'
  do i = 1, node
    write(50,'(3e15.6)')force(1,i), force(2,i), force(3,i)
  end do

  write(50,'(a)')'VECTORS unit_nom float'
  do i = 1, node
    write(50,'(3e15.6)')vec(1,i), vec(2,i), vec(3,i)
  end do

close(50)


end subroutine
! ---------------------------------------------------------------------
subroutine output_surf (node, nelem, xy, nc, melem)
! ---------------------------------------------------------------------
implicit none
integer, intent(in) :: node, nelem
integer, intent(in) :: nc(3,nelem), melem(node)
double precision, intent(in) :: xy(3,node)

integer :: i

open(50, file = 'out_surf.vtk')
  write(50,'(a)')'# vtk DataFile Version 2.0'
  write(50,'(a)')'Unstructured Grid Example'
  write(50,'(a)')'ASCII'
  write(50,'(a)')'     '

  write(50,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(50,'(a,i10,a)')'POINTS',node,' float'
  do i = 1, node
    write(50,'(3e15.6)')xy(1,i), xy(2,i), xy(3,i)
  end do

  write(50,'(a,2i10)')'CELLS   ',nelem, (nelem)*4
  do i = 1, nelem
    write(50,'(5i7)')3,nc(1,i)-1, nc(2,i)-1, nc(3,i)-1
  end do

  write(50,'(a,i10)')'CELL_TYPES    ',nelem
  do i = 1, nelem
    write(50,'(i7)')5
  end do

  write(50,'(a)')'     '
  write(50,'(a,i10)')'POINT_DATA   ',node
  write(50,'(a)')'SCALARS melem float'
  write(50,'(a)')'LOOKUP_TABLE default'
  do i = 1, node
    write(50,'(E15.5)') dble(melem(i))
  end do
close(50)


end subroutine

