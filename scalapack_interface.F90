module scalapack_interface

   use iso_fortran_env, only: sp=>real32, dp=>real64
   
   implicit none

   integer, parameter :: block_cyclic_2d = 1, dlen_ = 9, dtype_ = 1, &
      ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ = 6, &
      rsrc_ = 7, csrc_ = 8, lld_ = 9 


   interface
      integer function indxl2g(indxloc, nb, iproc, isrcproc, nprocs)
         integer, intent(in) :: indxloc, nb, iproc, isrcproc, nprocs
      end function indxl2g

      integer function indxg2l(indxglob, nb, iproc, isrcproc, nprocs)
      integer, intent(in) :: indxglob, nb, iproc, isrcproc, nprocs
      end function indxg2l
      
      integer function indxg2p(indxglob, nb, iproc, isrcproc, nprocs)
         integer, intent(in) :: indxglob, nb, iproc, isrcproc, nprocs
      end function indxg2p

      integer function numroc(n, size_block, myproc, srcproc, nprocs)
         integer, intent(in) :: n, size_block, myproc, srcproc, nprocs
      end function numroc

      subroutine blacs_get(ictxt, what, val)
         integer, intent(in) :: ictxt, what
         integer, intent(out) :: val
      end subroutine blacs_get

      subroutine blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)
         integer, intent(in) :: ictxt, nprow, npcol
         integer, intent(out) :: myrow, mycol
      end subroutine blacs_gridinfo

      subroutine blacs_gridinit(ictxt, order, nprow, npcol)
         integer, intent(inout) :: ictxt
         character, intent(in) :: order
         integer, intent(in) :: nprow, npcol
      end subroutine blacs_gridinit

      subroutine blacs_pinfo(id, nproc)
         integer, intent(out) :: id, nproc
      end subroutine blacs_pinfo

      subroutine blacs_barrier(ictxt, scope)
         integer, intent(in) :: ictxt
         character, intent(in) :: scope
      end subroutine blacs_barrier
      
      subroutine blacs_gridexit(ictxt)
         integer, intent(in) :: ictxt
      end subroutine blacs_gridexit

      subroutine blacs_exit(cont)
         integer, intent(in) :: cont
      end subroutine blacs_exit

      subroutine descinit(desc, mm, nn, mb, nb, irsrc, icsrc, ictxt, lld, info)
         integer, intent(out) :: desc(*)
         integer, intent(in) :: mm, nn, mb, nb, irsrc, icsrc, ictxt, lld
         integer, intent(out) :: info
      end subroutine descinit

      subroutine pzgemm(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
            & bb, ib, jb, descb, beta, cc, ic, jc, descc)
         import
         character, intent(in) :: transa, transb
         integer, intent(in) :: mm, nn, kk
         integer, intent(in) :: desca(*), ia, ja
         integer, intent(in) :: descb(*), ib, jb
         integer, intent(in) :: descc(*), ic, jc
         complex(kind=dp), intent(in) :: alpha, beta
         complex(kind=dp), intent(in) :: aa(*), bb(*)
         complex(kind=dp), intent(inout) :: cc(*)
      end subroutine pzgemm
      
      subroutine pdgemm(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
            & bb, ib, jb, descb, beta, cc, ic, jc, descc)
         import
         character, intent(in) :: transa, transb
         integer, intent(in) :: mm, nn, kk
         integer, intent(in) :: desca(*), ia, ja
         integer, intent(in) :: descb(*), ib, jb
         integer, intent(in) :: descc(*), ic, jc
         real(kind=dp), intent(in) :: alpha, beta
         real(kind=dp), intent(in) :: aa(*), bb(*)
         real(kind=dp), intent(inout) :: cc(*)
      end subroutine pdgemm

      subroutine pcgemm(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
            & bb, ib, jb, descb, beta, cc, ic, jc, descc)
         import
         character, intent(in) :: transa, transb
         integer, intent(in) :: mm, nn, kk
         integer, intent(in) :: desca(*), ia, ja
         integer, intent(in) :: descb(*), ib, jb
         integer, intent(in) :: descc(*), ic, jc
         complex(kind=sp), intent(in) :: alpha, beta
         complex(kind=sp), intent(in) :: aa(*), bb(*)
         complex(kind=sp), intent(inout) :: cc(*)
      end subroutine pcgemm
      
      subroutine psgemm(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
            & bb, ib, jb, descb, beta, cc, ic, jc, descc)
         import
         character, intent(in) :: transa, transb
         integer, intent(in) :: mm, nn, kk
         integer, intent(in) :: desca(*), ia, ja
         integer, intent(in) :: descb(*), ib, jb
         integer, intent(in) :: descc(*), ic, jc
         real(kind=sp), intent(in) :: alpha, beta
         real(kind=sp), intent(in) :: aa(*), bb(*)
         real(kind=sp), intent(inout) :: cc(*)
      end subroutine psgemm

      subroutine psger(mm, nn, alpha, xx, ix, jx, descx, incx, &
            yy, iy, jy, descy, incy, aa, ia, ja, desca)
         import
         integer, intent(in) :: mm, nn
         integer, intent(in) :: descx(*), ix, jx, incx
         integer, intent(in) :: descy(*), iy, jy, incy
         integer, intent(in) :: desca(*), ia, ja
         real(kind=sp), intent(in) :: alpha
         real(kind=sp), intent(in) :: xx(*)
         real(kind=sp), intent(in) :: yy(*)
         real(kind=sp), intent(inout) :: aa(*)
      end subroutine psger

      subroutine pdger(mm, nn, alpha, xx, ix, jx, descx, incx, &
            yy, iy, jy, descy, incy, aa, ia, ja, desca)
         import
         integer, intent(in) :: mm, nn
         integer, intent(in) :: descx(*), ix, jx, incx
         integer, intent(in) :: descy(*), iy, jy, incy
         integer, intent(in) :: desca(*), ia, ja
         real(kind=dp), intent(in) :: alpha
         real(kind=dp), intent(in) :: xx(*)
         real(kind=dp), intent(in) :: yy(*)
         real(kind=dp), intent(inout) :: aa(*)
      end subroutine pdger

      subroutine pcgeru(mm, nn, alpha, xx, ix, jx, descx, incx, &
            yy, iy, jy, descy, incy, aa, ia, ja, desca)
         import
         integer, intent(in) :: mm, nn
         integer, intent(in) :: descx(*), ix, jx, incx
         integer, intent(in) :: descy(*), iy, jy, incy
         integer, intent(in) :: desca(*), ia, ja
         complex(kind=sp), intent(in) :: alpha
         complex(kind=sp), intent(in) :: xx(*)
         complex(kind=sp), intent(in) :: yy(*)
         complex(kind=sp), intent(inout) :: aa(*)
      end subroutine pcgeru

      subroutine pzgeru(mm, nn, alpha, xx, ix, jx, descx, incx, &
            yy, iy, jy, descy, incy, aa, ia, ja, desca)
         import
         integer, intent(in) :: mm, nn
         integer, intent(in) :: descx(*), ix, jx, incx
         integer, intent(in) :: descy(*), iy, jy, incy
         integer, intent(in) :: desca(*), ia, ja
         complex(kind=dp), intent(in) :: alpha
         complex(kind=dp), intent(in) :: xx(*)
         complex(kind=dp), intent(in) :: yy(*)
         complex(kind=dp), intent(inout) :: aa(*)
      end subroutine pzgeru

   end interface

end module scalapack_interface
