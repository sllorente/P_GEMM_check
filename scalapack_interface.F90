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
         complex(kind=dp), intent(in) :: alpha
         integer, intent(in) :: desca(*)
         complex(kind=dp), intent(in) :: aa(desca(lld_), *)
         integer, intent(in) :: ia, ja
         integer, intent(in) :: descb(*)
         complex(kind=dp), intent(in) :: bb(descb(lld_), *)
         integer, intent(in) :: ib, jb
         complex(kind=dp), intent(in) :: beta
         integer, intent(in) :: descc(*)
         complex(kind=dp), intent(inout) :: cc(descb(lld_), *)
         integer, intent(in) :: ic, jc
      end subroutine pzgemm
      
      subroutine pdgemm(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
            & bb, ib, jb, descb, beta, cc, ic, jc, descc)
         import
         character, intent(in) :: transa, transb
         integer, intent(in) :: mm, nn, kk
         real(kind=dp), intent(in) :: alpha
         integer, intent(in) :: desca(*)
         real(kind=dp), intent(in) :: aa(desca(lld_), *)
         integer, intent(in) :: ia, ja
         integer, intent(in) :: descb(*)
         real(kind=dp), intent(in) :: bb(descb(lld_), *)
         integer, intent(in) :: ib, jb
         real(kind=dp), intent(in) :: beta
         integer, intent(in) :: descc(*)
         real(kind=dp), intent(inout) :: cc(descb(lld_), *)
         integer, intent(in) :: ic, jc
      end subroutine pdgemm

      subroutine pcgemm(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
            & bb, ib, jb, descb, beta, cc, ic, jc, descc)
         import
         character, intent(in) :: transa, transb
         integer, intent(in) :: mm, nn, kk
         complex(kind=sp), intent(in) :: alpha
         integer, intent(in) :: desca(*)
         complex(kind=sp), intent(in) :: aa(desca(lld_), *)
         integer, intent(in) :: ia, ja
         integer, intent(in) :: descb(*)
         complex(kind=sp), intent(in) :: bb(descb(lld_), *)
         integer, intent(in) :: ib, jb
         complex(kind=sp), intent(in) :: beta
         integer, intent(in) :: descc(*)
         complex(kind=sp), intent(inout) :: cc(descb(lld_), *)
         integer, intent(in) :: ic, jc
      end subroutine pcgemm
      
      subroutine psgemm(transa, transb, mm, nn, kk, alpha, aa, ia, ja, desca, &
            & bb, ib, jb, descb, beta, cc, ic, jc, descc)
         import
         character, intent(in) :: transa, transb
         integer, intent(in) :: mm, nn, kk
         real(kind=sp), intent(in) :: alpha
         integer, intent(in) :: desca(*)
         real(kind=sp), intent(in) :: aa(desca(lld_), *)
         integer, intent(in) :: ia, ja
         integer, intent(in) :: descb(*)
         real(kind=sp), intent(in) :: bb(descb(lld_), *)
         integer, intent(in) :: ib, jb
         real(kind=sp), intent(in) :: beta
         integer, intent(in) :: descc(*)
         real(kind=sp), intent(inout) :: cc(descb(lld_), *)
         integer, intent(in) :: ic, jc
      end subroutine psgemm
   end interface

end module scalapack_interface
