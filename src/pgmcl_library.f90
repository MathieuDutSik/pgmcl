#include "pgmcl_functions.h"
MODULE pgmcl_library
      USE the_mpi
      implicit none
      integer MyRankGlobal, MyRankLocal
      TYPE T_grid
         integer IsFE
         integer eta_rho, xi_rho
         REAL*8, dimension(:,:), pointer :: LON_fd, LAT_fd
         integer, dimension(:,:), pointer :: MSK_fd
         integer nbNode, nbTrig
         REAL*8, dimension(:), pointer :: LON_fe, LAT_fe
         integer, dimension(:,:), pointer :: ListTrig
      END TYPE T_grid
      TYPE T_sparse
         integer numelt, nRows, nColumns
         integer, dimension(:), pointer :: rows, columns
         real, dimension(:), pointer :: weights
      END TYPE T_sparse
      TYPE T_sparse_splitIO
         integer nbProc, nbNeedTot
         integer nProcInput, nProcOutput
         integer TheType
         integer comm
         integer, dimension(:), pointer :: ListProc, ListNbNeed
         integer, dimension(:,:), pointer :: ListIdx
         integer, dimension(:), pointer :: ListGlobalRankOutput
         integer, dimension(:), pointer :: ListGlobalRankInput
         type(T_sparse) :: RestrictedSparseMat
      END TYPE T_sparse_splitIO
      TYPE T_sparse_splitIO_as
         integer nbProc, nbNeedTot
         integer TheType
         integer Ntimes
         integer comm
         integer, dimension(:), pointer :: ListProc
         integer, dimension(:), pointer :: ListGlobalRankOutput
         integer, dimension(:), pointer :: ListGlobalRankInput
         integer, dimension(:), pointer :: as_send_rqst
         integer, dimension(:), pointer :: as_recv_rqst
         integer, dimension(:,:), pointer :: as_send_stat
         integer, dimension(:,:), pointer :: as_recv_stat
         integer, dimension(:), pointer :: as_send_type
         integer, dimension(:), pointer :: as_recv_type
#ifdef DEBUG
         integer, dimension(:), pointer :: as_nbsend
         integer, dimension(:), pointer :: as_nbrecv
#endif
         type(T_sparse) :: RestrictedSparseMat
      END TYPE T_sparse_splitIO_as
      TYPE T_listMapGlobal
        integer eColor
        integer NbLocalProc
        logical IsAssigned
        integer, dimension(:), pointer :: MapLocalToGlobal
      END TYPE T_listMapGlobal
      TYPE T_localNodes
        integer comm
        integer NProc
        integer maxColor
        integer, dimension(:), pointer :: ListLocalRank
        integer, dimension(:), pointer :: ListColor
        integer, dimension(:), pointer :: ListFirstRank
        integer, dimension(:), pointer :: NbCompNodes
        type(T_listMapGlobal), dimension(:), pointer :: ListMapGlobal
        integer rank_in_comm
        integer MyLocalRank
        integer MyColor
      END TYPE T_localNodes
      TYPE T_node_partition
        integer Nnode
        integer Nthread
        integer, dimension(:,:), pointer :: TheMatrix
      END TYPE T_node_partition
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GetString (TheNb, eStr)
      implicit none
      character*4, intent(out) :: eStr
      integer, intent(in) :: TheNb
      integer STAT_VALUE
      character*1 eStr1
      character*2 eStr2
      character*3 eStr3
      IF (TheNb.le.9) THEN
         WRITE (FMT=10, UNIT=eStr1, IOSTAT=STAT_VALUE) TheNb
         eStr='000' // eStr1
      ELSE IF (TheNb.le.99) THEN
         WRITE (FMT=20, UNIT=eStr2, IOSTAT=STAT_VALUE) TheNb
         eStr='00' // eStr2
      ELSE IF (TheNb.le.999) THEN
         WRITE (FMT=30, UNIT=eStr3, IOSTAT=STAT_VALUE) TheNb
         eStr='0' // eStr3
      ELSE
         WRITE (FMT=40, UNIT=eStr, IOSTAT=STAT_VALUE) TheNb
      END IF
 10   FORMAT (i1)
 20   FORMAT (i2)
 30   FORMAT (i3)
 40   FORMAT (i4)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TEXT_ReadSparseMatrix(sMat, FileName)
      implicit none
      type(T_sparse), intent(inout) :: sMat
      character(len=40), intent(in) :: FileName
      !
      integer numelements
      integer n
      integer nRows, nColumns
      integer mdev, MyError, MyRank
      mdev = 745
      open(mdev, file=TRIM(FileName), status="old")
      read(mdev,*) nRows, nColumns, numelements
      sMat % nRows=nRows
      sMat % nColumns=nColumns
      sMat % numelt=numelements
      allocate(sMat % rows(numelements))
      allocate(sMat % columns(numelements))
      allocate(sMat % weights(numelements))
      do n=1, numelements
         read(mdev,*) sMat % rows(n), sMat % columns(n), sMat % weights(n)
      end do
      close(mdev)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_GRID_ARRAY_FD_r8(eta_rho, xi_rho, LON, LAT, MSK, eGrid)
      implicit none
      integer, intent(in) :: eta_rho, xi_rho
      real*8, intent(in) :: LON(eta_rho, xi_rho)
      real*8, intent(in) :: LAT(eta_rho, xi_rho)
      integer, intent(in) :: MSK(eta_rho, xi_rho)
      type(T_grid), intent(inout) :: eGrid
      eGrid % IsFE = 0
      eGrid % eta_rho=eta_rho
      eGrid % xi_rho=xi_rho
      allocate(eGrid % LON_fd(eta_rho, xi_rho))
      allocate(eGrid % LAT_fd(eta_rho, xi_rho))
      allocate(eGrid % MSK_fd(eta_rho, xi_rho))
      eGrid % LON_fd=LON
      eGrid % LAT_fd=LAT
      eGrid % MSK_fd=MSK
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RETRIEVE_GRID_ARRAY_FD_r8(eta_rho, xi_rho, LON, LAT, MSK, eGrid)
      integer, intent(inout) :: eta_rho, xi_rho
      real*8, allocatable, intent(inout) :: LON(:,:)
      real*8, allocatable, intent(inout) :: LAT(:,:)
      integer, allocatable, intent(inout) :: MSK(:,:)
      type(T_grid), intent(in) :: eGrid
      IF (eGrid % IsFE .eq. 1) THEN
        Print *, 'Error on the grid'
        STOP
      END IF
      eta_rho=eGrid % eta_rho
      xi_rho =eGrid % xi_rho
      allocate(LON(eta_rho, xi_rho))
      allocate(LAT(eta_rho, xi_rho))
      allocate(MSK(eta_rho, xi_rho))
      LON = eGrid % LON_fd
      LAT = eGrid % LAT_fd
      MSK = eGrid % MSK_fd
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PRINT_GRID_INFO(IndexOutput, eGrid)
      implicit none
      integer, intent(in) :: IndexOutput
      type(T_grid), intent(inout) :: eGrid
      WRITE(IndexOutput, *) 'IsFE=', eGrid % IsFE
      IF (eGrid %IsFE .eq. 1) THEN
        WRITE(IndexOutput, *) 'LON(min/max)=', minval(eGrid % LON_fe), maxval(eGrid % LON_fe)
        WRITE(IndexOutput, *) 'LAT(min/max)=', minval(eGrid % LAT_fe), maxval(eGrid % LAT_fe)
      ELSE
        WRITE(IndexOutput, *) 'LON(min/max)=', minval(eGrid % LON_fd), maxval(eGrid % LON_fd)
        WRITE(IndexOutput, *) 'LAT(min/max)=', minval(eGrid % LAT_fd), maxval(eGrid % LAT_fd)
      END IF
      FLUSH(IndexOutput)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_GRID_ARRAY_FD_r4(eta_rho, xi_rho, LON, LAT, MSK, eGrid)
      implicit none
      integer, intent(in) :: eta_rho, xi_rho
      real, intent(in) :: LON(eta_rho, xi_rho)
      real, intent(in) :: LAT(eta_rho, xi_rho)
      integer, intent(in) :: MSK(eta_rho, xi_rho)
      type(T_grid), intent(inout) :: eGrid
      real*8 :: LONr8(eta_rho, xi_rho), LATr8(eta_rho, xi_rho)
      LONr8=DBLE(LON)
      LATr8=DBLE(LAT)
      CALL GET_GRID_ARRAY_FD_r8(eta_rho, xi_rho, LONr8, LATr8, MSK, eGrid)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_GRID_ARRAY_FE_r8(nbNode, nbTrig, LON, LAT, ListTrig, eGrid)
      implicit none
      integer, intent(in) :: nbNode, nbTrig
      real*8, intent(in) :: LON(nbNode)
      real*8, intent(in) :: LAT(nbNode)
      integer, intent(in) :: ListTrig(3, nbTrig)
      type(T_grid), intent(inout) :: eGrid
      eGrid % IsFE = 1
      eGrid % nbNode=nbNode
      eGrid % nbTrig=nbTrig
      allocate(eGrid % LON_fe(nbNode))
      allocate(eGrid % LAT_fe(nbNode))
      allocate(eGrid % ListTrig(3, nbTrig))
      eGrid % LON_fe=LON
      eGrid % LAT_fe=LAT
      eGrid % ListTrig=ListTrig
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOCATE_GRID_ARRAY(eGrid)
      implicit none
      type(T_grid), intent(inout) :: eGrid
      IF (eGrid % IsFE .eq. 1) THEN
        deallocate(eGrid % LON_fe)
        deallocate(eGrid % LAT_fe)
        deallocate(eGrid % ListTrig)
      ELSE
        deallocate(eGrid % LON_fd)
        deallocate(eGrid % LAT_fd)
        deallocate(eGrid % MSK_fd)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_GRID_ARRAY_FE_r4(nbNode, nbTrig, LON, LAT, ListTrig, eGrid)
      implicit none
      integer, intent(in) :: nbNode, nbTrig
      real, intent(in) :: LON(nbNode)
      real, intent(in) :: LAT(nbNode)
      integer, intent(in) :: ListTrig(3, nbTrig)
      type(T_grid), intent(inout) :: eGrid
      real*8 :: LONr8(nbNode), LATr8(nbNode)
      LONr8=DBLE(LON)
      LATr8=DBLE(LAT)
      CALL GET_GRID_ARRAY_FE_r8(nbNode, nbTrig, LONr8, LATr8, ListTrig, eGrid)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GENERIC_NETCDF_ERROR(CallFct, idx, iret)
      USE NETCDF
      implicit none
      integer, intent(in) :: iret, idx
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        WRITE(*,*) TRIM(CallFct), ' -', idx, '-', TRIM(CHRERR)
        STOP 'DEBUG the netcdf'
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NETCDF_ReadSparseMatrix(sMat, nbKey, ListKey, FileName)
      USE NETCDF
      implicit none
      type(T_sparse), intent(inout) :: sMat
      integer, intent(inout) :: nbKey
      integer, intent(inout), allocatable :: ListKey(:)
      character(len=40), intent(in) :: FileName
      !
      integer nRows, nColumns, numelements
      integer iret, ncid, var_id
      integer nbkey_dims, nrows_dims, ncolumns_dims, numelements_dims
      integer ielt
      character (len = *), parameter :: CallFct="NETCDF_ReadSparseMatric"
      Print *, 'Reading FileName=', TRIM(FileName)
      iret=nf90_open(TRIM(FileName), NF90_NOWRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      iret=nf90_inq_dimid(ncid,'nbKey',nbkey_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      iret=nf90_inquire_dimension(ncid, nbkey_dims, len=nbKey)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      !
      iret=nf90_inq_dimid(ncid,'nRows', nrows_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      iret=nf90_inquire_dimension(ncid, nrows_dims, len=nRows)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      !
      iret=nf90_inq_dimid(ncid,'nColumns', ncolumns_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      iret=nf90_inquire_dimension(ncid, ncolumns_dims, len=nColumns)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      !
      iret=nf90_inq_dimid(ncid,'numelements', numelements_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      iret=nf90_inquire_dimension(ncid, numelements_dims, len=numelements)
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
      !
      allocate(ListKey(nbKey))
      iret=nf90_inq_varid(ncid, "ListKey", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
      iret=nf90_get_var(ncid, var_id, ListKey, start=(/1/), count=(/nbKey/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
      !
      sMat % nRows=nRows
      sMat % nColumns=nColumns
      sMat % numelt=numelements
      allocate(sMat % rows(numelements))
      allocate(sMat % columns(numelements))
      allocate(sMat % weights(numelements))
      !
      iret=nf90_inq_varid(ncid, "rows", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
      iret=nf90_get_var(ncid, var_id, sMat % rows, start=(/1/), count=(/numelements/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      !
      iret=nf90_inq_varid(ncid, "columns", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      iret=nf90_get_var(ncid, var_id, sMat % columns, start=(/1/), count=(/numelements/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
      !
      iret=nf90_inq_varid(ncid, "weights", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
      iret=nf90_get_var(ncid, var_id, sMat % weights, start=(/1/), count=(/numelements/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)
      !
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 18, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NETCDF_WriteSparseMatrix(sMat, nbKey, ListKey, FileName)
      USE NETCDF
      implicit none
      type(T_sparse), intent(in) :: sMat
      integer, intent(in) :: nbKey
      integer, intent(in) :: ListKey(:)
      character(len=40), intent(in) :: FileName
      integer nRows, nColumns, numelements
      integer iret, ncid, var_id
      integer nbkey_dims, nrows_dims, ncolumns_dims, numelements_dims
      character (len = *), parameter :: CallFct="NETCDF_WriteSparseMatric"
      nRows=sMat % nRows
      nColumns=sMat % nColumns
      numelements=sMat % numelt
      iret=nf90_create(FileName, NF90_CLOBBER, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      iret = nf90_def_dim(ncid, 'nbKey', nbKey, nbkey_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      iret = nf90_def_dim(ncid, 'nRows', nRows, nrows_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      iret = nf90_def_dim(ncid, 'nColumns', nColumns, ncolumns_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      iret = nf90_def_dim(ncid, 'numelements', numelements, numelements_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      !
      iret = nf90_def_var(ncid,'ListKey', NF90_INT,(/ nbkey_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      iret = nf90_def_var(ncid,'rows', NF90_INT,(/ numelements_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      iret = nf90_def_var(ncid,'columns', NF90_INT,(/ numelements_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      iret = nf90_def_var(ncid,'weights', NF90_REAL,(/ numelements_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
      !
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
     
      iret=nf90_open(FileName, NF90_WRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
      iret=nf90_inq_varid(ncid, "ListKey", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
      iret=nf90_put_var(ncid,var_id,ListKey,start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      iret=nf90_inq_varid(ncid, "rows", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      iret=nf90_put_var(ncid,var_id,sMat % rows,start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
      iret=nf90_inq_varid(ncid, "columns", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
      iret=nf90_put_var(ncid,var_id,sMat % columns,start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)
      iret=nf90_inq_varid(ncid, "weights", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 18, iret)
      iret=nf90_put_var(ncid,var_id,sMat % weights,start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DeallocSparseMatrix(sMat)
      implicit none
      type(T_sparse), intent(inout) :: sMat
      deallocate(sMat % rows)
      deallocate(sMat % columns)
      deallocate(sMat % weights)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TwoPiNormalization(TheAng)
      implicit none
      REAL*8, intent(inout) :: TheAng
      !
      REAL*8 :: ThePi
      ThePi=3.141592653589792
      IF (TheAng < -ThePi) THEN
        TheAng=TheAng + 2*ThePi
      ENDIF
      IF (TheAng > ThePi) THEN
        TheAng=TheAng - 2*ThePi
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MySign(TheVal, TheSign)
      implicit none
      REAL*8, intent(in) :: TheVal
      REAL*8, intent(out) :: TheSign
      IF (TheVal > 0) THEN
        TheSign=1
      ELSEIF (TheVal < 0) THEN
        TheSign=-1
      ELSE
        TheSign=0
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CreateAngleMatrix(eta_rho, xi_rho, ANG_rho, LON_rho, LAT_rho)
      implicit none
      integer, intent(in) :: eta_rho, xi_rho
      REAL*8, DIMENSION(eta_rho, xi_rho), intent(in) :: LON_rho, LAT_rho
      REAL*8, DIMENSION(eta_rho, xi_rho), intent(out) :: ANG_rho
      !
      integer eta_u, xi_u, iEta, iXi
      real*8, allocatable :: LONrad_u(:,:)
      real*8, allocatable :: LATrad_u(:,:)
      real*8, allocatable :: azim(:,:)
      real*8 :: eAzim, fAzim, dlam, eFact1, eFact2
      real*8 :: signAzim, signDlam, ThePi, DegTwoRad
      real*8 :: eLon, eLat, phi1, phi2, xlam1, xlam2
      real*8 :: TPSI2, cta12
      eta_u=eta_rho
      xi_u=xi_rho-1
      allocate(LONrad_u(eta_u,xi_u), LATrad_u(eta_u,xi_u), azim(eta_u,xi_u-1))
      ThePi=3.141592653589792
      DegTwoRad=ThePi/180
      DO iEta=1,eta_u
        DO iXi=1,xi_u
          eLon=(LON_rho(iEta,iXi)+LON_rho(iEta,iXi+1))*0.5
          eLat=(LAT_rho(iEta,iXi)+LAT_rho(iEta,iXi+1))*0.5
          LONrad_u(iEta,iXi)=eLon*DegTwoRad
          LATrad_u(iEta,iXi)=eLat*DegTwoRad
        END DO
      END DO
      DO iEta=1,eta_u
        DO iXi=1,xi_u-1
          phi1=LATrad_u(iEta,iXi)
          xlam1=LONrad_u(iEta,iXi)
          phi2=LATrad_u(iEta,iXi+1)
          xlam2=LONrad_u(iEta,iXi+1)
          TPSI2=TAN(phi2)
          dlam=xlam2-xlam1
          CALL TwoPiNormalization(dlam)
          cta12=(cos(phi1)*TPSI2 - sin(phi1)*cos(dlam))/sin(dlam)
          eAzim=ATAN(1./cta12)
          CALL MySign(eAzim, signAzim)
          CALL MySign(dlam, signDlam)
          IF (signDlam.ne.signAzim) THEN
            eFact2=1
          ELSE
            eFact2=0
          END IF
          eFact1=-signAzim
          fAzim=eAzim+ThePi*eFact1*eFact2
          azim(iEta,iXi)=fAzim
        END DO
      END DO
      DO iEta=1,eta_u
        DO iXi=2,xi_u
          ANG_rho(iEta,iXi)=ThePi*0.5 - azim(iEta,iXi-1)
        END DO
      END DO
      DO iEta=1,eta_u
        ANG_rho(iEta,1)=ANG_rho(iEta,2)
        ANG_rho(iEta,xi_rho)=ANG_rho(iEta,xi_u)
      END DO
      deallocate(LONrad_u, LATrad_u, azim)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CreateAngleMatrix_v(eta_rho, xi_rho, ANG_rho, LON_rho, LAT_rho)
      implicit none
      integer, intent(in) :: eta_rho, xi_rho
      REAL*8, DIMENSION(eta_rho, xi_rho), intent(in) :: LON_rho, LAT_rho
      REAL*8, DIMENSION(eta_rho, xi_rho), intent(out) :: ANG_rho
      !
      integer eta_v, xi_v, iEta, iXi
      real*8, allocatable :: LONrad_v(:,:)
      real*8, allocatable :: LATrad_v(:,:)
      real*8, allocatable :: azim(:,:)
      real*8 :: eAzim, fAzim, dlam, eFact1, eFact2
      real*8 :: signAzim, signDlam, ThePi, DegTwoRad
      real*8 :: eLon, eLat, phi1, phi2, xlam1, xlam2
      real*8 :: TPSI2, cta12
      eta_v=eta_rho-1
      xi_v=xi_rho
      allocate(LONrad_v(eta_v,xi_v))
      allocate(LATrad_v(eta_v,xi_v))
      allocate(azim(eta_v-1,xi_v))
      ThePi=3.141592653589792
      DegTwoRad=ThePi/180
      DO iEta=1,eta_v
        DO iXi=1,xi_v
          eLon=(LON_rho(iEta,iXi)+LON_rho(iEta+1,iXi))*0.5
          eLat=(LAT_rho(iEta,iXi)+LAT_rho(iEta+1,iXi))*0.5
          LONrad_v(iEta,iXi)=eLon*DegTwoRad
          LATrad_v(iEta,iXi)=eLat*DegTwoRad
        END DO
      END DO
      DO iEta=1,eta_v-1
        DO iXi=1,xi_v
          phi1=LATrad_v(iEta,iXi)
          xlam1=LONrad_v(iEta,iXi)
          phi2=LATrad_v(iEta+1,iXi)
          xlam2=LONrad_v(iEta+1,iXi)
          TPSI2=TAN(phi2)
          dlam=xlam2-xlam1
          CALL TwoPiNormalization(dlam)
          cta12=(cos(phi1)*TPSI2 - sin(phi1)*cos(dlam))/sin(dlam)
          eAzim=ATAN(1./cta12)
          CALL MySign(eAzim, signAzim)
          CALL MySign(dlam, signDlam)
          IF (signDlam.ne.signAzim) THEN
            eFact2=1
          ELSE
            eFact2=0
          END IF
          eFact1=-signAzim
          fAzim=eAzim+ThePi*eFact1*eFact2
          azim(iEta,iXi)=fAzim
        END DO
      END DO
      DO iEta=2,eta_v
        DO iXi=1,xi_v
          ANG_rho(iEta,iXi)=ThePi*0.5 - azim(iEta-1,iXi)
        END DO
      END DO
      DO iXi=1,xi_v
        ANG_rho(1,iXi)=ANG_rho(2,iXi)
        ANG_rho(eta_rho,iXi)=ANG_rho(eta_v,iXi)
      END DO
      deallocate(LONrad_v, LATrad_v, azim)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ALLOCATE_node_partition(Nnode, Nthread, MatrixBelonging)
      implicit none
      integer, intent(in) :: Nnode, Nthread
      type(T_node_partition), intent(inout) :: MatrixBelonging
      MatrixBelonging % Nnode=Nnode
      MatrixBelonging % Nthread=Nthread
      allocate(MatrixBelonging % TheMatrix(Nnode, Nthread))
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOCATE_node_partition(MatrixBelonging)
      implicit none
      type(T_node_partition), intent(inout) :: MatrixBelonging
      deallocate(MatrixBelonging % TheMatrix)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_GetSystemOutputSide(ArrLocal, eColorInput, eColorOutput, MatrixBelongingInput, MatrixBelongingOutput, MySparseMatrix, TheArr)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: eColorInput, eColorOutput
      TYPE(T_node_partition), intent(in) :: MatrixBelongingInput
      TYPE(T_node_partition), intent(in) :: MatrixBelongingOutput
      type(T_sparse), intent(in) :: MySparseMatrix
      type(T_sparse_splitIO), intent(inout) :: TheArr
      !
      integer, allocatable :: eSetNeedInp(:)
      integer nInput, nOutput, idxs, numelt
      integer iCol, iRow, ielt, idx, iRowRed, iColRed
      integer iProcIn, nbProc, iProc
      integer, dimension(:), pointer :: est_rows, est_columns
      real, dimension(:), pointer :: est_weights
      integer, dimension(:), pointer :: ListIndexCol, ListIndexRow
      integer, dimension(:), pointer :: nbNeeded
      integer, dimension(:), pointer :: est_idxlocout(:)
      integer MaxNbNeed, idxlocout, numeltRed, idx_loc, idxlocin
      real eWeight
      integer nProcInput, nProcOutput, iProcOut
      integer nOutputRed, nInputRed
!      Print *, 'MPI_INTERP_GetOutputSide, step 1, MRG=', MyRankGlobal
      nProcInput =ArrLocal % NbCompNodes(eColorInput)
      nProcOutput=ArrLocal % NbCompNodes(eColorOutput)
      TheArr % nProcInput=nProcInput
      TheArr % nProcOutput=nProcOutput
      iProcOut=ArrLocal % MyLocalRank
      allocate(TheArr % ListGlobalRankOutput(nProcOutput))
      allocate(TheArr % ListGlobalRankInput(nProcInput))
      DO iProc=1,nProcInput
        TheArr % ListGlobalRankInput(iProc)=GetGlobalRank(ArrLocal, eColorInput, iProc)
      END DO
      DO iProc=1,nProcOutput
        TheArr % ListGlobalRankOutput(iProc)=GetGlobalRank(ArrLocal, eColorOutput, iProc)
      END DO
!      Print *, 'MPI_INTERP_GetOutputSide, step 2, MRG=', MyRankGlobal
      TheArr % comm=ArrLocal % comm
      nInput=MySparseMatrix%nColumns
      nOutput=MySparseMatrix%nRows
      ALLOCATE(eSetNeedInp(nInput))
      eSetNeedInp=0
      idxs=0
      numelt=MySparseMatrix%numelt
      allocate(est_rows(numelt))
      allocate(est_columns(numelt))
      allocate(est_weights(numelt))
      allocate(est_idxlocout(numelt))
      DO ielt=1,numelt
        iCol=MySparseMatrix%columns(ielt)
        iRow=MySparseMatrix%rows(ielt)
        eWeight=MySparseMatrix%weights(ielt)
        idxlocout=MatrixBelongingOutput % TheMatrix(iRow,iProcOut+1)
        IF (idxlocout .gt. 0) THEN
          eSetNeedInp(iCol)=1
          idxs=idxs+1
          est_rows(idxs)=iRow
          est_columns(idxs)=iCol
          est_weights(idxs)=eWeight
          est_idxlocout(idxs)=idxlocout
        END IF
      END DO
      numeltRed=idxs
      allocate(ListIndexCol(nInput))
      idx=0
      DO iCol=1,nInput
        IF (eSetNeedInp(iCol).eq.1) THEN
          idx=idx+1
          ListIndexCol(iCol)=idx
        END IF
      END DO
      nInputRed=idx  ! size of data needed for output
      TheArr % nbNeedTot=nInputRed
      allocate(ListIndexRow(nOutput))
      idx=0
      DO iRow=1,nOutput
        idx_loc=MatrixBelongingOutput % TheMatrix(iRow,iProcOut+1)
        IF (idx_loc .gt. 0) THEN
          idx=idx+1
          ListIndexRow(iRow)=idx
        END IF
      END DO
      nOutputRed=idx  ! this is the size of the data output.
      TheArr % RestrictedSparseMat % numelt=numeltRed
      TheArr % RestrictedSparseMat % nRows=nOutputRed
      TheArr % RestrictedSparseMat % nColumns=nInputRed
      ALLOCATE(TheArr % RestrictedSparseMat % rows(numeltRed))
      ALLOCATE(TheArr % RestrictedSparseMat % columns(numeltRed))
      ALLOCATE(TheArr % RestrictedSparseMat % weights(numeltRed))
      TheArr % TheType=1
      DO idx=1,numeltRed
!        iRow=est_rows(idx)
!        iRowRed=ListIndexRow(iRow)
        idxlocout=est_idxlocout(idx)
        TheArr % RestrictedSparseMat % rows(idx)=idxlocout
        !
        iCol=est_columns(idx)
        iColRed=ListIndexCol(iCol)
        TheArr % RestrictedSparseMat % columns(idx)=iColRed
        !
        TheArr % RestrictedSparseMat % weights(idx)=est_weights(idx)
      END DO
      deallocate(est_rows, est_columns, est_weights, est_idxlocout)
      ALLOCATE(nbNeeded(nProcInput))
      nbNeeded=0
      DO iCol=1,nInput
        IF (eSetNeedInp(iCol).eq.1) THEN
          DO iProcIn=1,nProcInput
            idxlocin=MatrixBelongingInput % TheMatrix(iCol,iProcIn)
            IF (idxlocin .gt. 0) THEN
               nbNeeded(iProcIn)=nbNeeded(iProcIn)+1
            END IF
          END DO
        END IF
      END DO
      nbProc=0
      MaxNbNeed=0
      DO iProcIn=1,nProcInput
        IF (nbNeeded(iProcIn).gt.0) THEN
          nbProc=nbProc+1
        END IF
        IF (nbNeeded(iProcIn).gt.MaxNbNeed) THEN
          MaxNbNeed=nbNeeded(iProcIn)
        END IF
      END DO
      TheArr % nbProc=nbProc
      ALLOCATE(TheArr % ListProc(nbProc))
      ALLOCATE(TheArr % ListNbNeed(nbProc))
      iProc=0
      DO iProcIn=1,nProcInput
        IF (nbNeeded(iProcIn).gt.0) THEN
          iProc=iProc+1
          TheArr % ListProc(iProc)=iProcIn
          TheArr % ListNbNeed(iProc)=nbNeeded(iProcIn)
        END IF
      END DO
      ALLOCATE(TheArr % ListIdx(MaxNbNeed,nbProc))
      TheArr % ListIdx=0
      DO iProc=1,nbProc
        iProcIn=TheArr % ListProc(iProc)
        idx=0
        DO iCol=1,nInput
          IF (eSetNeedInp(iCol).eq.1) THEN
            idxlocin=MatrixBelongingInput % TheMatrix(iCol,iProcIn)
            IF (idxlocin .gt. 0) THEN
              idx=idx+1
              iColRed=ListIndexCol(iCol)
              TheArr % ListIdx(idx,iProc)=iColRed
            END IF
          END IF
        END DO
      END DO
      DEALLOCATE(nbNeeded, eSetNeedInp, ListIndexRow, ListIndexCol)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_GetSystemInputSide(ArrLocal, eColorInput, eColorOutput, MatrixBelongingInput, MatrixBelongingOutput, MySparseMatrix, TheArr)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: eColorInput, eColorOutput
      type(T_node_partition), intent(in) :: MatrixBelongingInput
      type(T_node_partition), intent(in) :: MatrixBelongingOutput
      type(T_sparse), intent(in) :: MySparseMatrix
      type(T_sparse_splitIO), intent(inout) :: TheArr
      !
      integer numelt, nInput, nOutput, ielt, iCol, iRow
      integer iProcOut, nb, i, iProc, idxlocin
      real eWeight
      integer SumInput, SumOutput
      integer, allocatable :: eSetNeedInp(:)
      integer, allocatable :: nbNeededInp(:)
      integer, allocatable :: ListIndexCol(:)
      integer, allocatable :: ListIndexLocIn(:)
      integer idx, iColRed, nbProc, MaxNbNeed
      integer nProcInput, nProcOutput, iProcIn
      integer nInputRed, idxlocout
      iProcIn=ArrLocal % MyLocalRank
      nProcInput =ArrLocal % NbCompNodes(eColorInput)
      nProcOutput=ArrLocal % NbCompNodes(eColorOutput)
      TheArr % nProcInput=nProcInput
      TheArr % nProcOutput=nProcOutput
      allocate(TheArr % ListGlobalRankOutput(nProcOutput))
      allocate(TheArr % ListGlobalRankInput(nProcInput))
      DO iProc=1,nProcInput
        TheArr % ListGlobalRankInput(iProc)=GetGlobalRank(ArrLocal, eColorInput, iProc)
      END DO
      DO iProc=1,nProcOutput
        TheArr % ListGlobalRankOutput(iProc)=GetGlobalRank(ArrLocal, eColorOutput, iProc)
      END DO
      TheArr % comm=ArrLocal % comm
      numelt=MySparseMatrix % numelt
      nInput=MySparseMatrix % nColumns
      nOutput=MySparseMatrix % nRows
      ALLOCATE(nbNeededInp(nProcOutput))
      allocate(eSetNeedInp(nInput))
      DO iProcOut=1,nProcOutput
        eSetNeedInp=0
        DO ielt=1,numelt
          iCol=MySparseMatrix % columns(ielt)
          iRow=MySparseMatrix % rows(ielt)
          IF (MatrixBelongingInput % TheMatrix(iCol,iProcIn+1) .gt. 0) THEN
            IF (MatrixBelongingOutput % TheMatrix(iRow,iProcOut) .gt. 0) THEN
              eSetNeedInp(iCol)=1
            END IF
          END IF
        END DO
        nb=0
        DO iCol=1,nInput
          IF (eSetNeedInp(iCol).eq.1) THEN
            nb=nb+1
          END IF
        END DO
        nbNeededInp(iProcOut)=nb
      END DO
      allocate(ListIndexCol(nInput), ListIndexLocIn(nInput))
      idx=0
      DO iCol=1,nInput
        idxlocin=MatrixBelongingInput % TheMatrix(iCol,iProcIn+1)
        IF (idxlocin .gt. 0) THEN
          idx=idx+1
          ListIndexCol(iCol)=idx
          ListIndexLocIn(iCol)=idxlocin
        END IF
      END DO
      nInputRed=idx
      TheArr % RestrictedSparseMat % nColumns=nInputRed
      nbProc=0
      MaxNbNeed=0
      DO iProcOut=1,nProcOutput
        IF (nbNeededInp(iProcOut).gt.0) THEN
          nbProc=nbProc+1
        END IF
        IF (nbNeededInp(iProcOut).gt.MaxNbNeed) THEN
          MaxNbNeed=nbNeededInp(iProcOut)
        END IF
      END DO
      TheArr % nbProc=nbProc
      TheArr % TheType=0
      ALLOCATE(TheArr % ListProc(nbProc))
      ALLOCATE(TheArr % ListNbNeed(nbProc))
      ALLOCATE(TheArr % ListIdx(MaxNbNeed,nbProc))
      iProc=0
      DO iProcOut=1,nProcOutput
        IF (nbNeededInp(iProcOut).gt.0) THEN
          iProc=iProc+1
          TheArr % ListProc(iProc)=iProcOut
          TheArr % ListNbNeed(iProc)=nbNeededInp(iProcOut)
        END IF
      END DO
      TheArr % ListIdx=0
      DO iProc=1,nbProc
        iProcOut=TheArr % ListProc(iProc)
        eSetNeedInp=0
        DO ielt=1,numelt
          iCol=MySparseMatrix % columns(ielt)
          iRow=MySparseMatrix % rows(ielt)
          idxlocin=MatrixBelongingInput % TheMatrix(iCol,iProcIn+1)
          idxlocout=MatrixBelongingOutput % TheMatrix(iRow,iProcOut)
          IF (idxlocin .gt. 0) THEN
            IF (idxlocout .gt. 0) THEN
              eSetNeedInp(iCol)=1
            END IF
          END IF
        END DO
        idx=0
        DO iCol=1,nInput
          IF (eSetNeedInp(iCol).eq.1) THEN
            idx=idx+1
!            iColRed=ListIndexCol(iCol)
            iColRed=ListIndexLocIn(iCol)
            TheArr % ListIdx(idx,iProc)=iColRed
          END IF
        END DO
      END DO
      deallocate(eSetNeedInp, nbNeededInp, ListIndexCol, ListIndexLocIn)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_GetAsyncOutput_r8(TheArr, Ntimes, ArrAsync)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: Ntimes
      type(T_sparse_splitIO_as), intent(out) :: ArrAsync
      integer nbProc, numelt, nRows, nColumns, iProc
      integer nProcInput, nProcOutput
      integer nbNeed, i, ierr
      integer, allocatable :: dspl_recv(:)
      nbProc=TheArr%nbProc

      ArrAsync % nbProc = nbProc
      ALLOCATE(ArrAsync % ListProc(nbProc))
      ArrAsync % ListProc=TheArr % ListProc

      nProcInput =TheArr % nProcInput
      nProcOutput=TheArr % nProcOutput
      allocate(ArrAsync % ListGlobalRankOutput(nProcOutput))
      allocate(ArrAsync % ListGlobalRankInput(nProcInput))
      ArrAsync % ListGlobalRankOutput=TheArr % ListGlobalRankOutput
      ArrAsync % ListGlobalRankInput=TheArr % ListGlobalRankInput

      ArrAsync%comm=TheArr%comm
      ArrAsync % nbNeedTot=TheArr % nbNeedTot
      ArrAsync % TheType=1
      ArrAsync % Ntimes = Ntimes

      numelt   = TheArr % RestrictedSparseMat % numelt
      nRows    = TheArr % RestrictedSparseMat % nRows
      nColumns = TheArr % RestrictedSparseMat % nColumns
      ALLOCATE(ArrAsync % RestrictedSparseMat % rows(numelt))
      ALLOCATE(ArrAsync % RestrictedSparseMat % columns(numelt))
      ALLOCATE(ArrAsync % RestrictedSparseMat % weights(numelt))
      ArrAsync % RestrictedSparseMat % numelt=numelt
      ArrAsync % RestrictedSparseMat % nRows=nRows
      ArrAsync % RestrictedSparseMat % nColumns=nColumns
      ArrAsync % RestrictedSparseMat % rows=TheArr % RestrictedSparseMat % rows
      ArrAsync % RestrictedSparseMat % columns=TheArr % RestrictedSparseMat % columns
      ArrAsync % RestrictedSparseMat % weights=TheArr % RestrictedSparseMat % weights

      allocate(ArrAsync % as_recv_rqst(nbProc))
      allocate(ArrAsync % as_recv_stat(MPI_STATUS_SIZE,nbProc))
      allocate(ArrAsync % as_recv_type(nbProc))
#ifdef DEBUG
      allocate(ArrAsync % as_nbrecv(nbProc))
#endif
      DO iProc=1,nbProc
        nbNeed=TheArr % ListNbNeed(iProc)
        allocate(dspl_recv(nbNeed))
        DO i=1,nbNeed
          dspl_recv(i)=Ntimes*(TheArr%ListIdx(i,iProc)-1)
        END DO
#ifdef DEBUG
        ArrAsync % as_nbrecv(iProc)=nbNeed
#endif
        call mpi_type_create_indexed_block(nbNeed,Ntimes,dspl_recv,MPI_DOUBLE_PRECISION, ArrAsync % as_recv_type(iProc), ierr)
        call mpi_type_commit(ArrAsync % as_recv_type(iProc), ierr)
        deallocate(dspl_recv)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_GetAsyncInput_r8(TheArr, Ntimes, ArrAsync)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: Ntimes
      type(T_sparse_splitIO_as), intent(out) :: ArrAsync
      integer nbProc, numelt, nRows, nColumns, iProc
      integer nProcInput, nProcOutput
      integer nbNeed, i, ierr
      integer, allocatable :: dspl_send(:)
      nbProc=TheArr%nbProc

      ArrAsync % nbProc = nbProc
      ALLOCATE(ArrAsync % ListProc(nbProc))
      ArrAsync % ListProc=TheArr % ListProc

      nProcInput =TheArr % nProcInput
      nProcOutput=TheArr % nProcOutput
      allocate(ArrAsync % ListGlobalRankOutput(nProcOutput))
      allocate(ArrAsync % ListGlobalRankInput(nProcInput))
      ArrAsync % ListGlobalRankOutput=TheArr % ListGlobalRankOutput
      ArrAsync % ListGlobalRankInput=TheArr % ListGlobalRankInput

      ArrAsync%comm=TheArr%comm
      ArrAsync % nbNeedTot=TheArr % nbNeedTot
      ArrAsync % TheType=0
      ArrAsync % Ntimes = Ntimes


      allocate(ArrAsync % as_send_rqst(nbProc))
      allocate(ArrAsync % as_send_stat(MPI_STATUS_SIZE,nbProc))
      allocate(ArrAsync % as_send_type(nbProc))
#ifdef DEBUG
      allocate(ArrAsync % as_nbsend(nbProc))
#endif
      DO iProc=1,nbProc
        nbNeed=TheArr % ListNbNeed(iProc)
        allocate(dspl_send(nbNeed))
        DO i=1,nbNeed
          dspl_send(i)=Ntimes*(TheArr%ListIdx(i,iProc)-1)
        END DO
#ifdef DEBUG
        ArrAsync % as_nbsend(iProc)=nbNeed
#endif
        call mpi_type_create_indexed_block(nbNeed,Ntimes,dspl_send,MPI_DOUBLE_PRECISION, ArrAsync % as_send_type(iProc), ierr)
        call mpi_type_commit(ArrAsync % as_send_type(iProc), ierr)
        deallocate(dspl_send)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_GetAsyncOutput_r4(TheArr, Ntimes, ArrAsync)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: Ntimes
      type(T_sparse_splitIO_as), intent(out) :: ArrAsync
      integer nbProc, numelt, nRows, nColumns, iProc
      integer nProcInput, nProcOutput
      integer nbNeed, i, ierr
      integer, allocatable :: dspl_recv(:)
      nbProc=TheArr%nbProc

      ArrAsync % nbProc = nbProc
      ALLOCATE(ArrAsync % ListProc(nbProc))
      ArrAsync % ListProc=TheArr % ListProc

      nProcInput =TheArr % nProcInput
      nProcOutput=TheArr % nProcOutput
      allocate(ArrAsync % ListGlobalRankOutput(nProcOutput))
      allocate(ArrAsync % ListGlobalRankInput(nProcInput))
      ArrAsync % ListGlobalRankOutput=TheArr % ListGlobalRankOutput
      ArrAsync % ListGlobalRankInput=TheArr % ListGlobalRankInput

      ArrAsync%comm=TheArr%comm
      ArrAsync % nbNeedTot=TheArr % nbNeedTot
      ArrAsync % TheType=1
      ArrAsync % Ntimes = Ntimes

      numelt   = TheArr % RestrictedSparseMat % numelt
      nRows    = TheArr % RestrictedSparseMat % nRows
      nColumns = TheArr % RestrictedSparseMat % nColumns
      ALLOCATE(ArrAsync % RestrictedSparseMat % rows(numelt))
      ALLOCATE(ArrAsync % RestrictedSparseMat % columns(numelt))
      ALLOCATE(ArrAsync % RestrictedSparseMat % weights(numelt))
      ArrAsync % RestrictedSparseMat % numelt=numelt
      ArrAsync % RestrictedSparseMat % nRows=nRows
      ArrAsync % RestrictedSparseMat % nColumns=nColumns
      ArrAsync % RestrictedSparseMat % rows=TheArr % RestrictedSparseMat % rows
      ArrAsync % RestrictedSparseMat % columns=TheArr % RestrictedSparseMat % columns
      ArrAsync % RestrictedSparseMat % weights=TheArr % RestrictedSparseMat % weights

      allocate(ArrAsync % as_recv_rqst(nbProc))
      allocate(ArrAsync % as_recv_stat(MPI_STATUS_SIZE,nbProc))
      allocate(ArrAsync % as_recv_type(nbProc))
#ifdef DEBUG
      allocate(ArrAsync % as_nbrecv(nbProc))
#endif
      DO iProc=1,nbProc
        nbNeed=TheArr % ListNbNeed(iProc)
        allocate(dspl_recv(nbNeed))
        DO i=1,nbNeed
          dspl_recv(i)=Ntimes*(TheArr%ListIdx(i,iProc)-1)
        END DO
#ifdef DEBUG
        ArrAsync % as_nbrecv(iProc)=nbNeed
#endif
        call mpi_type_create_indexed_block(nbNeed,Ntimes,dspl_recv,MPI_REAL, ArrAsync % as_recv_type(iProc), ierr)
        call mpi_type_commit(ArrAsync % as_recv_type(iProc), ierr)
        deallocate(dspl_recv)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_GetAsyncInput_r4(TheArr, Ntimes, ArrAsync)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: Ntimes
      type(T_sparse_splitIO_as), intent(out) :: ArrAsync
      integer nbProc, numelt, nRows, nColumns, iProc
      integer nProcInput, nProcOutput
      integer nbNeed, i, ierr
      integer, allocatable :: dspl_send(:)
      nbProc=TheArr%nbProc

      ArrAsync % nbProc = nbProc
      ALLOCATE(ArrAsync % ListProc(nbProc))
      ArrAsync % ListProc=TheArr % ListProc

      nProcInput =TheArr % nProcInput
      nProcOutput=TheArr % nProcOutput
      allocate(ArrAsync % ListGlobalRankOutput(nProcOutput))
      allocate(ArrAsync % ListGlobalRankInput(nProcInput))
      ArrAsync % ListGlobalRankOutput=TheArr % ListGlobalRankOutput
      ArrAsync % ListGlobalRankInput=TheArr % ListGlobalRankInput

      ArrAsync%comm=TheArr%comm
      ArrAsync % nbNeedTot=TheArr % nbNeedTot
      ArrAsync % TheType=0
      ArrAsync % Ntimes = Ntimes


      allocate(ArrAsync % as_send_rqst(nbProc))
      allocate(ArrAsync % as_send_stat(MPI_STATUS_SIZE,nbProc))
      allocate(ArrAsync % as_send_type(nbProc))
#ifdef DEBUG
      allocate(ArrAsync % as_nbsend(nbProc))
#endif
      DO iProc=1,nbProc
        nbNeed=TheArr % ListNbNeed(iProc)
        allocate(dspl_send(nbNeed))
        DO i=1,nbNeed
          dspl_send(i)=Ntimes*(TheArr%ListIdx(i,iProc)-1)
        END DO
#ifdef DEBUG
        ArrAsync % as_nbsend(iProc)=nbNeed
#endif
        call mpi_type_create_indexed_block(nbNeed,Ntimes,dspl_send,MPI_REAL, ArrAsync % as_send_type(iProc), ierr)
        call mpi_type_commit(ArrAsync % as_send_type(iProc), ierr)
        deallocate(dspl_send)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_Dummy_Recv(ArrLocal, eTag, TheId)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: TheId, eTag
      integer status(MPI_STATUS_SIZE)
      integer rbuf_int(1)
      integer ierr, eOrig
      eOrig=ArrLocal % ListFirstRank(TheId)
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'MPI_INTERP_Dummy_Recv Bef eTag=', eTag, 'eOrig=', eOrig
      FLUSH(MyRankGlobal+740)
#endif
      CALL MPI_RECV(rbuf_int,1,MPI_INTEGER, eOrig, eTag, ArrLocal%comm, status, ierr)
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'MPI_INTERP_Dummy_Recv Aft eTag=', eTag, 'eOrig=', eOrig
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_Dummy_Send(ArrLocal, eTag, TheId)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: TheId, eTag
      integer ierr, dest, iProc, NbCompNode, MyLocalRank
      integer rbuf_int(1)
      MyLocalRank=ArrLocal % MyLocalRank
      rbuf_int(1)=3
      IF (MyLocalRank == 0) THEN
        NbCompNode=ArrLocal % NbCompNodes(TheId)
        DO iProc=1,NbCompNode
          dest=GetGlobalRank(ArrLocal, TheId, iProc)
#ifdef DEBUG
          WRITE(MyRankGlobal+740,*) 'MPI_INTERP_Dummy_Send Bef eTag=', eTag, 'dest=', dest
          FLUSH(MyRankGlobal+740)
#endif
          CALL MPI_SEND(rbuf_int,1,MPI_INTEGER, dest, eTag, ArrLocal%comm, ierr)
#ifdef DEBUG
          WRITE(MyRankGlobal+740,*) 'MPI_INTERP_Dummy_Send Aft eTag=', eTag, 'dest=', dest
          FLUSH(MyRankGlobal+740)
#endif
        END DO
      ENDIF
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_DummyB_Recv(ArrLocal, eTag, TheId)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: TheId, eTag
      integer status(MPI_STATUS_SIZE)
      integer rbuf_int(1)
      integer ierr, eOrig, iProc, NbCompNode
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'MPI_INTERP_Dummy_Recv Bef eTag=', eTag, 'eOrig=', eOrig
      FLUSH(MyRankGlobal+740)
#endif
      NbCompNode=ArrLocal % NbCompNodes(TheId)
      DO iProc=1,NbCompNode
        eOrig=GetGlobalRank(ArrLocal, TheId, iProc)
        CALL MPI_RECV(rbuf_int,1,MPI_INTEGER, eOrig, eTag, ArrLocal%comm, status, ierr)
      END DO
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'MPI_INTERP_Dummy_Recv Aft eTag=', eTag, 'eOrig=', eOrig
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_DummyB_Send(ArrLocal, eTag, TheId)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: TheId, eTag
      integer ierr, dest, iProc, NbCompNode, MyLocalRank
      integer rbuf_int(1)
      MyLocalRank=ArrLocal % MyLocalRank
      rbuf_int(1)=3
      NbCompNode=ArrLocal % NbCompNodes(TheId)
      DO iProc=1,NbCompNode
        dest=GetGlobalRank(ArrLocal, TheId, iProc)
#ifdef DEBUG
        WRITE(MyRankGlobal+740,*) 'MPI_INTERP_DummyB_Send Bef eTag=', eTag, 'dest=', dest
        FLUSH(MyRankGlobal+740)
#endif
        CALL MPI_SEND(rbuf_int,1,MPI_INTEGER, dest, eTag, ArrLocal%comm, ierr)
#ifdef DEBUG
        WRITE(MyRankGlobal+740,*) 'MPI_INTERP_DummyB_Send Aft eTag=', eTag, 'dest=', dest
        FLUSH(MyRankGlobal+740)
#endif
      ENDDO
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_RECV_r4(TheArr, eTag, TheFieldOut)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: eTag
      real, intent(out) :: TheFieldOut(:)
      !
      real :: eField(TheArr % nbNeedTot)
      integer :: status(MPI_STATUS_SIZE)
      integer ierror, nInput, nOutput, iProc, iProcIn, eSize
      real, allocatable :: rbuf_real_r4(:)
      integer eRankGlobal, iColRed, iRowRed, ielt, numelt, i
      real eWeight
#ifdef DEBUG
      integer :: rbuf_int(3)
      WRITE(MyRankGlobal+740,*) 'Beginning of MPI_INTERP_recv_r4, eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      nInput=TheArr%RestrictedSparseMat%nColumns
      nOutput=TheArr%RestrictedSparseMat%nRows
      DO iProc=1,TheArr%nbProc
        iProcIn=TheArr%ListProc(iProc)
        eSize=TheArr % ListNbNeed(iProc)
        ALLOCATE(rbuf_real_r4(eSize))
        eRankGlobal=TheArr%ListGlobalRankInput(iProcIn)
#ifdef DEBUG
        CALL MPI_RECV(rbuf_int,3,MPI_INTEGER,eRankGlobal, 144+eTag, TheArr%comm, status, ierror)
        IF ((rbuf_int(1).ne.4).or.(rbuf_int(2).ne.eSize).or. (rbuf_int(3).ne.eTag)) THEN
          Print *, 'RECV_r4, eTag=', eTag
          Print *, 'Arithmetic rbuf_int(1)=', rbuf_int(1), ' vs 4'
          Print *, '      size rbuf_int(2)=', rbuf_int(2), ' vs ', eSize
          Print *, '       tag rbuf_int(3)=', rbuf_int(3), ' vs ', eTag
          STOP
        END IF
#endif
        CALL MPI_RECV(rbuf_real_r4,eSize,MPI_REAL,eRankGlobal,eTag, TheArr%comm, status, ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
          Print *, 'Error in MPI_INTERP_RECV_r4'
          STOP
        END IF
        DO i=1,eSize
          iColRed=TheArr%ListIdx(i,iProc)
          eField(iColRed)=rbuf_real_r4(i)
        END DO
        DEALLOCATE(rbuf_real_r4)
      END DO
      TheFieldOut=0
      numelt=TheArr%RestrictedSparseMat%numelt
      DO ielt=1,numelt
        iColRed=TheArr%RestrictedSparseMat%columns(ielt)
        iRowRed=TheArr%RestrictedSparseMat%rows(ielt)
        eWeight=TheArr%RestrictedSparseMat%weights(ielt)
        TheFieldOut(iRowRed)=TheFieldOut(iRowRed)+ eWeight*eField(iColRed)
      END DO
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'Ending of MPI_INTERP_recv_r4'
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_RECV_r8(TheArr, eTag, TheFieldOut)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: eTag
      real*8, intent(out) :: TheFieldOut(:)
      !
      real*8 :: eField(TheArr % nbNeedTot)
      integer :: status(MPI_STATUS_SIZE)
      integer ierror, nInput, nOutput, iProc, iProcIn, eSize
      real*8, allocatable :: rbuf_real_r8(:)
      integer eRankGlobal, iColRed, iRowRed, ielt, numelt, i
      real eWeight
#ifdef DEBUG
      integer :: rbuf_int(3)
      WRITE(MyRankGlobal+740,*) 'MPI_INTERP_RECV_r8 eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      nInput=TheArr % RestrictedSparseMat % nColumns
      nOutput=TheArr % RestrictedSparseMat % nRows
      DO iProc=1,TheArr % nbProc
        iProcIn=TheArr % ListProc(iProc)
        eSize=TheArr % ListNbNeed(iProc)
        eRankGlobal=TheArr%ListGlobalRankInput(iProcIn)
#ifdef DEBUG
        WRITE(MyRankGlobal+740,*) 'iProc=',iProc,'iProcIn=',iProcIn
        WRITE(MyRankGlobal+740,*) 'eSize=', eSize, 'eRankGlobal=', eRankGlobal
        FLUSH(MyRankGlobal+740)
#endif
        ALLOCATE(rbuf_real_r8(eSize))
#ifdef DEBUG
        CALL MPI_RECV(rbuf_int,3,MPI_INTEGER,eRankGlobal, 144+eTag, TheArr%comm, status, ierror)
        IF ((rbuf_int(1).ne.8).or.(rbuf_int(2).ne.eSize).or.(rbuf_int(3).ne.eTag)) THEN
          Print *, 'RECV_r8, eTag=', eTag
          Print *, 'Arithmetic rbuf_int(1)=', rbuf_int(1), ' vs 8'
          Print *, '      size rbuf_int(2)=', rbuf_int(2), ' vs ', eSize
          Print *, '       tag rbuf_int(3)=', rbuf_int(3), ' vs ', eTag
          STOP
        END IF
        WRITE(MyRankGlobal+740,*) 'After recv 144 in MPI_INTERP_recv_r8'
        FLUSH(MyRankGlobal+740)
#endif
        CALL MPI_RECV(rbuf_real_r8,eSize,MPI_DOUBLE_PRECISION, eRankGlobal, eTag, TheArr%comm, status, ierror)
#ifdef DEBUG
        WRITE(MyRankGlobal+740,*) 'After recv in MPI_INTERP_recv_r8 eTag=', eTag
        FLUSH(MyRankGlobal+740)
#endif
        IF (ierror.ne.MPI_SUCCESS) THEN
          Print *, 'Error in MPI_INTERP_RECV_r8'
          STOP
        END IF
        DO i=1,eSize
          iColRed=TheArr%ListIdx(i,iProc)
          eField(iColRed)=rbuf_real_r8(i)
        END DO
        DEALLOCATE(rbuf_real_r8)
      END DO
      TheFieldOut=0
      numelt=TheArr%RestrictedSparseMat%numelt
      DO ielt=1,numelt
        iColRed=TheArr % RestrictedSparseMat % columns(ielt)
        iRowRed=TheArr % RestrictedSparseMat % rows(ielt)
        eWeight=TheArr % RestrictedSparseMat % weights(ielt)
        TheFieldOut(iRowRed)=TheFieldOut(iRowRed)+DBLE(eWeight)*eField(iColRed)
      END DO
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'End MPI_INTERP_recv_r8 eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_SEND_r4(TheArr, eTag, TheFieldIn)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: eTag
      real, intent(in) :: TheFieldIn(:)
      integer ierr, nbProc, iProc, eSize
      integer iProcOut, idx, iColRed, eRankGlobal
      real, allocatable :: rbuf_real_r4(:)
#ifdef DEBUG
      integer :: rbuf_int(3)
      WRITE(MyRankGlobal+740,*) 'MPI_INTERP_SEND_r4 eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      nbProc=TheArr%nbProc
      DO iProc=1,nbProc
        iProcOut=TheArr % ListProc(iProc)
        eSize=TheArr % ListNbNeed(iProc)
        allocate(rbuf_real_r4(eSize))
        DO idx=1,eSize
          iColRed=TheArr%ListIdx(idx,iProc)
          rbuf_real_r4(idx)=TheFieldIn(iColRed)
        END DO
        eRankGlobal=TheArr%ListGlobalRankOutput(iProcOut)
#ifdef DEBUG
        rbuf_int(1)=4
        rbuf_int(2)=eSize
        rbuf_int(3)=eTag
        CALL MPI_SEND(rbuf_int,3,MPI_INTEGER, eRankGlobal, 144+eTag, TheArr%comm, ierr)
#endif
        CALL MPI_SEND(rbuf_real_r4,eSize,MPI_REAL, eRankGlobal, eTag, TheArr%comm, ierr)
        IF (ierr.ne.MPI_SUCCESS) THEN
           Print *, 'Error in MPI_INTERP_SEND_r4'
           STOP
        END IF
        DEALLOCATE(rbuf_real_r4)
      END DO
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'Ending of MPI_INTERP_SEND_r4'
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_SEND_r8(TheArr, eTag, TheFieldIn)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: eTag
      real*8, intent(in) :: TheFieldIn(:)
      integer ierr, iProc, eSize
      integer iProcOut, idx, iColRed, eRankGlobal
      real*8, allocatable :: rbuf_real_r8(:)
#ifdef DEBUG
      integer :: rbuf_int(3)
      WRITE(MyRankGlobal+740,*) 'MPI_Interp_send_r8 eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      DO iProc=1,TheArr % nbProc
        iProcOut=TheArr % ListProc(iProc)
        eSize=TheArr % ListNbNeed(iProc)
        eRankGlobal=TheArr%ListGlobalRankOutput(iProcOut)
#ifdef DEBUG
        WRITE(MyRankGlobal+740,*) 'iProc=',iProc,'iProcOut=',iProcOut
        FLUSH(MyRankGlobal+740)
        WRITE(MyRankGlobal+740,*) 'eSize=', eSize, 'eRankGlobal=', eRankGlobal
        FLUSH(MyRankGlobal+740)
#endif
        allocate(rbuf_real_r8(eSize))
        DO idx=1,eSize
          iColRed=TheArr%ListIdx(idx,iProc)
          rbuf_real_r8(idx)=TheFieldIn(iColRed)
        END DO
#ifdef DEBUG
        rbuf_int(1)=8
        rbuf_int(2)=eSize
        rbuf_int(3)=eTag
        CALL MPI_SEND(rbuf_int,3,MPI_INTEGER, eRankGlobal, 144+eTag, TheArr%comm, ierr)
        WRITE(MyRankGlobal+740,*) 'After send 144 in MPI_INTERP_send_r8'
        FLUSH(MyRankGlobal+740)
#endif
        CALL MPI_SEND(rbuf_real_r8,eSize,MPI_DOUBLE_PRECISION, eRankGlobal, eTag, TheArr%comm, ierr)
        IF (ierr .ne. MPI_SUCCESS) THEN
           Print *, 'Error in MPI_INTERP_SEND_r8'
           STOP
        END IF
#ifdef DEBUG
        WRITE(MyRankGlobal+740,*) 'After send in MPI_INTERP_send_r8, eTag=', eTag
        FLUSH(MyRankGlobal+740)
#endif
        DEALLOCATE(rbuf_real_r8)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_RECV_3D_r4(TheArr, eTag, Ntimes, TheFieldOut)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: eTag, Ntimes
      real, intent(out) :: TheFieldOut(Ntimes,TheArr % RestrictedSparseMat % nRows)
      !
      real :: eField(Ntimes, TheArr % nbNeedTot)
      integer :: status(MPI_STATUS_SIZE)
      integer ierror, nInput, nOutput, iProc, iProcIn, eSize, iTime
      real, allocatable :: rbuf_real_r4(:)
      integer eRankGlobal, iColRed, iRowRed, ielt, numelt, i, idx
      real eWeight
#ifdef DEBUG
      integer :: rbuf_int(4)
#endif
      nInput=TheArr%RestrictedSparseMat%nColumns
      nOutput=TheArr%RestrictedSparseMat%nRows
      DO iProc=1,TheArr % nbProc
        iProcIn=TheArr % ListProc(iProc)
        eSize=TheArr % ListNbNeed(iProc)
        ALLOCATE(rbuf_real_r4(eSize))
        eRankGlobal=TheArr%ListGlobalRankInput(iProcIn)
#ifdef DEBUG
        CALL MPI_RECV(rbuf_int,4,MPI_INTEGER,eRankGlobal, 94+eTag, TheArr%comm, status, ierror)
        IF ((rbuf_int(1).ne.4).or.(rbuf_int(2).ne.eSize).or.            &
     &      (rbuf_int(3).ne.Ntimes).or.(rbuf_int(4).ne.eTag)) THEN
          Print *, 'RECV_3D_r4, eTag=', eTag
          Print *, 'Arithmetic rbuf_int(1)=', rbuf_int(1), ' vs 4'
          Print *, '      size rbuf_int(2)=', rbuf_int(2), ' vs ', eSize
          Print *, '    Ntimes rbuf_int(3)=', rbuf_int(3), ' vs ', Ntimes
          Print *, '       tag rbuf_int(4)=', rbuf_int(4), ' vs ', eTag
          STOP
        END IF
#endif
        CALL MPI_RECV(rbuf_real_r4,eSize*Ntimes,MPI_REAL,eRankGlobal, eTag, TheArr%comm, status, ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
          Print *, 'Error in MPI_INTERP_RECV_3D_r4'
          STOP
        END IF
        idx=0
        DO i=1,eSize
          iColRed=TheArr%ListIdx(i,iProc)
          DO iTime=1,Ntimes
            idx=idx+1
            eField(iTime,iColRed)=rbuf_real_r4(idx)
          END DO
        END DO
        DEALLOCATE(rbuf_real_r4)
      END DO
      TheFieldOut=0
      numelt=TheArr%RestrictedSparseMat%numelt
      DO ielt=1,numelt
        iColRed=TheArr%RestrictedSparseMat%columns(ielt)
        iRowRed=TheArr%RestrictedSparseMat%rows(ielt)
        eWeight=TheArr%RestrictedSparseMat%weights(ielt)
        TheFieldOut(:,iRowRed)=TheFieldOut(:,iRowRed)+eWeight*eField(:,iColRed)
      END DO
      END SUBROUTINE
!**********************************************************************
!* Some MPI_BARRIER are just not reliable. This construction makes    *
!* that every process receive and send to every other process         *
!**********************************************************************
      SUBROUTINE MYOWN_MPI_BARRIER(comm)
      IMPLICIT NONE
      integer, intent(in) :: comm
      integer eInt(1)
      integer iRank, jRank, eTag, nproc, myrank, ierr
      integer status(MPI_STATUS_SIZE)
      CALL MPI_COMM_SIZE(comm, nproc, ierr)
      CALL MPI_COMM_rank(comm, myrank, ierr)
      eInt(1)=4
      DO iRank=0,nproc-1
        eTag=137 + iRank
        IF (myrank .eq. iRank) THEN
          DO jRank=0,nproc-1
            IF (iRank.ne. jRank) THEN
              CALL MPI_SEND(eInt,1,MPI_INT,jRank,eTag,comm,ierr)
            END IF
          END DO
        ELSE
          CALL MPI_RECV(eInt, 1, MPI_INT,iRank,eTag,comm,status,ierr)
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_RECV_3D_r8(TheArr, eTag, Ntimes, TheFieldOut)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: eTag, Ntimes
      real*8, intent(out) :: TheFieldOut(Ntimes,TheArr % RestrictedSparseMat % nRows)
      !
      real*8 :: eField(Ntimes,TheArr % nbNeedTot)
      integer :: status(MPI_STATUS_SIZE)
      integer ierror, nInput, nOutput, iProc, iProcIn, eSize, idx
      real*8, allocatable :: rbuf_real_r8(:)
      integer eRankGlobal, iColRed, iRowRed, ielt, numelt, i, iTime
      real eWeight
#ifdef DEBUG
      integer :: rbuf_int(4)
#endif
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'B MPI_Interp_recv_3D_r8 eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      nInput=TheArr % RestrictedSparseMat % nColumns
      nOutput=TheArr % RestrictedSparseMat % nRows
      DO iProc=1,TheArr % nbProc
        iProcIn=TheArr % ListProc(iProc)
        eSize=TheArr % ListNbNeed(iProc)
        eRankGlobal=TheArr%ListGlobalRankInput(iProcIn)
#ifdef DEBUG
        WRITE(MyRankGlobal+740,*) 'iProc=',iProc,'iProcIn=',iProcIn
        FLUSH(MyRankGlobal+740)
        WRITE(MyRankGlobal+740,*) 'eSize=', eSize, 'eRankGlobal=', eRankGlobal
        FLUSH(MyRankGlobal+740)
#endif
        ALLOCATE(rbuf_real_r8(eSize*Ntimes))
#ifdef DEBUG
        CALL MPI_RECV(rbuf_int,4,MPI_INTEGER,eRankGlobal, 104+eTag, TheArr%comm, status, ierror)
        IF ((rbuf_int(1).ne.8).or.(rbuf_int(2).ne.eSize).or.            &
     &     (rbuf_int(3).ne.Ntimes).or.(rbuf_int(4).ne.eTag)) THEN
          Print *, 'RECV_3D_r8, eTag=', eTag
          Print *, 'Arithmetic rbuf_int(1)=', rbuf_int(1), ' vs 8'
          Print *, '      size rbuf_int(2)=', rbuf_int(2), ' vs ', eSize
          Print *, '    Ntimes rbuf_int(3)=', rbuf_int(3), ' vs ', Ntimes
          Print *, '       tag rbuf_int(4)=', rbuf_int(4), ' vs ', eTag
          STOP
        END IF
#endif
        CALL MPI_RECV(rbuf_real_r8,eSize*Ntimes,MPI_DOUBLE_PRECISION,eRankGlobal, eTag, TheArr%comm, status, ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
          Print *, 'Error in MPI_INTERP_RECV_3D_r8'
          STOP
        END IF
        idx=0
        DO i=1,eSize
          iColRed=TheArr%ListIdx(i,iProc)
          DO iTime=1,Ntimes
            idx=idx+1
            eField(iTime,iColRed)=rbuf_real_r8(idx)
          END DO
        END DO
        DEALLOCATE(rbuf_real_r8)
      END DO
      TheFieldOut=0
      numelt=TheArr%RestrictedSparseMat%numelt
      DO ielt=1,numelt
        iColRed=TheArr%RestrictedSparseMat%columns(ielt)
        iRowRed=TheArr%RestrictedSparseMat%rows(ielt)
        eWeight=TheArr%RestrictedSparseMat%weights(ielt)
        TheFieldOut(:,iRowRed)=TheFieldOut(:,iRowRed)+DBLE(eWeight)*eField(:,iColRed)
      END DO
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'E MPI_Interp_recv_3D_r8 eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_ARECV_3D_r8(TheArr, eTag, TheFieldOut)
      implicit none
      type(T_sparse_splitIO_as), intent(in) :: TheArr
      integer, intent(in) :: eTag
      real*8, intent(out) :: TheFieldOut(TheArr%Ntimes, TheArr % RestrictedSparseMat % nRows)
      !
      real*8 :: eField(TheArr % Ntimes, TheArr % nbNeedTot)
      integer :: status(MPI_STATUS_SIZE)
      integer ierror, nInput, nOutput, iProc, iProcIn, idx
      integer eRankGlobal, iColRed, iRowRed, ielt, numelt, i, iTime
      real eWeight
#ifdef DEBUG
      real*8 MinFieldIn(TheArr% Ntimes)
      real*8 MaxFieldIn(TheArr% Ntimes)
      integer :: rbuf_int(4), myrank
#endif
      nInput=TheArr % RestrictedSparseMat % nColumns
      nOutput=TheArr % RestrictedSparseMat % nRows
#ifdef DEBUG
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
      WRITE(MyRankGlobal+740,*) 'B MPI_INTERP_ARECV_3D_r8 ',  myrank, ' eTag=', eTag
      WRITE(MyRankGlobal+740,*) 'Ntimes=', TheArr % Ntimes
      WRITE(MyRankGlobal+740,*) 'nbProc=', TheArr % nbProc
      FLUSH(MyRankGlobal+740)
      DO iProc=1,TheArr % nbProc
        iProcIn=TheArr % ListProc(iProc)
        eRankGlobal=TheArr%ListGlobalRankInput(iProcIn)
        WRITE(MyRankGlobal+740,*) 'iProc=', iProc, ' iProcIn=', iProcIn, ' eRankGlobal=', eRankGlobal
        FLUSH(MyRankGlobal+740)
        CALL MPI_RECV(rbuf_int,4,MPI_INTEGER,eRankGlobal, 204+eTag,TheArr%comm, status, ierror)
        IF (ierror /= 0) THEN
          WRITE(MyRankGlobal+740,*) 'Error in MPI_RECV'
          FLUSH(MyRankGlobal+740)
        END IF
        IF ((rbuf_int(1).ne.8).or.(rbuf_int(2).ne.TheArr % Ntimes).or.(rbuf_int(3).ne.eTag).or.(rbuf_int(4).ne.TheArr % as_nbrecv(iProc))) THEN
          Print *, 'RECV_3D_r8, eTag=', eTag
          Print *, 'Arithmetic rbuf_int(1)=', rbuf_int(1), ' vs 8'
          Print *, '    Ntimes rbuf_int(3)=', rbuf_int(2), ' vs ', TheArr%Ntimes
          Print *, '       tag rbuf_int(3)=', rbuf_int(3), ' vs ', eTag
          Print *, '      size rbuf_int(4)=', rbuf_int(3), ' vs ', TheArr % as_nbrecv(iProc)
          STOP
        END IF
        WRITE(MyRankGlobal+740,*) 'After the tests'
        FLUSH(MyRankGlobal+740)
      END DO
      WRITE(MyRankGlobal+740,*) 'After MPI_RECV myrank=',  myrank, ' eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      DO iProc=1,TheArr % nbProc
        iProcIn=TheArr % ListProc(iProc)
        eRankGlobal=TheArr%ListGlobalRankInput(iProcIn)
        CALL MPI_IRECV(eField,1,TheArr%as_recv_type(iProc),eRankGlobal, eTag, TheArr%comm, TheArr%as_recv_rqst(iProc),ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
          Print *, 'Error in MPI_INTERP_ARECV_3D_r8'
          STOP
        END IF
      END DO
      IF (TheArr % nbProc .gt. 0) THEN
        call mpi_waitall(TheArr%nbProc, TheArr%as_recv_rqst, TheArr%as_recv_stat, ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
          Print *, 'Error in MPI_INTERP_ARECV_3D_r8'
          STOP
        END IF
      ENDIF
      TheFieldOut=0
      numelt=TheArr%RestrictedSparseMat%numelt
#ifdef DEBUG
      MinFieldIn = 1.0/TINY(1.)
      MaxFieldIn = -MinFieldIn
#endif
      DO ielt=1,numelt
        iColRed=TheArr%RestrictedSparseMat%columns(ielt)
        iRowRed=TheArr%RestrictedSparseMat%rows(ielt)
        eWeight=TheArr%RestrictedSparseMat%weights(ielt)
        TheFieldOut(:,iRowRed)=TheFieldOut(:,iRowRed)+DBLE(eWeight)*eField(:,iColRed)
#ifdef DEBUG
        MinFieldIn = MIN(MinFieldIn, eField(:,iColRed))
        MaxFieldIn = MAX(MaxFieldIn, eField(:,iColRed))
#endif
      END DO
#ifdef DEBUG
      DO iTime=1,TheArr% Ntimes
        WRITE(MyRankGlobal+740,*) 'i, Field(min/max)=', iTime, MinFieldIn(iTime), MaxFieldIn(iTime)
      END DO
      WRITE(MyRankGlobal+740,*) 'Weight(min/max)=', minval(TheArr%RestrictedSparseMat%weights), maxval(TheArr%RestrictedSparseMat%weights)
      WRITE(MyRankGlobal+740,*) 'E MPI_INTERP_ARECV_3D_r8 ', myrank, ' eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_ARECV_3D_r4(TheArr, eTag, TheFieldOut)
      implicit none
      type(T_sparse_splitIO_as), intent(in) :: TheArr
      integer, intent(in) :: eTag
      real, intent(out) :: TheFieldOut(TheArr%Ntimes, TheArr % RestrictedSparseMat % nRows)
      !
      real :: eField(TheArr % Ntimes, TheArr % nbNeedTot)
      integer :: status(MPI_STATUS_SIZE)
      integer ierror, nInput, nOutput, iProc, iProcIn, idx
      integer eRankGlobal, iColRed, iRowRed, ielt, numelt, i, iTime
      integer ierr
      real eWeight
#ifdef DEBUG
      integer :: rbuf_int(4), myrank
#endif
      nInput=TheArr % RestrictedSparseMat % nColumns
      nOutput=TheArr % RestrictedSparseMat % nRows
#ifdef DEBUG
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
      WRITE(MyRankGlobal+740,*) 'B MPI_INTERP_ARECV_3D_r4 ',  myrank, ' eTag=', eTag
      WRITE(MyRankGlobal+740,*) 'Ntimes=', TheArr % Ntimes
      WRITE(MyRankGlobal+740,*) 'nbProc=', TheArr % nbProc
      FLUSH(MyRankGlobal+740)
      DO iProc=1,TheArr % nbProc
        iProcIn=TheArr % ListProc(iProc)
        eRankGlobal=TheArr%ListGlobalRankInput(iProcIn)
        WRITE(MyRankGlobal+740,*) 'iProc=', iProc, ' iProcIn=', iProcIn, ' eRankGlobal=', eRankGlobal
        FLUSH(MyRankGlobal+740)
        CALL MPI_RECV(rbuf_int,4,MPI_INTEGER,eRankGlobal, 144+eTag,TheArr%comm, status, ierror)
        WRITE(MyRankGlobal+740,*) 'After MPI_RECV'
        FLUSH(MyRankGlobal+740)
        IF (ierror /= 0) THEN
          WRITE(MyRankGlobal+740,*) 'Error in MPI_RECV'
          FLUSH(MyRankGlobal+740)
          STOP
        END IF
        IF ((rbuf_int(1).ne.4).or.(rbuf_int(2).ne.TheArr % Ntimes).or.(rbuf_int(3).ne.eTag).or.(rbuf_int(4).ne.TheArr % as_nbrecv(iProc))) THEN
          Print *, 'RECV_3D_r4, eTag=', eTag
          Print *, 'Arithmetic rbuf_int(1)=', rbuf_int(1), ' vs 4'
          Print *, '    Ntimes rbuf_int(3)=', rbuf_int(2), ' vs ', TheArr%Ntimes
          Print *, '       tag rbuf_int(3)=', rbuf_int(3), ' vs ', eTag
          Print *, '      size rbuf_int(4)=', rbuf_int(3), ' vs ', TheArr % as_nbrecv(iProc)
          STOP
        END IF
        WRITE(MyRankGlobal+740,*) 'Before MPI_SEND'
        FLUSH(MyRankGlobal+740)
        CALL MPI_SEND(rbuf_int,4,MPI_INTEGER, eRankGlobal, 144+eTag,TheArr%comm, ierr)
        WRITE(MyRankGlobal+740,*) 'After MPI_SEND'
        WRITE(MyRankGlobal+740,*) 'After the tests'
        FLUSH(MyRankGlobal+740)
      END DO
      WRITE(MyRankGlobal+740,*) 'After MPI_RECV myrank=',  myrank, ' eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      DO iProc=1,TheArr % nbProc
        iProcIn=TheArr % ListProc(iProc)
        eRankGlobal=TheArr%ListGlobalRankInput(iProcIn)
        CALL MPI_IRECV(eField,1,TheArr%as_recv_type(iProc),eRankGlobal, eTag, TheArr%comm, TheArr%as_recv_rqst(iProc),ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
          Print *, 'Error in MPI_INTERP_ARECV_3D_r4'
          STOP
        END IF
      END DO
      IF (TheArr % nbProc .gt. 0) THEN
        call mpi_waitall(TheArr%nbProc, TheArr%as_recv_rqst, TheArr%as_recv_stat, ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
          Print *, 'Error in MPI_INTERP_ARECV_3D_r4'
          STOP
        END IF
      ENDIF
      TheFieldOut=0
      numelt=TheArr%RestrictedSparseMat%numelt
      DO ielt=1,numelt
        iColRed=TheArr%RestrictedSparseMat%columns(ielt)
        iRowRed=TheArr%RestrictedSparseMat%rows(ielt)
        eWeight=TheArr%RestrictedSparseMat%weights(ielt)
        TheFieldOut(:,iRowRed)=TheFieldOut(:,iRowRed)+DBLE(eWeight)*eField(:,iColRed)
      END DO
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'E MPI_INTERP_ARECV_3D_r4 ', myrank, ' eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_SEND_3D_r4(TheArr, eTag, Ntimes, TheFieldIn)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: eTag, Ntimes
      real, intent(in) :: TheFieldIn(Ntimes,TheArr % RestrictedSparseMat % nColumns)
      integer ierr, iProc, eSize
      integer iProcOut, idx, i, iColRed, eRankGlobal, iTime
      real, allocatable :: rbuf_real_r4(:)
#ifdef DEBUG
      integer :: rbuf_int(4)
#endif
      DO iProc=1,TheArr % nbProc
        iProcOut=TheArr % ListProc(iProc)
        eSize=TheArr % ListNbNeed(iProc)
        allocate(rbuf_real_r4(eSize*Ntimes))
        idx=0
        DO i=1,eSize
          iColRed=TheArr%ListIdx(i,iProc)
          DO iTime=1,Ntimes
            idx=idx+1
            rbuf_real_r4(idx)=TheFieldIn(iTime,iColRed)
          END DO
        END DO
        eRankGlobal=TheArr%ListGlobalRankOutput(iProcOut)
#ifdef DEBUG
        rbuf_int(1)=4
        rbuf_int(2)=eSize
        rbuf_int(3)=Ntimes
        rbuf_int(4)=eTag
        CALL MPI_SEND(rbuf_int,4,MPI_INTEGER, eRankGlobal, 94+eTag,TheArr%comm, ierr)
#endif
        CALL MPI_SEND(rbuf_real_r4,eSize*Ntimes,MPI_REAL, eRankGlobal, eTag, TheArr%comm, ierr)
        IF (ierr .ne. MPI_SUCCESS) THEN
           Print *, 'Error in MPI_INTERP_SEND_3D_r4'
           STOP
        END IF
        DEALLOCATE(rbuf_real_r4)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_SEND_3D_r8(TheArr, eTag, Ntimes, TheFieldIn)
      implicit none
      type(T_sparse_splitIO), intent(in) :: TheArr
      integer, intent(in) :: eTag, Ntimes
      real*8, intent(in) :: TheFieldIn(Ntimes,TheArr % RestrictedSparseMat % nColumns)
      integer ierr, iProc, eSize
      integer iProcOut, idx, i, iColRed, eRankGlobal, iTime
      real*8, allocatable :: rbuf_real_r8(:)
#ifdef DEBUG
      integer :: rbuf_int(4)
      WRITE(MyRankGlobal+740,*) 'B MPI_Interp_send_3d_r8 eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      DO iProc=1,TheArr % nbProc
        iProcOut=TheArr % ListProc(iProc)
        eSize=TheArr % ListNbNeed(iProc)
        eRankGlobal=TheArr%ListGlobalRankOutput(iProcOut)
#ifdef DEBUG
        WRITE(MyRankGlobal+740,*) 'iProc=',iProc,'iProcOut=',iProcOut
        FLUSH(MyRankGlobal+740)
        WRITE(MyRankGlobal+740,*) 'eSize=', eSize, 'eRankGlobal=', eRankGlobal
        FLUSH(MyRankGlobal+740)
#endif
        allocate(rbuf_real_r8(eSize*Ntimes))
        idx=0
        DO i=1,eSize
          iColRed=TheArr%ListIdx(i,iProc)
          DO iTime=1,Ntimes
            idx=idx+1
            rbuf_real_r8(idx)=TheFieldIn(iTime,iColRed)
          END DO
        END DO
#ifdef DEBUG
        rbuf_int(1)=8
        rbuf_int(2)=eSize
        rbuf_int(3)=Ntimes
        rbuf_int(4)=eTag
        CALL MPI_SEND(rbuf_int,4,MPI_INTEGER, eRankGlobal, 104+eTag, TheArr%comm, ierr)
#endif
        CALL MPI_SEND(rbuf_real_r8,eSize*Ntimes,MPI_DOUBLE_PRECISION, eRankGlobal, eTag, TheArr%comm, ierr)
        IF (ierr .ne. MPI_SUCCESS) THEN
           Print *, 'Error in MPI_INTERP_SEND_3D_r8'
           STOP
        END IF
        DEALLOCATE(rbuf_real_r8)
      END DO
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'B MPI_Interp_send_3d_r8 eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_ASEND_3D_r8(TheArr, eTag, TheFieldIn)
      implicit none
      type(T_sparse_splitIO_as), intent(in) :: TheArr
      integer, intent(in) :: eTag
      real*8, intent(in) :: TheFieldIn(TheArr % Ntimes,TheArr % RestrictedSparseMat % nColumns)
      integer ierror, iProc, eSize
      integer iProcOut, eRankGlobal
#ifdef DEBUG
      integer rbuf_int(4), myrank
      WRITE(MyRankGlobal+740,*) 'TheFieldIn(min/max)=', minval(TheFieldIn), maxval(TheFieldIn)
      FLUSH(MyRankGlobal+740)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
      WRITE(MyRankGlobal+740,*) 'B MPI_INTERP_ASEND_3D_r8 ', myrank, ' eTag=', eTag
      WRITE(MyRankGlobal+740,*) 'Ntimes=', TheArr % Ntimes
      FLUSH(MyRankGlobal+740)
      DO iProc=1,TheArr % nbProc
        iProcOut=TheArr % ListProc(iProc)
        eRankGlobal=TheArr%ListGlobalRankOutput(iProcOut)
        WRITE(MyRankGlobal+740,*) 'iProc=', iProc, ' iProcOut=', iProcOut, ' eRankGlobal=', eRankGlobal
        FLUSH(MyRankGlobal+740)
        rbuf_int(1)=8
        rbuf_int(2)=TheArr % Ntimes
        rbuf_int(3)=eTag
        rbuf_int(4)=TheArr % as_nbsend(iProc)
        CALL MPI_SEND(rbuf_int,4,MPI_INTEGER, eRankGlobal, 204+eTag,TheArr%comm, ierror)
        IF (ierror /= 0) THEN
          WRITE(MyRankGlobal+740,*) 'Error in MPI_SEND'
          FLUSH(MyRankGlobal+740)
        END IF
      END DO
      WRITE(MyRankGlobal+740,*) 'After MPI_SEND myrank=', myrank
      FLUSH(MyRankGlobal+740)
#endif
      DO iProc=1,TheArr % nbProc
        iProcOut=TheArr % ListProc(iProc)
        eRankGlobal=TheArr%ListGlobalRankOutput(iProcOut)
        CALL MPI_ISEND(TheFieldIn,1,TheArr%as_send_type(iProc), eRankGlobal, eTag, TheArr%comm, TheArr%as_send_rqst(iProc),ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
           Print *, 'Error in MPI_INTERP_ASEND_3D_r8'
           STOP
        END IF
      END DO
      IF (TheArr % nbProc .gt. 0) THEN
        call mpi_waitall(TheArr%nbProc,TheArr % as_send_rqst, TheArr%as_send_stat, ierror)
      END IF
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'E MPI_INTERP_ASEND_3D_r8, ', myrank, ' eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_INTERP_ASEND_3D_r4(TheArr, eTag, TheFieldIn)
      implicit none
      type(T_sparse_splitIO_as), intent(in) :: TheArr
      integer, intent(in) :: eTag
      real, intent(in) :: TheFieldIn(TheArr % Ntimes,TheArr % RestrictedSparseMat % nColumns)
      integer ierror, iProc, eSize
      integer iProcOut, eRankGlobal
      integer :: status(MPI_STATUS_SIZE)
#ifdef DEBUG
      integer rbuf_int(4), myrank
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
      WRITE(MyRankGlobal+740,*) 'B MPI_INTERP_ASEND_3D_r4 ', myrank, ' eTag=', eTag
      WRITE(MyRankGlobal+740,*) 'Ntimes=', TheArr % Ntimes
      FLUSH(MyRankGlobal+740)
      DO iProc=1,TheArr % nbProc
        iProcOut=TheArr % ListProc(iProc)
        eRankGlobal=TheArr%ListGlobalRankOutput(iProcOut)
        WRITE(MyRankGlobal+740,*) 'iProc=', iProc, ' iProcOut=', iProcOut, ' eRankGlobal=', eRankGlobal
        FLUSH(MyRankGlobal+740)
        rbuf_int(1)=4
        rbuf_int(2)=TheArr % Ntimes
        rbuf_int(3)=eTag
        rbuf_int(4)=TheArr % as_nbsend(iProc)
        CALL MPI_SEND(rbuf_int,4,MPI_INTEGER, eRankGlobal, 144+eTag,TheArr%comm, ierror)
        WRITE(MyRankGlobal+740,*) 'After the MPI_SEND'
        FLUSH(MyRankGlobal+740)
        IF (ierror /= 0) THEN
          WRITE(MyRankGlobal+740,*) 'Error in MPI_SEND'
          FLUSH(MyRankGlobal+740)
          STOP
        END IF
        CALL MPI_RECV(rbuf_int,4,MPI_INTEGER,eRankGlobal, 144+eTag,TheArr%comm, status, ierror)
        IF (ierror /= 0) THEN
          WRITE(MyRankGlobal+740,*) 'Error in MPI_RECV'
          FLUSH(MyRankGlobal+740)
          STOP
        END IF
        WRITE(MyRankGlobal+740,*) 'After the SEND/RECV of checks'
        FLUSH(MyRankGlobal+740)
      END DO
      WRITE(MyRankGlobal+740,*) 'After MPI_SEND myrank=', myrank
      FLUSH(MyRankGlobal+740)
#endif
      DO iProc=1,TheArr % nbProc
        iProcOut=TheArr % ListProc(iProc)
        eRankGlobal=TheArr%ListGlobalRankOutput(iProcOut)
        CALL MPI_ISEND(TheFieldIn,1,TheArr%as_send_type(iProc), eRankGlobal, eTag, TheArr%comm, TheArr%as_send_rqst(iProc),ierror)
        IF (ierror.ne.MPI_SUCCESS) THEN
           Print *, 'Error in MPI_INTERP_ASEND_3D_r4'
           STOP
        END IF
      END DO
      IF (TheArr % nbProc .gt. 0) THEN
        call mpi_waitall(TheArr%nbProc,TheArr % as_send_rqst, TheArr%as_send_stat, ierror)
      END IF
#ifdef DEBUG
      WRITE(MyRankGlobal+740,*) 'E MPI_INTERP_ASEND_3D_r4, ', myrank, ' eTag=', eTag
      FLUSH(MyRankGlobal+740)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOCATE_Arr(TheArr)
      implicit none
      type(T_sparse_splitIO), intent(inout) :: TheArr
      IF (TheArr % TheType.eq.1) THEN
        DEALLOCATE(TheArr % RestrictedSparseMat % rows)
        DEALLOCATE(TheArr % RestrictedSparseMat % columns)
        DEALLOCATE(TheArr % RestrictedSparseMat % weights)
      END IF
      DEALLOCATE(TheArr % ListProc)
      DEALLOCATE(TheArr % ListNbNeed)
      DEALLOCATE(TheArr % ListIdx)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOCATE_ArrAsync(TheArr)
      implicit none
      type(T_sparse_splitIO_as), intent(inout) :: TheArr
      IF (TheArr % TheType.eq.1) THEN
        DEALLOCATE(TheArr % RestrictedSparseMat % rows)
        DEALLOCATE(TheArr % RestrictedSparseMat % columns)
        DEALLOCATE(TheArr % RestrictedSparseMat % weights)
        DEALLOCATE(TheArr % as_recv_rqst)
        DEALLOCATE(TheArr % as_recv_type)
        DEALLOCATE(TheArr % as_recv_stat)
      ELSE
        DEALLOCATE(TheArr % as_send_rqst)
        DEALLOCATE(TheArr % as_send_type)
        DEALLOCATE(TheArr % as_send_stat)
      END IF
      DEALLOCATE(TheArr % ListProc)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SetComputationalNodes(ArrLocal, MyNbCompNode, eDestColor)
      type(T_localNodes), intent(inout) :: ArrLocal
      integer, intent(in) :: MyNbCompNode
      integer, intent(in) :: eDestColor
      integer :: rbuf_recv(1)
      integer :: rbuf_send(1)
      integer status(MPI_STATUS_SIZE)
      integer ierr, iProc, dest, eOrig, NbDestNodes, nbCompNode
      integer MyColor, MyLocalRank, MyGlobalRank
      MyColor=ArrLocal % MyColor
      MyLocalRank=ArrLocal % MyLocalRank
      MyGlobalRank=ArrLocal % rank_in_comm
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'SetCompNode: MyColor=', MyColor
      WRITE(740+MyRankGlobal,*) 'SetCompNode: MyLocalRank=', MyLocalRank
      WRITE(740+MyRankGlobal,*) 'SetCompNode: MyGlobalRank=', MyGlobalRank
      WRITE(740+MyRankGlobal,*) 'eDestColor=', eDestColor
      FLUSH(740+MyRankGlobal)
#endif
      IF (MyLocalRank.eq.0) THEN
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'Case 1'
        FLUSH(740+MyRankGlobal)
#endif
        rbuf_send(1)=MyNbCompNode
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'SetCompNode: We have rbuf_send'
        FLUSH(740+MyRankGlobal)
#endif
        dest=ArrLocal % ListFirstRank(eDestColor)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'SetCompNode: We have dest'
        WRITE(740+MyRankGlobal,*) 'MyNbCompNode=', MyNbCompNode
        WRITE(740+MyRankGlobal,*) 'dest=', dest
        FLUSH(740+MyRankGlobal)
#endif
        CALL MPI_SENDRECV(rbuf_send, 1, MPI_INTEGER, dest, 100, rbuf_recv, 1, MPI_INTEGER, dest, 100, ArrLocal%comm, status, ierr)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After MPI_SENDRECV'
        FLUSH(740+MyRankGlobal)
#endif
        NbDestNodes=rbuf_recv(1)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After NbDestNodes assignation'
        FLUSH(740+MyRankGlobal)
#endif
        rbuf_send(1)=NbDestNodes
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After rbuf_send assignation'
        WRITE(740+MyRankGlobal,*) 'NbDestNodes=', NbDestNodes
        FLUSH(740+MyRankGlobal)
#endif
        DO iProc=2,MyNbCompNode
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'Before GetGlobalRank'
          FLUSH(740+MyRankGlobal)
#endif
          dest=GetGlobalRank(ArrLocal, MyColor, iProc)
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After GetGlobalRank'
          FLUSH(740+MyRankGlobal)
#endif
          CALL MPI_SEND(rbuf_send,1,MPI_INTEGER, dest, 216, ArrLocal%comm, ierr)
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After MPI_SEND'
          FLUSH(740+MyRankGlobal)
#endif
        END DO
      ELSE
        eOrig=ArrLocal % ListFirstRank(MyColor)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After obtain eOrig'
        FLUSH(740+MyRankGlobal)
        WRITE(740+MyRankGlobal,*) 'eOrig=', eOrig
        FLUSH(740+MyRankGlobal)
#endif
        CALL MPI_RECV(rbuf_recv,1,MPI_INTEGER, eOrig, 216, ArrLocal%comm, status, ierr)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After MPI_RECV'
        FLUSH(740+MyRankGlobal)
#endif
        NbDestNodes=rbuf_recv(1)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After getting NbDestNodes'
        FLUSH(740+MyRankGlobal)
#endif
      ENDIF
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'SetCompNode: After the Main operation'
      FLUSH(740+MyRankGlobal)
#endif
      nbCompNode=ArrLocal % NbCompNodes(eDestColor)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) '1: nbCompNode=', nbCompNode
      FLUSH(740+MyRankGlobal)
#endif
      IF (nbCompNode.eq.-1) THEN
        ArrLocal % NbCompNodes(eDestColor)=NbDestNodes
      ELSE
        IF (nbCompNode.ne.NbDestNodes) THEN
          Print *, 'Inconsistency in number of computational nodes 1'
          STOP
        ENDIF
      ENDIF
      nbCompNode=ArrLocal % NbCompNodes(MyColor)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) '2: nbCompNode=', nbCompNode
      FLUSH(740+MyRankGlobal)
#endif
      IF (nbCompNode.eq.-1) THEN
        ArrLocal % NbCompNodes(MyColor)=MyNbCompNode
      ELSE
        IF (nbCompNode.ne.MyNbCompNode) THEN
          WRITE(740+MyRankGlobal,*) 'nbCompNode=', nbCompNode
          WRITE(740+MyRankGlobal,*) 'MyNbCompNode=', MyNbCompNode
          FLUSH(740+MyRankGlobal)
          Print *, 'nbCompNode=', nbCompNode
          Print *, 'MyNbCompNode=', MyNbCompNode
          Print *, 'Inconsistency in number of computational nodes 2'
          STOP
        ENDIF
      ENDIF
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Exiting SetComputationalNodes'
      FLUSH(740+MyRankGlobal)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GetAllLocalRanks(comm, MyLocalRank, MyColor, ArrLocal)
      implicit none
      integer, intent(in) :: comm
      integer, intent(in) :: MyLocalRank, MyColor
      type(T_localNodes), intent(inout) :: ArrLocal
      integer :: rbuf_int(2)
      integer :: status(MPI_STATUS_SIZE)
      integer MyError, ierror, ierr, iProc, MRG
      integer NProc, iColor, maxColor, NbLocalProc, NbTotalProc
      integer eRankLocal
      logical IsAssigned
      CALL mpi_comm_size (comm, NProc, MyError)
      CALL mpi_comm_rank (comm, MRG, ierror)
      MRG=MRG+1
      ArrLocal % comm=comm
      ArrLocal % NProc=NProc
      ArrLocal % rank_in_comm=MRG-1
      allocate(ArrLocal % ListLocalRank(NProc))
      allocate(ArrLocal % ListColor(NProc))
      IF (MRG .eq. 1) THEN
        ArrLocal % ListLocalRank(1)=MyLocalRank
        ArrLocal % ListColor(1)=MyColor
        DO iProc=2,NProc
          CALL MPI_RECV(rbuf_int,2,MPI_INTEGER, iProc-1, 39, comm, status, ierror)
          ArrLocal % ListLocalRank(iProc)=rbuf_int(1)
          ArrLocal % ListColor(iProc)=rbuf_int(2)
        END DO
        DO iProc=2,NProc
          CALL MPI_SEND(ArrLocal % ListLocalRank,NProc,MPI_INTEGER, iProc-1, 42, comm, ierr)
          CALL MPI_SEND(ArrLocal % ListColor    ,NProc,MPI_INTEGER, iProc-1, 43, comm, ierr)
        END DO
      ELSE
        rbuf_int(1)=MyLocalRank
        rbuf_int(2)=MyColor
        CALL MPI_SEND(rbuf_int,2,MPI_INTEGER, 0, 39, comm, ierr)
        CALL MPI_RECV(ArrLocal % ListLocalRank,NProc,MPI_INTEGER, 0, 42, comm, status, ierror)
        CALL MPI_RECV(ArrLocal % ListColor    ,NProc,MPI_INTEGER, 0, 43, comm, status, ierror)
      END IF
      maxColor=maxval(ArrLocal % ListColor)
      ArrLocal % maxColor=maxColor
      allocate(ArrLocal % ListMapGlobal(maxColor))
      allocate(ArrLocal % ListFirstRank(maxColor))
      NbTotalProc=0
      DO iColor=1,maxcolor
        NbLocalProc=0
        DO iProc=1,NProc
          IF (ArrLocal % ListColor(iProc).eq.iColor) THEN
            NbLocalProc=NbLocalProc+1
          END IF
        END DO
        NbTotalProc=NbTotalProc+NbLocalProc
        ArrLocal % ListMapGlobal(iColor) % eColor=iColor;
        ArrLocal % ListMapGlobal(iColor) % NbLocalProc=NbLocalProc

        IF (NbLocalProc .gt. 0) THEN
          allocate(ArrLocal % ListMapGlobal(iColor) % MapLocalToGlobal(NbLocalProc))
          DO iProc=1,NProc
            IF (ArrLocal % ListColor(iProc).eq.iColor) THEN
              eRankLocal=ArrLocal % ListLocalRank(iProc)
              ArrLocal % ListMapGlobal(iColor) % MapLocalToGlobal(eRankLocal+1)=iProc-1
            END IF
          END DO
          ArrLocal % ListFirstRank(iColor)=ArrLocal % ListMapGlobal(iColor) % MapLocalToGlobal(1)
          IsAssigned=.TRUE.
        ELSE
          IsAssigned=.FALSE.
        END IF
        ArrLocal % ListMapGlobal(iColor) % IsAssigned=IsAssigned
      END DO
      allocate(ArrLocal % NbCompNodes(maxColor))
      ArrLocal % NbCompNodes=-1
      ArrLocal % MyLocalRank=MyLocalRank
      ArrLocal % MyColor=MyColor
#ifdef DEBUG
      IF (NbTotalProc.ne.NProc) THEN
        Print *, 'NbTotalProc=', NbTotalProc, 'NProc=', NProc
        Print *, 'Clear inconsistency here'
        STOP
      ENDIF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DeallocateArrLocal(ArrLocal)
      implicit none
      type(T_localNodes), intent(inout) :: ArrLocal
      integer maxColor, iColor
      maxColor=ArrLocal % maxColor
      deallocate(ArrLocal % ListLocalRank)
      deallocate(ArrLocal % ListColor)
      deallocate(ArrLocal % ListFirstRank)
      deallocate(ArrLocal % NbCompNodes)
      DO iColor=1,maxColor
        deallocate(ArrLocal % ListMapGlobal(iColor) % MapLocalToGlobal)
      END DO
      deallocate(ArrLocal % ListMapGlobal)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION GetGlobalRank(ArrLocal, eColor, eRankLocal)
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: eColor, eRankLocal
      integer :: iRank, iProc
      integer GetGlobalRank
      GetGlobalRank=ArrLocal % ListMapGlobal(eColor) % MapLocalToGlobal(eRankLocal)
      END FUNCTION GetGlobalRank
      !
      SUBROUTINE MPI_SEND_SPARSE_MATRIX(sMat, comm, myrankrecv, eTag)
      implicit none
      type(T_sparse), intent(in) :: sMat
      integer, intent(in) :: comm, myrankrecv, eTag
      !
      integer, allocatable :: rbuf_int(:)
      real*8, allocatable :: rbuf_real(:)
      integer ierr
      integer idx, i, numelt
      allocate(rbuf_int(3))
      numelt=sMat % numelt
      rbuf_int(1)=sMat % nRows
      rbuf_int(2)=sMat % nColumns
      rbuf_int(3)=sMat % numelt
      CALL MPI_SEND(rbuf_int,3,MPI_INTEGER, myrankrecv, eTag+1300, comm, ierr)
      deallocate(rbuf_int)
      allocate(rbuf_real(numelt))
      DO i=1,numelt
        rbuf_real(i)=sMat%weights(i)
      END DO
      CALL MPI_SEND(rbuf_real, numelt, MPI_DOUBLE_PRECISION, myrankrecv, eTag+1301, comm, ierr)
      deallocate(rbuf_real)
      allocate(rbuf_int(2*numelt))
      idx=0
      DO i=1,numelt
        idx=idx+1
        rbuf_int(idx)=sMat%rows(i)
      END DO
      DO i=1,numelt
        idx=idx+1
        rbuf_int(idx)=sMat%columns(i)
      END DO
      CALL MPI_SEND(rbuf_int,2*numelt,MPI_INTEGER, myrankrecv, eTag+1302, comm, ierr)
      deallocate(rbuf_int)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MPI_RECV_SPARSE_MATRIX(sMat, comm, myranksend, eTag)
      implicit none
      type(T_sparse), intent(inout) :: sMat
      integer, intent(in) :: comm, myranksend, eTag
      !
      integer, allocatable :: rbuf_int(:)
      real*8, allocatable :: rbuf_real(:)
      integer :: status(MPI_STATUS_SIZE)
      integer nRows, nColumns, numelt
      integer idx, i
      integer ierror
      allocate(rbuf_int(3))
      CALL MPI_RECV(rbuf_int,3,MPI_INTEGER, myranksend, eTag+1300, comm, status, ierror)
      nRows=rbuf_int(1)
      nColumns=rbuf_int(2)
      numelt=rbuf_int(3)
      deallocate(rbuf_int)
      sMat % nRows=nRows
      sMat % nColumns=nColumns
      sMat % numelt=numelt
      allocate(sMat % weights(numelt))
      allocate(sMat % rows(numelt))
      allocate(sMat % columns(numelt))
      allocate(rbuf_real(numelt))
      CALL MPI_RECV(rbuf_real,numelt,MPI_DOUBLE_PRECISION, myranksend, eTag+1301, comm, status, ierror)
      DO i=1,numelt
        sMat % weights(i)=rbuf_real(i)
      END DO
      deallocate(rbuf_real)
      !
      allocate(rbuf_int(2*numelt))
      CALL MPI_RECV(rbuf_int,2*numelt,MPI_INTEGER, myranksend, eTag+1302, comm, status, ierror)
      idx=0
      DO i=1,numelt
        idx=idx+1
        sMat % rows(i)=rbuf_int(idx)
      END DO
      DO i=1,numelt
        idx=idx+1
        sMat % columns(i)=rbuf_int(idx)
      END DO
      deallocate(rbuf_int)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE M2M_send_sparseMatrix(ArrLocal, TargetColor, sMat)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      type(T_sparse), intent(in) :: sMat
      integer, intent(in) :: TargetColor
      integer :: NnodesTarget, iNode, iProc, eTag
      NnodesTarget=ArrLocal % NbCompNodes(TargetColor)
      eTag=0
      DO iNode=1,NnodesTarget
        iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
        CALL MPI_SEND_SPARSE_MATRIX(sMat, ArrLocal % comm, iProc, eTag)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE M2M_recv_sparseMatrix(ArrLocal, eColorSend, sMat)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      type(T_sparse), intent(inout) :: sMat
      integer, intent(in) :: eColorSend
      integer :: RankFirstSend, eTag
      RankFirstSend=ArrLocal % ListFirstRank(eColorSend)
      eTag=0
      CALL MPI_RECV_SPARSE_MATRIX(sMat, ArrLocal % comm, RankFirstSend, eTag)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE M2M_send_grid(ArrLocal, TargetColor, eGrid)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: TargetColor
      type(T_grid) :: eGrid
      !
      integer idx, i, j, ierr, NnodesTarget
      integer, allocatable :: rbuf_int(:)
      real*8, allocatable :: rbuf_real_r8(:)
      integer ierror, iNode, iProc, Ntot
      integer Nnode_S, Ntrig_S
      integer xSize_S, ySize_S
      allocate(rbuf_int(1))
      rbuf_int(1)=eGrid % IsFE
      NnodesTarget=ArrLocal % NbCompNodes(TargetColor)
      DO iNode=1,NnodesTarget
        iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
        CALL MPI_SEND(rbuf_int,1,MPI_INTEGER, iProc, 2048, ArrLocal % comm, ierr)
      END DO
      deallocate(rbuf_int)
      !
      IF (eGrid % IsFE .eq. 1) THEN
        Nnode_S=eGrid % nbNode
        Ntrig_S=eGrid % nbTrig
        allocate(rbuf_int(2))
        rbuf_int(1)=Nnode_S
        rbuf_int(2)=Ntrig_S
        NnodesTarget=ArrLocal % NbCompNodes(TargetColor)
        DO iNode=1,NnodesTarget
          iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
          CALL MPI_SEND(rbuf_int,2,MPI_INTEGER, iProc, 1000, ArrLocal % comm, ierr)
        END DO
        deallocate(rbuf_int)
        !
        allocate(rbuf_real_r8(2*Nnode_S))
        idx=0
        DO i=1,Nnode_S
          idx=idx+1
          rbuf_real_r8(idx)=eGrid % LON_fe(i)
          idx=idx+1
          rbuf_real_r8(idx)=eGrid % LAT_fe(i)
        END DO
        DO iNode=1,NnodesTarget
          iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
          CALL MPI_SEND(rbuf_real_r8,2*Nnode_S,MPI_DOUBLE_PRECISION,iProc, 1001, ArrLocal % comm, ierr)
        END DO
        deallocate(rbuf_real_r8)
        allocate(rbuf_int(Ntrig_S*3))
        idx=0
        DO i=1,Ntrig_S
          DO j=1,3
            idx=idx+1
            rbuf_int(idx)=eGrid % ListTrig(j,i)
          END DO
        END DO
        DO iNode=1,NnodesTarget
          iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
          CALL MPI_SEND(rbuf_int,3*Ntrig_S,MPI_INTEGER, iProc, 1002, ArrLocal % comm, ierr)
        END DO
        deallocate(rbuf_int)
      ELSE
        xSize_S=eGrid % eta_rho
        ySize_S=eGrid % xi_rho
        allocate(rbuf_int(2))
        rbuf_int(1)=xSize_S
        rbuf_int(2)=ySize_S
        NnodesTarget=ArrLocal % NbCompNodes(TargetColor)
        DO iNode=1,NnodesTarget
          iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
          CALL MPI_SEND(rbuf_int,2,MPI_INTEGER, iProc, 2000, ArrLocal % comm, ierr)
        END DO
        deallocate(rbuf_int)
        !
        allocate(rbuf_real_r8(2*xSize_S*ySize_S))
        idx=0
        DO i=1,xSize_S
          DO j=1,ySize_S
            idx=idx+1
            rbuf_real_r8(idx)=eGrid % LON_fd(i,j)
            idx=idx+1
            rbuf_real_r8(idx)=eGrid % LAT_fd(i,j)
          END DO
        END DO
        DO iNode=1,NnodesTarget
          iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
          CALL MPI_SEND(rbuf_real_r8,2*xSize_S*ySize_S, MPI_DOUBLE_PRECISION, iProc, 2001,ArrLocal % comm,ierr)
        END DO
        deallocate(rbuf_real_r8)
        allocate(rbuf_int(xSize_S*ySize_S))
        idx=0
        DO i=1,xSize_S
          DO j=1,ySize_S
            idx=idx+1
            rbuf_int(idx)=eGRid % MSK_fd(i,j)
          END DO
        END DO
        DO iNode=1,NnodesTarget
          iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
          CALL MPI_SEND(rbuf_int,xSize_S*ySize_S,MPI_INTEGER, iProc, 2002, ArrLocal % comm, ierr)
        END DO
        deallocate(rbuf_int)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE M2M_recv_grid(ArrLocal, eColorSend, eGrid)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: eColorSend
      type(T_grid), intent(inout) :: eGrid
      !
      integer idx, i, j
      integer :: status(MPI_STATUS_SIZE)
      integer, allocatable :: rbuf_int(:)
      real*8, allocatable :: rbuf_real_r8(:)
      integer ierr, Ntot, RankFirstSend
      integer Nnode_R, Ntrig_R
      integer xSize_R, ySize_R
      RankFirstSend=ArrLocal % ListFirstRank(eColorSend)
      allocate(rbuf_int(1))
      CALL MPI_RECV(rbuf_int,1,MPI_INTEGER, RankFirstSend, 2048, ArrLocal % comm, status, ierr)
      eGrid % IsFE = rbuf_int(1)
      deallocate(rbuf_int)
      !
      IF (eGrid % IsFE .eq. 1) THEN
        allocate(rbuf_int(2))
        CALL MPI_RECV(rbuf_int,2,MPI_INTEGER, RankFirstSend, 1000, ArrLocal%comm, status, ierr)
        Nnode_R=rbuf_int(1)
        Ntrig_R=rbuf_int(2)
        deallocate(rbuf_int)
        eGrid % nbNode = Nnode_R
        eGrid % nbTrig = Ntrig_R
        !
        allocate(rbuf_real_r8(2*Nnode_R))
        CALL MPI_RECV(rbuf_real_r8,2*Nnode_R,MPI_DOUBLE_PRECISION, RankFirstSend, 1001, ArrLocal%comm, status, ierr)
        allocate(eGrid % LON_fe(Nnode_R), eGrid % LAT_fe(Nnode_R))
        idx=0
        DO i=1,Nnode_R
          idx=idx+1
          eGrid % LON_fe(i)=rbuf_real_r8(idx)
          idx=idx+1
          eGrid % LAT_fe(i)=rbuf_real_r8(idx)
        END DO
        deallocate(rbuf_real_r8)
        !
        allocate(rbuf_int(3*Ntrig_R))
        allocate(eGrid % ListTrig(3,Ntrig_R))
        CALL MPI_RECV(rbuf_int,3*Ntrig_R,MPI_INTEGER, RankFirstSend, 1002, ArrLocal%comm, status, ierr)
        idx=0
        DO i=1,Ntrig_R
          DO j=1,3
            idx=idx+1
            eGrid % ListTrig(j,i)=rbuf_int(idx)
          END DO
        END DO
        deallocate(rbuf_int)
      ELSE
        allocate(rbuf_int(2))
        CALL MPI_RECV(rbuf_int,2,MPI_INTEGER, RankFirstSend, 2000, ArrLocal % comm, status, ierr)
        xSize_R=rbuf_int(1)
        ySize_R=rbuf_int(2)
        deallocate(rbuf_int)
        eGrid % eta_rho = xSize_R
        eGrid % xi_rho  = ySize_R
        !
        allocate(rbuf_real_r8(2*xSize_R*ySize_R))
        CALL MPI_RECV(rbuf_real_r8,2*xSize_R*ySize_R,MPI_DOUBLE_PRECISION, RankFirstSend, 2001, ArrLocal%comm, status, ierr)
        allocate(eGrid % LON_fd(xSize_R, ySize_R))
        allocate(eGrid % LAT_fd(xSize_R, ySize_R))
        idx=0
        DO i=1,xSize_R
          DO j=1,ySize_R
            idx=idx+1
            eGrid % LON_fd(i,j)=rbuf_real_r8(idx)
            idx=idx+1
            eGrid % LAT_fd(i,j)=rbuf_real_r8(idx)
          END DO
        END DO
        deallocate(rbuf_real_r8)
        !
        allocate(rbuf_int(xSize_R*ySize_R))
        CALL MPI_RECV(rbuf_int,xSize_R*ySize_R,MPI_INTEGER, RankFirstSend, 2002, ArrLocal%comm, status, ierr)
        allocate(eGrid % MSK_fd(xSize_R, ySize_R))
        idx=0
        DO i=1,xSize_R
          DO j=1,ySize_R
            idx=idx+1
            eGrid % MSK_fd(i,j)=rbuf_int(idx)
          END DO
        END DO
        deallocate(rbuf_int)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE M2M_send_node_partition(ArrLocal, TargetColor, MatrixBelongingS)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: TargetColor
      type(T_node_partition), intent(in) :: MatrixBelongingS
      !
      integer idx, i, j, ierr, NnodesTarget
      integer, allocatable :: rbuf_int(:)
      integer ierror, iNode, iProc
      integer Ntot, NnodesS
      Ntot = MatrixBelongingS % Nnode
      NnodesS = MatrixBelongingS % Nthread
      allocate(rbuf_int(2))
      rbuf_int(1)=Ntot
      rbuf_int(2)=NnodesS
      NnodesTarget=ArrLocal % NbCompNodes(TargetColor)
      DO iNode=1,NnodesTarget
        iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
        CALL MPI_SEND(rbuf_int,2,MPI_INTEGER, iProc, 945, ArrLocal % comm, ierr)
      END DO
      deallocate(rbuf_int)

      allocate(rbuf_int(Ntot*NnodesS))
      idx=0
      DO i=1,Ntot
        DO j=1,NnodesS
          idx=idx+1
          rbuf_int(idx)=MatrixBelongingS % TheMatrix(i,j)
        END DO
      END DO
      DO iNode=1,NnodesTarget
        iProc=GetGlobalRank(ArrLocal, TargetColor, iNode)
        CALL MPI_SEND(rbuf_int,Ntot*NnodesS,MPI_INTEGER, iProc, 946, ArrLocal % comm, ierr)
      END DO
      deallocate(rbuf_int)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE M2M_recv_node_partition(ArrLocal, eColorSend, MatrixBelongingR)
      implicit none
      type(T_localNodes), intent(in) :: ArrLocal
      integer, intent(in) :: eColorSend
      type(T_node_partition), intent(inout) :: MatrixBelongingR
      !
      integer Ntot, NnodesSend
      integer idx, i, j
      integer, allocatable :: rbuf_int(:)
      integer :: status(MPI_STATUS_SIZE)
      integer ierr, RankFirstSend
      RankFirstSend=ArrLocal % ListFirstRank(eColorSend)
      allocate(rbuf_int(2))
      CALL MPI_RECV(rbuf_int,2,MPI_INTEGER, RankFirstSend, 945, ArrLocal%comm, status, ierr)
      Ntot=rbuf_int(1)
      NnodesSend=rbuf_int(2)
      deallocate(rbuf_int)
      MatrixBelongingR % Nnode=Ntot
      MatrixBelongingR % Nthread=NnodesSend
      IF (NnodesSend.ne.ArrLocal % NbCompNodes(eColorSend)) THEN
        Print *, 'Inconsistency between different declarations'
        STOP
      ENDIF
      !
      allocate(rbuf_int(Ntot*NnodesSend))
      CALL MPI_RECV(rbuf_int,Ntot*NnodesSend,MPI_INTEGER, RankFirstSend, 946, ArrLocal%comm, status, ierr)
      allocate(MatrixBelongingR % TheMatrix(Ntot, NnodesSend))
      idx=0
      DO i=1,Ntot
        DO j=1,NnodesSend
          idx=idx+1
          MatrixBelongingR % TheMatrix(i,j)=rbuf_int(idx)
        END DO
      END DO
      deallocate(rbuf_int)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE pgmcl_library
