#include "pgmcl_functions.h"
MODULE pgmcl_interp
      USE pgmcl_library
      implicit none
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRID_INTERP_GetKey(aKey, bKey, eGrid)
      implicit none
      integer, intent(out) :: aKey, bKey
      type(T_grid), intent(in) :: eGrid
      IF (eGrid % IsFE .eq. 1) THEN
        aKey=eGrid % nbNode
        bKey=eGrid % nbTrig
      ELSE
        aKey=eGrid % eta_rho
        bKey=eGrid % xi_rho
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRID_PRINT_KEY_INFO(eId, eGrid)
      implicit none
      integer, intent(in) :: eId
      type(T_grid), intent(in) :: eGrid
      real*8 :: MinLon, MaxLon
      real*8 :: MinLat, MaxLat
      IF (eGrid % IsFE .eq. 1) THEN
        WRITE(eId,*) 'unstructured mesh'
        MinLon=minval(eGrid % LON_fe)
        MaxLon=maxval(eGrid % LON_fe)
        MinLat=minval(eGrid % LAT_fe)
        MaxLat=maxval(eGrid % LAT_fe)
        WRITE(eId,*) 'nbNode=', eGrid % nbNode
        WRITE(eId,*) 'nbTrig=', eGrid % nbTrig
      ELSE
        WRITE(eId,*) 'finite difference mesh'
        MinLon=minval(eGrid % LON_fd)
        MaxLon=maxval(eGrid % LON_fd)
        MinLat=minval(eGrid % LAT_fd)
        MaxLat=maxval(eGrid % LAT_fd)
        WRITE(eId,*) 'eta_rho=', eGrid % eta_rho
        WRITE(eId,*) ' xi_rho=', eGrid % xi_rho
      END IF
      WRITE(eId, *) 'Lon(min/max)=', MinLon, MaxLon
      WRITE(eId, *) 'Lat(min/max)=', MinLat, MaxLat
      FLUSH(eId)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SAVE_CreateInterpolationSparseMatrix_Parall(           &
     &    FileSave, sMat, DoNearest,                                    &
     &    eGrid_input, eGrid_output,                                    &
     &    comm, MatrixBelongOutput)
      implicit none
      character(len=*), intent(in) :: FileSave
      type(T_sparse), intent(inout) :: sMat
      logical, intent(in) :: DoNearest
      type(T_grid), intent(in) :: eGrid_input, eGrid_output
      integer, intent(in) :: comm
      type(T_node_partition), intent(in) :: MatrixBelongOutput
      !
      integer myrank, ierr
      character(len=256) :: FileSave1
      logical file_exists, NeedComputation
      integer nbKey, MyNbKey, i
      integer, allocatable :: MyListKey(:), ListKey(:)
      integer aKey_input, bKey_input
      integer aKey_output, bKey_output
      MyNbKey=4
      allocate(MyListKey(MyNbKey))
      CALL GRID_INTERP_GetKey(aKey_input, bKey_input, eGrid_input)
      CALL GRID_INTERP_GetKey(aKey_output, bKey_output, eGrid_output)
      MyListKey(1)=aKey_input
      MyListKey(2)=bKey_input
      MyListKey(3)=aKey_output
      MyListKey(4)=bKey_output
      !
      WRITE (FileSave1,10) TRIM(FileSave)
  10  FORMAT (a,'.nc')
      !
      INQUIRE(FILE=FileSave1, EXIST=file_exists)
      NeedComputation=.FALSE.
      IF (file_exists .eqv. .FALSE.) THEN
        NeedComputation=.TRUE.
      ELSE
        CALL NETCDF_ReadSparseMatrix(sMat, nbKey, ListKey, FileSave1)
        IF (nbKey.ne.MyNbKey) THEN
          NeedComputation=.TRUE.
        ELSE
          DO i=1,nbKey
            IF (ListKey(i).ne.MyListKey(i)) THEN
              NeedComputation=.TRUE.
            END IF
          ENDDO
        ENDIF
        deallocate(ListKey)
      ENDIF
      IF (NeedComputation) THEN
        CALL CreateInterpolationSparseMatrix_Parall(                    &
     &    sMat, DoNEarest, eGrid_input, eGrid_output,                   &
     &    comm, MatrixBelongOutput)
        CALL mpi_comm_rank(comm, myrank, ierr)
        IF (myrank .eq. 0) THEN
          CALL NETCDF_WriteSparseMatrix(sMat, MyNbKey, MyListKey, FileSave1)
          Print *, 'R8_FD_R8_FD Write interpol. to ', TRIM(FileSave1)
        END IF
      ELSE
        Print *, 'Using interpolation info from ', TRIM(FileSave1)
      ENDIF
      deallocate(MyListKey)
      END SUBROUTINE
!**********************************************************************
!* The routine 
!* Each node compute its block of the interpolation matrix            *
!**********************************************************************
      SUBROUTINE PRINT_TO_FILE_SPARSEMATRIX(id,                          &
     &    eGrid_input, eGrid_output, sMat)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: id
      type(T_sparse), intent(in) :: sMat
      type(T_grid), intent(in) :: eGrid_input
      type(T_grid), intent(in) :: eGrid_output
      integer ielt
      WRITE(id,*) 'numelt=', sMat % numelt
      WRITE(id,*) 'nRows=', sMat % nRows
      WRITE(id,*) 'nColumns=', sMat % nColumns
      DO ielt=1,sMat % numelt
        WRITE(id,*) 'ielt=', ielt
        WRITE(id,*) '  c/r/w=', sMat % columns(ielt), sMat % rows(ielt), sMat % weights(ielt)
      END DO
      FLUSH(id)
      END SUBROUTINE
!**********************************************************************
!* The routine below is fully parallelized                            *
!* Each node compute its block of the interpolation matrix            *
!**********************************************************************
      SUBROUTINE CreateInterpolationSparseMatrix_Parall(                 &
     &    sMat, DoNearest, eGrid_input, eGrid_output,                    &
     &    comm, MatrixBelongOutput)
      implicit none
      type(T_sparse), intent(inout) :: sMat
      logical, intent(in) :: DoNearest
      type(T_grid), intent(in) :: eGrid_input
      type(T_grid), intent(in) :: eGrid_output
      integer, intent(in) :: comm
      type(T_node_partition), intent(in) :: MatrixBelongOutput
      !
      integer nbPoint, nbPointB, iPoint
      real*8, allocatable :: LLon(:), LLat(:)
      integer, allocatable :: ListIdxOut(:)
      integer iout, jout, nx_output, ny_output
      integer eMsk, ielt, idxout, numelt
      real*8 eLon, eLat
      integer myrank, nbproc, nbNode, nbproc_total
      integer status(MPI_STATUS_SIZE)
      integer nbPointLoc, numeltLoc
      integer i, eTag, iNode, iProc, idx
      integer eRow, fRow, eCol
      integer Totalnumelt, eVal
      integer, allocatable :: rbuf_int(:)
      integer, allocatable :: ListNbPoint(:), Listnumelt(:)
      integer, allocatable :: IJ_to_I(:), IJ_to_J(:)
      integer, allocatable :: ListIdxOutLoc(:)
      real*8, allocatable :: ListWeightLoc(:), rbuf_real(:)
      real*8 :: eWei
      integer :: NbPointProc(MatrixBelongOutput % Nthread)
      integer :: AssignWorkDone(MatrixBelongOutput % Nnode)
      type(T_sparse) :: sMatPartial
      integer istatus, ierr
      CALL mpi_comm_rank(comm, myrank, ierr)
      CALL mpi_comm_size(comm, nbproc_total, ierr)
      nbproc=MatrixBelongOutput % Nthread
      nbNode=MatrixBelongOutput % Nnode
      IF (myrank .ge. nbproc) THEN
        Print *, 'This processor has nothing to do here'
        Print *, 'Please leave'
        STOP
      END IF
      IF (nbproc_total .lt. nbproc) THEN
        Print *, 'Error in nbproc_total. Please debug'
        STOP
      END IF
      AssignWorkDone=0
      DO iProc=1,nbproc
        DO iNode=1,nbNode
          eVal=MatrixBelongOutput % TheMatrix(iNode,iProc)
          IF (eVal .gt. 0) THEN
            IF (AssignWorkDone(iNode) .eq. 0) THEN
              AssignWorkDone(iNode)=iProc
            END IF
          END IF
        END DO
      END DO
      NbPointProc=0
      DO iNode=1,nbNode
        iProc=AssignWorkDone(iNode)
        IF (iProc .gt. 0) THEN
          NbPointProc(iProc)=NbPointProc(iProc) + 1
        END IF
      END DO
      nbPoint=NbPointProc(myrank+1)
      ALLOCATE(LLon(nbPoint), LLat(nbPoint), ListIdxOut(nbPoint))
      IF (eGrid_output % IsFE .eq. 0) THEN
        nx_output=eGrid_output % eta_rho
        ny_output=eGrid_output % xi_rho
        allocate(IJ_to_I(nx_output*ny_output), IJ_to_J(nx_output*ny_output))
        DO iout=1,nx_output
          DO jout=1,ny_output
            idxout=iout+nx_output*(jout-1)
            IJ_to_I(idxout)=iout
            IJ_to_J(idxout)=jout
          END DO
        END DO
        iPoint=0
        DO iNode=1,nbNode
          IF (AssignWorkDone(iNode) .eq. myrank+1) THEN
            iPoint=iPoint+1
            iout=IJ_to_I(iNode)
            jout=IJ_to_J(iNode)
            LLon(iPoint)=eGrid_output % LON_fd(iout,jout)
            LLat(iPoint)=eGrid_output % LAT_fd(iout,jout)
            ListIdxOut(iPoint)=iNode
          END IF
        END DO
        deallocate(IJ_to_I, IJ_to_J)
      ELSE
        iPoint=0
        DO iNode=1,nbNode
          IF (AssignWorkDone(iNode) .eq. myrank+1) THEN
            iPoint=iPoint+1
            LLon(iPoint)=eGrid_output % LON_fe(iNode)
            LLat(iPoint)=eGrid_output % LAT_fe(iNode)
            ListIdxOut(iPoint)=iNode
          END IF
        END DO
      END IF
      IF (eGrid_input % IsFE .eq. 0) THEN
        CALL CreateInterpolationSparseMatrix_r8_FD_2_Points(             &
     &    sMatPartial, DoNearest, eGrid_input, nbPoint, LLon, LLat)
      ELSE
        CALL CreateInterpolationSparseMatrix_r8_FE_2_Points(             &
     &    sMatPartial, DoNearest, eGrid_input, nbPoint, LLon, LLat)
      END IF
      deallocate(LLon, LLat)
      eTag=240
      IF (myrank .eq. 0) THEN
        allocate(ListNbPoint(nbproc), Listnumelt(nbproc))
        ListNbPoint(1)=nbPoint
        Listnumelt(1)=sMatPartial % numelt
        Totalnumelt=sMatPartial % numelt
        allocate(rbuf_int(2))
        DO iProc=2,nbproc
          CALL MPI_RECV(rbuf_int,2,MPI_INTEGER, iProc-1, 237, comm, status, ierr)
          nbPointLoc=rbuf_int(1)
          numeltLoc=rbuf_int(2)
          ListNbPoint(iProc)=nbPointLoc
          Listnumelt(iProc)=numeltLoc
          Totalnumelt=Totalnumelt + numeltLoc
        END DO
        deallocate(rbuf_int)
        sMat % numelt=Totalnumelt
        sMat % nRows=nbNode
        sMat % nColumns=sMatPartial % nColumns
        allocate(sMat % rows(Totalnumelt))
        allocate(sMat % columns(Totalnumelt))
        allocate(sMat % weights(Totalnumelt))
        i=0
        DO ielt=1,sMatPartial % numelt
          i=i+1
          eRow=sMatPartial % rows(ielt)
          fRow=ListIdxOut(eRow)
          sMat % rows(i)=fRow
          sMat % columns(i)=sMatPartial % columns(ielt)
          sMat % weights(i)=sMatPartial % weights(ielt)
        END DO
        DO iProc=2,nbproc
          nbPointLoc=ListNbPoint(iProc)
          numeltLoc=Listnumelt(iProc)
          nbPointB=NbPointProc(iProc)
          IF (nbPointB .ne. ListNbPoint(iProc)) THEN
            Print *, 'Inconsistency here between nbPointB and ListNbPoint(iProc)'
            STOP
          END IF
          allocate(ListIdxOutLoc(nbPointB), rbuf_int(2*numeltLoc), ListWeightLoc(numeltLoc))
          iPoint=0
          DO iNode=1,nbNode
            IF (AssignWorkDone(iNode) .eq. iProc) THEN
              iPoint=iPoint+1
              ListIdxOutLoc(iPoint)=iNode
            END IF
          END DO
          CALL MPI_RECV(rbuf_int,2*numeltLoc,MPI_INTEGER, iProc-1, 239, comm, status, ierr)
          CALL MPI_RECV(ListWeightLoc,numeltLoc,MPI_DOUBLE_PRECISION, iProc-1, 241, comm, status, ierr)
          DO ielt=1,numeltLoc
            eRow=rbuf_int(ielt)
            fRow=ListIdxOutLoc(eRow)
            eCol=rbuf_int(ielt + numeltLoc)
            eWei=ListWeightLoc(ielt)
            i=i+1
            sMat % rows(i)=fRow
            sMat % columns(i)=eCol
            sMat % weights(i)=eWei
          END DO
          deallocate(rbuf_int, ListWeightLoc, ListIdxOutLoc)
        END DO
        deallocate(ListNbPoint, Listnumelt)
        !
        DO iProc=2,nbproc
          CALL MPI_SEND_SPARSE_MATRIX(sMat, comm, iProc-1, eTag)
        END DO
      ELSE
        numelt=sMatPartial % numelt
        allocate(rbuf_int(2))
        rbuf_int(1)=nbPoint
        rbuf_int(2)=numelt
        CALL MPI_SEND(rbuf_int,2,MPI_INTEGER, 0, 237, comm, istatus)
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(2*numelt))
        idx=0
        DO i=1,numelt
          idx=idx+1
          rbuf_int(idx)=sMatPartial % rows(i)
        END DO
        DO i=1,numelt
          idx=idx+1
          rbuf_int(idx)=sMatPartial % columns(i)
        END DO
        CALL MPI_SEND(rbuf_int,2*numelt,MPI_INTEGER, 0, 239, comm, istatus)
        deallocate(rbuf_int)
        allocate(rbuf_real(numelt))
        DO i=1,numelt
          rbuf_real(i)=sMatPartial % weights(i)
        END DO
        CALL MPI_SEND(rbuf_real,numelt,MPI_DOUBLE_PRECISION, 0, 241, comm, istatus)
        deallocate(rbuf_real)
        CALL MPI_RECV_SPARSE_MATRIX(sMat, comm, 0, eTag)
      END IF
      CALL DeallocSparseMatrix(sMatPartial)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CreateInterpolationSparseMatrix_r8_FD_2_Points(         &
     &    sMat, DoNearest, eGrid_input, nbPoint, LLon, LLat)
      implicit none
      type(T_sparse), intent(inout) :: sMat
      logical, intent(in) :: DoNearest
      type(T_grid), intent(in) :: eGrid_input
      integer, intent(in) :: nbPoint
      real*8, intent(in) :: LLon(nbPoint), LLat(nbPoint)
      !
      integer upperest_numelements, numelements
      integer, dimension(:), pointer :: est_rows, est_columns
      integer, dimension(:), pointer :: ListIIN, ListJIN
      real*8, dimension(:), pointer :: ListDet, ListSign
      real*8 eQuot1, eQuot2, eQuot3
      real*8 eDet1, eDet2, eDet3, eDet, eSign
      real*8 eProd1, eProd2, eProd3, eDiff
      real*8, dimension(:,:), pointer :: ListTrigLON, ListTrigLAT
      real*8, dimension(:), pointer :: est_weights
      integer, dimension(4) :: ListI, ListJ
      real, dimension(4) :: ListV
      integer NumberDirect, iTrig, nbTrig_input, eRes
      real*8 eLon1, eLon2, eLon3, eLat1, eLat2, eLat3
      real*8 eLon11, eLat11, eLon12, eLat12, eLon21, eLat21, eLon22, eLat22
      integer idx
      integer iin, jin, iout, jout, idxin, idxmat
      real*8 eLonO, eLatO, eLonI, eLatI
      real*8 MinDist, eDist, eDiffDet
      real*8 MinLonInput, MaxLonInput, MinLatInput, MaxLatInput
      real*8 MinLonOutput, MaxLonOutput, MinLatOutput, MaxLatOutput
      logical WeFind
      integer is, iinSel, jinSel, i, nbWet
      integer nbWetOut, nbWetIn, nbNearest, nbFind
      real*8 eWeight, SumWeight, eCoeffLon, eCoeffLat
      integer MaxRowsIdx, MaxColumnsIdx
      integer nbContained, iPoint
      integer nRows, nColumns, iNode
      logical DoVerbose
      real*8 tolWeight
      real*8 eDist1, eDist2, eDist3, eCritDist
      real*8 eLon, eLat
      integer ix, iy, MyError, eRank
      character(len=40) :: FileSave
      character(len=3) :: eStr
      integer iTrigSelect
      real MinAbsValue, eValue
      integer nx_input, ny_input
      integer eMsk, eMsk1, eMsk2, eMsk3, eMsk4
      nx_input = eGrid_input % eta_rho
      ny_input = eGrid_input % xi_rho
      !
      DoVerbose=.TRUE.
      tolWeight=0.2
      upperest_numelements=4*nbPoint
      allocate(est_rows(upperest_numelements))
      allocate(est_columns(upperest_numelements))
      allocate(est_weights(upperest_numelements))
      idxmat=0
      nbNearest=0
      nbFind=0
      nbWetIn=0
      MinLonInput=909
      MaxLonInput=-900
      MinLatInput=900
      MaxLatInput=-900
      DO jin=1,ny_input
        DO iin=1,nx_input
          eMsk=eGrid_input % MSK_fd(iin,jin)
          IF (eMsk == 1) THEN
            nbWetIn=nbWetIn+1
            eLon=eGrid_input % LON_fd(iin,jin)
            eLat=eGrid_input % LAT_fd(iin,jin)
            IF (MinLonInput .gt. eLon) THEN
               MinLonInput=eLon
            END IF
            IF (MaxLonInput.lt.eLon) THEN
               MaxLonInput=eLon
            END IF
            IF (MinLatInput .gt. eLat) THEN
               MinLatInput=eLat
            END IF
            IF (MaxLatInput .lt. eLat) THEN
               MaxLatInput=eLat
            END IF
          END IF
        END DO
      END DO
      MinLonOutput=minval(LLon)
      MaxLonOutput=maxval(LLon)
      MinLatOutput=minval(LLat)
      MaxLatOutput=maxval(LLat)
      nbTrig_input=0
      DO jin=2,ny_input
        DO iin=2,nx_input
          eMsk1=eGrid_input % MSK_fd(iin-1,jin-1)
          eMsk2=eGrid_input % MSK_fd(iin  ,jin-1)
          eMsk3=eGrid_input % MSK_fd(iin-1,jin  )
          eMsk4=eGrid_input % MSK_fd(iin  ,jin)
          nbWet=eMsk1 + eMsk2 + eMsk3 + eMsk4
          IF (nbWet.gt.0) THEN
            nbTrig_input=nbTrig_input+2
          END IF
        END DO
      END DO
      ALLOCATE(ListIIN(2*nbTrig_input))
      ALLOCATE(ListJIN(2*nbTrig_input))
      ALLOCATE(ListTrigLON(3,2*nbTrig_input))
      ALLOCATE(ListTrigLAT(3,2*nbTrig_input))
      idx=0
      DO iin=2,nx_input
        DO jin=2,ny_input
          eMsk1=eGrid_input % MSK_fd(iin-1,jin-1)
          eMsk2=eGrid_input % MSK_fd(iin  ,jin-1)
          eMsk3=eGrid_input % MSK_fd(iin-1,jin  )
          eMsk4=eGrid_input % MSK_fd(iin  ,jin)
          nbWet=eMsk1 + eMsk2 + eMsk3 + eMsk4
          IF (nbWet.gt.0) THEN
            eLon11=eGrid_input % LON_fd(iin-1,jin-1)
            eLat11=eGrid_input % LAT_fd(iin-1,jin-1)
            eLon12=eGrid_input % LON_fd(iin-1,jin)
            eLat12=eGrid_input % LAT_fd(iin-1,jin)
            eLon21=eGrid_input % LON_fd(iin,jin-1)
            eLat21=eGrid_input % LAT_fd(iin,jin-1)
            eLon22=eGrid_input % LON_fd(iin,jin)
            eLat22=eGrid_input % LAT_fd(iin,jin)
            idx=idx+1
            ListIIN(idx)=iin
            ListJIN(idx)=jin
            ListTrigLON(1,idx)=eLon11
            ListTrigLON(2,idx)=eLon12
            ListTrigLON(3,idx)=eLon21
            ListTrigLAT(1,idx)=eLat11
            ListTrigLAT(2,idx)=eLat12
            ListTrigLAT(3,idx)=eLat21
            idx=idx+1
            ListIIN(idx)=iin
            ListJIN(idx)=jin
            ListTrigLON(1,idx)=eLon22
            ListTrigLON(2,idx)=eLon12
            ListTrigLON(3,idx)=eLon21
            ListTrigLAT(1,idx)=eLat22
            ListTrigLAT(2,idx)=eLat12
            ListTrigLAT(3,idx)=eLat21
          END IF
        END DO
      END DO
      ALLOCATE(ListDet(nbTrig_input), ListSign(nbTrig_input))
      DO iTrig=1,nbTrig_input
        eLon1=ListTrigLON(1,iTrig)
        eLon2=ListTrigLON(2,iTrig)
        eLon3=ListTrigLON(3,iTrig)
        eLat1=ListTrigLAT(1,iTrig)
        eLat2=ListTrigLAT(2,iTrig)
        eLat3=ListTrigLAT(3,iTrig)
        eDet=(eLon3-eLon1)*(eLat2-eLat1) - (eLat3-eLat1)*(eLon2-eLon1)
        ListDet(iTrig)=eDet
        IF (eDet.gt.0) THEN
          eSign=1
        ELSE
          eSign=-1
        END IF
        ListSign(iTrig)=eSign
      ENDDO
      NumberDirect=0
      eCritDist=0.01
      DO iPoint=1,nbPoint
        eLonO=LLon(iPoint)
        eLatO=LLat(iPoint)
        WeFind=.FALSE.
        nbContained=0
        DO iTrig=1,nbTrig_input
          IF (WeFind.eqv..FALSE.) THEN
            eLon1=ListTrigLON(1,iTrig)
            eLon2=ListTrigLON(2,iTrig)
            eLon3=ListTrigLON(3,iTrig)
            eLat1=ListTrigLAT(1,iTrig)
            eLat2=ListTrigLAT(2,iTrig)
            eLat3=ListTrigLAT(3,iTrig)
            eDet1=(eLonO-eLon2)*(eLat3-eLat2) - (eLatO-eLat2)*(eLon3-eLon2)
            eDet2=(eLonO-eLon3)*(eLat1-eLat3) - (eLatO-eLat3)*(eLon1-eLon3)
            eDet3=(eLonO-eLon1)*(eLat2-eLat1) - (eLatO-eLat1)*(eLon2-eLon1)
            eDet=ListDet(iTrig)
            eSign=ListSign(iTrig)
            eProd1=eDet1*eSign
            eProd2=eDet2*eSign
            eProd3=eDet3*eSign
            IF ((eProd1.ge.0).and.(eProd2.ge.0).and.(eProd3.ge.0)) THEN
              nbContained=nbContained+1
              WeFind=.TRUE.
              eDist1=SQRT((eLonO - eLon1)**2 + (eLatO - eLat1)**2)
              eDist2=SQRT((eLonO - eLon2)**2 + (eLatO - eLat2)**2)
              eDist3=SQRT((eLonO - eLon3)**2 + (eLatO - eLat3)**2)
              eQuot1=eDet1/eDet
              eQuot2=eDet2/eDet
              eQuot3=eDet3/eDet
              eRes=MOD(iTrig,2)
              ListV(1)=eQuot1
              ListV(2)=eQuot2
              ListV(3)=eQuot3
              iin=ListIIN(iTrig)
              jin=ListJIN(iTrig)
              IF (eRes.eq.1) THEN
                ListI(1)=iin-1
                ListI(2)=iin-1
                ListI(3)=iin
                ListJ(1)=jin-1
                ListJ(2)=jin
                ListJ(3)=jin-1
              ELSE
                ListI(1)=iin
                ListI(2)=iin-1
                ListI(3)=iin
                ListJ(1)=jin
                ListJ(2)=jin
                ListJ(3)=jin-1
              ENDIF
            ENDIF
          ENDIF
        END DO
        nbWet=0
        SumWeight=0
        IF (WeFind) THEN
          nbFind=nbFind+1
          DO is=1,3
            iin=ListI(is)
            jin=ListJ(is)
            eWeight=ListV(is)
            eMsk=eGrid_input % MSK_fd(iin,jin)
            IF (eMsk .eq. 1) THEN
              nbWet=nbWet+1
              SumWeight=SumWeight+eWeight
            END IF
          END DO
        END IF
        IF (WeFind.and.(nbWet.gt.0).and.(SumWeight.ge.tolWeight)) THEN
          NumberDirect=NumberDirect+1
          DO is=1,3
            iin=ListI(is)
            jin=ListJ(is)
            eWeight=ListV(is)
            eMsk=eGrid_input % MSK_fd(iin,jin)
            IF (eMsk .eq. 1) THEN
              idxin=iin+nx_input*(jin-1)
              idxmat=idxmat+1
              est_rows(idxmat)=iPoint
              est_columns(idxmat)=idxin
              est_weights(idxmat)=eWeight/SumWeight
            END IF
          END DO
        ELSE
          IF (DoNearest) THEN
            nbNearest=nbNearest+1
            MinDist=1000000000
            WeFind=.FALSE.
            DO iin=1,nx_input
              DO jin=1,ny_input
                eMsk=eGrid_input % MSK_fd(iin,jin)
                IF (eMsk .eq. 1) THEN
                  eLonI=eGrid_input % LON_fd(iin,jin)
                  eLatI=eGrid_input % LAT_fd(iin,jin)
                  eDist=sqrt((eLonO-eLonI)**2+(eLatO-eLatI)**2)
                  IF (MinDist.gt.eDist) THEN
                    iinSel=iin
                    jinSel=jin
                    MinDist=eDist
                    WeFind=.TRUE.
                  END IF
                END IF
              END DO
            END DO
            IF (WeFind.eqv..FALSE.) THEN
              Print *, 'We have a problem since we did not find'
              Print *, 'Leaving CreateInterpolationSparseMatrix_r8_FD_2_Points'
              Print *, 'any point'
              Print *, 'MinLonInput=', MinLonInput
              Print *, 'MaxLonInput=', MaxLonInput
              Print *, 'MinLatInput=', MinLatInput
              Print *, 'MaxLatInput=', MaxLatInput
              Print *, 'MinLonOutput=', MinLonOutput
              Print *, 'MaxLonOutput=', MaxLonOutput
              Print *, 'MinLatOutput=', MinLatOutput
              Print *, 'MaxLatOutput=', MaxLatOutput
              STOP
            END IF
            idxin=iinSel+nx_input*(jinSel-1)
            idxmat=idxmat+1
            est_rows(idxmat)=iPoint
            est_columns(idxmat)=idxin
            est_weights(idxmat)=1
          END IF
        END IF
      END DO
      numelements=idxmat
      nRows=nbPoint
      nColumns=nx_input*ny_input
      sMat % nRows=nRows
      sMat % nColumns=nColumns
      sMat % numelt=numelements
      allocate(sMat % rows(numelements))
      allocate(sMat % columns(numelements))
      allocate(sMat % weights(numelements))
      MaxRowsIdx=0
      MaxColumnsIdx=0
      DO i=1,numelements
         sMat % rows(i)=est_rows(i)
         sMat % columns(i)=est_columns(i)
         sMat % weights(i)=est_weights(i)
         IF (MaxRowsIdx.le.est_rows(i)) THEN
            MaxRowsIdx=est_rows(i)
         END IF
         IF (MaxColumnsIdx.le.est_columns(i)) THEN
            MaxColumnsIdx=est_columns(i)
         END IF
      END DO
      deallocate(est_rows, est_columns, est_weights)
      deallocate(ListIIN, ListJIN, ListDet, ListSign)
      deallocate(ListTrigLON, ListTrigLAT)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CreateInterpolationSparseMatrix_r8_FE_2_Points(        &
     &    sMat, DoNearest,                                              &
     &    eGrid_input, nbPoint, LLon, LLat)
      implicit none
      type(T_sparse), intent(inout) :: sMat
      logical, intent(in) :: DoNearest
      type(T_grid), intent(in) :: eGrid_input
      integer, intent(in) :: nbPoint
      real*8, intent(in) :: LLon(nbPoint), LLat(nbPoint)
      !
      integer upperest_numelements, numelements
      integer, dimension(:), pointer :: est_rows, est_columns
      real*8, dimension(:), pointer :: ListDet, ListSign
      real*8 eQuot1, eQuot2, eQuot3
      real*8 eDet1, eDet2, eDet3, eDet, eSign
      real*8 eProd1, eProd2, eProd3, eDiff
      real*8, dimension(:), pointer :: est_weights
      integer, dimension(4) :: ListI, ListJ
      real*8, dimension(4) :: ListV
      integer NumberDirect, iTrig, nbTrig, eRes
      real*8 eLon1, eLon2, eLon3, eLat1, eLat2, eLat3
      real*8 eLon11, eLat11, eLon12, eLat12, eLon21, eLat21
      real*8 eLon22, eLat22
      integer idx
      integer iin, jin, iout, jout, idxout, idxin, idxmat
      real*8 eLonO, eLatO, eLonI, eLatI
      real*8 MinDist, eDist
      real*8 MinLonInput, MaxLonInput, MinLatInput, MaxLatInput
      real*8 MinLonOutput, MaxLonOutput, MinLatOutput, MaxLatOutput
      logical WeFind
      integer is, iinSel, jinSel, i, nbWet
      integer nbNearest, nbFind
      real*8 eWeight, SumWeight, eCoeffLon, eCoeffLat
      integer MaxRowsIdx, MaxColumnsIdx
      integer nRows, nColumns, iNode, iPoint
      integer iNodeSel, iNode1, iNode2, iNode3
      logical DoVerbose
      real*8 tolWeight, eLon, eLat
      integer nbNode_input, nbTrig_input
      integer eMsk
      nbNode_input=eGrid_input % nbNode
      nbTrig_input=eGrid_input % nbTrig
      !
      DoVerbose=.TRUE.
      tolWeight=0.2
      upperest_numelements=4*nbPoint
      allocate(est_rows(upperest_numelements))
      allocate(est_columns(upperest_numelements))
      allocate(est_weights(upperest_numelements))
      idxmat=0
      nbNearest=0
      nbFind=0
      MinLonInput=minval(eGrid_input % LON_fe)
      MaxLonInput=maxval(eGrid_input % LON_fe)
      MinLatInput=minval(eGrid_input % LAT_fe)
      MaxLatInput=maxval(eGrid_input % LAT_fe)
      MinLonOutput=minval(LLon)
      MaxLonOutput=maxval(LLon)
      MinLatOutput=minval(LLat)
      MaxLatOutput=maxval(LLat)
      ALLOCATE(ListDet(nbTrig_input), ListSign(nbTrig_input))
      DO iTrig=1,nbTrig_input
        iNode1=eGrid_input % ListTrig(1,iTrig)
        iNode2=eGrid_input % ListTrig(2,iTrig)
        iNode3=eGrid_input % ListTrig(3,iTrig)
        eLon1=eGrid_input % LON_fe(iNode1)
        eLon2=eGrid_input % LON_fe(iNode2)
        eLon3=eGrid_input % LON_fe(iNode3)
        eLat1=eGrid_input % LAT_fe(iNode1)
        eLat2=eGrid_input % LAT_fe(iNode2)
        eLat3=eGrid_input % LAT_fe(iNode3)
        eDet=(eLon3-eLon1)*(eLat2-eLat1) - (eLat3-eLat1)*(eLon2-eLon1)
        ListDet(iTrig)=eDet
        IF (eDet.gt.0) THEN
          eSign=1
        ELSE
          eSign=-1
        END IF
        ListSign(iTrig)=eSign
      ENDDO
      NumberDirect=0
      DO iPoint=1,nbPoint
        eLonO=LLon(iPoint)
        eLatO=LLat(iPoint)
        WeFind=.FALSE.
        DO iTrig=1,nbTrig_input
          IF (WeFind.eqv..FALSE.) THEN
            iNode1=eGrid_input % ListTrig(1,iTrig)
            iNode2=eGrid_input % ListTrig(2,iTrig)
            iNode3=eGrid_input % ListTrig(3,iTrig)
            eLon1=eGrid_input % LON_fe(iNode1)
            eLon2=eGrid_input % LON_fe(iNode2)
            eLon3=eGrid_input % LON_fe(iNode3)
            eLat1=eGrid_input % LAT_fe(iNode1)
            eLat2=eGrid_input % LAT_fe(iNode2)
            eLat3=eGrid_input % LAT_fe(iNode3)
            eDet3=(eLonO-eLon1)*(eLat2-eLat1) -                         &
     &            (eLatO-eLat1)*(eLon2-eLon1)
            eDet1=(eLonO-eLon2)*(eLat3-eLat2) -                         &
     &            (eLatO-eLat2)*(eLon3-eLon2)
            eDet2=(eLonO-eLon3)*(eLat1-eLat3) -                         &
     &            (eLatO-eLat3)*(eLon1-eLon3)
            eDet=ListDet(iTrig)
            eSign=ListSign(iTrig)
            eProd1=eDet1*eSign
            eProd2=eDet2*eSign
            eProd3=eDet3*eSign
            IF ((eProd1.ge.0).and.(eProd2.ge.0).and.(eProd3.ge.0)) THEN
              WeFind=.TRUE.
              eQuot1=eDet1/eDet
              eQuot2=eDet2/eDet
              eQuot3=eDet3/eDet
              eRes=MOD(iTrig,2)
              ListV(1)=eQuot1
              ListV(2)=eQuot2
              ListV(3)=eQuot3
              ListI(1)=iNode1
              ListI(2)=iNode2
              ListI(3)=iNode3
            ENDIF
          ENDIF
        END DO
        IF (WeFind) THEN
          NumberDirect=NumberDirect+1
          DO is=1,3
            iNode=ListI(is)
            eWeight=ListV(is)
            idxmat=idxmat+1
            est_rows(idxmat)=iPoint
            est_columns(idxmat)=iNode
            est_weights(idxmat)=eWeight
          END DO
        ELSE
          IF (DoNearest) THEN
            nbNearest=nbNearest+1
            MinDist=1000000000
            WeFind=.FALSE.
            DO iNode=1,nbNode_input
              eLonI=eGrid_input % LON_fe(iNode)
              eLatI=eGrid_input % LAT_fe(iNode)
              eDist=sqrt((eLonO-eLonI)**2+(eLatO-eLatI)**2)
              IF (MinDist.gt.eDist) THEN
                iNodeSel=iNode
                MinDist=eDist
                WeFind=.TRUE.
              END IF
            END DO
            IF (WeFind.eqv..FALSE.) THEN
              Print *, 'We have a problem since we did not find'
              Print *, 'Leaving CreateInterpolationSparseMatrix_r8_FE_2_Points'
              Print *, 'any point'
              Print *, 'MinLonInput=', MinLonInput
              Print *, 'MaxLonInput=', MaxLonInput
              Print *, 'MinLatInput=', MinLatInput
              Print *, 'MaxLatInput=', MaxLatInput
              Print *, 'MinLonOutput=', MinLonOutput
              Print *, 'MaxLonOutput=', MaxLonOutput
              Print *, 'MinLatOutput=', MinLatOutput
              Print *, 'MaxLatOutput=', MaxLatOutput
              STOP
            END IF
            idxmat=idxmat+1
            est_rows(idxmat)=iPoint
            est_columns(idxmat)=iNodeSel
            est_weights(idxmat)=1
          END IF
        END IF
      END DO
      numelements=idxmat
      nRows=nbPoint
      nColumns=nbNode_input
      sMat % nRows=nRows
      sMat % nColumns=nColumns
      sMat % numelt=numelements
      allocate(sMat % rows(numelements))
      allocate(sMat % columns(numelements))
      allocate(sMat % weights(numelements))
      MaxRowsIdx=0
      MaxColumnsIdx=0
      DO i=1,numelements
         sMat % rows(i)=est_rows(i)
         sMat % columns(i)=est_columns(i)
         sMat % weights(i)=est_weights(i)
         IF (MaxRowsIdx.le.est_rows(i)) THEN
            MaxRowsIdx=est_rows(i)
         END IF
         IF (MaxColumnsIdx.le.est_columns(i)) THEN
            MaxColumnsIdx=est_columns(i)
         END IF
      END DO
      deallocate(est_rows, est_columns, est_weights)
      deallocate(ListDet, ListSign)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE Kernel_DetermineReachedLand(                           &
     &    MSK_attained, sMat, eGrid_output)
      implicit none
      integer, allocatable, intent(inout) :: MSK_attained(:,:)
      type(T_grid) :: eGrid_input, eGrid_output
      TYPE(T_sparse) :: sMat
      integer nx_output, ny_output, ielt, numelt
      integer res, eIdxOut, iout, jout
      nx_output=eGrid_output % eta_rho
      ny_output=eGrid_output % xi_rho
      allocate(MSK_attained(nx_output, ny_output))
      MSK_attained=0
      numelt=sMat % numelt
      DO ielt=1,numelt
        eIdxOut=sMat % rows(ielt)
        !idxout=iout+nx_output*(jout-1)
        res=MODULO(eIdxOut, nx_output)
        IF (res.eq.0) THEN
          iout=nx_output
          jout=eIdxOut/nx_output
        ELSE
          iout=res
          jout=1+(eIdxOut-res)/nx_output
        END IF
        MSK_attained(iout, jout)=1
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE pgmcl_interp
