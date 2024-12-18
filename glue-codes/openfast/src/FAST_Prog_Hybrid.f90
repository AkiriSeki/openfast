!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of FAST.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
PROGRAM FAST
! This program models 2- or 3-bladed turbines of a standard configuration.
!
! noted compilation switches:
!   OUTPUT_ADDEDMASS        (outputs a file called "<RootName>.AddedMass" that contains HydroDyn's added-mass matrix.
!   OUTPUT_JACOBIAN
!   FPE_TRAP_ENABLED        (use with gfortran when checking for floating point exceptions)
!   DOUBLE_PRECISION        (compile in double precision)
!.................................................................................................


USE FAST_Subs   ! all of the ModuleName and ModuleName_types modules are inherited from FAST_Subs
                       
IMPLICIT  NONE
   
! Local parameters:
REAL(DbKi),             PARAMETER     :: t_initial = 0.0_DbKi                    ! Initial time
INTEGER(IntKi),         PARAMETER     :: NumTurbines = 1
   
! Other/Misc variables
TYPE(FAST_TurbineType)                :: Turbine(NumTurbines)                    ! Data for each turbine instance

INTEGER(IntKi)                        :: i_turb                                  ! current turbine number
INTEGER(IntKi)                        :: n_t_global                              ! simulation time step, loop counter for global (FAST) simulation
INTEGER(IntKi)                        :: ErrStat                                 ! Error status
CHARACTER(1024)                       :: ErrMsg                                  ! Error message

! data for restart:
CHARACTER(1000)                       :: InputFile                               ! String to hold the intput file name
CHARACTER(1024)                       :: CheckpointRoot                          ! Rootname of the checkpoint file
CHARACTER(20)                         :: FlagArg                                 ! flag argument from command line
INTEGER(IntKi)                        :: Restart_step                            ! step to start on (for restart) 

! Varuables for socket communication @ AS
!INTEGER,                PARAMETER     :: NumInputs_c = 12                        ! Number of Inputs from OpenFresco (This may be removed in the future)
!INTEGER,                PARAMETER     :: NumFixedInputs = 12
INTEGER,                PARAMETER      :: dataSize = 256                        ! maximum data size sent/receive to/from OpenFresco
! Currently, assuming only one socket communication so belowe is commented out
!INTEGER,                PARAMETER      :: numSockIDs = 32                       ! maximum number of sockect communication)
INTEGER,                PARAMETER      :: sizeInt = 4                           ! data size of integer
INTEGER,                PARAMETER      :: sizeDouble = 8                        ! data size of double
INTEGER,                PARAMETER      :: port = 8090
INTEGER                                :: sizeMachineInet
! Currently, assuming only one socket communication so below is commented out
!INTEGER(numSocketIDs)                  :: socketIDs                             ! container of socketIDs
! In fortran 90, to my understanding, decalaring variable in PROGRAM makes the variable global variable. Fortran 77 does not, so it need common.
INTEGER                                :: socketID
! Below is fortran 77 format
!integer*4 socketID
! common    //socketID ! AS: check this
INTEGER                                         :: stat                              ! status for error
INTEGER,                   DIMENSION(11)        :: iData                             ! Array of Integer, store number of DOF of receiving/sending sinals
! Something wrong with this
REAL(8),                   DIMENSION(dataSize)  :: sData                             ! Array of sending (command) signal
REAL(8),                   DIMENSION(dataSize)  :: rData                             ! Array of receiving (feedback) signal
! Currently, there is no iteration between FAST_SetExternalInputs_Hybrid() and FAST_SetExternalInputs_Hybrid()
! So, if no need to check if time advances or not.
! REAL                                    :: timePast                            ! to cmopared current time step and previsou time step

! Below is not used for this progeam, will be cleaned up later.
!save socketIDs
!save timePast
!data socketIDs /numSockIDs*0/
!data timePast /0.0/
      
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! determine if this is a restart from checkpoint
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALL NWTC_Init() ! open console for writing
ProgName = 'OpenFAST'
InputFile = ""
CheckpointRoot = ""

CALL CheckArgs( InputFile, Flag=FlagArg, Arg2=CheckpointRoot )

IF ( TRIM(FlagArg) == 'H' ) THEN ! Exit after help prompt
    CALL NormStop()

ELSE IF ( TRIM(FlagArg) == 'RESTART' ) THEN ! Restart from checkpoint file
    CALL FAST_RestoreFromCheckpoint_Tary(t_initial, Restart_step, Turbine, CheckpointRoot, ErrStat, ErrMsg  )
    CALL CheckError( ErrStat, ErrMsg, 'during restore from checkpoint'  )            
ELSE
    Restart_step = 0
      
    DO i_turb = 1,NumTurbines
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! initialization
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
        CALL FAST_InitializeAll_T( t_initial, i_turb, Turbine(i_turb), ErrStat, ErrMsg )     ! bjj: we need to get the input files for each turbine (not necessarily the same one)
        CALL CheckError( ErrStat, ErrMsg, 'during module initialization' )
                        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! loose coupling
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
        !...............................................................................................................................
        ! Initialization: (calculate outputs based on states at t=t_initial as well as guesses of inputs and constraint states)
        !...............................................................................................................................     
        CALL FAST_Solution0_T( Turbine(i_turb), ErrStat, ErrMsg )
        CALL CheckError( ErrStat, ErrMsg, 'during simulation initialization'  )
      
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! linearization (bjj: we want to call FAST_Linearize_T whenever WriteOutputToFile is called, but I'll put it at the driver level for now)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! if we need to do linarization analysis at t=0, do it at this operating point 
        CALL FAST_Linearize_T(t_initial, 0, Turbine(i_turb), ErrStat, ErrMsg)
        CALL CheckError( ErrStat, ErrMsg  )
        
    END DO
END IF
   

      
!...............................................................................................................................
! Time Stepping:
!...............................................................................................................................         
   
DO n_t_global = Restart_step, Turbine(1)%p_FAST%n_TMax_m1 
     
    ! bjj: we have to make sure the n_TMax_m1 and n_ChkptTime are the same for all turbines or have some different logic here
        
    ! write checkpoint file if requested
    IF (mod(n_t_global, Turbine(1)%p_FAST%n_ChkptTime) == 0 .AND. Restart_step /= n_t_global) then
        CheckpointRoot = TRIM(Turbine(1)%p_FAST%OutFileRoot)//'.'//TRIM(Num2LStr(n_t_global))
         
        CALL FAST_CreateCheckpoint_Tary(t_initial, n_t_global, Turbine, CheckpointRoot, ErrStat, ErrMsg)
        IF(ErrStat >= AbortErrLev .and. AbortErrLev >= ErrID_Severe) THEN
            ErrStat = MIN(ErrStat,ErrID_Severe) ! We don't need to stop simulation execution on this error
            ErrMsg = TRIM(ErrMsg)//Newline//'WARNING: Checkpoint file could not be generated. Simulation continuing.'
        END IF
            CALL CheckError( ErrStat, ErrMsg  )
    END IF

    ! this takes data from n_t_global and gets values at n_t_global + 1
    DO i_turb = 1,NumTurbines
        ! ##########################################################################################
        ! Receive motion from OpenFresco  @ AS 
        ! ##########################################################################################
        CALL FAST_SetExternalInputs_Hybrid(i_turb, Turbine(i_turb)%m_FAST)
        
        CALL FAST_Solution_T( t_initial, n_t_global, Turbine(i_turb), ErrStat, ErrMsg )
        CALL CheckError( ErrStat, ErrMsg  )
                                   
            
        ! if we need to do linarization analysis, do it at this operating point (which is now n_t_global + 1) 
        ! put this at the end of the loop so that we can output linearization analysis at last OP if desired
        CALL FAST_Linearize_T(t_initial, n_t_global+1, Turbine(i_turb), ErrStat, ErrMsg)
        CALL CheckError( ErrStat, ErrMsg  )
        
        ! ##########################################################################################
        ! Send force to OpenFresco @ AS
        ! ##########################################################################################
        CALL FillOutputAry_T_Hybrid(Turbine(i_turb))
    END DO
      
END DO ! n_t_global
  
  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Write simulation times and stop
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DO i_turb = 1,NumTurbines
    CALL ExitThisProgram_T( Turbine(i_turb), ErrID_None, i_turb==NumTurbines )
END DO
   

    CONTAINS
    !...............................................................................................................................
    ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
    !...............................................................................................................................
    SUBROUTINE CheckError(ErrID,Msg,ErrLocMsg)
   
        ! Passed arguments
        INTEGER(IntKi), INTENT(IN)           :: ErrID       ! The error identifier (ErrStat)
        CHARACTER(*),   INTENT(IN)           :: Msg         ! The error message (ErrMsg)
        CHARACTER(*),   INTENT(IN), OPTIONAL :: ErrLocMsg   ! an optional message describing the location of the error

        CHARACTER(1024)                      :: SimMsg      
        integer(IntKi)                       :: i_turb2
      
      
        IF ( ErrID /= ErrID_None ) THEN
            CALL WrScr( NewLine//TRIM(Msg)//NewLine )
         
            IF ( ErrID >= AbortErrLev ) THEN
            
                IF (PRESENT(ErrLocMsg)) THEN
                   SimMsg = ErrLocMsg
                ELSE
                   SimMsg = 'at simulation time '//TRIM(Num2LStr(Turbine(1)%m_FAST%t_global))//' of '//TRIM(Num2LStr(Turbine(1)%p_FAST%TMax))//' seconds'
                END IF   
                    DO i_turb2 = 1,NumTurbines
                    CALL ExitThisProgram_T( Turbine(i_turb2), ErrID, i_turb2==NumTurbines, SimMsg )
                END DO
                        
             END IF
         
      END IF

    END SUBROUTINE CheckError
    
    !...............................................................................................................................
    ! This subroutine gets motion, which is recorded by DAQ system and sent to OpenFresco, from OpenFresco. Then it assigns the 
    ! motion to the turbine model, m_FAST. 
    !...............................................................................................................................
    SUBROUTINE FAST_SetExternalInputs_Hybrid(iTurb, m_FAST)
    
        USE, INTRINSIC :: ISO_C_Binding
        USE FAST_Types
        ! USE FAST_Data, only: NumFixedInputs
   
        IMPLICIT  NONE

        INTEGER(C_INT),         INTENT(IN   )           :: iTurb                  ! Turbine number 
        !INTEGER(C_INT),         INTENT(IN   )           :: NumInputs_c            ! May 
        !REAL(C_DOUBLE),         INTENT(IN   )           :: InputAry(NumInputs_c)  ! Inputs array from OpenFresco
        TYPE(FAST_MiscVarType), INTENT(INOUT)           :: m_FAST                 ! Miscellaneous variables
   
        INTEGER  :: num_twr_nodes
        INTEGER  :: i
        
        ! May remove this later -----------------------------------------------
        ! set the inputs from external code here...
        ! transfer inputs from Simulink to FAST
        !IF ( NumInputs_c < NumFixedInputs ) RETURN ! This is an error
        
        ! extract socketID
        ! (jtype = user-defined integer value n in element type VUn)
        ! if (jtype .le. numSockIDs) then
        !    socketID = socketIDs(jtype)
        ! else
        !    write(*,*) 'ERROR - Only ',numSockIDs,' genericClient_exp ',
        !*              'elements supported: consider increasing ',
        !*              'numSockIDs parameter in genericClient_exp.for'
        !    call xplb_exit
        ! endif
        ! setup connection with SimAppSiteServer
        ! if (socketID .eq. 0 .and. time(iTotalTime) .gt. 0.0) then
        !IF (socketID .eq. 0 .and. n_t_global .gt. 0.0) THEN
        ! ---------------------------------------------------------------------
        IF (socketID .eq. 0) THEN
            sizeMachineInet = 9+1                               ! length of Internet address including char(0) at the end
            CALL setupconnectionclient(port, &
                                       '127.0.0.1'//char(0), &
                                       sizeMachineInet, &
                                       socketID)
            IF (socketID .le. 0) THEN
             !   call xplb_exit AS:terminate
            ENDIF
             ! socketIDs(jtype) = socketID: AS: not saveing
         
            ! set the data size for the experimental site
            ! sizeCtrl(disp)
            iData(1)  = 0
            ! sizeCtrl(vel)
            iData(2)  = 0
            ! sizeCtrl(accel)
            iData(3)  = 0
            ! sizeCtrl(force)
            iData(4)  = 6
            ! sizeCtrl(time)
            iData(5)  = 1
            ! sizeDaq(disp)
            iData(6)  = 6
            ! sizeDaq(vel)
            iData(7)  = 6
            ! sizeDaq(accel)
            iData(8)  = 0
            ! sizeDaq(force)
            iData(9)  = 0! ndofel
            ! sizeDaq(time)
            iData(10) = 1
            ! dataSize
            iData(11) = dataSize
            CALL senddata(socketID, sizeInt, &
                          iData, 11, stat)
        ENDIF
      
        ! Ask OpenFresc to send measured motion
        sData(1) = 6                                ! OF_RemoteTest_getDaqResponse = 6
        CALL senddata(socketID, sizeDouble,&
                      sData, dataSize, stat)
        ! Receive measured motion from OpenFresco
        CALL recvdata(socketID, sizeDouble,&
                      rData, dataSize, stat)
        
        ! Remove this later --------------
        ! DO i = 1, 12! getting motions
        !    InputAry(i) = rData(i)
        !ENDDO
        ! --------------------------------
        
        ! Assign received motion to wind turbine model
        !IF ( NumInputs_c == NumFixedInputs ) THEN  ! Default for hybrid model use: ElastoDyn inputs
        m_FAST%ExternInput%PtfmSurge    = rData(1)
        m_FAST%ExternInput%PtfmSway     = rData(2)
        m_FAST%ExternInput%PtfmHeave    = rData(3)
        m_FAST%ExternInput%PtfmRoll     = rData(4)
        m_FAST%ExternInput%PtfmPitch    = rData(5)
        m_FAST%ExternInput%PtfmYaw      = rData(6)
        m_FAST%ExternInput%PtfmSurgeVel = rData(7)
        m_FAST%ExternInput%PtfmSwayVel  = rData(8)
        m_FAST%ExternInput%PtfmHeaveVel = rData(9)
        m_FAST%ExternInput%PtfmRollVel  = rData(10)
        m_FAST%ExternInput%PtfmPitchVel = rData(11)
        m_FAST%ExternInput%PtfmYawVel   = rData(12)
        !ENDIF
            
        ! May remove this later ---------------------------------------------------------------------------------------------------------------------
        ! Some other modular configuration is being used for Simulink (@mcd: these functions comprised the traditional use of Simulink with OpenFAST)
        !IF ( NumInputs_c > NumFixedInputs ) THEN
        !    IF ( NumInputs_c == NumFixedInputs + 8 ) THEN  ! ED + SrvD inputs, no IfW inputs
        !        m_FAST%ExternInput%GenTrq      = InputAry(13)
        !        m_FAST%ExternInput%ElecPwr     = InputAry(14)
        !        m_FAST%ExternInput%YawPosCom   = InputAry(15)
        !        m_FAST%ExternInput%YawRateCom  = InputAry(16)
        !        m_FAST%ExternInput%BlPitchCom  = InputAry(17:19)
        !        m_FAST%ExternInput%HSSBrFrac   = InputAry(20)  
        !    ELSEIF ( Numinputs_c == NumFixedInputs + 11 ) THEN  ! SrvD + IfW + ED inputs
        !        m_FAST%ExternInput%GenTrq      = InputAry(13)
        !        m_FAST%ExternInput%ElecPwr     = InputAry(14)
        !        m_FAST%ExternInput%YawPosCom   = InputAry(15)
        !        m_FAST%ExternInput%YawRateCom  = InputAry(16)
        !        m_FAST%ExternInput%BlPitchCom  = InputAry(17:19)
        !        m_FAST%ExternInput%HSSBrFrac   = InputAry(20)  
        !        m_FAST%ExternInput%LidarFocus  = InputAry(21:23)
        !    ENDIF
        !  ENDIF
        ! -------------------------------------------------------------------------------------------------------------------------------------------
        
    END SUBROUTINE FAST_SetExternalInputs_Hybrid
   
    !...............................................................................................................................
    ! This subroutine gets all the simulated force of the turbine at t = n+1 time step. It also commits AeroDyn force and time at 
    ! previous time step, t = n,  before assigning the obtained AeroDyn force and time to sData array.
    !...............................................................................................................................
    SUBROUTINE FillOutputAry_T_Hybrid(Turbine)
    
        TYPE(FAST_TurbineType), INTENT(IN   )   :: Turbine                  ! All data for one instance of a turbine
        REAL(ReKi)                :: Outputs(81)               ! Single array of output
        ! figure out what Reki means (CHECK type.c file in registory)
        INTEGER,                PARAMETER       :: First_idx_AeroDyn = 67 ! I do not know why but this is the index
        INTEGER,                PARAMETER       :: idx_Azimuth = 4 
        INTEGER                                 :: i
        REAL(8),                PARAMETER       :: pi = 4 *ATAN(1.0_8)
        REAL(8)                                :: AzimuthPi
   
        ! Fill Outputs array with simulated force
        CALL FillOutputAry(Turbine%p_FAST, Turbine%y_FAST, Turbine%IfW%y%WriteOutput, Turbine%OpFM%y%WriteOutput, &
                           Turbine%ED%Output(1)%WriteOutput, Turbine%AD%y%WriteOutput, Turbine%SrvD%y%WriteOutput, &
                           Turbine%HD%y%WriteOutput, Turbine%SD%y%WriteOutput, Turbine%ExtPtfm%y%WriteOutput, Turbine%MAP%y%WriteOutput, &
                           Turbine%FEAM%y%WriteOutput, Turbine%MD%y%WriteOutput, Turbine%Orca%y%WriteOutput, &
                           Turbine%IceF%y%WriteOutput, Turbine%IceD%y, Turbine%BD%y, Outputs)
    
    
        ! Commit states variables before updating sData array
        sData(1) = 5                                            ! OF_RemoteTest_commitState = 5
        CALL senddata(socketID, sizeDouble,&
                      sData, dataSize, stat)
        
        ! get Azmith in degree and convert it in rad 
        AzimuthPi = Outputs(idx_Azimuth)*pi/180
        
        ! Extract AeroDyn force from OutPuts and assign it to sData array
        sData(1) = 3                                    ! OF_RemoteTest_setTrialResponse = 3
        
        ! Fill from sData(2) to (7) with RtAeroFxh, RtAeroFyh, RtAeroFzh, RtAeroMxh, RtAeroMyh, and RtAeroMzh.
        sData(2) = Outputs(First_idx_AeroDyn + 1) 
        sData(3) = Outputs(First_idx_AeroDyn + 2)*COS(AzimuthPI) - Outputs(First_idx_AeroDyn + 3)*SIN(AzimuthPI)
        sData(4) = Outputs(First_idx_AeroDyn + 2)*SIN(AzimuthPI) + Outputs(First_idx_AeroDyn + 3)*COS(AzimuthPI)
        sData(5) = Outputs(First_idx_AeroDyn + 4)
        sData(6) = Outputs(First_idx_AeroDyn + 5)*COS(AzimuthPI) - Outputs(First_idx_AeroDyn + 6)*SIN(AzimuthPI)
        sData(7) = Outputs(First_idx_AeroDyn + 5)*SIN(AzimuthPI) + Outputs(First_idx_AeroDyn + 6)*COS(AzimuthPI)
    
        
        ! Assign numerical time to sData array
        sData(8) = n_t_global                     ! Fill sData(8)->(9) with current time
        
        
        ! Send AeroDyn force and time to OpenFresco
        CALL senddata(socketID, sizeDouble,&
                      sData, dataSize, stat)
                        
    END SUBROUTINE FillOutputAry_T_Hybrid
    
!...............................................................................................................................
END PROGRAM FAST
!=======================================================================

    
