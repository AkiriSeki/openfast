!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2018  Envision Energy USA, LTD
!
!    This file is part of AeroDyn.
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
module AeroDyn_Model_Subs
   
   use AeroDyn_Driver_Types   
   use AeroDyn
   use VersionInfo

   implicit none   
   
   ! state array indices
   INTEGER(IntKi), PARAMETER :: current_state           = 1     ! index for current states (i.e. t_global)
   INTEGER(IntKi), PARAMETER :: predicted_state         = 2     ! index for predicted states (i.e. t_global_next)
   
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDyn_driver', '', '' )  ! The version number of this program.
                                                    
   contains

!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_Init(DvrData,errStat,errMsg )

   type(Dvr_SimData),            intent(  out) :: DvrData       ! driver data
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Dvr_Init'

   CHARACTER(1000)                             :: inputFile     ! String to hold the file name.
   CHARACTER(200)                              :: git_commit    ! String containing the current git commit hash

   TYPE(ProgDesc), PARAMETER                   :: version   = ProgDesc( 'AeroDyn Driver', '', '' )  ! The version number of this program.

   ErrStat = ErrID_None
   ErrMsg  = ""


   DvrData%OutFileData%unOutFile   = -1
   
   CALL NWTC_Init()
      ! Display the copyright notice
   CALL DispCopyrightLicense( version )   
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetNVD( version )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )

   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, ErrStat2 )
   IF (LEN_TRIM(InputFile) == 0) THEN ! no input file was specified
      call SetErrStat(ErrID_Fatal, 'The required input file was not specified on the command line.', ErrStat, ErrMsg, RoutineName) 
      
         !bjj:  if people have compiled themselves, they should be able to figure out the file name, right?         
      IF (BITS_IN_ADDR==32) THEN
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver_Win32.exe' )
      ELSEIF( BITS_IN_ADDR == 64) THEN
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver_x64.exe' )
      ELSE
         CALL NWTC_DisplaySyntax( InputFile, 'AeroDyn_Driver.exe' )
      END IF
         
      return
   END IF        
         
      ! Read the AeroDyn driver input file
   call Dvr_ReadInputFile(inputFile, DvrData, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 
      if (errStat >= AbortErrLev) return
      
      ! validate the inputs
   call ValidateInputs(DvrData, errStat2, errMsg2)      
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 
      
end subroutine Dvr_Init 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Init_AeroDyn(iCase, DvrData, AD, PhysData, dt, Phys_HubFile, Phys_TwrFile, errStat, errMsg)

   integer(IntKi),                 intent(in   ) :: iCase          ! driver case
   type(Dvr_SimData),              intent(inout) :: DvrData        ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(AeroDyn_Data),             intent(inout) :: AD             ! AeroDyn data
   type(AD_InputType),             intent(inout) :: PhysData       ! Physical model data
   real(DbKi),                     intent(inout) :: dt             ! interval
   character(*),                   intent(  out) :: Phys_HubFile   ! Name of file containing current physical hub data
   character(*),                   intent(  out) :: Phys_TwrFile   ! Name of file containing current physical tower data
      
   integer(IntKi)                , intent(  out) :: errStat        ! Status of error message
   character(*)                  , intent(  out) :: errMsg         ! Error message if ErrStat /= ErrID_None

      ! locals
   real(reKi)                                    :: theta(3)
   integer(IntKi)                                :: j, k   
   integer(IntKi)                                :: errStat2       ! local status of error message
   character(ErrMsgLen)                          :: errMsg2        ! local error message if ErrStat /= ErrID_None
   character(*), parameter                       :: RoutineName = 'Init_AeroDyn'
                                                  
   ! local data                                
   type(AD_InitInputType)                        :: InitInData     ! Input data for initialization
   type(AD_InitOutputType)                       :: InitOutData    ! Output data from initialization
      
      
   errStat = ErrID_None
   errMsg  = ''
   
   InitInData%InputFile      = DvrData%AD_InputFile
   InitInData%NumBlades      = DvrData%numBlades
   InitInData%RootName       = DvrData%outFileData%Root
   InitInData%Gravity        = 9.80665_ReKi                
   
      ! set initialization data:
   call AllocAry( InitInData%BladeRootPosition, 3, InitInData%NumBlades, 'BladeRootPosition', errStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( InitInData%BladeRootOrientation, 3, 3, InitInData%NumBlades, 'BladeRootOrientation', errStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
      
   InitInData%HubPosition = (/ DvrData%Overhang * cos(DvrData%shftTilt), 0.0_ReKi, DvrData%HubHt /)
   theta(1) = 0.0_ReKi
   theta(2) = -DvrData%shftTilt
   theta(3) = 0.0_ReKi
   InitInData%HubOrientation = EulerConstruct( theta ) ! @mcd: I think this is fine if we're assuming shftTilt=0
   
   ! @mcd: these variables might be better suited within a type, but I'll leave it like this for now.
   Phys_HubFile = DvrData%Phys_HubFile
   Phys_TwrFile = DvrData%Phys_TwrFile
     
   
   do k=1,InitInData%numBlades
                     
      theta(1) = (k-1)*TwoPi/real(InitInData%numBlades,ReKi)
      theta(2) = DvrData%precone
      theta(3) = 0.0_ReKi
      InitInData%BladeRootOrientation(:,:,k) = matmul( EulerConstruct( theta ), InitInData%HubOrientation )
                  
      InitInData%BladeRootPosition(:,k)   = InitInData%HubPosition + DvrData%hubRad * InitInData%BladeRootOrientation(3,:,k)      
      
   end do

   ! @mcd: Dustin edited the input args for x, xd, z, and OtherState to be at STATE_CURR. I think it is related to the temporary time stepping methods they use, but keep that in mind for debugging.
   call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
   call PhysMod_Init(AD%p, InitInData, InitOutData, PhysData, ErrStat2, ErrMsg2)
   
   ! @mcd: Dustin also added calls to AD_Copy(Cont/Disc/Constr/Other)State here, I think due to the "temporary steps forward" thing for the loose integration coupling w/ Proteus.
   !       I don't think I need to add this for our purposes.
      
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if   
         
   do j = 2, numInp
      call AD_CopyInput (AD%u(1),  AD%u(j),  MESH_NEWCOPY, errStat2, errMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end do
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
   
      ! we know exact values, so we're going to initialize inputs this way (instead of using the input guesses from AD_Init)
      ! @mcd: since this is just setting up the framework/initial extrapolation for the inputs, just do normal Set_AD_Inputs instead of the modified versions.
      ! (the physical motions for the first time step will be read in during the time marching scheme)
   AD%InputTime = -999
   DO j = 1-numInp, 0
        call Set_AD_Inputs(iCase,j,DvrData,AD,errStat2,errMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END DO              
      
      ! move AD initOut data to AD Driver
   call move_alloc( InitOutData%WriteOutputHdr, DvrData%OutFileData%WriteOutputHdr )
   call move_alloc( InitOutData%WriteOutputUnt, DvrData%OutFileData%WriteOutputUnt )   
     
   DvrData%OutFileData%AD_ver = InitOutData%ver
   
contains
   subroutine cleanup()
      call AD_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )   
      call AD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )      
   end subroutine cleanup
   
end subroutine Init_AeroDyn

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes PhysData meshes and variables for use during the simulation.
subroutine PhysMod_Init(p, InitInp, InitOut, PhysData, ErrStat, ErrMsg)
    type(AD_ParameterType)  ,  intent(in   )  :: p                  ! Parameters
    type(AD_InitInputType)  ,  intent(in   )  :: InitInp            ! Input data for AD initialization routine
    type(AD_InitOutputType) ,  intent(in   )  :: InitOut            ! Output data from initialization
    type(AD_InputType)      ,  intent(  out)  :: PhysData           ! Physical model data
    integer(IntKi)          ,  intent(  out)  :: errStat            ! Error status of the operation
    character(*)            ,  intent(  out)  :: errMsg             ! Error message if ErrStat /= ErrID_None

    ! Local variables
    integer                                   :: j, k               ! loop variables
    real(ReKi)                                :: position(3)        ! node reference position
    integer(intKi)                            :: ErrStat2           ! Temporary error status
    character(ErrMsgLen)                      :: ErrMsg2            ! Temporary error message
    character(*), parameter                   :: RoutineName = 'PhysMod_Init'
    
    ! Initialize variables for this routine
    errStat = ErrID_None
    errMsg = ""
    
    ! Meshes for motion inputs
    
        !.............
        ! tower
        !.............
    if (p%NumTwrNds > 0) then
        call MeshCreate ( BlankMesh         = PhysData%TowerMotion  &
                         ,IOS               = COMPONENT_INPUT       &
                         ,Nnodes            = p%NumTwrNds           &
                         ,ErrStat           = ErrStat2              &
                         ,ErrMess           = ErrMsg2               &
                         ,Orientation       = .true.                &
                         ,TranslationDisp   = .true.                &
                         ,TranslationVel    = .true.                &
                        )
            call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
            
        if (errStat >= AbortErrLev) return

    ! set node initial position/orientation
    position = 0.0_ReKi
    do j=1,p%NumTwrNds 
        position(3) = InitOut%TwrElev(j) ! @mcd: TwrElev in AD input file should match sensor locations on physical model
         
        call MeshPositionNode(PhysData%TowerMotion, j, position, errStat2, errMsg2)  ! orientation is identity by default
        call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
    end do !j
         
    ! create line2 elements
    do k=1,p%NumTwrNds-1
        call MeshConstructElement(PhysData%TowerMotion, ELEMENT_LINE2, errStat2, errMsg2, p1=k, p2=k+1)
        call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
    end do !j
            
    call MeshCommit(PhysData%TowerMotion, errStat2, errMsg2 )
        call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
    if (errStat >= AbortErrLev) return

      
    PhysData%TowerMotion%Orientation     = PhysData%TowerMotion%RefOrientation
    PhysData%TowerMotion%TranslationDisp = 0.0_R8Ki
    PhysData%TowerMotion%TranslationVel  = 0.0_ReKi
      
    end if
   
        !................
        ! hub
        !................
   
    call MeshCreate (BlankMesh         = PhysData%HubMotion    &
                    ,IOS                = COMPONENT_INPUT       &
                    ,Nnodes             = 1                     &
                    ,ErrStat            = ErrStat2              &
                    ,ErrMess            = ErrMsg2               &
                    ,Orientation        = .true.                &
                    ,TranslationDisp    = .true.                &
                    ,RotationVel        = .true.                &
                    )
        call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

    if (errStat >= AbortErrLev) return
                     
    call MeshPositionNode(PhysData%HubMotion, 1, InitInp%HubPosition, errStat2, errMsg2, InitInp%HubOrientation)
        call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
    call MeshConstructElement(PhysData%HubMotion, ELEMENT_POINT, errStat2, errMsg2, p1=1 )
        call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
    call MeshCommit(PhysData%HubMotion, errStat2, errMsg2 )
        call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
    if (errStat >= AbortErrLev) return

         
    PhysData%HubMotion%Orientation     = PhysData%HubMotion%RefOrientation
    PhysData%HubMotion%TranslationDisp = 0.0_R8Ki ! @mcd: originally was 0.0_R8Ki, but changed it to match what Set_AD_Inputs initially sets. Keep in mind for debugging.
    PhysData%HubMotion%RotationVel     = 0.0_ReKi
    
    end subroutine PhysMod_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads the present data from the physical model and formats it into meshes that match the standard for AeroDyn.
subroutine PhysMod_Get_Physical_Motions(PhysData, HubUn, TwrUn)
   type(AD_InputType)          , intent(inout) :: PhysData      ! Physical model data
   integer                     , intent(in   ) :: HubUn         ! Unit for the input hub file
   integer                     , intent(in   ) :: TwrUn         ! Unit for the input tower file
   ! local variables
   character(*), parameter                     :: RoutineName = 'PhysMod_Get_Physical_Motions'
   integer                                     :: j, k
   
   ! Map data to AeroDyn-readable type
   ! The hub data is in format Position(1x3), Orientation(3x3), TranslationDisp(1x3), RotationVel(1x3)

   read(HubUn, *) PhysData%HubMotion%Position(:,1)
   do j = 1, 3
       read(HubUn, *) PhysData%HubMotion%Orientation(:,j,1)
   end do
   read(HubUn, *) PhysData%HubMotion%TranslationDisp(:,1)
   read(HubUn, *) PhysData%HubMotion%RotationVel(:,1)
   close(HubUn, status='DELETE')
   
   ! @mcd: I'm doing this by tower node for now, but I doubt we will have enough sensors to analyze many point along the tower, so this will almost certainly change.
   !       Not to mention the eventual partial integration with ElastoDyn.
   ! The tower data is in format Position(1x3), Orientation(3x3), TranslationDisp(1x3), TranslationVel(1x3), working node by node

   do j = 1, PhysData%TowerMotion%NNodes
       read(TwrUn, *) PhysData%TowerMotion%Position(:,j)
       do k = 1, 3
           read(TwrUn, *) PhysData%TowerMotion%Orientation(:,k,j)
       end do
       read(TwrUn, *) PhysData%TowerMotion%TranslationDisp(:,j)
       read(TwrUn, *) PhysData%TowerMotion%TranslationVel(:,j)
   end do
   close(TwrUn, status='DELETE')
   
   end subroutine PhysMod_Get_Physical_Motions
!----------------------------------------------------------------------------------------------------------------------------------
!> this routine returns time=(nt-1) * DvrData%Cases(iCase)%dT, and cycles values in the input array AD%InputTime and AD%u.
!! it then sets the inputs for nt * DvrData%Cases(iCase)%dT, which are index values 1 in the arrays.
subroutine Set_AD_Inputs(iCase,nt,DvrData,AD,errStat,errMsg)

   integer(IntKi)              , intent(in   ) :: iCase         ! case number 
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   
   type(Dvr_SimData),            intent(inout) :: DvrData       ! Driver data 
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Set_AD_Inputs'

   integer(intKi)                              :: j             ! loop counter for nodes
   integer(intKi)                              :: k             ! loop counter for blades

   real(ReKi)                                  :: z             ! height (m)
   !real(ReKi)                                  :: angle
   real(R8Ki)                                  :: theta(3)
   real(R8Ki)                                  :: position(3)
   real(R8Ki)                                  :: orientation(3,3)
   real(R8Ki)                                  :: rotateMat(3,3)
   
   
   errStat = ErrID_None
   errMsg  = ""
   
      ! note that this initialization is a little different than the general algorithm in FAST because here
      ! we can get exact values, so we are going to ignore initial guesses and not extrapolate
      
   !................
   ! shift previous calculations:
   !................
   do j = numInp-1,1,-1
      call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      AD%InputTime(j+1) = AD%InputTime(j)
   end do
   AD%inputTime(1) = nt * DvrData%Cases(iCase)%dT
         
   !................
   ! calculate new values
   !................
   
      ! Tower motions:
      do j=1,AD%u(1)%TowerMotion%nnodes
         AD%u(1)%TowerMotion%Orientation(  :,:,j) = AD%u(1)%TowerMotion%RefOrientation(:,:,j) ! identity
         AD%u(1)%TowerMotion%TranslationDisp(:,j) = 0.0_ReKi
         AD%u(1)%TowerMotion%TranslationVel( :,j) = 0.0_ReKi
      end do !j=nnodes
      
      ! Hub motions:
      theta(1) = 0.0_ReKi
      theta(2) = 0.0_ReKi
      theta(3) = DvrData%Cases(iCase)%Yaw
      orientation = EulerConstruct(theta)
            
      AD%u(1)%HubMotion%TranslationDisp(:,1) = matmul( AD%u(1)%HubMotion%Position(:,1), orientation ) - AD%u(1)%HubMotion%Position(:,1) ! = matmul( transpose(orientation) - eye(3), AD%u(1)%HubMotion%Position(:,1) )
      
      theta(1) = AD%inputTime(1) * DvrData%Cases(iCase)%RotSpeed
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi
      AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( AD%u(1)%HubMotion%RefOrientation(:,:,1), orientation )
      orientation = EulerConstruct( theta )      
      AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%HubMotion%Orientation(  :,:,1) )
      
      AD%u(1)%HubMotion%RotationVel(    :,1) = AD%u(1)%HubMotion%Orientation(1,:,1) * DvrData%Cases(iCase)%RotSpeed
                  
      ! Blade root motions:
      do k=1,DvrData%numBlades         
         theta(1) = (k-1)*TwoPi/real(DvrData%numBlades,ReKi)
         theta(2) =  DvrData%precone
         theta(3) = -DvrData%Cases(iCase)%pitch
         orientation = EulerConstruct(theta)
         
         AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%HubMotion%Orientation(  :,:,1) )
         
      end do !k=numBlades
            
      ! Blade motions:
      do k=1,DvrData%numBlades
         rotateMat = transpose( AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) )
         rotateMat = matmul( rotateMat, AD%u(1)%BladeRootMotion(k)%RefOrientation(  :,:,1) ) 
         orientation = transpose(rotateMat)
         
         rotateMat(1,1) = rotateMat(1,1) - 1.0_ReKi
         rotateMat(2,2) = rotateMat(2,2) - 1.0_ReKi
         rotateMat(3,3) = rotateMat(3,3) - 1.0_ReKi
                  
         do j=1,AD%u(1)%BladeMotion(k)%nnodes        
            position = AD%u(1)%BladeMotion(k)%Position(:,j) - AD%u(1)%HubMotion%Position(:,1) 
            AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) = AD%u(1)%HubMotion%TranslationDisp(:,1) + matmul( rotateMat, position )
            
            AD%u(1)%BladeMotion(k)%Orientation(  :,:,j) = matmul( AD%u(1)%BladeMotion(k)%RefOrientation(:,:,j), orientation )
            
            
            position =  AD%u(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) &
                      - AD%u(1)%HubMotion%Position(:,1) - AD%u(1)%HubMotion%TranslationDisp(:,1)
            AD%u(1)%BladeMotion(k)%TranslationVel( :,j) = cross_product( AD%u(1)%HubMotion%RotationVel(:,1), position )

         end do !j=nnodes
                                    
      end do !k=numBlades       
      
      ! Inflow wind velocities:
      ! InflowOnBlade
      do k=1,DvrData%numBlades
         do j=1,AD%u(1)%BladeMotion(k)%nnodes
            z = AD%u(1)%BladeMotion(k)%Position(3,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(3,j)
            AD%u(1)%InflowOnBlade(1,j,k) = GetU(  DvrData%Cases(iCase)%WndSpeed, DvrData%HubHt, DvrData%Cases(iCase)%ShearExp, z )
            AD%u(1)%InflowOnBlade(2,j,k) = 0.0_ReKi !V
            AD%u(1)%InflowOnBlade(3,j,k) = 0.0_ReKi !W      
         end do !j=nnodes
      end do !k=numBlades
      
      !InflowOnTower
      do j=1,AD%u(1)%TowerMotion%nnodes
         z = AD%u(1)%TowerMotion%Position(3,j) + AD%u(1)%TowerMotion%TranslationDisp(3,j)
         AD%u(1)%InflowOnTower(1,j) = GetU(  DvrData%Cases(iCase)%WndSpeed, DvrData%HubHt, DvrData%Cases(iCase)%ShearExp, z )
         AD%u(1)%InflowOnTower(2,j) = 0.0_ReKi !V
         AD%u(1)%InflowOnTower(3,j) = 0.0_ReKi !W         
      end do !j=nnodes
                     
end subroutine Set_AD_Inputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine cycles InputTime values and sets all the AeroDyn motion inputs (no wind inflow values).
!!  For the time being, wind inflows are handled by Set_AD_Inflows.
subroutine Set_AD_Motion_Inputs_NoIfW(iCase,nt,DvrData,AD,PhysData,errStat,errMsg)
   integer(IntKi)              , intent(in   ) :: iCase         ! case number 
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   
   type(Dvr_SimData)           , intent(inout) :: DvrData       ! Driver data 
   type(AeroDyn_Data)          , intent(inout) :: AD            ! AeroDyn data 
   type(AD_InputType)          , intent(inout) :: PhysData      ! Physical model data
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Set_AD_Motion_Inputs_NoIfW'

   integer(intKi)                              :: j             ! loop counter for nodes
   integer(intKi)                              :: k             ! loop counter for blades

   real(ReKi)                                  :: z             ! height (m)
   real(R8Ki)                                  :: theta(3)
   real(R8Ki)                                  :: position(3)
   real(R8Ki)                                  :: orientation(3,3)
   real(R8Ki)                                  :: rotateMat(3,3)
   
   
   ! note that this initialization is a little different than the general algorithm in FAST because here
   ! we can get exact values, so we are going to ignore initial guesses and not extrapolate
      
   !................
   ! shift previous calculations:
   !................
   ! @mcd: this sets it so AD%InputTime(2) = (nt-1) * dT, while AD%InputTime(1) = nt * dT (yes, this seems completely backwards but it's what they did)
   do j = numInp-1,1,-1
      call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      AD%InputTime(j+1) = AD%InputTime(j)
   end do

   AD%inputTime(1) = nt * DvrData%Cases(iCase)%dT
         
   !.........................................
   ! Set the inputs from the physical model:
   !.........................................
      ! @mcd: I will eventually replace much of this section with material similar to AD_InputSolve_NoIfW from the FAST solution once I start incorporating other modules
   
      ! Tower motions:
      do j=1,AD%u(1)%TowerMotion%nnodes
         AD%u(1)%TowerMotion%Position(       :,j) = PhysData%TowerMotion%Position(:,j)
         AD%u(1)%TowerMotion%Orientation(  :,:,j) = PhysData%TowerMotion%Orientation(:,:,j)  ! this will likely need to be updated at some point, since the physical model orientation may not capture yaw well.
      end do !j=nnodes
      
      ! Hub motions:
      ! @mcd: Position can be fed in directly, but orientation will still be dependent on controls since this won't be on the physical model.
      theta(1) = 0.0_ReKi
      theta(2) = 0.0_ReKi
      theta(3) = DvrData%Cases(iCase)%Yaw
      orientation = EulerConstruct(theta)

      AD%u(1)%HubMotion%Position(:,1) = PhysData%HubMotion%Position(:,1)
      AD%u(1)%HubMotion%Orientation(:,:,1) = PhysData%HubMotion%Orientation(:,:,1)
      AD%u(1)%HubMotion%TranslationDisp(:,1) = matmul(AD%u(1)%HubMotion%Position(:,1), orientation) - AD%u(1)%HubMotion%Position(:,1)    
      
      theta(1) = AD%inputTime(1) * DvrData%Cases(iCase)%RotSpeed ! this will need to be updated once ServoDyn is incorporated
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi
      AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( AD%u(1)%HubMotion%RefOrientation(:,:,1), orientation )
      orientation = EulerConstruct( theta )      
      AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%HubMotion%Orientation(  :,:,1) )
      
      AD%u(1)%HubMotion%RotationVel(    :,1) = AD%u(1)%HubMotion%Orientation(1,:,1) * DvrData%Cases(iCase)%RotSpeed
                  
      ! Blade root motions:
      ! @mcd: blades will not be physical, so I'm not touching this until we incorporate ServoDyn
      do k=1,DvrData%numBlades         
         theta(1) = (k-1)*TwoPi/real(DvrData%numBlades,ReKi)
         theta(2) =  DvrData%precone
         theta(3) = -DvrData%Cases(iCase)%pitch
         orientation = EulerConstruct(theta)
         
         AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%HubMotion%Orientation(  :,:,1) )
         
      end do !k=numBlades
            
      ! Blade motions:
      ! @mcd: same deal as blade roots.
      do k=1,DvrData%numBlades
         rotateMat = transpose( AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) )
         rotateMat = matmul( rotateMat, AD%u(1)%BladeRootMotion(k)%RefOrientation(  :,:,1) ) 
         orientation = transpose(rotateMat)
         
         rotateMat(1,1) = rotateMat(1,1) - 1.0_ReKi
         rotateMat(2,2) = rotateMat(2,2) - 1.0_ReKi
         rotateMat(3,3) = rotateMat(3,3) - 1.0_ReKi
                  
         do j=1,AD%u(1)%BladeMotion(k)%nnodes        
            position = AD%u(1)%BladeMotion(k)%Position(:,j) - AD%u(1)%HubMotion%Position(:,1) 
            AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) = AD%u(1)%HubMotion%TranslationDisp(:,1) + matmul( rotateMat, position )
            
            AD%u(1)%BladeMotion(k)%Orientation(  :,:,j) = matmul( AD%u(1)%BladeMotion(k)%RefOrientation(:,:,j), orientation )
            
            
            position =  AD%u(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) &
                      - AD%u(1)%HubMotion%Position(:,1) - AD%u(1)%HubMotion%TranslationDisp(:,1)
            AD%u(1)%BladeMotion(k)%TranslationVel( :,j) = cross_product( AD%u(1)%HubMotion%RotationVel(:,1), position )

         end do !j=nnodes
                                    
      end do !k=numBlades      
   
end subroutine Set_AD_Motion_Inputs_NoIfW

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine cycles inflow values in AD%u. This will soon be replaced with AD_InputSolve_IfW from the FAST solution once InflowWind is coupled to the numerical model.
subroutine Set_AD_Inflows(iCase,nt,DvrData,AD,errStat,errMsg)

   integer(IntKi)              , intent(in   ) :: iCase         ! case number 
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   
   type(Dvr_SimData),            intent(inout) :: DvrData       ! Driver data 
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'Set_AD_Inflows'

   integer(intKi)                              :: j             ! loop counter for nodes
   integer(intKi)                              :: k             ! loop counter for blades

   real(ReKi)                                  :: z             ! height (m)
   real                                     :: x1,x2
   real(R8Ki)                                  :: theta(3)
   real(R8Ki)                                  :: position(3)
   real(R8Ki)                                  :: orientation(3,3)
   real(R8Ki)                                  :: rotateMat(3,3)
   
   
   errStat = ErrID_None
   errMsg  = ""
    
      
      ! Inflow wind velocities:
      ! InflowOnBlade
      do k=1,DvrData%numBlades
         do j=1,AD%u(1)%BladeMotion(k)%nnodes
            z = AD%u(1)%BladeMotion(k)%Position(3,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(3,j)
            x1 = AD%u(1)%BladeMotion(k)%Position(3,j)
            x2 = AD%u(1)%BladeMotion(k)%TranslationDisp(3,j)
            AD%u(1)%InflowOnBlade(1,j,k) = GetU(  DvrData%Cases(iCase)%WndSpeed, DvrData%HubHt, DvrData%Cases(iCase)%ShearExp, z )
            AD%u(1)%InflowOnBlade(2,j,k) = 0.0_ReKi !V
            AD%u(1)%InflowOnBlade(3,j,k) = 0.0_ReKi !W      
         end do !j=nnodes
      end do !k=numBlades
      
      !InflowOnTower
      do j=1,AD%u(1)%TowerMotion%nnodes
         z = AD%u(1)%TowerMotion%Position(3,j) + AD%u(1)%TowerMotion%TranslationDisp(3,j)
         AD%u(1)%InflowOnTower(1,j) = GetU(  DvrData%Cases(iCase)%WndSpeed, DvrData%HubHt, DvrData%Cases(iCase)%ShearExp, z )
         AD%u(1)%InflowOnTower(2,j) = 0.0_ReKi !V
         AD%u(1)%InflowOnTower(3,j) = 0.0_ReKi !W         
      end do !j=nnodes
                     
end subroutine Set_AD_Inflows
!----------------------------------------------------------------------------------------------------------------------------------
function GetU( URef, ZRef, PLExp, z ) result (U)
   real(ReKi), intent(in) :: URef
   real(ReKi), intent(in) :: ZRef
   real(ReKi), intent(in) :: PLExp
   real(ReKi), intent(in) :: z
   real(ReKi)             :: U
   
   U = URef*(z/ZRef)**PLExp

end function GetU
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_ReadInputFile(fileName, DvrData, errStat, errMsg )
   ! This routine opens the gets the data from the input files.

   character(*),                  intent( in    )   :: fileName
   type(Dvr_SimData),             intent(   out )   :: DvrData
   integer,                       intent(   out )   :: errStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: errMsg               ! Error message if errStat /= ErrID_None
   

      ! Local variables
   character(1024)              :: PriPath
   character(1024)              :: inpVersion                               ! String containing the input-version information.
   character(1024)              :: line                                     ! String containing a line of input.
   integer                      :: unIn, unEc
   integer                      :: ICase
   integer                      :: Sttus
   character( 11)               :: DateNow                                  ! Date shortly after the start of execution.
   character(  8)               :: TimeNow                                  ! Time of day shortly after the start of execution.
   
   integer, parameter           :: NumCols = 7                              ! number of columns to be read from the input file
   real(DbKi)                   :: InpCase(NumCols)                         ! Temporary array to hold combined-case input parameters. (note that we store in double precision so the time is read correctly)
   logical                      :: TabDel      
   logical                      :: echo   

   INTEGER(IntKi)               :: ErrStat2                                 ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                  ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'Dvr_ReadInputFile'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   UnIn = -1
   UnEc = -1
   
   ! Open the input file
   CALL GetPath( fileName, PriPath )     ! Input files will be relative to the path where the primary input file is located.

   call GetNewUnit( unIn )   
   call OpenFInpFile( unIn, fileName, errStat2, ErrMsg2 )
   call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   if ( errStat >= AbortErrLev ) then
      call cleanup()
      return
   end if

   
   call WrScr( 'Opening input file:  '//fileName )

      ! Skip a line, read the run title information.

   CALL ReadStr( UnIn, fileName, inpVersion, 'inpVersion', 'File Header: (line 1)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL ReadStr( UnIn, fileName, DvrData%OutFileData%runTitle, 'runTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   call WrScr1 ( ' '//DvrData%OutFileData%runTitle )
   
      ! Read in the title line for the input-configuration subsection.
   CALL ReadStr( UnIn, fileName, line, 'line', 'File Header: (line 3)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! See if we should echo the output.     
   call ReadVar ( unIn, fileName, echo, 'Echo', 'Echo Input', errStat2, errMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if ( echo )  then
         ! Get date and time.
      dateNow = CurDate()
      timeNow = CurTime()
      call GetNewUnit( unEc ) 
      call getroot(fileName,DvrData%OutFileData%Root)      
      call  OpenFOutFile ( unEc, trim( DvrData%OutFileData%Root )//'.ech', errStat2, errMsg2 )
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
         if ( errStat >= AbortErrLev ) then
            call cleanup()
            return
         end if
      
      write (unEc,'(A)')      'Echo of Input File:'
      write (unEc,'(A)')      ' "'//fileName//'"'
      write (unEc,'(A)')      'Generated on: '//trim( dateNow )//' at '//trim( timeNow )//'.'
      write (unEc,'(A)')      inpVersion
      write (unEc,'(A)')      DvrData%OutFileData%runTitle
      write (unEc,'(A)')      line
      write (unEc,Ec_LgFrmt)  echo, 'Echo', 'Echo input parameters to "rootname.ech"?'
   end if


      ! Read the rest of input-configuration section.
      
   call ReadVar ( unIn, fileName, DvrData%AD_InputFile,   'AD_InputFile',   'Name of the AeroDyn input file', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
   IF ( PathIsRelative( DvrData%AD_InputFile ) ) DvrData%AD_InputFile = TRIM(PriPath)//TRIM(DvrData%AD_InputFile)

   ! @mcd: Added in physical model output file locations.
   call ReadVar ( unIn, fileName, DvrData%Phys_HubFile,   'Phys_HubFile',   'Name of file containing current physical hub data', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
   IF ( PathIsRelative( DvrData%Phys_HubFile ) ) DvrData%Phys_HubFile = TRIM(PriPath)//TRIM(DvrData%Phys_HubFile)
   
   call ReadVar ( unIn, fileName, DvrData%Phys_TwrFile,   'Phys_TwrFile',   'Name of file containing current physical tower data', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
   IF ( PathIsRelative( DvrData%Phys_TwrFile ) ) DvrData%Phys_TwrFile = TRIM(PriPath)//TRIM(DvrData%Phys_TwrFile)
   
      ! Read the turbine-data section.

   call ReadCom ( unIn, fileName, 'the turbine-data subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%numBlades,'NumBlades','Number of blades', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%HubRad,   'HubRad',   'Hub radius (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%HubHt,    'HubHt',    'Hub height (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%Overhang, 'Overhang',  'Overhang (m)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%ShftTilt, 'ShftTilt',  'Shaft tilt (deg)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      DvrData%ShftTilt = DvrData%ShftTilt*D2R
   call ReadVar ( unIn, fileName, DvrData%precone, 'Precone',  'Precone (deg)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      DvrData%precone = DvrData%precone*D2R
            
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if           


      ! Read the I/O-configuration section.

   call ReadCom ( unIn, fileName, 'the I/O-configuration subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar ( unIn, fileName, DvrData%OutFileData%Root, 'OutFileRoot', 'Root name for any output files', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   if (len_trim(DvrData%OutFileData%Root) == 0) then
      call getroot(fileName,DvrData%OutFileData%Root)
   end if
   
   call ReadVar ( unIn, fileName, TabDel,   'TabDel',   'Make output tab-delimited (fixed-width otherwise)?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
      if (TabDel) then
         DvrData%OutFileData%delim = TAB
      else
         DvrData%OutFileData%delim = " "
      end if
               
      ! OutFmt - Format used for text tabular output (except time).  Resulting field should be 10 characters. (-):
   call ReadVar( UnIn, fileName, DvrData%OutFileData%OutFmt, "OutFmt", "Format used for text tabular output (except time).  Resulting field should be 10 characters. (-)", ErrStat2, ErrMsg2, UnEc)
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName ) !bjj: this is a global variable in NWTC_Library            
   call ReadVar ( unIn, fileName, Beep,  'Beep',     'Beep on exit?', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName ) !bjj: this is a global variable in NWTC_Library
      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if


      ! Read the combined-case section.

   call ReadCom  ( unIn, fileName, 'the combined-case subtitle', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadVar  ( unIn, fileName, DvrData%NumCases, 'NumCases', 'Number of cases to run', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadCom  ( unIn, fileName, 'the combined-case-block header (names)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
   call ReadCom  ( unIn, fileName, 'the combined-case-block header (units)', errStat2, errMsg2, UnEc )
      call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )

      if ( errStat >= AbortErrLev ) then
         call cleanup()
         return
      end if
      
   if ( DvrData%NumCases < 1 )  then
      call setErrStat( ErrID_Fatal,'Variable "NumCases" must be > 0.' ,errstat,errmsg,routinename)
      call cleanup()
      return
   end if
   
   allocate ( DvrData%Cases(DvrData%NumCases) , STAT=Sttus )
   if ( Sttus /= 0 )  then
      call setErrStat( ErrID_Fatal,'Error allocating memory for the Cases array.',errstat,errmsg,routinename)
      call cleanup()
      return
   end if

   do ICase=1,DvrData%NumCases

      call ReadAry ( unIn, fileName, InpCase,  NumCols, 'InpCase',  'parameters for Case #' &
                     //trim( Int2LStr( ICase ) )//'.', errStat2, errMsg2, UnEc )
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
            
      DvrData%Cases(iCase)%WndSpeed        = InpCase( 1)
      DvrData%Cases(ICase)%ShearExp        = InpCase( 2)
      DvrData%Cases(ICase)%RotSpeed        = InpCase( 3)*RPM2RPS
      DvrData%Cases(ICase)%Pitch           = InpCase( 4)*D2R
      DvrData%Cases(ICase)%Yaw             = InpCase( 5)*D2R
      DvrData%Cases(iCase)%dT              = InpCase( 6)
      DvrData%Cases(iCase)%Tmax            = InpCase( 7)
               
   end do ! ICase
   
   call cleanup ( )


   RETURN
contains
   subroutine cleanup()
      if (UnIn>0) close(UnIn)
      if (UnEc>0) close(UnEc)
   end subroutine cleanup
end subroutine Dvr_ReadInputFile
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ValidateInputs(DvrData, errStat, errMsg)

   type(Dvr_SimData),             intent(in)    :: DvrData
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None

   ! local variables:
   integer(intKi)                               :: i
   integer(intKi)                               :: FmtWidth          ! number of characters in string produced by DvrData%OutFmt
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ValidateInputs'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! Turbine Data:
   if ( DvrData%numBlades < 1 ) call SetErrStat( ErrID_Fatal, "There must be at least 1 blade (numBlades).", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%numBlades > 3 ) call SetErrStat( ErrID_Fatal, "There can be no more than 3 blades (numBlades).", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%HubRad < 0.0_ReKi .or. EqualRealNos(DvrData%HubRad, 0.0_ReKi) ) call SetErrStat( ErrID_Fatal, "HubRad must be a positive number.", ErrStat, ErrMsg, RoutineName)
   if ( DvrData%HubHt < DvrData%HubRad ) call SetErrStat( ErrID_Fatal, "HubHt must be at least HubRad.", ErrStat, ErrMsg, RoutineName)
   
      
      ! I-O Settings:
      ! Check that DvrData%OutFileData%OutFmt is a valid format specifier and will fit over the column headings
   call ChkRealFmtStr( DvrData%OutFileData%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if ( FmtWidth /= ChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column width of '// &
      TRIM(Num2LStr(FmtWidth))//' instead of '//TRIM(Num2LStr(ChanLen))//' characters.', ErrStat, ErrMsg, RoutineName )

      ! Combined-Case Analysis:
   do i=1,DvrData%NumCases
   
      if (DvrData%Cases(i)%DT < epsilon(0.0_ReKi) ) call SetErrStat(ErrID_Fatal,'dT must be larger than 0 in case '//trim(num2lstr(i))//'.',ErrStat, ErrMsg,RoutineName)
      if (DvrData%Cases(i)%TMax < DvrData%Cases(i)%DT ) call SetErrStat(ErrID_Fatal,'TMax must be larger than dT in case '//trim(num2lstr(i))//'.',ErrStat, ErrMsg,RoutineName)
      
   end do
   
   
   
end subroutine ValidateInputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputLine(OutFileData, t, output, OutUn, errStat, errMsg)

   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(Dvr_OutputFile)   ,  intent(in   )   :: OutFileData
   real(ReKi)             ,  intent(in   )   :: output(:)            ! Rootname for the output file
   integer                ,  intent(in   )   :: OutUn                ! Unit for the output hybrid file.
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
      
   ! Local variables.

   character(200)                   :: frmt                                      ! A string to hold a format specifier
   character(15)                    :: tmpStr                                    ! temporary string to print the time output as text
   integer :: numOuts
   
   errStat = ErrID_None
   errMsg  = ''
   numOuts = size(output,1)
   frmt = '"'//OutFileData%delim//'"'//trim(OutFileData%outFmt)      ! format for array elements from individual modules
   
      ! time
   write( tmpStr, '(F15.4)' ) t
   call WrFileNR( OutFileData%unOutFile, tmpStr )
   call WrNumAryFileNR ( OutFileData%unOutFile, output,  frmt, errStat, errMsg ) ! @mcd: normal output file
   call WrNumAryFileNR(OutUn, output, frmt, errStat, errMsg) ! @mcd: hybrid interface output file
   if ( errStat >= AbortErrLev ) return
   
     ! write a new line (advance to the next line)
   write (OutFileData%unOutFile,'()')
      
end subroutine Dvr_WriteOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_InitializeOutputFile( iCase, CaseData, OutFileData, errStat, errMsg)
      type(Dvr_OutputFile),     intent(inout)   :: OutFileData 
      
      integer(IntKi)         ,  intent(in   )   :: iCase                ! case number (to write in file description line and use for file name)
      type(Dvr_Case),           intent(in   )   :: CaseData
      
      integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
      character(*)           ,  intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None

         ! locals
      integer(IntKi)                            ::  i      
      integer(IntKi)                            :: numOuts
      
      
      
      call GetNewUnit( OutFileData%unOutFile, ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) then
            OutFileData%unOutFile = -1
            return
         end if
         

      call OpenFOutFile ( OutFileData%unOutFile, trim(outFileData%Root)//'.'//trim(num2lstr(iCase))//'.out', ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) return
         
      write (OutFileData%unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(version))
      write (OutFileData%unOutFile,'(1X,A)') trim(GetNVD(OutFileData%AD_ver))
      write (OutFileData%unOutFile,'()' )    !print a blank line
     ! write (OutFileData%unOutFile,'(A,11(1x,A,"=",ES11.4e2,1x,A))'   ) 'Case '//trim(num2lstr(iCase))//':' &
      write (OutFileData%unOutFile,'(A,11(1x,A,"=",A,1x,A))'   ) 'Case '//trim(num2lstr(iCase))//':' &
         , 'WndSpeed', trim(num2lstr(CaseData%WndSpeed)), 'm/s;' &
         , 'ShearExp', trim(num2lstr(CaseData%ShearExp)), ';' &
         , 'RotSpeed', trim(num2lstr(CaseData%RotSpeed*RPS2RPM)),'rpm;' &
         , 'Pitch',    trim(num2lstr(CaseData%Pitch*R2D)), 'deg;' &
         , 'Yaw',      trim(num2lstr(CaseData%Yaw*R2D)), 'deg;' &
         , 'dT',       trim(num2lstr(CaseData%dT)), 's;' &
         , 'Tmax',     trim(num2lstr(CaseData%Tmax)),'s'
      
      write (OutFileData%unOutFile,'()' )    !print a blank line
              

      numOuts = size(OutFileData%WriteOutputHdr)
         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................

      call WrFileNR ( OutFileData%unOutFile, '     Time           ' )

      do i=1,NumOuts
         call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputHdr(i) )
      end do ! i

      write (OutFileData%unOutFile,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................

      call WrFileNR ( OutFileData%unOutFile, '      (s)           ' )

      do i=1,NumOuts
         call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputUnt(i) )
      end do ! i

      write (OutFileData%unOutFile,'()')      
      

      
end subroutine Dvr_InitializeOutputFile

end module AeroDyn_Model_Subs
