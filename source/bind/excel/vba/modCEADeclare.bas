Attribute VB_Name = "modCEADeclare"
Option Explicit

#If Win64 Then
    Private Const CEA_EXCEL_DLL_NAME As String = "cea_excel.dll"
    Private Const LOAD_WITH_ALTERED_SEARCH_PATH As Long = &H8

    Private Declare PtrSafe Function LoadLibraryExW Lib "kernel32" ( _
        ByVal lpLibFileName As LongPtr, _
        ByVal hFile As LongPtr, _
        ByVal dwFlags As Long) As LongPtr
    Private Declare PtrSafe Function SetDllDirectoryW Lib "kernel32" ( _
        ByVal lpPathName As LongPtr) As Long
    Private Declare PtrSafe Function GetLastError Lib "kernel32" () As Long

    Public Declare PtrSafe Function cea_excel_version Lib "cea_excel.dll" ( _
        ByVal versionBuffer As String, _
        ByVal versionLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_test_add Lib "cea_excel.dll" ( _
        ByVal a As Long, _
        ByVal b As Long) As Long
    Public Declare PtrSafe Function cea_excel_last_error Lib "cea_excel.dll" ( _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_to_si Lib "cea_excel.dll" ( _
        ByVal value As Double, _
        ByVal units As String, _
        ByRef siValue As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_weights_from_of Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef baseAmounts As Double, _
        ByRef roles As Long, _
        ByRef bases As Long, _
        ByVal nReactants As Long, _
        ByVal ofRatio As Double, _
        ByRef weights As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_of_from_equivalence Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef baseAmounts As Double, _
        ByRef roles As Long, _
        ByRef bases As Long, _
        ByVal nReactants As Long, _
        ByVal equivalence As Double, _
        ByRef ofRatio As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_of_from_phi Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef baseAmounts As Double, _
        ByRef roles As Long, _
        ByRef bases As Long, _
        ByVal nReactants As Long, _
        ByVal phi As Double, _
        ByRef ofRatio As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_equivalence_from_of Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef baseAmounts As Double, _
        ByRef roles As Long, _
        ByRef bases As Long, _
        ByVal nReactants As Long, _
        ByVal ofRatio As Double, _
        ByRef equivalence As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_phi_from_of Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef baseAmounts As Double, _
        ByRef roles As Long, _
        ByRef bases As Long, _
        ByVal nReactants As Long, _
        ByVal ofRatio As Double, _
        ByRef phi As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_weights_from_moles Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef moles As Double, _
        ByVal nReactants As Long, _
        ByRef weights As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_moles_from_weights Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef weights As Double, _
        ByVal nReactants As Long, _
        ByRef moles As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_per_mole_from_per_weight Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef perWeight As Double, _
        ByVal nReactants As Long, _
        ByRef perMole As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_per_weight_from_per_mole Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef perMole As Double, _
        ByVal nReactants As Long, _
        ByRef perWeight As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_calc_thermo Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, _
        ByRef weights As Double, _
        ByVal nReactants As Long, _
        ByVal propertyType As Long, _
        ByRef temperatures As Double, _
        ByVal nTemperatures As Long, _
        ByVal pressure As Double, _
        ByVal usePressure As Long, _
        ByRef value As Double, _
        ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_eq_solve Lib "cea_excel.dll" ( _
        ByVal eqType As Long, ByVal reactantNames As String, _
        ByRef baseAmounts As Double, ByRef roles As Long, ByRef bases As Long, _
        ByVal nReactants As Long, ByRef explicitWeights As Double, _
        ByVal explicitWeightsLen As Long, ByVal amountMode As Long, _
        ByVal amountValue As Double, ByRef reactantTemperatures As Double, _
        ByVal nReactantTemperatures As Long, ByVal state1 As Double, _
        ByVal state2 As Double, ByVal onlyNames As String, ByVal omitNames As String, _
        ByVal insertNames As String, ByVal propertyNames As String, _
        ByVal speciesNames As String, ByVal transport As Long, ByVal ions As Long, _
        ByVal trace As Double, ByRef values As Double, ByVal valuesCap As Long, _
        ByVal headers As String, ByVal headersLen As Long, ByRef nValues As Long, _
        ByRef converged As Long, ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_rocket_solve Lib "cea_excel.dll" ( _
        ByVal finiteArea As Long, ByVal reactantNames As String, _
        ByRef baseAmounts As Double, ByRef roles As Long, ByRef bases As Long, _
        ByVal nReactants As Long, ByRef explicitWeights As Double, _
        ByVal explicitWeightsLen As Long, ByVal amountMode As Long, _
        ByVal amountValue As Double, ByRef reactantTemperatures As Double, _
        ByVal nReactantTemperatures As Long, ByVal pc As Double, _
        ByRef piP As Double, ByVal nPiP As Long, ByRef subar As Double, _
        ByVal nSubar As Long, ByRef supar As Double, ByVal nSupar As Long, _
        ByVal nFrz As Long, ByVal hcOrTc As Double, ByVal useHc As Long, _
        ByVal mdotOrAcat As Double, ByVal useMdot As Long, ByVal tcEst As Double, _
        ByVal useTcEst As Long, ByVal omitNames As String, ByVal insertNames As String, _
        ByVal propertyNames As String, ByVal speciesNames As String, _
        ByVal transport As Long, ByVal ions As Long, ByVal trace As Double, _
        ByRef values As Double, ByVal valuesCap As Long, ByVal headers As String, _
        ByVal headersLen As Long, ByRef nValues As Long, ByRef converged As Long, _
        ByVal messageBuffer As String, ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_shock_solve Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, ByRef baseAmounts As Double, _
        ByRef roles As Long, ByRef bases As Long, ByVal nReactants As Long, _
        ByRef explicitWeights As Double, ByVal explicitWeightsLen As Long, _
        ByVal amountMode As Long, ByVal amountValue As Double, ByVal T0 As Double, _
        ByVal p0 As Double, ByVal mach1OrU1 As Double, ByVal useMach As Long, _
        ByVal reflected As Long, ByVal incidentFrozen As Long, _
        ByVal reflectedFrozen As Long, ByVal omitNames As String, _
        ByVal insertNames As String, ByVal propertyNames As String, _
        ByVal speciesNames As String, ByVal transport As Long, ByVal ions As Long, _
        ByVal trace As Double, ByRef values As Double, ByVal valuesCap As Long, _
        ByVal headers As String, ByVal headersLen As Long, ByRef nValues As Long, _
        ByRef converged As Long, ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_detonation_solve Lib "cea_excel.dll" ( _
        ByVal reactantNames As String, ByRef baseAmounts As Double, _
        ByRef roles As Long, ByRef bases As Long, ByVal nReactants As Long, _
        ByRef explicitWeights As Double, ByVal explicitWeightsLen As Long, _
        ByVal amountMode As Long, ByVal amountValue As Double, ByVal T1 As Double, _
        ByVal p1 As Double, ByVal frozen As Long, ByVal omitNames As String, _
        ByVal insertNames As String, ByVal propertyNames As String, _
        ByVal speciesNames As String, ByVal transport As Long, ByVal ions As Long, _
        ByVal trace As Double, ByRef values As Double, ByVal valuesCap As Long, _
        ByVal headers As String, ByVal headersLen As Long, ByRef nValues As Long, _
        ByRef converged As Long, ByVal messageBuffer As String, _
        ByVal messageLen As Long) As Long

    Private mCEALibraryHandle As LongPtr
    Private mCEALibraryPath As String
    Private mCEALibraryLoadError As String

    Public Function EnsureCEALibraryLoaded(Optional ByRef errorMessage As String) As Boolean
        Dim candidatePaths As Variant
        Dim candidatePath As Variant
        Dim dllFolder As String
        Dim details As String
        Dim loadErrorCode As Long
        Dim setDirectoryErrorCode As Long
        Dim workbookFolder As String

        If mCEALibraryHandle <> 0 Then
            EnsureCEALibraryLoaded = True
            errorMessage = vbNullString
            Exit Function
        End If

        workbookFolder = ThisWorkbook.Path
        If Len(workbookFolder) = 0 Then
            mCEALibraryLoadError = "Save the workbook before loading " & _
                                   CEA_EXCEL_DLL_NAME & ". The loader resolves " & _
                                   "the native DLL relative to the workbook folder."
            errorMessage = mCEALibraryLoadError
            EnsureCEALibraryLoaded = False
            Exit Function
        End If

        candidatePaths = Array( _
            CEAJoinPath(workbookFolder, CEA_EXCEL_DLL_NAME), _
            CEAJoinPath(workbookFolder, "lib\" & CEA_EXCEL_DLL_NAME), _
            CEAJoinPath(workbookFolder, "lib\Release\" & CEA_EXCEL_DLL_NAME), _
            CEAJoinPath(workbookFolder, "Release\" & CEA_EXCEL_DLL_NAME))

        For Each candidatePath In candidatePaths
            If Len(Dir$(CStr(candidatePath), vbNormal)) > 0 Then
                dllFolder = CEAFolderName(CStr(candidatePath))
                setDirectoryErrorCode = 0
                If SetDllDirectoryW(StrPtr(dllFolder)) = 0 Then
                    setDirectoryErrorCode = GetLastError()
                End If

                mCEALibraryHandle = LoadLibraryExW( _
                    StrPtr(CStr(candidatePath)), _
                    0, _
                    LOAD_WITH_ALTERED_SEARCH_PATH)
                If mCEALibraryHandle <> 0 Then
                    mCEALibraryPath = CStr(candidatePath)
                    mCEALibraryLoadError = vbNullString
                    errorMessage = vbNullString
                    EnsureCEALibraryLoaded = True
                    Exit Function
                End If

                loadErrorCode = GetLastError()
                details = details & vbCrLf & "  " & CStr(candidatePath) & _
                          " (LoadLibraryExW error " & CStr(loadErrorCode)
                If setDirectoryErrorCode <> 0 Then
                    details = details & "; SetDllDirectoryW error " & _
                              CStr(setDirectoryErrorCode)
                End If
                details = details & ")"
            Else
                details = details & vbCrLf & "  " & CStr(candidatePath) & _
                          " (file not found)"
            End If
        Next candidatePath

        mCEALibraryLoadError = "Unable to load " & CEA_EXCEL_DLL_NAME & "." & _
                               vbCrLf & vbCrLf & _
                               "Place " & CEA_EXCEL_DLL_NAME & " in the same " & _
                               "folder as this workbook, or in a lib or " & _
                               "lib\Release folder below it. Keep any native " & _
                               "DLL dependencies in that same folder." & _
                               vbCrLf & vbCrLf & "Checked:" & details
        errorMessage = mCEALibraryLoadError
        EnsureCEALibraryLoaded = False
    End Function

    Public Function CEALibraryPath() As String
        CEALibraryPath = mCEALibraryPath
    End Function

    Public Function CEALibraryLoadError() As String
        CEALibraryLoadError = mCEALibraryLoadError
    End Function

    Private Function CEAJoinPath(ByVal folderPath As String, ByVal childPath As String) As String
        If Right$(folderPath, 1) = "\" Then
            CEAJoinPath = folderPath & childPath
        Else
            CEAJoinPath = folderPath & "\" & childPath
        End If
    End Function

    Private Function CEAFolderName(ByVal filePath As String) As String
        Dim separatorPosition As Long

        separatorPosition = InStrRev(filePath, "\")
        If separatorPosition > 0 Then
            CEAFolderName = Left$(filePath, separatorPosition - 1)
        Else
            CEAFolderName = vbNullString
        End If
    End Function
#Else
    ' The Excel wrapper currently supports 64-bit Windows VBA only.
#End If
