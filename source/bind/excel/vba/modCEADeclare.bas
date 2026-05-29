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
