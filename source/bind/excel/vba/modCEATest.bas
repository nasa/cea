Attribute VB_Name = "modCEATest"
Option Explicit

Private Const README_SHEET_NAME As String = "README"
Private Const VERSION_LABEL_CELL As String = "A1"
Private Const VERSION_VALUE_CELL As String = "B2"
Private Const STATUS_VALUE_CELL As String = "B3"
Private Const ADD_STATUS_CELL As String = "B5"

Private Function CallCEAAdd(ByVal a As Long, ByVal b As Long, ByRef result As Long) As Boolean
#If Win64 Then
    Dim loadError As String

    On Error GoTo LoadError

    If Not EnsureCEALibraryLoaded(loadError) Then
        result = 0
        MsgBox loadError, vbCritical, "CEA Excel wrapper"
        Exit Function
    End If

    result = cea_excel_test_add(a, b)
    CallCEAAdd = True
    Exit Function

LoadError:
    result = 0
    MsgBox "Unable to load or call the native CEA Excel add smoke test." & _
           vbCrLf & vbCrLf & _
           "VBA error " & CStr(Err.Number) & ": " & Err.Description, _
           vbCritical, "CEA Excel wrapper"
#Else
    result = 0
    MsgBox "CEA Excel wrapper support is limited to 64-bit Windows VBA.", _
           vbExclamation, "CEA Excel wrapper"
#End If
End Function

Private Function CallCEAVersion(ByRef version As String, ByRef status As Long) As Boolean
#If Win64 Then
    Dim buffer As String
    Dim loadError As String
    Dim nulPos As Long

    On Error GoTo LoadError

    If Not EnsureCEALibraryLoaded(loadError) Then
        version = ""
        status = -1
        MsgBox loadError, vbCritical, "CEA Excel wrapper"
        Exit Function
    End If

    buffer = String$(256, vbNullChar)
    status = cea_excel_version(buffer, Len(buffer))

    nulPos = InStr(1, buffer, vbNullChar)
    If nulPos > 0 Then
        version = Left$(buffer, nulPos - 1)
    Else
        version = buffer
    End If

    CallCEAVersion = True
    Exit Function

LoadError:
    version = ""
    status = -1
    MsgBox "Unable to load or call the native CEA Excel library." & vbCrLf & _
           "The loader tried workbook-relative locations before this call. " & _
           "Confirm " & CEALibraryPath() & " and its native dependencies " & _
           "are available." & vbCrLf & vbCrLf & _
           "VBA error " & CStr(Err.Number) & ": " & Err.Description, _
           vbCritical, "CEA Excel wrapper"
#Else
    version = ""
    status = -1
    MsgBox "CEA Excel wrapper support is limited to 64-bit Windows VBA.", _
           vbExclamation, "CEA Excel wrapper"
#End If
End Function

Public Sub TestCEAVersion()
    Dim version As String
    Dim status As Long
    Dim messageText As String

    If Not CallCEAVersion(version, status) Then
        Exit Sub
    End If

    If status = 0 Then
        messageText = "CEA version: " & version
    Else
        messageText = "CEA Excel wrapper error " & CStr(status) & ": " & version
    End If

    MsgBox messageText, vbInformation, "CEA Excel wrapper"
End Sub

Public Sub TestCEAAdd()
    Dim result As Long

    If Not CallCEAAdd(2, 3, result) Then
        Exit Sub
    End If

    MsgBox "cea_excel_test_add(2, 3) returned " & CStr(result) & ".", _
           vbInformation, "CEA Excel wrapper"
End Sub

Public Sub WriteCEAVersionToReadme()
    Dim addResult As Long
    Dim version As String
    Dim status As Long
    Dim readmeSheet As Worksheet

    If Not CallCEAAdd(2, 3, addResult) Then
        Exit Sub
    End If

    If Not CallCEAVersion(version, status) Then
        Exit Sub
    End If

    On Error Resume Next
    Set readmeSheet = ThisWorkbook.Worksheets(README_SHEET_NAME)
    On Error GoTo 0

    If readmeSheet Is Nothing Then
        MsgBox "Workbook is missing the README sheet required for the smoke " & _
               "test.", vbExclamation, "CEA Excel wrapper"
        Exit Sub
    End If

    readmeSheet.Range(VERSION_LABEL_CELL).Value = "CEA Excel wrapper smoke test"
    readmeSheet.Range("A2").Value = "Detected CEA version:"
    readmeSheet.Range("A3").Value = "Status:"
    readmeSheet.Range("A5").Value = "cea_excel_test_add(2, 3):"
    readmeSheet.Range(VERSION_VALUE_CELL).Value = version
    readmeSheet.Range(STATUS_VALUE_CELL).Value = status
    readmeSheet.Range(ADD_STATUS_CELL).Value = addResult

    If status = 0 Then
        MsgBox "README sheet updated with CEA version " & version & ".", _
               vbInformation, "CEA Excel wrapper"
    Else
        MsgBox "README sheet updated, but the native wrapper returned status " & _
               CStr(status) & ".", vbExclamation, "CEA Excel wrapper"
    End If
End Sub
