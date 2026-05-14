Attribute VB_Name = "modCEADeclare"
Option Explicit

#If Win64 Then
    Public Declare PtrSafe Function cea_excel_version Lib "cea_excel.dll" ( _
        ByVal versionBuffer As String, _
        ByVal versionLen As Long) As Long
    Public Declare PtrSafe Function cea_excel_test_add Lib "cea_excel.dll" ( _
        ByVal a As Long, _
        ByVal b As Long) As Long
#Else
    ' The Excel wrapper currently supports 64-bit Windows VBA only.
#End If
