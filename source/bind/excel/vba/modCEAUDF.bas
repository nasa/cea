Attribute VB_Name = "modCEAUDF"
Option Explicit

Private Const CEA_TP As Long = 0
Private Const CEA_HP As Long = 1
Private Const CEA_SP As Long = 2
Private Const CEA_TV As Long = 3
Private Const CEA_UV As Long = 4
Private Const CEA_SV As Long = 5

Private Const CEA_EXCEL_SUCCESS As Long = 0
Private Const CEA_EXCEL_AMOUNT_WEIGHTS As Long = 0
Private Const CEA_EXCEL_AMOUNT_OF As Long = 1
Private Const CEA_EXCEL_AMOUNT_PHI As Long = 2
Private Const CEA_EXCEL_AMOUNT_R_EQ As Long = 3
Private Const CEA_EXCEL_AMOUNT_PCT_FUEL As Long = 4
Private Const CEA_EXCEL_ROLE_FUEL As Long = 1
Private Const CEA_EXCEL_ROLE_OXIDIZER As Long = 2
Private Const CEA_EXCEL_ROLE_INERT As Long = 3
Private Const CEA_EXCEL_BASIS_WEIGHT As Long = 1
Private Const CEA_EXCEL_BASIS_MOLE As Long = 2
Private Const CEA_EXCEL_BUFFER_LEN As Long = 4096
Private Const CEA_EXCEL_HEADER_LEN As Long = 16384
Private Const CEA_ENTHALPY As Long = 6
Private Const CEA_ENERGY As Long = 7
Private Const CEA_ENTROPY As Long = 8
Private Const CEA_GIBBS_ENERGY As Long = 9
Private Const CEA_FROZEN_CP As Long = 11
Private Const CEA_FROZEN_CV As Long = 12

Private Type CEAReactants
    NamesText As String
    BaseAmounts() As Double
    Roles() As Long
    Bases() As Long
    Count As Long
End Type

Private Function CEATrimNull(ByVal text As String) As String
    Dim nulPos As Long

    nulPos = InStr(1, text, vbNullChar)
    If nulPos > 0 Then
        CEATrimNull = Left$(text, nulPos - 1)
    Else
        CEATrimNull = text
    End If
End Function

Private Function CEAIsProvided(ByVal value As Variant) As Boolean
    If IsMissing(value) Then
        CEAIsProvided = False
    ElseIf IsEmpty(value) Then
        CEAIsProvided = False
    ElseIf IsObject(value) Then
        CEAIsProvided = Not value Is Nothing
    Else
        CEAIsProvided = True
    End If
End Function

Private Function CEABool(ByVal value As Variant, Optional ByVal defaultValue As Boolean = False) As Long
    If CEAIsProvided(value) Then
        CEABool = IIf(CBool(value), 1, 0)
    Else
        CEABool = IIf(defaultValue, 1, 0)
    End If
End Function

Private Function CEADouble(ByVal value As Variant, Optional ByVal defaultValue As Double = 0#) As Double
    If CEAIsProvided(value) Then
        If TypeName(value) = "Range" Then
            CEADouble = CDbl(value.Cells(1, 1).Value2)
        Else
            CEADouble = CDbl(value)
        End If
    Else
        CEADouble = defaultValue
    End If
End Function

Private Function CEAPackStrings(ByVal value As Variant) As String
    Dim arr As Variant
    Dim r As Long
    Dim c As Long
    Dim text As String

    If Not CEAIsProvided(value) Then
        CEAPackStrings = vbNullString
        Exit Function
    End If

    If TypeName(value) = "Range" Then
        arr = value.Value2
    Else
        arr = value
    End If

    If IsArray(arr) Then
        For r = LBound(arr, 1) To UBound(arr, 1)
            For c = LBound(arr, 2) To UBound(arr, 2)
                If Len(CStr(arr(r, c))) > 0 Then
                    If Len(text) > 0 Then text = text & vbLf
                    text = text & CStr(arr(r, c))
                End If
            Next c
        Next r
    ElseIf Len(CStr(arr)) > 0 Then
        text = CStr(arr)
    End If

    CEAPackStrings = text
End Function

Private Function CEAReadVector(ByVal value As Variant, ByRef out() As Double) As Long
    Dim arr As Variant
    Dim r As Long
    Dim c As Long
    Dim n As Long

    ReDim out(0 To 0)
    If Not CEAIsProvided(value) Then
        CEAReadVector = 0
        Exit Function
    End If

    If TypeName(value) = "Range" Then
        arr = value.Value2
    Else
        arr = value
    End If

    If IsArray(arr) Then
        ReDim out(0 To (UBound(arr, 1) - LBound(arr, 1) + 1) * _
                      (UBound(arr, 2) - LBound(arr, 2) + 1) - 1)
        For r = LBound(arr, 1) To UBound(arr, 1)
            For c = LBound(arr, 2) To UBound(arr, 2)
                out(n) = CDbl(arr(r, c))
                n = n + 1
            Next c
        Next r
    Else
        ReDim out(0 To 0)
        out(0) = CDbl(arr)
        n = 1
    End If

    CEAReadVector = n
End Function

Private Function CEAHeaderIndex(ByVal headers As Variant, ByVal name As String) As Long
    Dim c As Long
    Dim key As String

    For c = LBound(headers, 2) To UBound(headers, 2)
        key = LCase$(Trim$(CStr(headers(LBound(headers, 1), c))))
        If key = LCase$(name) Then
            CEAHeaderIndex = c
            Exit Function
        End If
    Next c
    CEAHeaderIndex = 0
End Function

Private Function CEARoleCode(ByVal text As String) As Long
    text = LCase$(Trim$(text))
    Select Case text
        Case "fuel", "fu"
            CEARoleCode = CEA_EXCEL_ROLE_FUEL
        Case "oxidizer", "oxidiser", "oxidant", "ox"
            CEARoleCode = CEA_EXCEL_ROLE_OXIDIZER
        Case "inert"
            CEARoleCode = CEA_EXCEL_ROLE_INERT
        Case Else
            Err.Raise vbObjectError + 1001, , "Unknown reactant role: " & text
    End Select
End Function

Private Function CEABasisCode(ByVal text As String) As Long
    text = LCase$(Trim$(text))
    Select Case text
        Case "weight", "weights", "mass", "wt"
            CEABasisCode = CEA_EXCEL_BASIS_WEIGHT
        Case "mole", "moles", "mol"
            CEABasisCode = CEA_EXCEL_BASIS_MOLE
        Case Else
            Err.Raise vbObjectError + 1002, , "Unknown reactant basis: " & text
    End Select
End Function

Private Function CEAReadReactants(ByVal tableRange As Range) As CEAReactants
    Dim data As Variant
    Dim nameCol As Long
    Dim roleCol As Long
    Dim amountCol As Long
    Dim basisCol As Long
    Dim r As Long
    Dim n As Long
    Dim result As CEAReactants

    data = tableRange.Value2
    nameCol = CEAHeaderIndex(data, "name")
    roleCol = CEAHeaderIndex(data, "role")
    amountCol = CEAHeaderIndex(data, "base_amount")
    basisCol = CEAHeaderIndex(data, "basis")
    If nameCol = 0 Or roleCol = 0 Or amountCol = 0 Or basisCol = 0 Then
        Err.Raise vbObjectError + 1003, , _
            "Reactant table must include name, role, base_amount, and basis headers."
    End If

    result.Count = UBound(data, 1) - LBound(data, 1)
    If result.Count <= 0 Then
        Err.Raise vbObjectError + 1004, , "Reactant table has no data rows."
    End If
    ReDim result.BaseAmounts(0 To result.Count - 1)
    ReDim result.Roles(0 To result.Count - 1)
    ReDim result.Bases(0 To result.Count - 1)

    For r = LBound(data, 1) + 1 To UBound(data, 1)
        If Len(result.NamesText) > 0 Then result.NamesText = result.NamesText & vbLf
        result.NamesText = result.NamesText & CStr(data(r, nameCol))
        result.Roles(n) = CEARoleCode(CStr(data(r, roleCol)))
        result.BaseAmounts(n) = CDbl(data(r, amountCol))
        result.Bases(n) = CEABasisCode(CStr(data(r, basisCol)))
        n = n + 1
    Next r

    CEAReadReactants = result
End Function

Private Function CEAResolveAmounts( _
    ByVal ofRatio As Variant, ByVal phi As Variant, ByVal rEq As Variant, _
    ByVal pctFuel As Variant, ByVal weightsArg As Variant, _
    ByRef explicitWeights() As Double, ByRef explicitWeightsLen As Long, _
    ByRef amountMode As Long, ByRef amountValue As Double) As Boolean

    Dim countModes As Long

    explicitWeightsLen = CEAReadVector(weightsArg, explicitWeights)
    If explicitWeightsLen > 0 Then countModes = countModes + 1
    If CEAIsProvided(ofRatio) Then countModes = countModes + 1
    If CEAIsProvided(phi) Then countModes = countModes + 1
    If CEAIsProvided(rEq) Then countModes = countModes + 1
    If CEAIsProvided(pctFuel) Then countModes = countModes + 1

    If countModes <> 1 Then
        Err.Raise vbObjectError + 1005, , _
            "Provide exactly one amount mode: weights, of_ratio, phi, r_eq, or pct_fuel."
    End If

    If explicitWeightsLen > 0 Then
        amountMode = CEA_EXCEL_AMOUNT_WEIGHTS
        amountValue = 0#
    ElseIf CEAIsProvided(ofRatio) Then
        amountMode = CEA_EXCEL_AMOUNT_OF
        amountValue = CEADouble(ofRatio)
    ElseIf CEAIsProvided(phi) Then
        amountMode = CEA_EXCEL_AMOUNT_PHI
        amountValue = CEADouble(phi)
    ElseIf CEAIsProvided(rEq) Then
        amountMode = CEA_EXCEL_AMOUNT_R_EQ
        amountValue = CEADouble(rEq)
    Else
        amountMode = CEA_EXCEL_AMOUNT_PCT_FUEL
        amountValue = CEADouble(pctFuel)
    End If

    CEAResolveAmounts = True
End Function

Private Function CEABuildRow( _
    ByVal status As Long, ByVal converged As Long, ByVal messageText As String, _
    ByRef values() As Double, ByVal nValues As Long, ByVal headersText As String, _
    ByVal includeHeaders As Boolean) As Variant

    Dim output As Variant
    Dim headers As Variant
    Dim c As Long

    headers = Split(CEATrimNull(headersText), vbTab)
    If includeHeaders Then
        ReDim output(1 To 2, 1 To nValues + 3)
        output(1, 1) = "status"
        output(1, 2) = "converged"
        output(1, 3) = "message"
        output(2, 1) = status
        output(2, 2) = converged <> 0
        output(2, 3) = messageText
        For c = 0 To nValues - 1
            If c <= UBound(headers) Then output(1, c + 4) = headers(c)
            output(2, c + 4) = values(c)
        Next c
    Else
        ReDim output(1 To 1, 1 To nValues + 3)
        output(1, 1) = status
        output(1, 2) = converged <> 0
        output(1, 3) = messageText
        For c = 0 To nValues - 1
            output(1, c + 4) = values(c)
        Next c
    End If

    CEABuildRow = output
End Function

Private Function CEAErrorRow(ByVal messageText As String) As Variant
    Dim output(1 To 1, 1 To 3) As Variant
    output(1, 1) = -1
    output(1, 2) = False
    output(1, 3) = messageText
    CEAErrorRow = output
End Function

Private Function CEAEqSolveImpl( _
    ByVal eqType As Long, ByVal reactants As Range, ByVal state1 As Double, _
    ByVal state2 As Double, ByVal ofRatio As Variant, ByVal phi As Variant, _
    ByVal rEq As Variant, ByVal pctFuel As Variant, ByVal weightsArg As Variant, _
    ByVal TReac As Variant, ByVal properties As Variant, ByVal species As Variant, _
    ByVal transport As Variant, ByVal ions As Variant, ByVal includeHeaders As Variant) As Variant

#If Win64 Then
    Dim loadError As String
    Dim reac As CEAReactants
    Dim explicitWeights() As Double
    Dim explicitLen As Long
    Dim amountMode As Long
    Dim amountValue As Double
    Dim temperatures() As Double
    Dim nTemps As Long
    Dim values(0 To CEA_EXCEL_BUFFER_LEN - 1) As Double
    Dim headers As String
    Dim msg As String
    Dim nValues As Long
    Dim converged As Long
    Dim status As Long

    On Error GoTo HandleError
    If Not EnsureCEALibraryLoaded(loadError) Then
        CEAEqSolveImpl = CEAErrorRow(loadError)
        Exit Function
    End If
    reac = CEAReadReactants(reactants)
    CEAResolveAmounts ofRatio, phi, rEq, pctFuel, weightsArg, _
        explicitWeights, explicitLen, amountMode, amountValue
    nTemps = CEAReadVector(TReac, temperatures)
    headers = String$(CEA_EXCEL_HEADER_LEN, vbNullChar)
    msg = String$(1024, vbNullChar)
    status = cea_excel_eq_solve(eqType, reac.NamesText, reac.BaseAmounts(0), _
        reac.Roles(0), reac.Bases(0), reac.Count, explicitWeights(0), _
        explicitLen, amountMode, amountValue, temperatures(0), nTemps, _
        state1, state2, vbNullString, vbNullString, vbNullString, _
        CEAPackStrings(properties), CEAPackStrings(species), CEABool(transport), _
        CEABool(ions), -1#, values(0), CEA_EXCEL_BUFFER_LEN, headers, _
        Len(headers), nValues, converged, msg, Len(msg))
    CEAEqSolveImpl = CEABuildRow(status, converged, CEATrimNull(msg), _
        values, nValues, headers, CEABool(includeHeaders) <> 0)
    Exit Function

HandleError:
    CEAEqSolveImpl = CEAErrorRow(Err.Description)
#Else
    CEAEqSolveImpl = CEAErrorRow("CEA Excel wrapper support is limited to 64-bit Windows VBA.")
#End If
End Function

Public Function CEA_TP_SOLVE(ByVal reactants As Range, ByVal T As Double, ByVal P As Double, _
    Optional ByVal of_ratio As Variant, Optional ByVal phi As Variant, _
    Optional ByVal r_eq As Variant, Optional ByVal pct_fuel As Variant, _
    Optional ByVal weights As Variant, Optional ByVal T_reac As Variant, _
    Optional ByVal properties As Variant, Optional ByVal species As Variant, _
    Optional ByVal transport As Variant, Optional ByVal ions As Variant, _
    Optional ByVal include_headers As Variant) As Variant
    CEA_TP_SOLVE = CEAEqSolveImpl(CEA_TP, reactants, T, P, of_ratio, phi, _
        r_eq, pct_fuel, weights, T_reac, properties, species, transport, ions, include_headers)
End Function

Public Function CEA_HP_SOLVE(ByVal reactants As Range, ByVal H As Double, ByVal P As Double, _
    Optional ByVal of_ratio As Variant, Optional ByVal phi As Variant, _
    Optional ByVal r_eq As Variant, Optional ByVal pct_fuel As Variant, _
    Optional ByVal weights As Variant, Optional ByVal T_reac As Variant, _
    Optional ByVal properties As Variant, Optional ByVal species As Variant, _
    Optional ByVal transport As Variant, Optional ByVal ions As Variant, _
    Optional ByVal include_headers As Variant) As Variant
    CEA_HP_SOLVE = CEAEqSolveImpl(CEA_HP, reactants, H, P, of_ratio, phi, _
        r_eq, pct_fuel, weights, T_reac, properties, species, transport, ions, include_headers)
End Function

Public Function CEA_SP_SOLVE(ByVal reactants As Range, ByVal S As Double, ByVal P As Double, _
    Optional ByVal of_ratio As Variant, Optional ByVal phi As Variant, _
    Optional ByVal r_eq As Variant, Optional ByVal pct_fuel As Variant, _
    Optional ByVal weights As Variant, Optional ByVal T_reac As Variant, _
    Optional ByVal properties As Variant, Optional ByVal species As Variant, _
    Optional ByVal transport As Variant, Optional ByVal ions As Variant, _
    Optional ByVal include_headers As Variant) As Variant
    CEA_SP_SOLVE = CEAEqSolveImpl(CEA_SP, reactants, S, P, of_ratio, phi, _
        r_eq, pct_fuel, weights, T_reac, properties, species, transport, ions, include_headers)
End Function

Public Function CEA_TV_SOLVE(ByVal reactants As Range, ByVal T As Double, ByVal V As Double, _
    Optional ByVal of_ratio As Variant, Optional ByVal phi As Variant, _
    Optional ByVal r_eq As Variant, Optional ByVal pct_fuel As Variant, _
    Optional ByVal weights As Variant, Optional ByVal T_reac As Variant, _
    Optional ByVal properties As Variant, Optional ByVal species As Variant, _
    Optional ByVal transport As Variant, Optional ByVal ions As Variant, _
    Optional ByVal include_headers As Variant) As Variant
    CEA_TV_SOLVE = CEAEqSolveImpl(CEA_TV, reactants, T, V, of_ratio, phi, _
        r_eq, pct_fuel, weights, T_reac, properties, species, transport, ions, include_headers)
End Function

Public Function CEA_UV_SOLVE(ByVal reactants As Range, ByVal U As Double, ByVal V As Double, _
    Optional ByVal of_ratio As Variant, Optional ByVal phi As Variant, _
    Optional ByVal r_eq As Variant, Optional ByVal pct_fuel As Variant, _
    Optional ByVal weights As Variant, Optional ByVal T_reac As Variant, _
    Optional ByVal properties As Variant, Optional ByVal species As Variant, _
    Optional ByVal transport As Variant, Optional ByVal ions As Variant, _
    Optional ByVal include_headers As Variant) As Variant
    CEA_UV_SOLVE = CEAEqSolveImpl(CEA_UV, reactants, U, V, of_ratio, phi, _
        r_eq, pct_fuel, weights, T_reac, properties, species, transport, ions, include_headers)
End Function

Public Function CEA_SV_SOLVE(ByVal reactants As Range, ByVal S As Double, ByVal V As Double, _
    Optional ByVal of_ratio As Variant, Optional ByVal phi As Variant, _
    Optional ByVal r_eq As Variant, Optional ByVal pct_fuel As Variant, _
    Optional ByVal weights As Variant, Optional ByVal T_reac As Variant, _
    Optional ByVal properties As Variant, Optional ByVal species As Variant, _
    Optional ByVal transport As Variant, Optional ByVal ions As Variant, _
    Optional ByVal include_headers As Variant) As Variant
    CEA_SV_SOLVE = CEAEqSolveImpl(CEA_SV, reactants, S, V, of_ratio, phi, _
        r_eq, pct_fuel, weights, T_reac, properties, species, transport, ions, include_headers)
End Function

Private Function CEARocketSolveImpl( _
    ByVal finiteArea As Long, ByVal reactants As Range, ByVal pc As Double, _
    ByVal piPArg As Variant, ByVal subarArg As Variant, ByVal suparArg As Variant, _
    ByVal ofRatio As Variant, ByVal phi As Variant, ByVal rEq As Variant, _
    ByVal pctFuel As Variant, ByVal weightsArg As Variant, ByVal TReac As Variant, _
    ByVal nFrzArg As Variant, ByVal hcArg As Variant, ByVal tcArg As Variant, _
    ByVal mdotArg As Variant, ByVal acAtArg As Variant, ByVal tcEstArg As Variant, _
    ByVal omitArg As Variant, ByVal insertArg As Variant, ByVal properties As Variant, _
    ByVal species As Variant, ByVal transport As Variant, ByVal ions As Variant, _
    ByVal includeHeaders As Variant) As Variant

#If Win64 Then
    Dim loadError As String
    Dim reac As CEAReactants
    Dim explicitWeights() As Double
    Dim explicitLen As Long
    Dim amountMode As Long
    Dim amountValue As Double
    Dim temperatures() As Double
    Dim nTemps As Long
    Dim piP() As Double
    Dim nPiP As Long
    Dim subar() As Double
    Dim nSubar As Long
    Dim supar() As Double
    Dim nSupar As Long
    Dim values(0 To CEA_EXCEL_BUFFER_LEN - 1) As Double
    Dim headers As String
    Dim msg As String
    Dim nValues As Long
    Dim converged As Long
    Dim status As Long
    Dim hcOrTc As Double
    Dim useHc As Long
    Dim mdotOrAcat As Double
    Dim useMdot As Long

    On Error GoTo HandleError
    If Not EnsureCEALibraryLoaded(loadError) Then
        CEARocketSolveImpl = CEAErrorRow(loadError)
        Exit Function
    End If
    reac = CEAReadReactants(reactants)
    CEAResolveAmounts ofRatio, phi, rEq, pctFuel, weightsArg, _
        explicitWeights, explicitLen, amountMode, amountValue
    nTemps = CEAReadVector(TReac, temperatures)
    nPiP = CEAReadVector(piPArg, piP)
    nSubar = CEAReadVector(subarArg, subar)
    nSupar = CEAReadVector(suparArg, supar)
    If CEAIsProvided(hcArg) Then
        hcOrTc = CEADouble(hcArg)
        useHc = 1
    ElseIf CEAIsProvided(tcArg) Then
        hcOrTc = CEADouble(tcArg)
        useHc = 0
    Else
        hcOrTc = 0#
        useHc = 0
    End If
    If CEAIsProvided(mdotArg) Then
        mdotOrAcat = CEADouble(mdotArg)
        useMdot = 1
    Else
        mdotOrAcat = CEADouble(acAtArg)
        useMdot = 0
    End If
    headers = String$(CEA_EXCEL_HEADER_LEN, vbNullChar)
    msg = String$(1024, vbNullChar)
    status = cea_excel_rocket_solve(finiteArea, reac.NamesText, reac.BaseAmounts(0), _
        reac.Roles(0), reac.Bases(0), reac.Count, explicitWeights(0), explicitLen, _
        amountMode, amountValue, temperatures(0), nTemps, pc, piP(0), nPiP, _
        subar(0), nSubar, supar(0), nSupar, CLng(CEADouble(nFrzArg, 0#)), _
        hcOrTc, useHc, mdotOrAcat, useMdot, CEADouble(tcEstArg, 0#), _
        CEABool(tcEstArg), CEAPackStrings(omitArg), CEAPackStrings(insertArg), _
        CEAPackStrings(properties), CEAPackStrings(species), CEABool(transport), _
        CEABool(ions), -1#, values(0), CEA_EXCEL_BUFFER_LEN, headers, _
        Len(headers), nValues, converged, msg, Len(msg))
    CEARocketSolveImpl = CEABuildRow(status, converged, CEATrimNull(msg), _
        values, nValues, headers, CEABool(includeHeaders) <> 0)
    Exit Function

HandleError:
    CEARocketSolveImpl = CEAErrorRow(Err.Description)
#Else
    CEARocketSolveImpl = CEAErrorRow("CEA Excel wrapper support is limited to 64-bit Windows VBA.")
#End If
End Function

Public Function CEA_ROCKET_IAC_SOLVE(ByVal reactants As Range, ByVal pc As Double, _
    Optional ByVal pi_p As Variant, Optional ByVal subar As Variant, _
    Optional ByVal supar As Variant, Optional ByVal of_ratio As Variant, _
    Optional ByVal phi As Variant, Optional ByVal r_eq As Variant, _
    Optional ByVal pct_fuel As Variant, Optional ByVal weights As Variant, _
    Optional ByVal T_reac As Variant, Optional ByVal n_frz As Variant, _
    Optional ByVal hc As Variant, Optional ByVal tc As Variant, _
    Optional ByVal tc_est As Variant, Optional ByVal omit As Variant, _
    Optional ByVal insert As Variant, Optional ByVal properties As Variant, _
    Optional ByVal species As Variant, Optional ByVal transport As Variant, _
    Optional ByVal ions As Variant, Optional ByVal include_headers As Variant) As Variant
    CEA_ROCKET_IAC_SOLVE = CEARocketSolveImpl(0, reactants, pc, pi_p, subar, _
        supar, of_ratio, phi, r_eq, pct_fuel, weights, T_reac, n_frz, hc, tc, _
        Empty, Empty, tc_est, omit, insert, properties, species, transport, ions, include_headers)
End Function

Public Function CEA_ROCKET_FAC_SOLVE(ByVal reactants As Range, ByVal pc As Double, _
    Optional ByVal pi_p As Variant, Optional ByVal subar As Variant, _
    Optional ByVal supar As Variant, Optional ByVal of_ratio As Variant, _
    Optional ByVal phi As Variant, Optional ByVal r_eq As Variant, _
    Optional ByVal pct_fuel As Variant, Optional ByVal weights As Variant, _
    Optional ByVal T_reac As Variant, Optional ByVal n_frz As Variant, _
    Optional ByVal hc As Variant, Optional ByVal tc As Variant, _
    Optional ByVal mdot As Variant, Optional ByVal ac_at As Variant, _
    Optional ByVal tc_est As Variant, Optional ByVal omit As Variant, _
    Optional ByVal insert As Variant, Optional ByVal properties As Variant, _
    Optional ByVal species As Variant, Optional ByVal transport As Variant, _
    Optional ByVal ions As Variant, Optional ByVal include_headers As Variant) As Variant
    CEA_ROCKET_FAC_SOLVE = CEARocketSolveImpl(1, reactants, pc, pi_p, subar, _
        supar, of_ratio, phi, r_eq, pct_fuel, weights, T_reac, n_frz, hc, tc, _
        mdot, ac_at, tc_est, omit, insert, properties, species, transport, ions, include_headers)
End Function

Public Function CEA_SHOCK_SOLVE(ByVal reactants As Range, ByVal T0 As Double, _
    ByVal p0 As Double, Optional ByVal u1 As Variant, Optional ByVal Mach1 As Variant, _
    Optional ByVal of_ratio As Variant, Optional ByVal phi As Variant, _
    Optional ByVal r_eq As Variant, Optional ByVal pct_fuel As Variant, _
    Optional ByVal weights As Variant, Optional ByVal reflected As Variant, _
    Optional ByVal incident_frozen As Variant, Optional ByVal reflected_frozen As Variant, _
    Optional ByVal omit As Variant, Optional ByVal insert As Variant, _
    Optional ByVal properties As Variant, Optional ByVal species As Variant, _
    Optional ByVal transport As Variant, Optional ByVal ions As Variant, _
    Optional ByVal include_headers As Variant) As Variant

#If Win64 Then
    Dim loadError As String
    Dim reac As CEAReactants
    Dim explicitWeights() As Double
    Dim explicitLen As Long
    Dim amountMode As Long
    Dim amountValue As Double
    Dim values(0 To CEA_EXCEL_BUFFER_LEN - 1) As Double
    Dim headers As String
    Dim msg As String
    Dim nValues As Long
    Dim converged As Long
    Dim status As Long
    Dim useMach As Long
    Dim speedArg As Double

    On Error GoTo HandleError
    If Not EnsureCEALibraryLoaded(loadError) Then
        CEA_SHOCK_SOLVE = CEAErrorRow(loadError)
        Exit Function
    End If
    If CEAIsProvided(Mach1) Then
        useMach = 1
        speedArg = CEADouble(Mach1)
    Else
        useMach = 0
        speedArg = CEADouble(u1)
    End If
    reac = CEAReadReactants(reactants)
    CEAResolveAmounts of_ratio, phi, r_eq, pct_fuel, weights, _
        explicitWeights, explicitLen, amountMode, amountValue
    headers = String$(CEA_EXCEL_HEADER_LEN, vbNullChar)
    msg = String$(1024, vbNullChar)
    status = cea_excel_shock_solve(reac.NamesText, reac.BaseAmounts(0), _
        reac.Roles(0), reac.Bases(0), reac.Count, explicitWeights(0), explicitLen, _
        amountMode, amountValue, T0, p0, speedArg, useMach, CEABool(reflected, True), _
        CEABool(incident_frozen), CEABool(reflected_frozen), CEAPackStrings(omit), _
        CEAPackStrings(insert), CEAPackStrings(properties), CEAPackStrings(species), _
        CEABool(transport), CEABool(ions), -1#, values(0), CEA_EXCEL_BUFFER_LEN, _
        headers, Len(headers), nValues, converged, msg, Len(msg))
    CEA_SHOCK_SOLVE = CEABuildRow(status, converged, CEATrimNull(msg), _
        values, nValues, headers, CEABool(include_headers) <> 0)
    Exit Function

HandleError:
    CEA_SHOCK_SOLVE = CEAErrorRow(Err.Description)
#Else
    CEA_SHOCK_SOLVE = CEAErrorRow("CEA Excel wrapper support is limited to 64-bit Windows VBA.")
#End If
End Function

Public Function CEA_DETONATION_SOLVE(ByVal reactants As Range, ByVal T1 As Double, _
    ByVal p1 As Double, Optional ByVal of_ratio As Variant, Optional ByVal phi As Variant, _
    Optional ByVal r_eq As Variant, Optional ByVal pct_fuel As Variant, _
    Optional ByVal weights As Variant, Optional ByVal frozen As Variant, _
    Optional ByVal omit As Variant, Optional ByVal insert As Variant, _
    Optional ByVal properties As Variant, Optional ByVal species As Variant, _
    Optional ByVal transport As Variant, Optional ByVal ions As Variant, _
    Optional ByVal include_headers As Variant) As Variant

#If Win64 Then
    Dim loadError As String
    Dim reac As CEAReactants
    Dim explicitWeights() As Double
    Dim explicitLen As Long
    Dim amountMode As Long
    Dim amountValue As Double
    Dim values(0 To CEA_EXCEL_BUFFER_LEN - 1) As Double
    Dim headers As String
    Dim msg As String
    Dim nValues As Long
    Dim converged As Long
    Dim status As Long

    On Error GoTo HandleError
    If Not EnsureCEALibraryLoaded(loadError) Then
        CEA_DETONATION_SOLVE = CEAErrorRow(loadError)
        Exit Function
    End If
    reac = CEAReadReactants(reactants)
    CEAResolveAmounts of_ratio, phi, r_eq, pct_fuel, weights, _
        explicitWeights, explicitLen, amountMode, amountValue
    headers = String$(CEA_EXCEL_HEADER_LEN, vbNullChar)
    msg = String$(1024, vbNullChar)
    status = cea_excel_detonation_solve(reac.NamesText, reac.BaseAmounts(0), _
        reac.Roles(0), reac.Bases(0), reac.Count, explicitWeights(0), explicitLen, _
        amountMode, amountValue, T1, p1, CEABool(frozen), CEAPackStrings(omit), _
        CEAPackStrings(insert), CEAPackStrings(properties), CEAPackStrings(species), _
        CEABool(transport), CEABool(ions), -1#, values(0), CEA_EXCEL_BUFFER_LEN, _
        headers, Len(headers), nValues, converged, msg, Len(msg))
    CEA_DETONATION_SOLVE = CEABuildRow(status, converged, CEATrimNull(msg), _
        values, nValues, headers, CEABool(include_headers) <> 0)
    Exit Function

HandleError:
    CEA_DETONATION_SOLVE = CEAErrorRow(Err.Description)
#Else
    CEA_DETONATION_SOLVE = CEAErrorRow("CEA Excel wrapper support is limited to 64-bit Windows VBA.")
#End If
End Function

Private Function CEAArrayRow(ByRef values() As Double, ByVal count As Long) As Variant
    Dim output As Variant
    Dim i As Long

    ReDim output(1 To 1, 1 To count)
    For i = 0 To count - 1
        output(1, i + 1) = values(i)
    Next i
    CEAArrayRow = output
End Function

Private Function CEAPropertyCode(ByVal propertyType As Variant) As Long
    Dim text As String

    If IsNumeric(propertyType) Then
        CEAPropertyCode = CLng(propertyType)
        Exit Function
    End If
    text = LCase$(Trim$(CStr(propertyType)))
    Select Case text
        Case "enthalpy", "h"
            CEAPropertyCode = CEA_ENTHALPY
        Case "energy", "internal_energy", "u"
            CEAPropertyCode = CEA_ENERGY
        Case "entropy", "s"
            CEAPropertyCode = CEA_ENTROPY
        Case "gibbs_energy", "gibbs", "g"
            CEAPropertyCode = CEA_GIBBS_ENERGY
        Case "cp", "cp_fr", "frozen_cp"
            CEAPropertyCode = CEA_FROZEN_CP
        Case "cv", "cv_fr", "frozen_cv"
            CEAPropertyCode = CEA_FROZEN_CV
        Case Else
            Err.Raise vbObjectError + 1006, , "Unknown thermo property: " & CStr(propertyType)
    End Select
End Function

Public Function CEA_VERSION() As Variant
#If Win64 Then
    Dim loadError As String
    Dim version As String
    Dim status As Long

    If Not EnsureCEALibraryLoaded(loadError) Then
        CEA_VERSION = CVErr(xlErrValue)
        Exit Function
    End If
    version = String$(256, vbNullChar)
    status = cea_excel_version(version, Len(version))
    If status = CEA_EXCEL_SUCCESS Then
        CEA_VERSION = CEATrimNull(version)
    Else
        CEA_VERSION = CVErr(xlErrValue)
    End If
#Else
    CEA_VERSION = CVErr(xlErrNA)
#End If
End Function

Public Function CEA_LAST_ERROR() As Variant
#If Win64 Then
    Dim loadError As String
    Dim msg As String

    If Not EnsureCEALibraryLoaded(loadError) Then
        CEA_LAST_ERROR = loadError
        Exit Function
    End If
    msg = String$(1024, vbNullChar)
    cea_excel_last_error msg, Len(msg)
    CEA_LAST_ERROR = CEATrimNull(msg)
#Else
    CEA_LAST_ERROR = "CEA Excel wrapper support is limited to 64-bit Windows VBA."
#End If
End Function

Public Function CEA_TO_SI(ByVal value As Double, ByVal units As String) As Variant
#If Win64 Then
    Dim loadError As String
    Dim msg As String
    Dim status As Long
    Dim siValue As Double

    If Not EnsureCEALibraryLoaded(loadError) Then
        CEA_TO_SI = CVErr(xlErrValue)
        Exit Function
    End If
    msg = String$(256, vbNullChar)
    status = cea_excel_to_si(value, units, siValue, msg, Len(msg))
    If status = CEA_EXCEL_SUCCESS Then
        CEA_TO_SI = siValue
    Else
        CEA_TO_SI = CVErr(xlErrValue)
    End If
#Else
    CEA_TO_SI = CVErr(xlErrNA)
#End If
End Function

Public Function CEA_PRESSURE_SI(ByVal value As Double, ByVal units As String) As Variant
    CEA_PRESSURE_SI = CEA_TO_SI(value, units)
End Function

Public Function CEA_TEMPERATURE_SI(ByVal value As Double, ByVal units As String) As Variant
    CEA_TEMPERATURE_SI = CEA_TO_SI(value, units)
End Function

Public Function CEA_ENERGY_SI(ByVal value As Double, ByVal units As String) As Variant
    CEA_ENERGY_SI = CEA_TO_SI(value, units)
End Function

Public Function CEA_DENSITY_SI(ByVal value As Double, ByVal units As String) As Variant
    CEA_DENSITY_SI = CEA_TO_SI(value, units)
End Function

Public Function CEA_VOLUME_SI(ByVal value As Double, ByVal units As String) As Variant
    CEA_VOLUME_SI = CEA_TO_SI(value, units)
End Function

Public Function CEA_WEIGHTS_FROM_OF(ByVal reactants As Range, ByVal of_ratio As Double) As Variant
#If Win64 Then
    Dim loadError As String
    Dim reac As CEAReactants
    Dim weights() As Double
    Dim msg As String
    Dim status As Long

    On Error GoTo HandleError
    If Not EnsureCEALibraryLoaded(loadError) Then
        CEA_WEIGHTS_FROM_OF = CEAErrorRow(loadError)
        Exit Function
    End If
    reac = CEAReadReactants(reactants)
    ReDim weights(0 To reac.Count - 1)
    msg = String$(1024, vbNullChar)
    status = cea_excel_weights_from_of(reac.NamesText, reac.BaseAmounts(0), _
        reac.Roles(0), reac.Bases(0), reac.Count, of_ratio, weights(0), msg, Len(msg))
    If status = CEA_EXCEL_SUCCESS Then
        CEA_WEIGHTS_FROM_OF = CEAArrayRow(weights, reac.Count)
    Else
        CEA_WEIGHTS_FROM_OF = CEAErrorRow(CEATrimNull(msg))
    End If
    Exit Function
HandleError:
    CEA_WEIGHTS_FROM_OF = CEAErrorRow(Err.Description)
#Else
    CEA_WEIGHTS_FROM_OF = CEAErrorRow("CEA Excel wrapper support is limited to 64-bit Windows VBA.")
#End If
End Function

Private Function CEAScalarRatioHelper( _
    ByVal reactants As Range, ByVal inputValue As Double, ByVal mode As String) As Variant
#If Win64 Then
    Dim loadError As String
    Dim reac As CEAReactants
    Dim msg As String
    Dim status As Long
    Dim result As Double

    On Error GoTo HandleError
    If Not EnsureCEALibraryLoaded(loadError) Then
        CEAScalarRatioHelper = CVErr(xlErrValue)
        Exit Function
    End If
    reac = CEAReadReactants(reactants)
    msg = String$(1024, vbNullChar)
    Select Case mode
        Case "of_from_equivalence"
            status = cea_excel_of_from_equivalence(reac.NamesText, reac.BaseAmounts(0), _
                reac.Roles(0), reac.Bases(0), reac.Count, inputValue, result, msg, Len(msg))
        Case "of_from_phi"
            status = cea_excel_of_from_phi(reac.NamesText, reac.BaseAmounts(0), _
                reac.Roles(0), reac.Bases(0), reac.Count, inputValue, result, msg, Len(msg))
        Case "equivalence_from_of"
            status = cea_excel_equivalence_from_of(reac.NamesText, reac.BaseAmounts(0), _
                reac.Roles(0), reac.Bases(0), reac.Count, inputValue, result, msg, Len(msg))
        Case "phi_from_of"
            status = cea_excel_phi_from_of(reac.NamesText, reac.BaseAmounts(0), _
                reac.Roles(0), reac.Bases(0), reac.Count, inputValue, result, msg, Len(msg))
    End Select
    If status = CEA_EXCEL_SUCCESS Then
        CEAScalarRatioHelper = result
    Else
        CEAScalarRatioHelper = CVErr(xlErrValue)
    End If
    Exit Function
HandleError:
    CEAScalarRatioHelper = CVErr(xlErrValue)
#Else
    CEAScalarRatioHelper = CVErr(xlErrNA)
#End If
End Function

Public Function CEA_OF_FROM_EQUIVALENCE(ByVal reactants As Range, ByVal equivalence As Double) As Variant
    CEA_OF_FROM_EQUIVALENCE = CEAScalarRatioHelper(reactants, equivalence, "of_from_equivalence")
End Function

Public Function CEA_OF_FROM_PHI(ByVal reactants As Range, ByVal phi As Double) As Variant
    CEA_OF_FROM_PHI = CEAScalarRatioHelper(reactants, phi, "of_from_phi")
End Function

Public Function CEA_EQUIVALENCE_FROM_OF(ByVal reactants As Range, ByVal of_ratio As Double) As Variant
    CEA_EQUIVALENCE_FROM_OF = CEAScalarRatioHelper(reactants, of_ratio, "equivalence_from_of")
End Function

Public Function CEA_PHI_FROM_OF(ByVal reactants As Range, ByVal of_ratio As Double) As Variant
    CEA_PHI_FROM_OF = CEAScalarRatioHelper(reactants, of_ratio, "phi_from_of")
End Function

Private Function CEAArrayConversion(ByVal reactants As Range, ByVal inputValues As Variant, _
    ByVal conversion As String) As Variant
#If Win64 Then
    Dim loadError As String
    Dim reac As CEAReactants
    Dim inputArray() As Double
    Dim nInput As Long
    Dim outputArray() As Double
    Dim msg As String
    Dim status As Long

    On Error GoTo HandleError
    If Not EnsureCEALibraryLoaded(loadError) Then
        CEAArrayConversion = CEAErrorRow(loadError)
        Exit Function
    End If
    reac = CEAReadReactants(reactants)
    nInput = CEAReadVector(inputValues, inputArray)
    If nInput <> reac.Count Then Err.Raise vbObjectError + 1007, , "Input length must match reactants."
    ReDim outputArray(0 To reac.Count - 1)
    msg = String$(1024, vbNullChar)
    Select Case conversion
        Case "weights_from_moles"
            status = cea_excel_weights_from_moles(reac.NamesText, inputArray(0), _
                reac.Count, outputArray(0), msg, Len(msg))
        Case "moles_from_weights"
            status = cea_excel_moles_from_weights(reac.NamesText, inputArray(0), _
                reac.Count, outputArray(0), msg, Len(msg))
        Case "per_mole_from_per_weight"
            status = cea_excel_per_mole_from_per_weight(reac.NamesText, inputArray(0), _
                reac.Count, outputArray(0), msg, Len(msg))
        Case "per_weight_from_per_mole"
            status = cea_excel_per_weight_from_per_mole(reac.NamesText, inputArray(0), _
                reac.Count, outputArray(0), msg, Len(msg))
    End Select
    If status = CEA_EXCEL_SUCCESS Then
        CEAArrayConversion = CEAArrayRow(outputArray, reac.Count)
    Else
        CEAArrayConversion = CEAErrorRow(CEATrimNull(msg))
    End If
    Exit Function
HandleError:
    CEAArrayConversion = CEAErrorRow(Err.Description)
#Else
    CEAArrayConversion = CEAErrorRow("CEA Excel wrapper support is limited to 64-bit Windows VBA.")
#End If
End Function

Public Function CEA_WEIGHTS_FROM_MOLES(ByVal reactants As Range, ByVal moles As Variant) As Variant
    CEA_WEIGHTS_FROM_MOLES = CEAArrayConversion(reactants, moles, "weights_from_moles")
End Function

Public Function CEA_MOLES_FROM_WEIGHTS(ByVal reactants As Range, ByVal weights As Variant) As Variant
    CEA_MOLES_FROM_WEIGHTS = CEAArrayConversion(reactants, weights, "moles_from_weights")
End Function

Public Function CEA_PER_MOLE_FROM_PER_WEIGHT(ByVal reactants As Range, ByVal per_weight As Variant) As Variant
    CEA_PER_MOLE_FROM_PER_WEIGHT = CEAArrayConversion(reactants, per_weight, "per_mole_from_per_weight")
End Function

Public Function CEA_PER_WEIGHT_FROM_PER_MOLE(ByVal reactants As Range, ByVal per_mole As Variant) As Variant
    CEA_PER_WEIGHT_FROM_PER_MOLE = CEAArrayConversion(reactants, per_mole, "per_weight_from_per_mole")
End Function

Public Function CEA_CALC_THERMO(ByVal reactants As Range, ByVal weights As Variant, _
    ByVal property_type As Variant, ByVal temperatures As Variant, _
    Optional ByVal pressure As Variant) As Variant
#If Win64 Then
    Dim loadError As String
    Dim reac As CEAReactants
    Dim weightArray() As Double
    Dim tempArray() As Double
    Dim nWeights As Long
    Dim nTemps As Long
    Dim msg As String
    Dim status As Long
    Dim result As Double

    On Error GoTo HandleError
    If Not EnsureCEALibraryLoaded(loadError) Then
        CEA_CALC_THERMO = CVErr(xlErrValue)
        Exit Function
    End If
    reac = CEAReadReactants(reactants)
    nWeights = CEAReadVector(weights, weightArray)
    nTemps = CEAReadVector(temperatures, tempArray)
    If nWeights <> reac.Count Then Err.Raise vbObjectError + 1008, , "Weights length must match reactants."
    If nTemps <> 1 And nTemps <> reac.Count Then _
        Err.Raise vbObjectError + 1009, , "Temperatures must be scalar or match reactants."
    msg = String$(1024, vbNullChar)
    status = cea_excel_calc_thermo(reac.NamesText, weightArray(0), reac.Count, _
        CEAPropertyCode(property_type), tempArray(0), nTemps, CEADouble(pressure, 0#), _
        CEABool(pressure), result, msg, Len(msg))
    If status = CEA_EXCEL_SUCCESS Then
        CEA_CALC_THERMO = result
    Else
        CEA_CALC_THERMO = CVErr(xlErrValue)
    End If
    Exit Function
HandleError:
    CEA_CALC_THERMO = CVErr(xlErrValue)
#Else
    CEA_CALC_THERMO = CVErr(xlErrNA)
#End If
End Function
