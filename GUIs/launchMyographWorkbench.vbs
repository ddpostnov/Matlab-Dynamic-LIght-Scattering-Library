' launchMyographWorkbench.vbs
' Console-free launcher for the Myograph workbench: runs the .bat hidden so no
' command window flashes, then MATLAB opens the GUI. Double-click this file, or
' point a desktop shortcut at it. Keep it next to launchMyographWorkbench.bat.
Dim sh, fso, scriptDir, bat
Set sh  = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
scriptDir = fso.GetParentFolderName(WScript.ScriptFullName)
bat = scriptDir & "\launchMyographWorkbench.bat"
If Not fso.FileExists(bat) Then
    MsgBox "Cannot find launchMyographWorkbench.bat next to this file.", 48, "Myograph workbench"
    WScript.Quit 1
End If
' 0 = hidden window, False = do not wait for it to finish
sh.Run """" & bat & """", 0, False
