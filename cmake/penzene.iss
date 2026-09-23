; Inno Setup script: iscc /DVersion=x.y.z /DSource=dist\penzene cmake\penzene.iss
[Setup]
AppName=Penzene
AppVersion={#Version}
AppPublisher=Penzene contributors
AppPublisherURL=https://github.com/JamesOBrien2/penzene
DefaultDirName={autopf}\Penzene
DefaultGroupName=Penzene
UninstallDisplayIcon={app}\penzene.exe
OutputDir=.
OutputBaseFilename=penzene-windows-x64-setup
LicenseFile=..\LICENSE
Compression=lzma2
SolidCompression=yes
ArchitecturesInstallIn64BitMode=x64compatible
ArchitecturesAllowed=x64compatible
PrivilegesRequiredOverridesAllowed=dialog
ChangesAssociations=yes

[Files]
Source: "{#Source}\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Icons]
Name: "{group}\Penzene"; Filename: "{app}\penzene.exe"
Name: "{autodesktop}\Penzene"; Filename: "{app}\penzene.exe"; Tasks: desktopicon

[Tasks]
Name: desktopicon; Description: "Create a desktop shortcut"; Flags: unchecked

[Registry]
Root: HKA; Subkey: "Software\Classes\.penz"; ValueType: string; ValueData: "Penzene.Document"; Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\Penzene.Document"; ValueType: string; ValueData: "Penzene document"; Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\Penzene.Document\DefaultIcon"; ValueType: string; ValueData: "{app}\penzene.exe,0"
Root: HKA; Subkey: "Software\Classes\Penzene.Document\shell\open\command"; ValueType: string; ValueData: """{app}\penzene.exe"" ""%1"""

[Run]
Filename: "{app}\penzene.exe"; Description: "Launch Penzene"; Flags: nowait postinstall skipifsilent
