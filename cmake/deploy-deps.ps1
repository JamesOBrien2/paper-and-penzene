# Copies every DLL that penzene.exe and its Qt plugins load from the build env
# into the package. windeployqt handles Qt itself; this picks up RDKit, Boost, ICU...
# Only follows DLLs that exist in -Search, so system DLLs are never crawled.
param([string]$Dest, [string]$Search)
$queue = [System.Collections.Generic.Queue[string]]::new()
Get-ChildItem $Dest -Recurse -Include *.exe, *.dll | ForEach-Object { $queue.Enqueue($_.FullName) }
$seen = @{}
while ($queue.Count) {
    $file = $queue.Dequeue()
    foreach ($line in (dumpbin /nologo /dependents $file)) {
        $name = $line.Trim()
        if ($name -notmatch '\.dll$' -or $seen[$name]) { continue }
        $seen[$name] = $true
        $src = Join-Path $Search $name
        if ((Test-Path $src) -and -not (Test-Path (Join-Path $Dest $name))) {
            Copy-Item $src $Dest
            $queue.Enqueue((Join-Path $Dest $name))
        }
    }
}
Write-Host "Bundled $((Get-ChildItem $Dest -Filter *.dll).Count) DLLs"
