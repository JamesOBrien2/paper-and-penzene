# Installing

Download the installer for your system from the
[latest release](https://github.com/JamesOBrien2/penzene/releases/latest).

| System | File |
|---|---|
| macOS, Apple silicon | `penzene-macos-arm64.dmg` |
| macOS, Intel | `penzene-macos-x86_64.dmg` |
| Windows | `penzene-windows-x64-setup.exe` (or the `.zip`, no install needed) |
| Linux | `penzene-linux-x86_64.AppImage` |

## Opening it the first time

:::{note}
Penzene's installers aren't signed with a paid Apple or Microsoft certificate, so each system
asks once whether to trust it.
:::

**macOS.** Drag Penzene to Applications and open it. If macOS says it "can't be opened" or can't be
checked, go to **System Settings → Privacy & Security**, scroll down to the message about Penzene,
and click **Open Anyway**. You only need to do this once. From a terminal, this does the same:

```sh
xattr -dr com.apple.quarantine /Applications/Penzene.app
```

**Windows.** If SmartScreen shows "Windows protected your PC", click **More info → Run anyway**.

**Linux.** Make the AppImage executable, then run it:

```sh
chmod +x penzene-linux-x86_64.AppImage
./penzene-linux-x86_64.AppImage
```

## Python package

```sh
pip install penzene
```

Python 3.12+ on macOS, Linux (glibc 2.17+) and Windows. See [Python](python.md).

## Building from source

See [Contributing](contributing.md).
