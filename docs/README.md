# pydock3 documentation

For users, the repository's [README](../README.md) is the starting point. These pages are for
maintainers and contributors: how the package works, and why it is built the way it is.

| Page | Read it to… |
|---|---|
| [architecture.md](architecture.md) | get the lay of the land: repository layout, packages, the command line, core abstractions |
| [blastermaster.md](blastermaster.md) | understand blastermaster: its job directory, pipeline steps, files, configuration, and how steps run |
| [native-programs.md](native-programs.md) | know what each bundled program (reduce, dms, sphgen, chemgrid, qnifft, solvmap, filt) does, where it comes from, and what was changed in it |
| [cross-platform.md](cross-platform.md) | understand how blastermaster was made to run on Linux, macOS and Windows, on x86 and ARM, and why results can differ between machines |
| [build-and-ci.md](build-and-ci.md) | build from source, understand the compiler flags, the wheels and the CI |
| [testing.md](testing.md) | run and read the tests, understand the control files and tolerances, update them |
| [other-tools.md](other-tools.md) | work on dockopt, retrodock, lsd and the job schedulers |
| [maintenance.md](maintenance.md) | make common changes (steps, programs, vendored code, CI), release, and see what's unfinished |

**Where to start.** New to the code: architecture → blastermaster → testing. Changing a program
or the build: native-programs → cross-platform → build-and-ci. Something differs on one
platform: cross-platform (floating-point section) → testing.
