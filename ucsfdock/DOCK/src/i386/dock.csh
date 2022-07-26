#!/bin/csh -f

# Changes to this file must be made in parallel with get_dock_files.csh
#
# Determine and execute the best version of dock for this architecture
#
# Currently just determines 64 or 32 bit.  

# In the future, we could use info from /proc/cpuinfo
#   to dispatch to an architecture specific version of dock

# determine base directory where this script lives
set base=$0:h
if ( "$base" == "$0" ) then
    set base="."
endif

set arch = `uname -p`

if ( $arch == 'x86_64') then
    # base preceded by $ is grepped by get_dock_files.csh
    $base/dock64
else
    # base preceded by $ is grepped by get_dock_files.csh
    $base/dock32
endif

# Changes to this file must be made in parallel with get_dock_files.csh
