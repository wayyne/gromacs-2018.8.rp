# csh/tcsh configuration file for Gromacs.
# First we remove previous Gromacs stuff from paths 
# by selecting everything else. 
# Removal is not 100% necessary, but very useful when we
# repeatedly switch between gmx versions in a shell.

# zero possibly unset vars to avoid warnings
if (! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH ""
if (! $?PKG_CONFIG_PATH) setenv PKG_CONFIG_PATH ""
if (! $?PATH) setenv PATH ""
if (! $?MANPATH) setenv MANPATH ""
if (! $?GMXLDLIB) setenv GMXLDLIB ""
if (! $?GMXBIN) setenv GMXBIN ""
if (! $?GMXMAN) setenv GMXMAN ""

# remove previous gromacs part from ld_library_path
set tmppath = ""
foreach i ( `echo $LD_LIBRARY_PATH | sed "s/:/ /g"` )
  if ( "$i" != "$GMXLDLIB" ) then
    if ("${tmppath}" == "") then
      set tmppath = "$i"
    else
      set tmppath = "${tmppath}:$i"
    endif
  endif
end
setenv LD_LIBRARY_PATH $tmppath

# remove previous gromacs part from PKG_CONFIG_PATH
set tmppath = ""
foreach i ( `echo $PKG_CONFIG_PATH | sed "s/:/ /g"` )
  if ( "$i" != "$GMXLDLIB/pkgconfig" ) then
    if ("${tmppath}" == "") then
      set tmppath = "$i"
    else
      set tmppath = "${tmppath}:$i"
    endif
  endif
end
setenv PKG_CONFIG_PATH $tmppath

# remove gromacs stuff from binary path
set tmppath = ""
foreach i ( `echo $PATH | sed "s/:/ /g"` )
  if ( "$i" != "$GMXBIN" ) then
    if ("${tmppath}" == "") then
      set tmppath = "$i"
    else
      set tmppath = "${tmppath}:$i"
    endif
  endif
end
setenv PATH $tmppath

# and remove stuff from manual path
set tmppath = ""
foreach i ( `echo $MANPATH | sed "s/:/ /g"` )
  if ( "$i" != "$GMXMAN" ) then 
    if ("${tmppath}" == "") then
      set tmppath = "$i"
    else
      set tmppath = "${tmppath}:$i"
    endif
  endif
end
setenv MANPATH $tmppath

##########################################################
# This is the real configuration part. We save the Gromacs
# things in separate vars, so we can remove them later.
# If you move gromacs, change the first line.
##########################################################
setenv GMXPREFIX /work_bgfs/g/gdayhoff/grayson/lab/gmx+/local
setenv GMXBIN ${GMXPREFIX}/bin
setenv GMXLDLIB ${GMXPREFIX}/lib64
setenv GMXMAN ${GMXPREFIX}/share/man
setenv GMXDATA ${GMXPREFIX}/share/gromacs
setenv GROMACS_DIR ${GMXPREFIX}

#make them begin with :
if ($?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH ":${LD_LIBRARY_PATH}"
if ($?PKG_CONFIG_PATH) setenv PKG_CONFIG_PATH ":${PKG_CONFIG_PATH}"

#path is never empty
setenv PATH ${GMXBIN}:${PATH}
setenv LD_LIBRARY_PATH ${GMXLDLIB}${LD_LIBRARY_PATH}
setenv PKG_CONFIG_PATH ${GMXLDLIB}/pkgconfig${PKG_CONFIG_PATH}
#debian/ubuntu needs a : at the end
setenv MANPATH ${GMXMAN}:${MANPATH}

setenv GMX_FONT	10x20

# Read completions if we understand it (i.e. have tcsh)
# Currently disabled, since the completions don't work with the new
# gmx binary with subcommands.
# Contributions to get the functionality back are welcome.
#if { complete >& /dev/null } then
#  if ( -f $GMXBIN/completion.csh ) source $GMXBIN/completion.csh
#endif
