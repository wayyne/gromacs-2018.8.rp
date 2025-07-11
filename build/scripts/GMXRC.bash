# sh/bash/zsh configuration file for Gromacs
# First we remove old gromacs stuff from the paths
# by selecting everything else.
# This is not 100% necessary, but very useful when we
# repeatedly switch between gmx versions in a shell.

#we make use of IFS, which needs shwordsplit in zsh
test -n "${ZSH_VERSION+set}" && setopt shwordsplit
old_IFS="$IFS"
IFS=":"

# First remove gromacs part of ld_library_path
tmppath=""
for i in $LD_LIBRARY_PATH; do
  if test "$i" != "$GMXLDLIB"; then
    tmppath="${tmppath}${tmppath:+:}${i}"
  fi
done
LD_LIBRARY_PATH=$tmppath

# remove gromacs part of PKG_CONFIG_PATH
tmppath=""
for i in $PKG_CONFIG_PATH; do
  if test "$i" != "$GMXLDLIB/pkgconfig"; then
    tmppath="${tmppath}${tmppath:+:}${i}"
  fi
done
PKG_CONFIG_PATH=$tmppath

# remove gromacs part of path
tmppath=""
for i in $PATH; do
  if test "$i" != "$GMXBIN"; then
    tmppath="${tmppath}${tmppath:+:}${i}"
  fi
done
PATH=$tmppath

# and remove the gmx part of manpath
tmppath=""
for i in $MANPATH; do
  if test "$i" != "$GMXMAN"; then
    tmppath="${tmppath}${tmppath:+:}${i}"
  fi
done
MANPATH=$tmppath

##########################################################
# This is the real configuration part. We save the Gromacs
# things in separate vars, so we can remove them later.
# If you move gromacs, change the first line.
##########################################################
GMXPREFIX=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/local
GMXBIN=${GMXPREFIX}/bin
GMXLDLIB=${GMXPREFIX}/lib64
GMXMAN=${GMXPREFIX}/share/man
GMXDATA=${GMXPREFIX}/share/gromacs
GROMACS_DIR=${GMXPREFIX}

LD_LIBRARY_PATH=${GMXLDLIB}${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH}
PKG_CONFIG_PATH=${GMXLDLIB}/pkgconfig${PKG_CONFIG_PATH:+:}${PKG_CONFIG_PATH}
PATH=${GMXBIN}${PATH:+:}${PATH}
#debian/ubuntu needs a : at the end
MANPATH=${GMXMAN}:${MANPATH}

# export should be separate, so /bin/sh understands it
export GMXBIN GMXLDLIB GMXMAN GMXDATA LD_LIBRARY_PATH PATH MANPATH
export PKG_CONFIG_PATH GROMACS_DIR

IFS="$old_IFS"
unset old_IFS

# read bash completions if understand how to use them
# and this shell supports extended globbing
if test -n "${BASH_VERSION+set}" && (complete) > /dev/null 2>&1; then
  if (shopt -s extglob) > /dev/null 2>&1; then
    shopt -s extglob
    if [ -f $GMXBIN/gmx-completion.bash ]; then
      source $GMXBIN/gmx-completion.bash
      for cfile in $GMXBIN/gmx-completion-*.bash ; do
        source $cfile
      done
    fi
  fi
elif test -n "${ZSH_VERSION+set}" > /dev/null 2>&1 ; then
  autoload bashcompinit
  if (bashcompinit) > /dev/null 2>&1; then
    bashcompinit
    if [ -f $GMXBIN/gmx-completion.bash ]; then
      source $GMXBIN/gmx-completion.bash
      for cfile in $GMXBIN/gmx-completion-*.bash ; do
        source $cfile
      done
    fi
  fi
fi
