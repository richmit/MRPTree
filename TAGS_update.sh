#!/bin/sh

if [ -e 'CMakeLists.txt' -a -e 'lib' -a -e 'examples-lib3d' ]; then
  if [ -e 'TAGS' ]; then
    if [ -e 'TAGS.old' ]; then
      mv TAGS.old TAGS.old.old
    fi
    mv TAGS TAGS.old
  fi
  etags lib/*.hpp examples-lib3d/*.hpp
  if [ -e 'TAGS' ]; then
    echo 'TAGS_update.sh: INFO: New tag file created: TAGS'
    if [ -e 'TAGS.old' ]; then
      echo 'TAGS_update.sh: INFO: Backup tag file: TAGS.old'
      echo 'TAGS_update.sh: INFO: Backup & new tag file diff:'
      diff -s TAGS.old TAGS 
    else
      echo 'TAGS_update.sh: INFO: No old TAG file to backup.'
    fi
  else
    echo 'TAGS_update.sh: ERROR: Unable to create TAGS file!'
  fi
else
  echo 'TAGS_update.sh: ERROR: Could not find files -- should be in root of repo'
fi
