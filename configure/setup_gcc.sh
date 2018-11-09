me="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
above="$( cd "$( dirname "$me" )" && pwd )"
echo $above

export SBNFITDIR=$above
export SBNFIT_LIBDIR=$me/build

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SBNFIT_LIBDIR

export CXX=/usr/bin/g++
export CC=/usr/bin/gcc

