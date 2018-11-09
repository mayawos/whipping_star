me="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
above="$( cd "$( dirname "$me" )" && pwd )"
echo $above

export SBNFITDIR=$above
export SBNFIT_LIBDIR=$me/build

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SBNFIT_LIBDIR

export CXX=/data/gpu/linux86-64/2018/bin/pgc++
export CC=/data/gpu/linux86-64/2018/bin/pgcc



