import os
import re
import fnmatch
import multiprocessing

LIBPATH = [
    # insert libs here
    os.path.join(GetOption('outdir'), 'lib'),
]

LIBS = [
    # add dependency names here as they should be linked against
]

COMMONFLAGS=[
    # add language/compiler agnostic flags here+
    '-fPIC',
    '-Wabi',
    '-msse4.2',
    '-Wno-unused-parameter',
    '-Wall',
    '-Wextra',
    '-Werror',
    '-I{}'.format(os.path.join(Dir('#').abspath, 'src/')),

    # phys units degree symbol
    '-Wno-invalid-source-encoding',
]

CPPPATH = [
    # add c++ paths here
]

CPPFLAGS = [
    # required C++ flags
    '-std=c++14',

    '-I{}'.format(os.path.join(Dir('#').abspath, 'src/ext/phys_units')),   # to satisfy the inter-including
]

GCC_CPPFLAGS = [
    # TODO: GCC 4.9 give ABI incompatibility warnings upgraded to errors  without this
    '-fabi-version=10',
]

CFLAGS = [
    # required C flags
    '-std=c14',
]

DBGFLAGS = [
    # debug flags
    '-g',
    '-O0',
]

NDBGFLAGS = [
    # non-debug build flags
    '-O3',
]

CXXCOMPILER='g++'
CCCOMPILER='gcc'

if GetOption('clang'):
    CXXCOMPILER='clang++' if not GetOption('travis') else 'clang++-3.8'
    CCCOMPILER='clang' if not GetOption('travis') else 'clang-3.8'



def append_flags(src, to):
    for i in src:
        if i not in to:
            to.append(i)

# what will actually be used to build
FLAGSET = COMMONFLAGS

if GetOption('dbg'):
    append_flags(DBGFLAGS, FLAGSET)
else:
    append_flags(NDBGFLAGS, FLAGSET)


##
## Helpers
##

def get_files(env, root, exts):
    #directory of SConstruct
    rootpath = Dir('#').abspath
    srcroot = os.path.join(rootpath, root)
    if not os.path.exists(srcroot):
        raise Exception('{} does not exist'.format(srcroot))

    src = []
    filt = re.compile('.*\.({})$'.format('|'.join(exts)))
    for root, dirs, files in os.walk(srcroot, followlinks=True):
        for filename in files:
	    if filt.match(filename):
            	src.append(os.path.join(root, filename)[len(srcroot)+1:])
    return src
AddMethod(Environment, get_files)


##
## Build Targets
##

builds = [
    # add jobs here
    'src/prog',
]
for b in builds:
    idir =  os.path.join(Dir('#').abspath, b)
    ldir = os.path.join(Dir('#').abspath, GetOption('outdir'), os.path.basename(b))

    CPPFLAGS.extend(['-I{}'.format(idir)])
    CFLAGS.extend(['-I{}'.format(idir)])
    LIBPATH.extend([ldir])


##
## Build Entry
##
SetOption('num_jobs', int(multiprocessing.cpu_count()))

CPPFLAGS.extend(COMMONFLAGS),
CFLAGS.extend(COMMONFLAGS),
global_env = Environment(
    LIBS=LIBS,
    LIBPATH=LIBPATH,
    CPPPATH=CPPPATH,
    CC=CCCOMPILER,
    CXX=CXXCOMPILER,
    CXXFLAGS=CPPFLAGS,
    CFLAGS=CFLAGS,

    CCCOMSTR       = "[ CC  ]\t$TARGET" if not GetOption('ugly') else None,
    CXXCOMSTR      = "[ CPP ]\t$TARGET" if not GetOption('ugly') else None,
    LINKCOMSTR     = "[ LD  ]\t$TARGET" if not GetOption('ugly') else None,
    ARCOMSTR       = "[ AR  ]\t$TARGET" if not GetOption('ugly') else None,
    RANLIBCOMSTR   = "[ LIB ]\t$TARGET" if not GetOption('ugly') else None,
)

Export('global_env')

build_targets = {}
for build in builds:
    out = '/'.join(build.split('/')[1:])
    build_targets[build] = SConscript(os.path.join(build, "SConscript"), variant_dir=out, duplicate=False)

#Depends(build_targets['src/prog'], build_targets['src/lib'])
