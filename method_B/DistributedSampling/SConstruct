#/******************************************************************************
# * SConstruct
# *
# * Source of the sampling routine
# ******************************************************************************
# * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
# *
# * This program is free software: you can redistribute it and/or modify it
# * under the terms of the GNU General Public License as published by the Free
# * Software Foundation, either version 3 of the License, or (at your option)
# * any later version.
# *
# * This program is distributed in the hope that it will be useful, but WITHOUT
# * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# * more details.
# *
# * You should have received a copy of the GNU General Public License along with
# * this program.  If not, see <http://www.gnu.org/licenses/>.
# *****************************************************************************/

# scons build file for the graph generator.
#
# You can build it in the following variants:
#
#   optimized            no debug symbols, no assertions, optimization.
#
#   scons variant=${variant} program=${program}
import os
import platform
import sys

# Get the current platform.
SYSTEM = platform.uname()[0]
HOST = platform.uname()[1]

# Get shortcut to $HOME.
HOME = os.environ['HOME']

def GetEnvironment():
    opts = Variables()
    opts.Add('variant', 'the variant to build', 'optimized')
    opts.Add('program', 'program or interface to compile', 'library, run_methodH, run_methodR, run_methodSH, run_experiments')

    env = Environment(options=opts, ENV=os.environ)
    if not env['variant'] in ['optimized', 'debug']:
        print 'Illegal value for variant: %s' % env['variant']
        sys.exit(1)

    if not env['program'] in ['library', 'run_methodH', 'run_methodR', 'run_methodSH', 'run_experiments']:
        print 'Illegal value for program: %s' % env['program']
        sys.exit(1)

    # Special configuration for 64 bit machines.
    if platform.architecture()[0] == '64bit':
        env.Append(CPPFLAGS=['-DPOINTER64=1'])

    return env

# Get the common environment.
env = GetEnvironment()

# env.Append(LIBPATH=['./extern/argtable2/lib'])
# env.Append(CPPPATH=['./extern/argtable2/include'])
# env.Append(CPPPATH=['./extern/stocc'])
# env.Append(CPPPATH=['./extern/mersenne'])
env.Append(LIBPATH=['./optimized'])
env.Append(CPPPATH=['./extern/stocc'])
env.Append(CPPPATH=['./extern/mersenn4'])
env.Append(CPPPATH=['./extern/dSFMT'])
env.Append(CPPPATH=['./lib'])
env.Append(CPPPATH=['./lib/io'])
env.Append(CPPPATH=['./lib/tools'])
env.Append(CPPPATH=['./lib/sampling'])
# env.Append(LIBPATH=['../extern/argtable2/lib'])
# env.Append(CPPPATH=['../extern/argtable2/include'])
# env.Append(CPPPATH=['../extern/stocc'])
# env.Append(CPPPATH=['../extern/mersenne'])
env.Append(LIBPATH=['../optimized'])
env.Append(CPPPATH=['../extern/stocc64'])
env.Append(CPPPATH=['../extern/mersenne64'])
env.Append(CPPPATH=['../extern/spooky'])
env.Append(CPPPATH=['../extern/city'])
env.Append(CPPPATH=['../extern/dSFMT'])
env.Append(CPPPATH=['../lib'])
env.Append(CPPPATH=['../lib/io'])
env.Append(CPPPATH=['../lib/tools'])
env.Append(CPPPATH=['../lib/sampling'])

conf = Configure(env)

if SYSTEM == 'Darwin':
    env.Append(CPPPATH=['/opt/local/include/','../include'])
    env.Append(LIBPATH=['/opt/local/lib/'])
    # homebrew related paths
    env.Append(LIBPATH=['/usr/local/lib/'])

env.Append(CXXFLAGS = '-fopenmp')

# Apply variant specific settings.
if env['variant'] == 'optimized':
    env.Append(CXXFLAGS = '-DNDEBUG -DHAVE_SSE2 -Wall -march=native -mavx2 -fno-stack-limit -O3 -g -std=c++0x -Wall -Wno-literal-suffix -msse4.2 -ffast-math')
    env.Append(CCFLAGS  = '-DNDEBUG -DHAVE_SSE2 -march=native -mavx2 -O3 -g -ffast-math')

elif env['variant'] == 'debug':
# A little bit more output on the console
    env.Append(CXXFLAGS = '-DDEBUG -Wall -fno-stack-limit -g -pg -std=c++0x')
    env.Append(CCFLAGS  = '-DDEBUG -g -pg -std=c++0x')
    env.Append(LINKFLAGS = '-pg')

else:
    env.Append(CXXFLAGS = '-DNDEBUG -Wall -funroll-loops -fno-stack-limit -O3')
    env.Append(CCFLAGS  = '-O3 -DNDEBUG -funroll-loops')

# Execute the SConscript.
SConscript('SConscript', exports=['env'],variant_dir=env['variant'], duplicate=False)
