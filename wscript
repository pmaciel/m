#!/usr/bin/env python


# the following two variables are used by the target "waf dist"
VERSION='0.4.0'
APPNAME='m'

# these variables are mandatory ('/' are converted automatically)
srcdir = '.'
blddir = 'build'


def set_options(ctx):

  # build environment, plugins and compilers/flags
  ctx.add_option('--build',   type='string', default='', dest='m_build',   help='build environment ([default|debug|optim])')
  ctx.add_option('--plugins', type='string', default='', dest='m_plugins', help='plugins directories (''plugins'' looks recursively into plugins directory)')
  ctx.add_option('--enable-cc', action='store_true', default=False, help='Check for C compiler')
  ctx.add_option('--enable-fc', action='store_true', default=False, help='Check for Fortran compiler')
  ctx.tool_options('compiler_cxx')
  ctx.tool_options('compiler_cc')
  ctx.tool_options('compiler_fortran')
  ctx.add_option('--flags',    type='string', default='', dest='m_flags',    help='Compilation flags (generic, like -fPIC,-fopenmp,...)')
  ctx.add_option('--cxxflags', type='string', default='', dest='m_cxxflags', help='Compilation flags (C++)')
  ctx.add_option('--ccflags',  type='string', default='', dest='m_ccflags',  help='Compilation flags (C)')
  ctx.add_option('--fcflags',  type='string', default='', dest='m_fcflags',  help='Compilation flags (Fortran)')


def configure(ctx):
  import Options

  # set build environment
  Options.options.m_flags = Options.options.m_flags.split(',')
  if    Options.options.m_build=='debug': Options.options.m_flags += ['-Wall','-O0','-g','-ggdb']
  elif  Options.options.m_build=='optim': Options.options.m_flags += ['-Wall','-O3']
  else: Options.options.m_build = 'default'

  # set C++ compiler/flags
  ctx.check_tool('compiler_cxx')
  ctx.env['CXXFLAGS'] = filter( None, ctx.env['CXXFLAGS']
    + Options.options.m_flags
    + Options.options.m_cxxflags.split(',') )

  # set C compiler/flags
  if Options.options.enable_cc:
    ctx.check_tool('compiler_cc')
    ctx.env['CCFLAGS'] = filter( None, ctx.env['CCFLAGS']
      + Options.options.m_flags
      + Options.options.m_ccflags.split(',') )

  # set Fortran compiler/flags (and run some checks)
  if Options.options.enable_fc:
    ctx.check_tool('compiler_fortran')
    ctx.env['FCFLAGS'] = filter( None, ctx.env['FCFLAGS']
      + Options.options.m_flags
      + Options.options.m_fcflags.split(',') )
    if not ctx.check_fortran():
      ctx.fatal('cannot compile a simple fortran program!')
    #ctx.check_fortran_clib()
    ##st, mangler = ctx.check_fortran_mangling()
    ##if not st:
    ##  ctx.fatal('cannot detect the mangling scheme')
    # check Fortran if dummy main routine is needed
    ctx.check_fortran_dummy_main()

  # set mkernel and plugins
  ctx.env.m_mlibraries = []
  ctx.env.m_mversion = VERSION
  ctx.env.m_plugins  = ['mkernel']
  if len(Options.options.m_plugins):  ctx.env.m_plugins.extend(Options.options.m_plugins.split(','))
  if 'plugins' in ctx.env.m_plugins:
    while 'plugins' in ctx.env.m_plugins: ctx.env.m_plugins.remove('plugins')
    def lookforwscripts(args,dir,files):
      if 'wscript' in files: ctx.env.m_plugins.append(dir)
    from os import path
    path.walk('plugins',lookforwscripts,None)

  # summary
  print 'Info: compilation flags (C++,     ',Options.options.m_build,'): ',ctx.env['CXXFLAGS']
  if Options.options.enable_cc: print 'Info: compilation flags (C,       ',Options.options.m_build,'): ',ctx.env['CCFLAGS']
  if Options.options.enable_fc: print 'Info: compilation flags (Fortran, ',Options.options.m_build,'): ',ctx.env['FCFLAGS']
  print 'Info: plugins: ',ctx.env.m_plugins
  ctx.sub_config(ctx.env.m_plugins)


def build(ctx):

  # plugins build:
  # - single installation directory
  # - prefix imp/exporting symbols qualifiers for dll's
  ctx.add_subdirs(ctx.env.m_plugins)
  for t in [] + ctx.all_task_gen:
    t.install_path = '${PREFIX}'
    if 'cshlib' in t.features.join(' ').split(' '):
      t.defines = t.defines.split(' ') + ['M_CSHLIB']

  ## nice target
  #ctx.new_task_gen(name='love', rule='echo not war?')

