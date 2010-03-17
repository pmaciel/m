#!/usr/bin/env python


# the following two variables are used by the target "waf dist"
VERSION='0.4.0'
APPNAME='m'

# these variables are mandatory ('/' are converted automatically)
top = '.'
out = 'build'


def set_options(ctx):

  # build environment, plugins and compilers/flags
  ctx.add_option('--build',   type='string', default='default', dest='m_build',   help='build environment(s) ([default|debug|optim],...)')
  ctx.add_option('--plugins', type='string', default='',        dest='m_plugins', help='plugins directories (''plugins'' looks recursively into plugins directory)')
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
  from os import path

  # set generic, C++, C and Fortran compiler/flags
  ctx.check_tool('compiler_cxx')
  ctx.env['CXXFLAGS'] = filter( None, ctx.env['CXXFLAGS']
    + Options.options.m_flags.split(',')
    + Options.options.m_cxxflags.split(',') )

  if Options.options.enable_cc:
    ctx.check_tool('compiler_cc')
    ctx.env['CCFLAGS'] = filter( None, ctx.env['CCFLAGS']
      + Options.options.m_flags.split(',')
      + Options.options.m_ccflags.split(',') )

  if Options.options.enable_fc:
    ctx.check_tool('compiler_fortran')
    ctx.env['FCFLAGS'] = filter( None, ctx.env['FCFLAGS']
      + Options.options.m_flags.split(',')
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
    path.walk('plugins',lookforwscripts,None)
  print 'Info: plugins: ',ctx.env.m_plugins
  ctx.sub_config(ctx.env.m_plugins)

  # set build environment(s)
  ctx.env.m_build = Options.options.m_build.split(',')
  if not 'default' in ctx.env.m_build: ctx.env.m_build += ['default']
  for b in ctx.env.m_build:
    if b!='default':
      env = ctx.all_envs['default'].copy()
      env.set_variant(b)
      ctx.set_env_name(b,env)
    ctx.setenv(b)
    ctx.env['PREFIX'] += path.sep + b
    if b=='debug':
      ctx.                              env.append_unique('CXXFLAGS',['-Wall','-O0','-g3'])
      if Options.options.enable_cc: ctx.env.append_unique('CCFLAGS', ['-Wall','-O0','-g3'])
      if Options.options.enable_fc: ctx.env.append_unique('FCFLAGS', ['-Wall','-O0','-g3'])
    elif b=='optim':
      ctx.                              env.append_unique('CXXFLAGS',['-Wall','-O3'])
      if Options.options.enable_cc: ctx.env.append_unique('CCFLAGS', ['-Wall','-O3'])
      if Options.options.enable_fc: ctx.env.append_unique('FCFLAGS', ['-Wall','-O3'])
    ctx.write_config_header('mkernel' + path.sep + 'mconfig.h')
    ctx.setenv('default')

  # summary
  for b in ctx.env.m_build:
    print                               'Info:',b,'compilation flags (C++):     ',ctx.all_envs[b]['CXXFLAGS']
    if Options.options.enable_cc: print 'Info:',b,'compilation flags (C):       ',ctx.all_envs[b]['CCFLAGS']
    if Options.options.enable_fc: print 'Info:',b,'compilation flags (Fortran): ',ctx.all_envs[b]['FCFLAGS']


def build(ctx):

  # plugins build:
  # - single installation directory
  # - prefix imp/exporting symbols qualifiers for dll's
  ctx.add_subdirs(ctx.env.m_plugins)
  tasks = [] + ctx.all_task_gen
  for b in ctx.env.m_build:
    for task in tasks:
      t = task if b=='default' else task.clone(b)
      t.name  = t.target + '_'+b
      t.rpath = t.install_path = ctx.env_of_name(b)['PREFIX']
      if 'uselib_local' in dir(t):  t.uselib_local = [l + '_'+b for l in t.uselib_local]
      # if 'cshlib' in t.features.split(' '): t.defines = t.defines.split(' ') + ['M_CSHLIB']

  ## nice target
  #ctx.new_task_gen(name='love', rule='echo not war?')

