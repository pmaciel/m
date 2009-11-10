#!/usr/bin/env python

# the following two variables are used by the target "waf dist"
VERSION='0.4.0'
APPNAME='m'

# these variables are mandatory ('/' are converted automatically)
srcdir = '.'
blddir = 'build'


def init():
  pass


def set_options(opt):

  # options provided by compiler_cxx, and some extra ones
  opt.tool_options('compiler_cxx')
  opt.add_option('--build',    type='string', default='',       dest='m_build',    help='build environment ([default|debug|release])')
  opt.add_option('--cxxflags', type='string', default='',       dest='m_cxxflags', help='C++ compilation flags (like -fopenmp)')
  opt.add_option('--modules',  type='string', default='msmurf', dest='m_modules',  help='modules directories (like smurf)')


def configure(ctx):
  import Options

  # set compiler_cxx options
  ctx.check_tool('compiler_cxx')

  # set kernel/modules and build environment
  ctx.env.m_mlibraries = []
  ctx.env.m_mversion = VERSION
  ctx.env.m_modules  = ['mkernel']
  if len(Options.options.m_modules):  ctx.env.m_modules.extend(Options.options.m_modules.split(','))
  if len(Options.options.m_cxxflags): ctx.env['CXXFLAGS'].extend(Options.options.m_cxxflags.split(','))
  print 'Info: modules directories: ',ctx.env.m_modules

  ctx.env.m_minikernel = Options.options.m_build=='minikernel'
  if    Options.options.m_build=='debug':    ctx.env['CXXFLAGS'] += ['-Wall','-O0','-g','-ggdb']
  elif  Options.options.m_build=='release':  ctx.env['CXXFLAGS'] += ['-Wall','-O3']
  else: Options.options.m_build = 'default'
  print 'Info: C++ compilation flags (',Options.options.m_build,'): ',ctx.env['CXXFLAGS']

  # modules configure
  ctx.sub_config(ctx.env.m_modules)


def build(ctx):

  # modules build (will install to the same directory)
  ctx.add_subdirs(ctx.env.m_modules)
  for t in [] + ctx.all_task_gen:
    t.install_path = '${PREFIX}'

  # nice target
  ctx.new_task_gen(name='love', rule='echo not war?')

