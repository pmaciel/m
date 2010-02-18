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
  opt.add_option('--build',    type='string', default='',     dest='m_build',    help='build environment ([default|debug|release])')
  opt.add_option('--cxxflags', type='string', default='',     dest='m_cxxflags', help='compilation flags (like -fopenmp)')
  opt.add_option('--plugins',  type='string', default='auto', dest='m_plugins',  help='plugins directories (''auto'' looks recursively into plugins directory)')


def configure(ctx):
  import Options

  # set compiler_cxx options and build environment
  ctx.check_tool('compiler_cxx')
  if len(Options.options.m_cxxflags):        ctx.env['CXXFLAGS'] += Options.options.m_cxxflags.split(',')
  if    Options.options.m_build=='debug':    ctx.env['CXXFLAGS'] += ['-Wall','-O0','-g','-ggdb']
  elif  Options.options.m_build=='release':  ctx.env['CXXFLAGS'] += ['-Wall','-O3']
  else: Options.options.m_build = 'default'
  print 'Info: compilation flags (',Options.options.m_build,'): ',ctx.env['CXXFLAGS']

  # set mkernel and plugins
  ctx.env.m_mlibraries = []
  ctx.env.m_mversion = VERSION
  ctx.env.m_plugins  = ['mkernel']

  if len(Options.options.m_plugins):  ctx.env.m_plugins.extend(Options.options.m_plugins.split(','))
  if 'auto' in ctx.env.m_plugins:
    while 'auto' in ctx.env.m_plugins: ctx.env.m_plugins.remove('auto')
    def lookforwscripts(args,dir,files):
      if 'wscript' in files: ctx.env.m_plugins.append(dir)
    from os import path
    path.walk('plugins',lookforwscripts,None)

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

