#!/usr/bin/env python


top = '.'
out = 'build'
VERSION='0.4.0'
APPNAME='m'
# missing MPI support: CC=mpicc, CFLAGS+=['-DMPI','-Dlowercase_']


def options(opt):

  # build type, plugins and compilers/flags
  opt.add_option('--build',   type='string', default='optim', dest='m_mbuild',   help='build type ([debug|optim|...])')
  opt.add_option('--plugins', type='string', default='',      dest='m_mplugins', help='plugins directories (''plugins'' looks recursively into plugins directory)')

  opt.load('compiler_cxx')  # mandatory compiler
  opt.load('compiler_c')    # optional compilers
  opt.load('compiler_fc')   # ...
  opt.add_option('--enable-c',  action='store_true', default=False, help='Check for C compiler')
  opt.add_option('--enable-fc', action='store_true', default=False, help='Check for Fortran compiler')
  opt.add_option('--flags', type='string', default='-fPIC,-fopenmp', dest='m_flags', help='Generic compilation flags, comma-separated (default: ''-fPIC,-fopenmp'')')
  opt.add_option('--flags-cxx', type='string', default='', dest='m_cxxflags', help='Compilation flags, comma-separated (C++)')
  opt.add_option('--flags-c',   type='string', default='', dest='m_cflags',   help='Compilation flags, comma-separated (C)')
  opt.add_option('--flags-fc',  type='string', default='', dest='m_fcflags',  help='Compilation flags, comma-separated (Fortran)')


def configure(conf):
  from Options import options
  from os import path

  # general things
  conf.env.m_mbuild   = options.m_mbuild if options.m_mbuild else 'user'
  conf.env.m_mversion = VERSION

  # libraries search in (flat) installation directory
  conf.env['RPATH'] = path.abspath(conf.options.prefix + path.sep + conf.env.m_mbuild)

  # set compiler flags
  if True:
    conf.load('compiler_cxx')
    conf.env.append_unique('CXXFLAGS', filter(None,conf.env['CXXFLAGS']
      + options.m_flags.split(',') + options.m_cxxflags.split(',') ))
  if options.enable_c:
    conf.load('compiler_c')
    conf.env.append_unique('CFLAGS',filter(None,conf.env['CFLAGS']
      + options.m_flags.split(',') + options.m_cflags.split(',') ))
  if options.enable_fc:
    conf.load('compiler_fc')
    conf.env.append_unique('FCFLAGS',filter(None,conf.env['FCFLAGS']
      + options.m_flags.split(',') + options.m_fcflags.split(',') ))

  if conf.env.m_mbuild=='debug':
    conf.                      env.append_unique('CXXFLAGS',['-Wall','-O0','-g3'])
    if options.enable_c:  conf.env.append_unique('CFLAGS',  ['-Wall','-O0','-g3'])
    if options.enable_fc: conf.env.append_unique('FCFLAGS', ['-Wall','-O0','-g3'])
  elif conf.env.m_mbuild=='devel':
    conf.                      env.append_unique('CXXFLAGS',['-Wall','-O0','-g3','-Wextra','-Wno-unused-parameter','-Wshadow','-Winline','-Wundef'])
    if options.enable_c:  conf.env.append_unique('CFLAGS',  ['-Wall','-O0','-g3','-Wextra','-Wno-unused-parameter','-Wshadow','-Winline','-Wundef'])
    if options.enable_fc: conf.env.append_unique('FCFLAGS', ['-Wall','-O0','-g3','-Wextra','-Wno-unused-parameter','-Wshadow','-Winline'])
  elif conf.env.m_mbuild=='optim':
    conf.                      env.append_unique('CXXFLAGS',['-Wall','-O3','-ffast-math'])
    if options.enable_c:  conf.env.append_unique('CFLAGS',  ['-Wall','-O3','-ffast-math'])
    if options.enable_fc: conf.env.append_unique('FCFLAGS', ['-Wall','-O3','-ffast-math'])

  print                       'I: build type: ',conf.env.m_mbuild
  print                       'I: compilation flags (C++):     ',conf.env['CXXFLAGS']
  if options.enable_c:  print 'I: compilation flags (C):       ',conf.env['CFLAGS']
  if options.enable_fc: print 'I: compilation flags (Fortran): ',conf.env['FCFLAGS']


  # set plugins & libraries
  # 1. add libraries from plugins directory, recursively)
  # 2. add mkernel library unconditionally, in the first place)
  # 3. set (unique plugins)
  # 4. configure plugins/libraries
  p = options.m_mplugins.split(',') if len(options.m_mplugins) else []
  if 'plugins' in p:
    while 'plugins' in p: p.remove('plugins')
    def lookforwscripts(args,dir,files):
      if 'wscript' in files: p.append(dir)
    path.walk('plugins',lookforwscripts,None)
  conf.env.m_mplugins = list(set(p))
  while 'mkernel' in conf.env.m_mplugins: conf.env.m_mplugins.remove('mkernel')
  conf.env.m_mplugins.insert(0,'mkernel')
  conf.env.m_mlibraries = []
  conf.recurse(conf.env.m_mplugins)
  print 'I: plugins:   ',conf.env.m_mplugins
  print 'I: libraries: ',conf.env.m_mlibraries


def build(bld):
  from waflib.TaskGen import feature, before

  # bld(name='love', rule='echo not war?')
  bld.recurse(bld.env.m_mplugins)

  # force install path to be the same as library search path
  @feature('*')
  @before('process_rule') 
  def process_install_path(self):
    if not getattr(self,'install_path',None):
      self.install_path = bld.env['RPATH']

