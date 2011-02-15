#!/usr/bin/env python


top = '.'
out = 'build'
VERSION='0.0.1'
APPNAME='plas'


# missing things:
# 1. CC = mpicc, CFLAGS += -DMPI -Dlowercase_
# 2. set configuration header options
# # (import/export symbols depending on compiler)
# # - prefix imp/exporting symbols qualifiers for dll's
# # - if 'cshlib': t.defines.extend(['M_CSHLIB'])
#  conf.define( 'M_SYMBOL_EXPORT',
#    '__attribute__ ((dllexport))' if conf.env['CXX_NAME']=='gcc'
#    else '__declspec(dllexport)'  if conf.env['CXX_NAME']=='msvc'
#    else '', quote=False )
#  conf.define( 'M_SYMBOL_IMPORT',
#    '__attribute__ ((dllimport))' if conf.env['CXX_NAME']=='gcc'
#    else '__declspec(dllimport)'  if conf.env['CXX_NAME']=='msvc'
#    else '', quote=False )
#  conf.write_config_header('mconfig.h')


def options(opt):

  # compiler flags
  opt.load('compiler_cxx')
  opt.add_option('--build',  type='string', default='optim', dest='m_mbuild',  help='Build type ([debug|optim|...])')
  opt.add_option('--flags',  type='string', default='',      dest='m_mflags',  help='Compilation flags ([-fPIC[,-fopenmp[,...]]])')


def configure(conf):
  from Options import options
  from os import path

  if not options.m_mbuild:
    options.m_mbuild = 'default'
  conf.env['RPATH'] = path.abspath(conf.options.prefix)

  conf.load('compiler_cxx')
  conf.env.append_unique('CXXFLAGS', filter(None,conf.env['CXXFLAGS'] + options.m_mflags.split(',')) )
  if   options.m_mbuild=='optim': conf.env.append_unique('CXXFLAGS',['-Wall','-fPIC','-O3','-ffast-math'])
  elif options.m_mbuild=='debug': conf.env.append_unique('CXXFLAGS',['-Wall','-fPIC','-O0','-g3'])
  elif options.m_mbuild=='devel': conf.env.append_unique('CXXFLAGS',[
    '-Wall','-fPIC','-O0','-g3',
    '-Wextra','-Wno-unused-parameter','-Wundef','-Wshadow','-Winline' ])

  conf.recurse('src')


def build(bld):
  # bld(name='love', rule='echo not war?')
  bld.recurse('src')

  # force install path to be the same as library search path
  from waflib.TaskGen import feature, before
  @feature('*')
  @before('process_rule') 
  def process_install_path(self):
    if not getattr(self,'install_path',None):
      self.install_path = bld.env['RPATH']

