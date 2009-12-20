#!/usr/bin/env python


def configure(ctx):
  ctx.env.m_mlibraries.append('ls_trilinos')
  ctx.env['CPPFLAGS_TRILINOS']  = ['-I/home/pedro/local/x86_64/local/trilinos/install/include']
  ctx.env['LIB_TRILINOS']       = ['stratimikos']  # apparently links with aztec and others...
# ctx.env['LIB_TRILINOS']       = ['aztecoo','epetra','ml']
  ctx.env['LIBPATH_TRILINOS']   = ['/home/pedro/local/x86_64/local/trilinos/install/lib']
  ctx.env['RPATH_TRILINOS']     = ['/home/pedro/local/x86_64/local/trilinos/install/lib']
  ctx.env['LINKFLAGS_TRILINOS'] = ['-fopenmp']


def build(ctx):
  ctx.new_task_gen(
    features = 'cxx cshlib',
    name     = 'ls_trilinos',
    target   = 'ls_trilinos',
    includes = ['.','../../mkernel'],
    source   = ['ls_trilinos_aztecoo.cpp','ls_trilinos_stratimikos.cpp'],
    uselib   = ['TRILINOS'] )

