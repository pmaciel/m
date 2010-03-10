#!/usr/bin/env python


# the following two variables are used by the target "waf dist"
VERSION='2.1'
APPNAME='aztec'

# these variables are mandatory ('/' are converted automatically)
top = '.'
out = 'build'


def set_options(ctx):
  ctx.tool_options('compiler_cc')
  ctx.tool_options('compiler_fortran')
  ctx.add_option('--flags',    type='string', default='', dest='m_flags',    help='Compilation flags (generic, like -fPIC,...)')
  ctx.add_option('--ccflags',  type='string', default='', dest='m_ccflags',  help='Compilation flags (C)')
  ctx.add_option('--fcflags',  type='string', default='', dest='m_fcflags',  help='Compilation flags (Fortran)')


def configure(ctx):
  import Options

  ctx.check_tool('compiler_cc')
  if len(Options.options.m_flags):   ctx.env['CCFLAGS']  += Options.options.m_flags.split(',')
  if len(Options.options.m_ccflags): ctx.env['CCFLAGS']  += Options.options.m_ccflags.split(',')

  ctx.check_tool('compiler_fortran')
  if len(Options.options.m_flags):   ctx.env['FCFLAGS']  += Options.options.m_flags.split(',')
  if len(Options.options.m_fcflags): ctx.env['FCFLAGS']  += Options.options.m_fcflags.split(',')
  if not ctx.check_fortran():
    ctx.fatal('cannot compile a simple fortran program!')
  #ctx.check_fortran_clib()
  ##st, mangler = ctx.check_fortran_mangling()
  ##if not st:
  ##  ctx.fatal('cannot detect the mangling scheme')
  # check Fortran if dummy main routine is needed
  ctx.check_fortran_dummy_main()


def build(ctx):

  # aztec static/dynamic library (standalone version, needs Lapack and BLAS)
  l = ctx(
    features = 'fortran cc cstaticlib',
    target   = 'aztec',
    defines  = ['append_','AZ_SERIAL'],
    lib      = ['lapack','blas','gfortran'],
    includes = ['.'],
    source   = [
      # aztec
      'az_bilu.c', 'az_cg.c', 'az_cgs.c', 'az_cgstab.c', 'az_check.c',
      'az_comm.c', 'az_converge.c', 'az_dd_overlap.c', 'az_dgemv2.c',
      'az_dgemv3.c', 'az_domain_decomp.c', 'az_fortran_wrap.c', 'az_scaling.c',
      'az_flop_cnt.c', 'az_gmres.c', 'az_gmresr.c', 'az_ilu_util.c',
      'az_ilut.c', 'az_interface.c', 'az_lu_y12.c', 'az_matrix_util.c',
      'az_matvec_mult.c', 'az_old_matvec_mult.c', 'az_pad_utils.c', 'az_poly.c',
      'az_precond.c', 'az_qmrcgs.c', 'az_reorder.f', 'az_rilu.c', 'az_solve.c',
      'az_sort.c', 'az_subdomain_solver.c', 'az_tools.c', 'az_util.c',
      'az_icc.c', 'az_fix_pt.c', 'md_wrap_scalar_c.c', 'md_timer_generic.c',
      # netlib Y12M
      'y12m/y12m.f', 'y12m/y12mae.f', 'y12m/y12maf.f', 'y12m/y12mbe.f',
      'y12m/y12mbf.f', 'y12m/y12mce.f', 'y12m/y12mcf.f', 'y12m/y12mde.f',
      'y12m/y12mdf.f', 'y12m/y12mfe.f', 'y12m/y12mge.f', 'y12m/y12mhe.f',
      'y12m/y12cck.f' ] )
  ctx.new_task_gen(
    features = 'fortran cc cshlib',
    ccflags  = ['-fPIC'],
    target   = l.target,
    defines  = l.defines,
    lib      = l.lib,
    includes = l.includes,
    source   = l.source )

