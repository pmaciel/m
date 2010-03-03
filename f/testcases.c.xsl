<?xml version="1.0" encoding="utf-8"?>
<html xsl:version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml">
 <head>
  <style type="text/css">
   table.simple td {
   	border-style: solid;
   	border-color: gray;
   	border-width: 1px;
   }
   td.h { background-color: #EEEEEE }
   td.H { background-color: #DDDDDD }
</style>
 </head>
 <body style="font-family:sans-serif;font-size:small">

<table class="simple">
 <thead style="background-color:#DDDDDD;font-weight:bold">
  <tr>
   <td>setup</td>
   <td>terms</td>
   <td>solver_nonlinear</td>
   <td>solver_linear</td>
   <td>vars</td>
   <td>bcs</td>
   <td>mitremassembler</td>
  </tr>
 </thead>
 <tbody>
  <xsl:for-each select="cdb/c">
   <tr><td class="h">label</td><td class="h" colspan="6"><xsl:value-of select="@label"/></td></tr>
   <tr><td class="h">comment</td><td colspan="6"><xsl:value-of select="comment"/></td></tr>
  <tr valign="top">

   <!-- setup -->
   <td><table><tbody>
    <tr><td class="h">dimensions </td><td><xsl:value-of select="setup/@dimensions" /></td></tr>
    <tr><td class="h">restart    </td><td><xsl:value-of select="setup/@restart"    /></td></tr>
    <tr><td class="h">file_input </td><td><xsl:value-of select="setup/@file_input" /></td></tr>
    <tr><td class="h">file_inlet </td><td><xsl:value-of select="setup/@file_inlet" /></td></tr>
    <tr><td class="h">file_output</td><td><xsl:value-of select="setup/@file_output"/></td></tr>
   </tbody></table></td>

   <!-- terms -->
   <td><table><tbody>
    <tr><td class="h">density         </td><td><xsl:value-of select="terms/@density"         /></td></tr>
    <tr><td class="h">kviscosity      </td><td><xsl:value-of select="terms/@kviscosity"      /></td></tr>
    <tr><td class="h">prandlt         </td><td><xsl:value-of select="terms/@prandlt"         /></td></tr>
    <tr><td class="h">gx              </td><td><xsl:value-of select="terms/@gx"              /></td></tr>
    <tr><td class="h">gy              </td><td><xsl:value-of select="terms/@gy"              /></td></tr>
    <tr><td class="h">gz              </td><td><xsl:value-of select="terms/@gz"              /></td></tr>
    <tr><td class="h">diffusion       </td><td><xsl:value-of select="terms/@diffusion"       /></td></tr>
    <tr><td class="h">temperature     </td><td><xsl:value-of select="terms/@temperature"     /></td></tr>
    <tr><td class="h">buoyancy        </td><td><xsl:value-of select="terms/@buoyancy"        /></td></tr>
    <tr><td class="h">vardensity      </td><td><xsl:value-of select="terms/@vardensity"      /></td></tr>
    <tr><td class="h">vardensity_value</td><td><xsl:value-of select="terms/@vardensity_value"/></td></tr>
    <tr><td class="h">scalarscheme    </td><td><xsl:value-of select="terms/@scalarscheme"    /></td></tr>
    <tr><td class="h">turbulence      </td><td><xsl:value-of select="terms/@turbulence"      /></td></tr>
   </tbody></table></td>

   <!-- solver_nonlinear -->
   <td><table><tbody>
    <tr><td class="h">method              </td><td><xsl:value-of select="solver_nonlinear/@method"              /></td></tr>
    <tr><td class="h">maxiterations       </td><td><xsl:value-of select="solver_nonlinear/@maxiterations"       /></td></tr>
    <tr><td class="h">convergence_variable</td><td><xsl:value-of select="solver_nonlinear/@convergence_variable"/></td></tr>
    <tr><td class="h">convergence_level   </td><td><xsl:value-of select="solver_nonlinear/@convergence_level"   /></td></tr>
    <tr><td class="h">newton_eps          </td><td><xsl:value-of select="solver_nonlinear/@newton_eps"          /></td></tr>
    <tr><td class="h">newton_switch       </td><td><xsl:value-of select="solver_nonlinear/@newton_switch"       /></td></tr>
    <tr><td class="h">turb_iterinit       </td><td><xsl:value-of select="solver_nonlinear/@turb_iterinit"       /></td></tr>
    <tr><td class="h">couple_t            </td><td><xsl:value-of select="solver_nonlinear/@couple_t"            /></td></tr>
    <tr><td class="h">couple_ke           </td><td><xsl:value-of select="solver_nonlinear/@couple_ke"           /></td></tr>
    <tr><td class="h">relax_puvw_type     </td><td><xsl:value-of select="solver_nonlinear/@relax_puvw_type"     /></td></tr>
    <tr><td class="h">relax_t_type        </td><td><xsl:value-of select="solver_nonlinear/@relax_t_type"        /></td></tr>
    <tr><td class="h">relax_ke_type       </td><td><xsl:value-of select="solver_nonlinear/@relax_ke_type"       /></td></tr>
    <tr><td class="h">relax_puvw_value    </td><td><xsl:value-of select="solver_nonlinear/@relax_puvw_value"    /></td></tr>
    <tr><td class="h">relax_t_value       </td><td><xsl:value-of select="solver_nonlinear/@relax_t_value"       /></td></tr>
    <tr><td class="h">relax_ke_value      </td><td><xsl:value-of select="solver_nonlinear/@relax_ke_value"      /></td></tr>
   </tbody></table></td>

   <!-- system_(coupled|scalar|turb) -->
   <td>
    <xsl:for-each select="system_coupled">
     <table><thead><tr><td class="h" colspan="2">system_coupled</td></tr></thead><tbody>
      <tr><td class="h">type      </td><td><xsl:value-of select="@type"      /></td></tr>
      <!-- aztec options -->
      <tr><td class="h">mtype     </td><td><xsl:value-of select="@mtype"     /></td></tr>
      <tr><td class="h">solver    </td><td><xsl:value-of select="@solver"    /></td></tr>
      <tr><td class="h">scaling   </td><td><xsl:value-of select="@scaling"   /></td></tr>
      <tr><td class="h">precond   </td><td><xsl:value-of select="@precond"   /></td></tr>
      <tr><td class="h">subdomain_solve</td><td><xsl:value-of select="@subdomain_solve"/></td></tr>
      <tr><td class="h">conv      </td><td><xsl:value-of select="@conv"      /></td></tr>
      <tr><td class="h">output    </td><td><xsl:value-of select="@output"    /></td></tr>
      <tr><td class="h">pre_calc  </td><td><xsl:value-of select="@pre_calc"  /></td></tr>
      <tr><td class="h">graph_fill</td><td><xsl:value-of select="@graph_fill"/></td></tr>
      <tr><td class="h">max_iter  </td><td><xsl:value-of select="@max_iter"  /></td></tr>
      <tr><td class="h">poly_ord  </td><td><xsl:value-of select="@poly_ord"  /></td></tr>
      <tr><td class="h">overlap   </td><td><xsl:value-of select="@overlap"   /></td></tr>
      <tr><td class="h">type_overlap   </td><td><xsl:value-of select="@type_overlap"   /></td></tr>
      <tr><td class="h">kspace    </td><td><xsl:value-of select="@kspace"    /></td></tr>
      <tr><td class="h">reorder   </td><td><xsl:value-of select="@reorder"   /></td></tr>
      <tr><td class="h">keep_info </td><td><xsl:value-of select="@keep_info" /></td></tr>
      <tr><td class="h">orthog    </td><td><xsl:value-of select="@orthog"    /></td></tr>
      <tr><td class="h">aux_vec   </td><td><xsl:value-of select="@aux_vec"   /></td></tr>
      <!-- aztec params -->
      <tr><td class="h">tol      </td><td><xsl:value-of select="@tol"      /></td></tr>
      <tr><td class="h">drop     </td><td><xsl:value-of select="@drop"     /></td></tr>
      <tr><td class="h">ilut_fill</td><td><xsl:value-of select="@ilut_fill"/></td></tr>
      <tr><td class="h">omega    </td><td><xsl:value-of select="@omega"    /></td></tr>
     </tbody></table>
    </xsl:for-each>
    <xsl:for-each select="system_scalar">
     <table><thead><tr><td class="h" colspan="2">system_scalar</td></tr></thead><tbody>
      <tr><td class="h">type      </td><td><xsl:value-of select="@type"      /></td></tr>
      <!-- aztec options -->
      <tr><td class="h">mtype     </td><td><xsl:value-of select="@mtype"     /></td></tr>
      <tr><td class="h">solver    </td><td><xsl:value-of select="@solver"    /></td></tr>
      <tr><td class="h">scaling   </td><td><xsl:value-of select="@scaling"   /></td></tr>
      <tr><td class="h">precond   </td><td><xsl:value-of select="@precond"   /></td></tr>
      <tr><td class="h">subdomain_solve</td><td><xsl:value-of select="@subdomain_solve"/></td></tr>
      <tr><td class="h">conv      </td><td><xsl:value-of select="@conv"      /></td></tr>
      <tr><td class="h">output    </td><td><xsl:value-of select="@output"    /></td></tr>
      <tr><td class="h">pre_calc  </td><td><xsl:value-of select="@pre_calc"  /></td></tr>
      <tr><td class="h">graph_fill</td><td><xsl:value-of select="@graph_fill"/></td></tr>
      <tr><td class="h">max_iter  </td><td><xsl:value-of select="@max_iter"  /></td></tr>
      <tr><td class="h">poly_ord  </td><td><xsl:value-of select="@poly_ord"  /></td></tr>
      <tr><td class="h">overlap   </td><td><xsl:value-of select="@overlap"   /></td></tr>
      <tr><td class="h">type_overlap   </td><td><xsl:value-of select="@type_overlap"   /></td></tr>
      <tr><td class="h">kspace    </td><td><xsl:value-of select="@kspace"    /></td></tr>
      <tr><td class="h">reorder   </td><td><xsl:value-of select="@reorder"   /></td></tr>
      <tr><td class="h">keep_info </td><td><xsl:value-of select="@keep_info" /></td></tr>
      <tr><td class="h">orthog    </td><td><xsl:value-of select="@orthog"    /></td></tr>
      <tr><td class="h">aux_vec   </td><td><xsl:value-of select="@aux_vec"   /></td></tr>
      <!-- aztec params -->
      <tr><td class="h">tol      </td><td><xsl:value-of select="@tol"      /></td></tr>
      <tr><td class="h">drop     </td><td><xsl:value-of select="@drop"     /></td></tr>
      <tr><td class="h">ilut_fill</td><td><xsl:value-of select="@ilut_fill"/></td></tr>
      <tr><td class="h">omega    </td><td><xsl:value-of select="@omega"    /></td></tr>
     </tbody></table>
    </xsl:for-each>
    <xsl:for-each select="system_turb">
     <table><thead><tr><td class="h" colspan="2">system_turb</td></tr></thead><tbody>
      <tr><td class="h">type      </td><td><xsl:value-of select="@type"      /></td></tr>
      <!-- aztec options -->
      <tr><td class="h">mtype     </td><td><xsl:value-of select="@mtype"     /></td></tr>
      <tr><td class="h">solver    </td><td><xsl:value-of select="@solver"    /></td></tr>
      <tr><td class="h">scaling   </td><td><xsl:value-of select="@scaling"   /></td></tr>
      <tr><td class="h">precond   </td><td><xsl:value-of select="@precond"   /></td></tr>
      <tr><td class="h">subdomain_solve</td><td><xsl:value-of select="@subdomain_solve"/></td></tr>
      <tr><td class="h">conv      </td><td><xsl:value-of select="@conv"      /></td></tr>
      <tr><td class="h">output    </td><td><xsl:value-of select="@output"    /></td></tr>
      <tr><td class="h">pre_calc  </td><td><xsl:value-of select="@pre_calc"  /></td></tr>
      <tr><td class="h">graph_fill</td><td><xsl:value-of select="@graph_fill"/></td></tr>
      <tr><td class="h">max_iter  </td><td><xsl:value-of select="@max_iter"  /></td></tr>
      <tr><td class="h">poly_ord  </td><td><xsl:value-of select="@poly_ord"  /></td></tr>
      <tr><td class="h">overlap   </td><td><xsl:value-of select="@overlap"   /></td></tr>
      <tr><td class="h">type_overlap   </td><td><xsl:value-of select="@type_overlap"   /></td></tr>
      <tr><td class="h">kspace    </td><td><xsl:value-of select="@kspace"    /></td></tr>
      <tr><td class="h">reorder   </td><td><xsl:value-of select="@reorder"   /></td></tr>
      <tr><td class="h">keep_info </td><td><xsl:value-of select="@keep_info" /></td></tr>
      <tr><td class="h">orthog    </td><td><xsl:value-of select="@orthog"    /></td></tr>
      <tr><td class="h">aux_vec   </td><td><xsl:value-of select="@aux_vec"   /></td></tr>
      <!-- aztec params -->
      <tr><td class="h">tol      </td><td><xsl:value-of select="@tol"      /></td></tr>
      <tr><td class="h">drop     </td><td><xsl:value-of select="@drop"     /></td></tr>
      <tr><td class="h">ilut_fill</td><td><xsl:value-of select="@ilut_fill"/></td></tr>
      <tr><td class="h">omega    </td><td><xsl:value-of select="@omega"    /></td></tr>
     </tbody></table>
    </xsl:for-each>
   </td>

   <!-- vars -->
   <td><table>
    <thead><tr>
     <td class="h">label</td>
     <td class="h">init</td>
    </tr></thead>
    <tbody><xsl:for-each select="vars/v"><tr>
     <td><xsl:value-of select="@label"/></td>
     <td><xsl:value-of select="@init"/></td>
    </tr></xsl:for-each></tbody>
   </table></td>

   <!-- bcs -->
   <td><table>
    <thead><tr>
     <td class="h">label</td>
     <td class="h">zone</td>
     <td class="h">type</td>
     <td class="h">option</td>
     <td class="h">values</td>
    </tr></thead>
    <tbody><xsl:for-each select="bcs/bc"><tr>
     <td><xsl:value-of select="@label"/></td>
     <td><xsl:value-of select="@zone"/></td>
     <td><xsl:value-of select="@type"/></td>
     <td><xsl:value-of select="@option"/></td>
     <td><xsl:value-of select="@values"/></td>
    </tr></xsl:for-each></tbody>
   </table></td>

   <!-- mitremassembler -->
   <td><table>
    <tr><td><xsl:for-each select="mitremassembler">
     <table><tbody>
      <tr><td class="H" colspan="3">mitremassembler</td></tr>
      <tr><td class="H" rowspan="2">MITReM</td>
          <td class="h">file </td><td><xsl:value-of select="MITReM/@file" /></td></tr>
      <tr><td class="h">label</td><td><xsl:value-of select="MITReM/@label"/></td></tr>
      <tr><td class="H" rowspan="9">Element<br/>Matrix<br/>Assembler</td>
          <td class="h">convectionScheme    </td><td><xsl:value-of select="ElementMatrixAssembler/@convectionScheme"    /></td></tr>
      <tr><td class="h">diffusionScheme     </td><td><xsl:value-of select="ElementMatrixAssembler/@diffusionScheme"     /></td></tr>
      <tr><td class="h">migrationScheme     </td><td><xsl:value-of select="ElementMatrixAssembler/@migrationScheme"     /></td></tr>
      <tr><td class="h">magneticScheme      </td><td><xsl:value-of select="ElementMatrixAssembler/@magneticScheme"      /></td></tr>
      <tr><td class="h">homReactionScheme   </td><td><xsl:value-of select="ElementMatrixAssembler/@homReactionScheme"   /></td></tr>
      <tr><td class="h">electrostaticsScheme</td><td><xsl:value-of select="ElementMatrixAssembler/@electrostaticsScheme"/></td></tr>
      <tr><td class="h">timeScheme          </td><td><xsl:value-of select="ElementMatrixAssembler/@timeScheme"          /></td></tr>
      <tr><td class="h">elecReactionScheme  </td><td><xsl:value-of select="ElementMatrixAssembler/@elecReactionScheme"  /></td></tr>
      <tr><td class="h">gasReactionScheme   </td><td><xsl:value-of select="ElementMatrixAssembler/@gasReactionScheme"   /></td></tr>
      <tr><td class="H" rowspan="6">ls</td>
          <td class="h">type           </td><td><xsl:value-of select="ls/@type"      /></td></tr>
      <!-- aztec options -->
      <tr><td class="h">mtype          </td><td><xsl:value-of select="ls/@mtype"     /></td></tr>
      <tr><td class="h">solver         </td><td><xsl:value-of select="ls/@solver"    /></td></tr>
      <tr><td class="h">precond        </td><td><xsl:value-of select="ls/@precond"   /></td></tr>
      <tr><td class="h">subdomain_solve</td><td><xsl:value-of select="ls/@subdomain_solve"/></td></tr>
      <tr><td class="h" colspan="2" align="center">(...)</td></tr><!-- skip parameters, otherwise it looks like the bible-->
      <tr><td class="H">iterinit</td><td class="h">value</td><td><xsl:value-of select="iterinit/@value"/></td></tr>
      <tr><td class="H">linrelx </td><td class="h">value</td><td><xsl:value-of select="linrelx/@value" /></td></tr>
      <xsl:for-each select="bcs/bc">
       <tr><td class="H" rowspan="6">bcs/bc</td>
           <td class="h">type          </td><td><xsl:value-of select="@type"          /></td></tr>
       <tr><td class="h">label         </td><td><xsl:value-of select="@label"         /></td></tr>
       <tr><td class="h">zone          </td><td><xsl:value-of select="@zone"          /></td></tr>
       <tr><td class="h">metalpotential</td><td><xsl:value-of select="@metalpotential"/></td></tr>
       <tr><td class="h">elecreaction  </td><td><xsl:value-of select="@elecreaction"  /></td></tr>
       <tr><td class="h">gasreaction   </td><td><xsl:value-of select="@gasreaction"   /></td></tr>
      </xsl:for-each>
     </tbody></table>
    </xsl:for-each></td></tr>
   </table></td>

  </tr></xsl:for-each>
 </tbody>
</table>

 </body>
</html>

