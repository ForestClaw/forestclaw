# Copyright 2015 Lorenz Hüdepohl
#
# This file is part of fdep and licensed under the MIT license
# see the file LICENSE for more information
#

define translate_name
$(subst /,_,$(subst -,_,$(subst .,_,$(subst ./,,$1))))
endef

_f90_verbose = $(_f90_verbose_$(V))
_f90_verbose_ = $(_f90_verbose_$(AM_DEFAULT_VERBOSITY))
_f90_verbose_0 = @echo "  $1";
_f90_only_verbose = $(_f90_only_verbose_$(V))
_f90_only_verbose_ = @
_f90_only_verbose_0 = @
_f90_targets = $(call translate_name,$(PROGRAMS) $(LTLIBRARIES))

FORTRAN_CPP ?= cpp -P -traditional -Wall -Werror

# $1 source files
#
# returns: file without any .F90 .f90 .F .f extension
define strip_fortran_ext
$(patsubst %.F90,%,$(patsubst %.f90,%,$(patsubst %.F,%,$(patsubst %.f,%,$1))))
endef

# $1 program
#
# returns:
#  '1' if object files for target $1 are prefixed due to 'per-target' flags,
#  '' (the empty string) otherwise. See the automake manual for 'per-target'
#  compilation
#
define is_per_target
$(if $(filter $(call strip_fortran_ext,$(firstword $(call fortran_sources,$1))),$(patsubst %.o,%,$(patsubst %.lo,%,$($1_OBJECTS)))),,1)
endef

# $1 top-level target name (i.e. an entry of _f90_targets)
#
# returns: all target source files matching *.F90 *.f90 *.F *.f
define fortran_sources
$(filter %.F90 %.f90 %.F %.f,$($1_SOURCES))
endef

# $1 top-level target name
#
# returns: the appropriate extension (i.e. 'o' for normal programs, '.lo' for libraries)
define object_extension
$(if $(filter $1,$(PROGRAMS)),o,o)
endef

# $1 source file
# $2 stem
# $3 program
# $4 kind of file ('use' or 'def')
define modinfo_name
$(dir $1)$(2)$(call strip_fortran_ext,$(notdir $1)).$4_mods_$(patsubst .,_,$3).$(call object_extension,$3)
endef

# $1 source_file
# $2 stem
# $3 program
define module_targets
$(eval _$(3)_use_mods += $(call modinfo_name,$1,$2,$3,use))
$(call modinfo_name,$1,$2,$3,use): $1 $(dir $1)$(am__dirstamp)
	$(call _f90_verbose,F90 USE  [$3] $$<)$(FORTRAN_CPP) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $($p_CPPFLAGS) $(CPPFLAGS) -o /dev/stdout $$< | \
		grep -i -o '^ *use [^ ,!:]*' | sed 's/^[[:space:]]*//;' | tr '[:upper:]' '[:lower:]' | sort -u > $$@

$(eval _$(3)_def_mods += $(call modinfo_name,$1,$2,$3,def))
$(call modinfo_name,$1,$2,$3,def): $1 $(dir $1)$(am__dirstamp)
	$(call _f90_verbose,F90 MOD  [$3] $$<)$(FORTRAN_CPP) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $($p_CPPFLAGS) $(CPPFLAGS) -o /dev/stdout $$< | \
		grep -i -o '^ *module [^!]*' | sed 's/^[[:space:]]*//;' | tr '[:upper:]' '[:lower:]' | grep -v "\<procedure\>\|\<intrinsic\>" > $$@ || true

endef
$(foreach p,$(_f90_targets),$(if $(call is_per_target,$p),$(foreach s,$(call fortran_sources,$p),$(eval $(call module_targets,$s,$p-,$p))),$(foreach s,$(call fortran_sources,$p),$(eval $(call module_targets,$s,,$p)))))

_f90_depdir=$(abs_builddir)/.fortran_dependencies
_f90_depfile = $(_f90_depdir)/dependencies.mk

# $1 target-name
define recursive_lib_deps
$(foreach l,$(call translate_name,$($1_LDADD) $($1_LIBADD)),$l $(call recursive_lib_deps,$l))
endef

define is_clean
$(if $(filter-out mostlyclean clean distclean maintainer-clean am--depfiles,$(MAKECMDGOALS)),0,1)
endef

define newline


endef

ifneq ($(call is_clean),1)
include $(_f90_depfile)
endif

# $1 program
define program_dependencies
	$(_f90_only_verbose){ $(foreach argument,$(_$p_use_mods) $(_$p_def_mods) $(foreach l,$(call recursive_lib_deps,$p),$(_$l_use_mods) $(_$l_def_mods)),echo $(argument); ) true; } | \
	$(top_srcdir)/fdep/fortran_dependencies.pl $p >> $@ || { rm $@; exit 1; }

endef

$(_f90_depfile): $(top_srcdir)/fdep/fortran_dependencies.pl $(foreach p,$(_f90_targets),$(_$p_use_mods) $(_$p_def_mods)) | $(foreach p,$(_f90_targets),$(_f90_depdir)/$p)
	$(call _f90_verbose,F90 DEPS $@)echo > $@;
	$(foreach p,$(_f90_targets),$(call program_dependencies,$p))

$(_f90_depdir):
	@mkdir $@

$(foreach p,$(_f90_targets),$(_f90_depdir)/$p): | $(_f90_depdir)
	@mkdir -p $@

CLEANFILES += $(foreach p,$(_f90_targets),$(_$p_def_mods) $(_$p_use_mods))
CLEANFILES += $(foreach p,$(_f90_targets),$(_f90_depdir)/$p/*)
CLEANFILES += $(_f90_depfile)
