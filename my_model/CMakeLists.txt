################################################################################
# Aerosol library helper function

include(model_config)

################################################################################
# Aerosol model build configuration

create_aerosol_model(NAME "aero oslo"      # String a host model will use to find your model
                     TYPE aero_oslo_t      # Derived type describing your model
                     F90_MODULE aero_oslo_optics  # Fortran module where your derived type is defined
                     LIBRARY aero_oslo_optics     # Name of the library your model will be built into
                     SOURCE aero_oslo_optics.F90 aero_oslo/cam_logfile.F90 aero_oslo/commondefinitions.F90 aero_oslo/opttab.F90 aero_oslo/opttab_lw.F90 aero_oslo/inputForInterpol.F90 aero_oslo/pmxsub_light.F90 aero_oslo/optinterpol.F90 aero_oslo/shr_kind_mod.F90 aero_oslo/get_path.F90) # Space-separated list of all your aerosol model source files

################################################################################
