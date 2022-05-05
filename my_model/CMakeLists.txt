################################################################################
# Aerosol library helper function

include(model_config)

################################################################################
# Aerosol model build configuration

create_aerosol_model(NAME "oslo aero"      # String a host model will use to find your model
                     TYPE oslo_aero_t      # Derived type describing your model
                     F90_MODULE oslo_aero_optics  # Fortran module where your derived type is defined
                     LIBRARY oslo_aero_optics     # Name of the library your model will be built into
                     SOURCE my_model.F90 cam_logfile.F90  ) # Space-separated list of all your aerosol model source files

################################################################################
