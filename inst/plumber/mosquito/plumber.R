library(plumber)

#* @apiTitle API for mosquito-only simulation

#* Setup global parameters
#* @get /config_global_parameters
#* @param path a file path to the global config JSON file
MicroMoB:::api_setup_global_parameters

#* Setup aquatic component parameters
#* @get /config_aqua_parameters
MicroMoB:::api_setup_aqua_parameters

#* Setup adult component parameters
#* @get /config_adult_parameters
MicroMoB:::api_setup_adult_parameters

#* Get global parameters
#* @get /parameters_global
#* @serializer json list(pretty = TRUE)
MicroMoB:::api_get_parameters_global

#* Get adult component parameters
#* @get /parameters_adult
#* @serializer json list(pretty = TRUE)
MicroMoB:::api_get_parameters_adult

#* Get aquatic component parameters
#* @get /parameters_aqua
#* @serializer json list(pretty = TRUE)
MicroMoB:::api_get_parameters_aqua

#* Setup model object
#* @get /setup_model
MicroMoB:::api_setup_model_object

#* Setup aquatic component
#* @get /setup_aqua
MicroMoB:::api_setup_aqua

#* Setup adult component
#* @get /setup_adult
MicroMoB:::api_setup_adult

#* Update (step) the aquatic component
#* @get /step_aqua
MicroMoB:::api_step_aqua

#* Update (step) the adult mosquito component
#* @get /step_adult
MicroMoB:::api_step_adult

#* Get output from the aquatic (immature) mosquito component
#* @get /output_aqua
#* @serializer csv
MicroMoB:::api_get_output_aqua

#* Get output from the adult mosquito component
#* @get /output_adult
#* @serializer csv
MicroMoB:::api_get_output_adult

#* Increase clock by one time step
#* @get /clock_tick
MicroMoB:::api_clock_tick

