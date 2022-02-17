library(plumber)

#* @apiTitle API for mosquito-only simulation

#* Return hello world
#* @get /hello
MicroMoB:::hello_world

#* Setup parameters for mosquito-only simulation
#* @get /config_mosquito
#* @param path a file path to the global config
MicroMoB:::put_config_mosquito

#* Setup model object and components for mosquito-only simulation
#* @get /config_model_object_mosquito
MicroMoB:::put_model_object_mosquito

#* Get parameters for adult mosquitoes for mosquito-only simulation
#* @get /parameters_adults
#* @serializer json list(pretty = TRUE)
MicroMoB:::get_parameters_adult_mosquito

#* Get parameters for aquatic (immature) mosquitoes for mosquito-only simulation
#* @get /parameters_aqua
#* @serializer json list(pretty = TRUE)
MicroMoB:::get_parameters_aqua_mosquito

#* Update (step) the aquatic (immature) mosquito component
#* @get /step_aqua
MicroMoB:::put_step_aqua

#* Update (step) the adult mosquito component
#* @get /step_adult
MicroMoB:::put_step_adult

#* Get output from the aquatic (immature) mosquito component
#* @get /output_aqua
#* @serializer csv
MicroMoB:::get_output_aqua

#* Get output from the adult mosquito component
#* @get /output_adult
#* @serializer csv
MicroMoB:::get_output_adult

#* Increase clock by one time step
#* @get /clock_tick
MicroMoB:::clock_tick
