# plumber API for mosquito-only simulation

#* Setup parameters for mosquito-only simulation
#* @put /config_mosquito
MicroMoB:::put_config_mosquito

#* Get parameters for adult mosquitoes for mosquito-only simulation
#* @get /parameters_adults
#* @serializer json list(pretty = TRUE)
MicroMoB:::get_parameters_adult_mosquito

#* Get parameters for aquatic (immature) mosquitoes for mosquito-only simulation
#* @get /parameters_aqua
#* @serializer json list(pretty = TRUE)
MicroMoB:::get_parameters_aqua_mosquito

#* Setup model object and components for mosquito-only simulation
#* @put /config_model_object_mosquito
MicroMoB:::put_model_object_mosquito
