# MultiSpeciesSpatialEcoModel.jl 

## Overview

This is an implementation of several models that represent individuals as a discrete units spatially located in a two dimensional grid 
The most simple spatially explicit model in the context of ecology, is the contact process which is a straithforward implementation of the logistic model in 2D @Gastner2009

## The multiple species contact model

An natural extension of the one species contact model is the multiple species contact model, which automatically adds competition for space to the system. Considering $n$ species $S_i$ the processes of birth and death could be represented by the following transitions:

$S_i + Ø \xrightarrow{\lambda_i} S_i + S_i$

$S_i \xrightarrow{\delta_i} Ø$





# References

@Martin2020
@Black2012 
@Neuhauser1998 

@Levine2017

@Keymer2000
