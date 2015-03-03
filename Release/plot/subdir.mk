################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../plot/gnufor2.f90 

OBJS += \
./plot/gnufor2.o 


# Each subdirectory must supply rules for building sources it contributes
plot/%.o: ../plot/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

plot/gnufor2.o: ../plot/gnufor2.f90


