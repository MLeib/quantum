################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F03_SRCS += \
../testrun.f03 

OBJS += \
./testrun.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f03
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

testrun.o: ../testrun.f03 LinAlg/IR_Precision.o mps/mpo_node_module.o


