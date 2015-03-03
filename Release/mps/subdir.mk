################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F03_SRCS += \
../mps/mpo_module.f03 \
../mps/mpo_node_module.f03 \
../mps/mps_class.f03 \
../mps/mps_node.f03 \
../mps/mps_parameters.f03 \
../mps/mps_test.f03 

OBJS += \
./mps/mpo_module.o \
./mps/mpo_node_module.o \
./mps/mps_class.o \
./mps/mps_node.o \
./mps/mps_parameters.o \
./mps/mps_test.o 


# Each subdirectory must supply rules for building sources it contributes
mps/%.o: ../mps/%.f03
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

mps/mpo_module.o: ../mps/mpo_module.f03 mps/mpo_node_module.o

mps/mpo_node_module.o: ../mps/mpo_node_module.f03 LinAlg/IR_Precision.o LinAlg/lapack_wrap.o

mps/mps_class.o: ../mps/mps_class.f03 LinAlg/IR_Precision.o mps/mps_node.o mps/mps_parameters.o

mps/mps_node.o: ../mps/mps_node.f03 LinAlg/IR_Precision.o LinAlg/lapack_wrap.o

mps/mps_parameters.o: ../mps/mps_parameters.f03

mps/mps_test.o: ../mps/mps_test.f03 LinAlg/IR_Precision.o mps/mps_class.o


