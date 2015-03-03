################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F03_SRCS += \
../LinAlg/LinAlg_test.f03 \
../LinAlg/lapack_wrap.f03 \
../LinAlg/random.f03 

F90_SRCS += \
../LinAlg/IR_Precision.f90 

OBJS += \
./LinAlg/IR_Precision.o \
./LinAlg/LinAlg_test.o \
./LinAlg/lapack_wrap.o \
./LinAlg/random.o 


# Each subdirectory must supply rules for building sources it contributes
LinAlg/IR_Precision.o: ../LinAlg/IR_Precision.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -cpp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

LinAlg/IR_Precision.o: ../LinAlg/IR_Precision.f90

LinAlg/%.o: ../LinAlg/%.f03
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

LinAlg/LinAlg_test.o: ../LinAlg/LinAlg_test.f03 LinAlg/IR_Precision.o LinAlg/lapack_wrap.o LinAlg/random.o plot/gnufor2.o

LinAlg/lapack_wrap.o: ../LinAlg/lapack_wrap.f03 LinAlg/IR_Precision.o

LinAlg/random.o: ../LinAlg/random.f03 LinAlg/IR_Precision.o


