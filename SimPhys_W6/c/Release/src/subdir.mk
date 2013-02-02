################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/IsingMatrix.cpp \
../src/SimPhys_W6_C.cpp 

OBJS += \
./src/IsingMatrix.o \
./src/SimPhys_W6_C.o 

CPP_DEPS += \
./src/IsingMatrix.d \
./src/SimPhys_W6_C.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Ofast -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


