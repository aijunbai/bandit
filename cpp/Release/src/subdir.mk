################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/agent.cpp \
../src/bandit.cpp \
../src/main.cpp \
../src/random.cpp 

CPP_DEPS += \
./src/agent.d \
./src/bandit.d \
./src/main.d \
./src/random.d 

OBJS += \
./src/agent.o \
./src/bandit.o \
./src/main.o \
./src/random.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


