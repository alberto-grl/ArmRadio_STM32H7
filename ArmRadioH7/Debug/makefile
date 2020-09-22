################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include Drivers/STM32H7xx_HAL_Driver/Src/subdir.mk
-include Core/Startup/subdir.mk
-include Core/Src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 
EXECUTABLES += \
ArmRadioH7.elf \

SIZE_OUTPUT += \
default.size.stdout \

OBJDUMP_LIST += \
ArmRadioH7.list \

OBJCOPY_BIN += \
ArmRadioH7.bin \


# All Target
all: ArmRadioH7.elf secondary-outputs

# Tool invocations
ArmRadioH7.elf: $(OBJS) $(USER_OBJS) C:\Users\alberto\STM32CubeIDE\workspace_1.3.0\ArmRadio_STM32H7\ArmRadioH7\STM32H743ZITX_FLASH-DTCM.ld
	arm-none-eabi-gcc -o "ArmRadioH7.elf" @"objects.list" $(USER_OBJS) $(LIBS) -mcpu=cortex-m7 -T"C:\Users\alberto\STM32CubeIDE\workspace_1.3.0\ArmRadio_STM32H7\ArmRadioH7\STM32H743ZITX_FLASH-DTCM.ld" --specs=nosys.specs -Wl,-Map="ArmRadioH7.map" -Wl,--gc-sections -static -L"C:\Users\alberto\STM32Cube\Repository\STM32Cube_FW_H7_V1.7.0\Drivers\CMSIS/Lib/GCC" --specs=nano.specs -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb -u _printf_float -Wl,--start-group -lc -lm -Wl,--end-group
	@echo 'Finished building target: $@'
	@echo ' '

default.size.stdout: $(EXECUTABLES)
	arm-none-eabi-size  $(EXECUTABLES)
	@echo 'Finished building: $@'
	@echo ' '

ArmRadioH7.list: $(EXECUTABLES)
	arm-none-eabi-objdump -h -S $(EXECUTABLES) > "ArmRadioH7.list"
	@echo 'Finished building: $@'
	@echo ' '

ArmRadioH7.bin: $(EXECUTABLES)
	arm-none-eabi-objcopy  -O binary $(EXECUTABLES) "ArmRadioH7.bin"
	@echo 'Finished building: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) *
	-@echo ' '

secondary-outputs: $(SIZE_OUTPUT) $(OBJDUMP_LIST) $(OBJCOPY_BIN)

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
