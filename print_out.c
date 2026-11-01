#include "print_out.h"

#define ITM_PORT         (*(volatile unsigned char *)(0xe0000000))
#define DEMCR_           (*(volatile unsigned long *)(0xE000EDFC))
#define TRCENA_          0X01000000

uint8_t fputcp(char ch)
{
	if (DEMCR_ & TRCENA_)
	{
		while (ITM_PORT == 0);
		ITM_PORT = ch;
	}
	else
	{
		return 0;
	}
	return 1;
}

uint8_t print_str(char *str)
{
	uint32_t i;
	i = 0;
	while (1)
	{
		if (*(str + i) == '\0')
		{
			return 1;
		}
		if (fputcp(*(str + i)) == 0)
		{
			return 0;
		}
		i++;
	}
}
