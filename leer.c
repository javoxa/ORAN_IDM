#include <stdio.h>
#include <stdlib.h>

int main()
{
	printf("Ejecuta en R un grafico\n");
	system("echo \"source(\'grafico.txt\')\"| r --vanilla");
	exit(0);
}
