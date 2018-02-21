#include <stdio.h>
int main()
{
	int i=0, j=1;

	while(i<8) {
		printf("i %d\n",i);
		if(j==1) {
			i++;
			continue;
		}
		
		i++;
	}

}

