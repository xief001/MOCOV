#ifndef _CORRECT_ERROES_
#define _CORRECT_ERRORS_

typedef struct com
{
	int flag;
	int len;
}com;

typedef struct rawbase
{
	char base;
	int new_stat;						//-1   --->   insertion
										//-2   --->   deletion
										//-3   --->   complex
}rawbase;


extern void correct_errors();



#endif