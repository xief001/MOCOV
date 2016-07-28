#ifndef _CORRECT_ERROES_
#define _CORRECT_ERRORS_
#include"main.h"
extern char **corrected_read ;
extern int *corrected_len;
extern int **insertion;
extern int **deletion;



typedef struct rawbase
{
	char base;
	int new_stat;						//-1   --->   insertion
										//-2   --->   deletion
										//-3   --->   complex
}rawbase;
extern void correct_errors();
#endif