#include"main.h"
#include"hash_table.h"
#include"dataset.h"
#include"find_errors.h"
#include"correct_errors.h"
#include<iostream>
#include<stdio.h>
#include<time.h>
#include<fstream>
#include<string>

using namespace std;




int main()
{ 
	while (1)
	{
		dataset();
		hash_table();
		//find_error();
		correct_errors();

		getchar();
	}
	
   

	

return 0;
}

