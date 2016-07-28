#ifndef _MAIN_H_
#define _MAIN_H_

#include"main.h"
#include"create_kmer_table.h"
#include"correct_errors.h"

#include<iostream>
#include<stdio.h>
#include<time.h>
#include<fstream>
#include<string>
#include<math.h>
#include<cstdlib>
#include<string.h>
#include<malloc.h>

//==================================================================================================
//read_count is defined as 200000 , it is used to define other arraies
//num_read is up to input data , it is used to control program
//==================================================================================================

extern int  read_count ;
extern int  len_read ;         //read length 1000-10000

extern int *length;                       //array for every long read's length
extern int *name_length;

extern char **sample;
extern char **readname;
extern int num_read;


extern string query_file_name;
extern string corrected_file_name;

extern char ** get_char_Array(int row, int line);

extern int short_k;//-s
extern int long_k;//-l
extern int medium_k;//-l

extern int len_index;//-i
extern int valid_value;//v
extern int len_k_mer_error;//e

#endif 