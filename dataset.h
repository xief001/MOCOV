#ifndef _DATASET_H_
#define _DATASET_H_



#include"main.h"

using namespace std;


/*
#define len_genome 10000      //genome size 10000000
#define num_read 500           //random read position 1000
#define len_read 1000            //read length 1K-10K
#define short_k 11				//k-mer length for short aligning
#define long_k 20                // k-mer length for checking correction
#define len_index 3              //index length for comparing k-mers in hash table with k-mers in read 
#define len_index_error 2        //hamming distance between k-mers' index in hash table with k-mers in read
#define len_k_mer_error 0        //hamming distance between k-mers in hash table with k-mers in read
#define invalid_value  3         //if match value < invalid_value ,the match is invalid
#define valid_value  3           //if match value > valid_value ,the match is valid
#define count_value 2
*/
/*
#define len_genome 4000000       //genome size 10000000
#define num_read 1000            //random read position 1000
#define len_read 10000           //read length 1000-10000
#define short_k 10               //k-mer length for short aligning
#define long_k 6                 //k-mer length for checking correction
#define len_index 3              //index length for comparing k-mers in hash table with k-mers in read 
#define len_index_error 3        //hamming distance between k-mers' index in hash table with k-mers in read
#define len_k_mer_error 1        //hamming distance between k-mers in hash table with k-mers in read
#define invalid_value  2         //if match value < invalid_value ,the match is invalid
#define valid_value  5           //if match value > valid_value ,the match is valid
*/
#define len_genome 1000          //genome size 1000
#define num_read 500              //random read position 50
#define len_read 100             //read length 10-100
#define short_k 10                //k-mer length for short aligning
#define long_k 20                 //k-mer length for checking correction
#define medium_k short_k+1
//#define medium_k short_k+short_k/4
#define len_index 3              //index length for comparing k-mers in hash table with k-mers in read 
#define len_index_error 2        //hamming distance between k-mers' index in hash table with k-mers in read
#define len_k_mer_error 0        //hamming distance between k-mers in hash table with k-mers in read
#define invalid_value  3         //if match value < invalid_value ,the match is invalid
#define valid_value  3           //if match value > valid_value ,the match is valid
#define count_value 2      
//#define medium_len short_k/2

/*
typedef struct k_mer
{
	string ch;                   
	int count_ch;                //the k-mer appearing times
}k_mer;
*/
extern char genome[len_genome+1];                  //genome
extern int position[num_read];                     //array for long read starting positions
extern int length[num_read];                       //array for every long read's length
extern char sample[num_read][int(len_read*2)];     //array for store every long read with error
extern int total_num_of_bases;                     //total num of bases in sample
extern int num_short_k_mer;                        //number of k-mers for short match


extern void dataset();

#endif