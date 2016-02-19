#ifndef _HASH_TABLE_
#define _HASH_TABLE_

#include<iostream>
#include<stdio.h>
#include<time.h>
#include<fstream>
#include<string>
#include<math.h>
//#include"main.h"
#include"dataset.h"
using namespace std;


#define LEN sizeof(struct k_mer_node)
#define OLD_KMER_NODE 1
#define NEW_KMER_NODE 0

#define max_len_head_hash 256
typedef struct times_and_flag
{
	int times;
	int flag;
}times_and_flag;

//extern int count_of_node[num_read][int(len_read*1.1)];//节点计数数组

//extern int get_index_num(char *ch);
//extern head_node_index head_of_hash_table[max_len_head_hash];//哈希表头节点数组
//extern void output_hash_table();
extern void hash_table();


#endif