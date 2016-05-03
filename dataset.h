#ifndef _DATASET_H_
#define _DATASET_H_


#include<iostream>
#include<stdio.h>
#include<time.h>
#include<fstream>
#include<string>
//#include"main.h"
using namespace std;

/*
#define len_genome 1000000     //基因组长度为10000000
#define num_read 5000           //随机起点个数为1000
#define len_read 1000            //read 的长度 1K-10K
#define short_k 7               //进行匹配的短k-mer长度
#define long_k 20                // 进行匹配的长k-mer长度
#define len_index 3              //进行匹配时的索引长度
#define len_index_error 3        //匹配k-mer时允许有的错误匹配数
#define len_array_k_mer 64         //k-mer存储数组的长度
#define len_k_mer_error 3        //匹配时容许k-mer错误
#define invalid_value  3            //无效匹配
#define valid_value  3            //有效匹配
*/

#define len_genome 4000000     //基因组长度为10000000
#define num_read 1000         //随机起点个数为1000
#define len_read 10000            //read 的长度 10-20
#define short_k 10              //进行匹配的短k-mer长度，
#define long_k 6               // 进行匹配的长k-mer长度
#define len_index 3              //进行匹配时的索引长度
#define len_index_error 3        //匹配k-mer时允许有的错误匹配数
#define len_array_k_mer 64         //k-mer存储数组的长度
#define len_k_mer_error 1        //匹配时容许k-mer错误
#define invalid_value  2            //无效匹配
#define valid_value  5//有效匹配

/*
#define len_genome 1000  //基因组长度为10000000
#define num_read 50         //随机起点个数为1000
#define len_read 100            //read 的长度 10-20
#define short_k 6              //进行匹配的短k-mer长度
#define long_k 6               // 进行匹配的长k-mer长度
#define len_index 3              //进行匹配时的索引长度
#define len_index_error 3        //匹配k-mer时允许有的错误匹配数
#define len_array_k_mer 64         //k-mer存储数组的长度
#define len_k_mer_error 0        //匹配时容许k-mer错误
#define invalid_value  2            //无效匹配
#define valid_value  2//有效匹配

*/
typedef struct k_mer
{
	//char* ch ;
	string ch;                   //short_k长度k-mer
	int count_ch;                //本k-mer出现次数
}k_mer;

extern char genome[len_genome+1];//基因组字符串
extern int position[num_read];//长read起始位置数组
extern int length[num_read];//1000条read的随机长度。1K-10K
extern char sample[num_read][int(len_read*2)];//存放长read，长度扩大为1.1倍，用来存放碱基插入错误
extern int total_num_of_bases;//生成的长read中所有碱基总数
extern k_mer array_k_mer[len_array_k_mer];
extern int num_short_k_mer; //标记短k-mer组中出现的k-mer的数量


extern void dataset();

#endif