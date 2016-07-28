#ifndef _CREATE_KMER_TABLE_
#define _CREATE_KMER_TABLE_

using namespace std;

#define LEN sizeof(k_mer_node)
#define OLD_KMER_NODE 1
#define NEW_KMER_NODE 0

typedef struct times_and_flag
{
	int times;
	int flag;
}times_and_flag;

typedef struct k_mer_node     
{
	int i_read;
	int j_position;
	struct k_mer_node *next;
}k_mer_node;

typedef struct head_node_index   
{
	string head_hash_table;
	struct k_mer_node *next;
}head_node_index;

typedef struct node_position                                            //the position of kmer's first appearance
{
	int i_position;
	int j_position;
}node_position;



extern node_position **every_short_kmer_times_position;          
extern node_position **every_long_kmer_times_position;          
extern node_position **every_medium_kmer_times_position;       
extern void free_Array_node_position(node_position **p, int row);

extern int **count_of_short_node;                           //array of counting every kmer (in read) times
extern int **count_of_long_node;                           //array of counting every kmer (in read) times
extern int **count_of_medium_node;                           //array of counting every kmer (in read) times
extern int ** get_Array_int(int row, int line);
extern void free_Array_int(int **p, int row);

extern int get_index_num(char *ch);                                             //check the first len_index bases of kmer len_index to get its index_num£¨AA---0,AT---1,AC---2,AG---3...£©


extern void hash_table();


#endif