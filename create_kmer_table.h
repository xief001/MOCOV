#ifndef _HASH_TABLE_
#define _HASH_TABLE_


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

typedef struct k_mer_node     
{
	int i_read;
	int j_position;
	struct k_mer_node *next;
}k_mer_node;

typedef struct head_node_index   
{
	char *head_hash_table;
	struct k_mer_node *next;
}head_node_index;

extern void hash_table();


#endif