#include"create_kmer_table.h"
#include"dataset.h"
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<malloc.h>

int index_num;
extern int length[num_read];

times_and_flag short_k_mer_times[num_read][int(len_read*2)];
times_and_flag long_k_mer_times[num_read][int(len_read*2)];
times_and_flag medium_k_mer_times[num_read][int(len_read*2)];

times_and_flag short_base_score[num_read][int(len_read*2)];
times_and_flag long_base_score[num_read][int(len_read*2)];
times_and_flag medium_base_score[num_read][int(len_read*2)];

int times[num_read][int(len_read*2)]={0};

head_node_index head_of_short_hash_table[max_len_head_hash];                  //hash table's head node
head_node_index head_of_long_hash_table[max_len_head_hash];                  //hash table's head node
head_node_index head_of_medium_hash_table[max_len_head_hash];                  //hash table's head node

int len_short_hash_table[max_len_head_hash]={0};                                  //number of hash table's head node
int len_long_hash_table[max_len_head_hash]={0};                                  //number of hash table's head node
int len_medium_hash_table[max_len_head_hash]={0};                                  //number of hash table's head node

int count_of_short_node[num_read][int(len_read*2)]={0};                           //array of counting every kmer (in read) times
int count_of_long_node[num_read][int(len_read*2)]={0};                           //array of counting every kmer (in read) times
int count_of_medium_node[num_read][int(len_read*2)]={0};                           //array of counting every kmer (in read) times

typedef struct node_position                                            //the position of kmer's first appearance
{
	int i_position;
	int j_position;
}node_position;

node_position every_short_kmer_times_position[num_read][len_read*2];          
node_position every_long_kmer_times_position[num_read][len_read*2];          
node_position every_medium_kmer_times_position[num_read][len_read*2];          

void init_k_mer_times()
{
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i];j++)
		{
			short_k_mer_times[i][j].times=0;
			short_k_mer_times[i][j].flag=0;
			long_k_mer_times[i][j].times=0;
			long_k_mer_times[i][j].flag=0;
			medium_k_mer_times[i][j].times=0;
			medium_k_mer_times[i][j].flag=0;
		}
	}
}

void init_every_kmer_times_position()
{
	for(int i = 0 ;i<num_read;i++)
	{
		for(int j= 0; j< length[i];j++)
		{
			every_short_kmer_times_position[i][j].i_position=0;
			every_short_kmer_times_position[i][j].j_position=0;
			every_long_kmer_times_position[i][j].i_position=0;
			every_long_kmer_times_position[i][j].j_position=0;
			every_medium_kmer_times_position[i][j].i_position=0;
			every_medium_kmer_times_position[i][j].j_position=0;
		}
	}
}

void init_short_hash_table()
{
	int real_len_head_hash;
	real_len_head_hash=int(pow(4,len_index));
	if(len_index==1)
	{
		head_of_short_hash_table[0].head_hash_table="A";head_of_short_hash_table[0].next=NULL;
		head_of_short_hash_table[1].head_hash_table="T";head_of_short_hash_table[1].next=NULL;
		head_of_short_hash_table[2].head_hash_table="C";head_of_short_hash_table[2].next=NULL;
		head_of_short_hash_table[3].head_hash_table="G";head_of_short_hash_table[3].next=NULL;
	}
	if(len_index==2)
	{
		head_of_short_hash_table[0].head_hash_table="AA";head_of_short_hash_table[0].next=NULL;//0
		head_of_short_hash_table[1].head_hash_table="AT";head_of_short_hash_table[1].next=NULL;//1
		head_of_short_hash_table[2].head_hash_table="AC";head_of_short_hash_table[2].next=NULL;//2
		head_of_short_hash_table[3].head_hash_table="AG";head_of_short_hash_table[3].next=NULL;//3
		head_of_short_hash_table[4].head_hash_table="TA";head_of_short_hash_table[4].next=NULL;//4
		head_of_short_hash_table[5].head_hash_table="TT";head_of_short_hash_table[5].next=NULL;//5
		head_of_short_hash_table[6].head_hash_table="TC";head_of_short_hash_table[6].next=NULL;//6
		head_of_short_hash_table[7].head_hash_table="TG";head_of_short_hash_table[7].next=NULL;//7
		head_of_short_hash_table[8].head_hash_table="CA";head_of_short_hash_table[8].next=NULL;//8
		head_of_short_hash_table[9].head_hash_table="CT";head_of_short_hash_table[9].next=NULL;//9
		head_of_short_hash_table[10].head_hash_table="CC";head_of_short_hash_table[10].next=NULL;//10
		head_of_short_hash_table[11].head_hash_table="CG";head_of_short_hash_table[11].next=NULL;//11
		head_of_short_hash_table[12].head_hash_table="GA";head_of_short_hash_table[12].next=NULL;//12
		head_of_short_hash_table[13].head_hash_table="GT";head_of_short_hash_table[13].next=NULL;//13
		head_of_short_hash_table[14].head_hash_table="GC";head_of_short_hash_table[14].next=NULL;//14
		head_of_short_hash_table[15].head_hash_table="GG";head_of_short_hash_table[15].next=NULL;//15
	}
	if(len_index==3)
	{
		head_of_short_hash_table[0].head_hash_table="AAA";head_of_short_hash_table[0].next=NULL;//0
		head_of_short_hash_table[1].head_hash_table="AAT";head_of_short_hash_table[1].next=NULL;//1
		head_of_short_hash_table[2].head_hash_table="AAC";head_of_short_hash_table[2].next=NULL;//2
		head_of_short_hash_table[3].head_hash_table="AAG";head_of_short_hash_table[3].next=NULL;//3
		head_of_short_hash_table[4].head_hash_table="ATA";head_of_short_hash_table[4].next=NULL;//4
		head_of_short_hash_table[5].head_hash_table="ATT";head_of_short_hash_table[5].next=NULL;//5
		head_of_short_hash_table[6].head_hash_table="ATC";head_of_short_hash_table[6].next=NULL;//6
		head_of_short_hash_table[7].head_hash_table="ATG";head_of_short_hash_table[7].next=NULL;//7
		head_of_short_hash_table[8].head_hash_table="ACA";head_of_short_hash_table[8].next=NULL;//8
		head_of_short_hash_table[9].head_hash_table="ACT";head_of_short_hash_table[9].next=NULL;//9
		head_of_short_hash_table[10].head_hash_table="ACC";head_of_short_hash_table[10].next=NULL;//10
		head_of_short_hash_table[11].head_hash_table="ACG";head_of_short_hash_table[11].next=NULL;//11
		head_of_short_hash_table[12].head_hash_table="AGA";head_of_short_hash_table[12].next=NULL;//12
		head_of_short_hash_table[13].head_hash_table="AGT";head_of_short_hash_table[13].next=NULL;//13
		head_of_short_hash_table[14].head_hash_table="AGC";head_of_short_hash_table[14].next=NULL;//14
		head_of_short_hash_table[15].head_hash_table="AGG";head_of_short_hash_table[15].next=NULL;//15

		head_of_short_hash_table[16].head_hash_table="TAA";head_of_short_hash_table[16].next=NULL;//16
		head_of_short_hash_table[17].head_hash_table="TAT";head_of_short_hash_table[17].next=NULL;//17
		head_of_short_hash_table[18].head_hash_table="TAC";head_of_short_hash_table[18].next=NULL;//18
		head_of_short_hash_table[19].head_hash_table="TAG";head_of_short_hash_table[19].next=NULL;//19
		head_of_short_hash_table[20].head_hash_table="TTA";head_of_short_hash_table[20].next=NULL;//20
		head_of_short_hash_table[21].head_hash_table="TTT";head_of_short_hash_table[21].next=NULL;//21
		head_of_short_hash_table[22].head_hash_table="TTC";head_of_short_hash_table[22].next=NULL;//22
		head_of_short_hash_table[23].head_hash_table="TTG";head_of_short_hash_table[23].next=NULL;//23
		head_of_short_hash_table[24].head_hash_table="TCA";head_of_short_hash_table[24].next=NULL;//24
		head_of_short_hash_table[25].head_hash_table="TCT";head_of_short_hash_table[25].next=NULL;//25
		head_of_short_hash_table[26].head_hash_table="TCC";head_of_short_hash_table[26].next=NULL;//26
		head_of_short_hash_table[27].head_hash_table="TCG";head_of_short_hash_table[27].next=NULL;//27
		head_of_short_hash_table[28].head_hash_table="TGA";head_of_short_hash_table[28].next=NULL;//28
		head_of_short_hash_table[29].head_hash_table="TGT";head_of_short_hash_table[29].next=NULL;//29
		head_of_short_hash_table[30].head_hash_table="TGC";head_of_short_hash_table[30].next=NULL;//30
		head_of_short_hash_table[31].head_hash_table="TGG";head_of_short_hash_table[31].next=NULL;//31

		head_of_short_hash_table[32].head_hash_table="CAA";head_of_short_hash_table[32].next=NULL;//32
		head_of_short_hash_table[33].head_hash_table="CAT";head_of_short_hash_table[33].next=NULL;//33
		head_of_short_hash_table[34].head_hash_table="CAC";head_of_short_hash_table[34].next=NULL;//34
		head_of_short_hash_table[35].head_hash_table="CAG";head_of_short_hash_table[35].next=NULL;//35
		head_of_short_hash_table[36].head_hash_table="CTA";head_of_short_hash_table[36].next=NULL;//36
		head_of_short_hash_table[37].head_hash_table="CTT";head_of_short_hash_table[37].next=NULL;//37
		head_of_short_hash_table[38].head_hash_table="CTC";head_of_short_hash_table[38].next=NULL;//38
		head_of_short_hash_table[39].head_hash_table="CTG";head_of_short_hash_table[39].next=NULL;//39
		head_of_short_hash_table[40].head_hash_table="CCA";head_of_short_hash_table[40].next=NULL;//40
		head_of_short_hash_table[41].head_hash_table="CCT";head_of_short_hash_table[41].next=NULL;//41
		head_of_short_hash_table[42].head_hash_table="CCC";head_of_short_hash_table[42].next=NULL;//42
		head_of_short_hash_table[43].head_hash_table="CCG";head_of_short_hash_table[43].next=NULL;//43
		head_of_short_hash_table[44].head_hash_table="CGA";head_of_short_hash_table[44].next=NULL;//44
		head_of_short_hash_table[45].head_hash_table="CGT";head_of_short_hash_table[45].next=NULL;//45
		head_of_short_hash_table[46].head_hash_table="CGC";head_of_short_hash_table[46].next=NULL;//46
		head_of_short_hash_table[47].head_hash_table="CGG";head_of_short_hash_table[47].next=NULL;//47

		head_of_short_hash_table[48].head_hash_table="GAA";head_of_short_hash_table[48].next=NULL;//48
		head_of_short_hash_table[49].head_hash_table="GAT";head_of_short_hash_table[49].next=NULL;//49
		head_of_short_hash_table[50].head_hash_table="GAC";head_of_short_hash_table[50].next=NULL;//50
		head_of_short_hash_table[51].head_hash_table="GAG";head_of_short_hash_table[51].next=NULL;//51
		head_of_short_hash_table[52].head_hash_table="GTA";head_of_short_hash_table[52].next=NULL;//52
		head_of_short_hash_table[53].head_hash_table="GTT";head_of_short_hash_table[53].next=NULL;//53
		head_of_short_hash_table[54].head_hash_table="GTC";head_of_short_hash_table[54].next=NULL;//54
		head_of_short_hash_table[55].head_hash_table="GTG";head_of_short_hash_table[55].next=NULL;//55
		head_of_short_hash_table[56].head_hash_table="GCA";head_of_short_hash_table[56].next=NULL;//56
		head_of_short_hash_table[57].head_hash_table="GCT";head_of_short_hash_table[57].next=NULL;//57
		head_of_short_hash_table[58].head_hash_table="GCC";head_of_short_hash_table[58].next=NULL;//58
		head_of_short_hash_table[59].head_hash_table="GCG";head_of_short_hash_table[59].next=NULL;//59
		head_of_short_hash_table[60].head_hash_table="GGA";head_of_short_hash_table[60].next=NULL;//60
		head_of_short_hash_table[61].head_hash_table="GGT";head_of_short_hash_table[61].next=NULL;//61
		head_of_short_hash_table[62].head_hash_table="GGC";head_of_short_hash_table[62].next=NULL;//62
		head_of_short_hash_table[63].head_hash_table="GGG";head_of_short_hash_table[63].next=NULL;//63
	}
	if(len_index==4)
	{


	}

}

void init_long_hash_table()
{
	int real_len_head_hash;
	real_len_head_hash=int(pow(4,len_index));
	if(len_index==1)
	{
		head_of_long_hash_table[0].head_hash_table="A";head_of_long_hash_table[0].next=NULL;
		head_of_long_hash_table[1].head_hash_table="T";head_of_long_hash_table[1].next=NULL;
		head_of_long_hash_table[2].head_hash_table="C";head_of_long_hash_table[2].next=NULL;
		head_of_long_hash_table[3].head_hash_table="G";head_of_long_hash_table[3].next=NULL;
	}
	if(len_index==2)
	{
		head_of_long_hash_table[0].head_hash_table="AA";head_of_long_hash_table[0].next=NULL;//0
		head_of_long_hash_table[1].head_hash_table="AT";head_of_long_hash_table[1].next=NULL;//1
		head_of_long_hash_table[2].head_hash_table="AC";head_of_long_hash_table[2].next=NULL;//2
		head_of_long_hash_table[3].head_hash_table="AG";head_of_long_hash_table[3].next=NULL;//3
		head_of_long_hash_table[4].head_hash_table="TA";head_of_long_hash_table[4].next=NULL;//4
		head_of_long_hash_table[5].head_hash_table="TT";head_of_long_hash_table[5].next=NULL;//5
		head_of_long_hash_table[6].head_hash_table="TC";head_of_long_hash_table[6].next=NULL;//6
		head_of_long_hash_table[7].head_hash_table="TG";head_of_long_hash_table[7].next=NULL;//7
		head_of_long_hash_table[8].head_hash_table="CA";head_of_long_hash_table[8].next=NULL;//8
		head_of_long_hash_table[9].head_hash_table="CT";head_of_long_hash_table[9].next=NULL;//9
		head_of_long_hash_table[10].head_hash_table="CC";head_of_long_hash_table[10].next=NULL;//10
		head_of_long_hash_table[11].head_hash_table="CG";head_of_long_hash_table[11].next=NULL;//11
		head_of_long_hash_table[12].head_hash_table="GA";head_of_long_hash_table[12].next=NULL;//12
		head_of_long_hash_table[13].head_hash_table="GT";head_of_long_hash_table[13].next=NULL;//13
		head_of_long_hash_table[14].head_hash_table="GC";head_of_long_hash_table[14].next=NULL;//14
		head_of_long_hash_table[15].head_hash_table="GG";head_of_long_hash_table[15].next=NULL;//15
	}
	if(len_index==3)
	{
		head_of_long_hash_table[0].head_hash_table="AAA";head_of_long_hash_table[0].next=NULL;//0
		head_of_long_hash_table[1].head_hash_table="AAT";head_of_long_hash_table[1].next=NULL;//1
		head_of_long_hash_table[2].head_hash_table="AAC";head_of_long_hash_table[2].next=NULL;//2
		head_of_long_hash_table[3].head_hash_table="AAG";head_of_long_hash_table[3].next=NULL;//3
		head_of_long_hash_table[4].head_hash_table="ATA";head_of_long_hash_table[4].next=NULL;//4
		head_of_long_hash_table[5].head_hash_table="ATT";head_of_long_hash_table[5].next=NULL;//5
		head_of_long_hash_table[6].head_hash_table="ATC";head_of_long_hash_table[6].next=NULL;//6
		head_of_long_hash_table[7].head_hash_table="ATG";head_of_long_hash_table[7].next=NULL;//7
		head_of_long_hash_table[8].head_hash_table="ACA";head_of_long_hash_table[8].next=NULL;//8
		head_of_long_hash_table[9].head_hash_table="ACT";head_of_long_hash_table[9].next=NULL;//9
		head_of_long_hash_table[10].head_hash_table="ACC";head_of_long_hash_table[10].next=NULL;//10
		head_of_long_hash_table[11].head_hash_table="ACG";head_of_long_hash_table[11].next=NULL;//11
		head_of_long_hash_table[12].head_hash_table="AGA";head_of_long_hash_table[12].next=NULL;//12
		head_of_long_hash_table[13].head_hash_table="AGT";head_of_long_hash_table[13].next=NULL;//13
		head_of_long_hash_table[14].head_hash_table="AGC";head_of_long_hash_table[14].next=NULL;//14
		head_of_long_hash_table[15].head_hash_table="AGG";head_of_long_hash_table[15].next=NULL;//15

		head_of_long_hash_table[16].head_hash_table="TAA";head_of_long_hash_table[16].next=NULL;//16
		head_of_long_hash_table[17].head_hash_table="TAT";head_of_long_hash_table[17].next=NULL;//17
		head_of_long_hash_table[18].head_hash_table="TAC";head_of_long_hash_table[18].next=NULL;//18
		head_of_long_hash_table[19].head_hash_table="TAG";head_of_long_hash_table[19].next=NULL;//19
		head_of_long_hash_table[20].head_hash_table="TTA";head_of_long_hash_table[20].next=NULL;//20
		head_of_long_hash_table[21].head_hash_table="TTT";head_of_long_hash_table[21].next=NULL;//21
		head_of_long_hash_table[22].head_hash_table="TTC";head_of_long_hash_table[22].next=NULL;//22
		head_of_long_hash_table[23].head_hash_table="TTG";head_of_long_hash_table[23].next=NULL;//23
		head_of_long_hash_table[24].head_hash_table="TCA";head_of_long_hash_table[24].next=NULL;//24
		head_of_long_hash_table[25].head_hash_table="TCT";head_of_long_hash_table[25].next=NULL;//25
		head_of_long_hash_table[26].head_hash_table="TCC";head_of_long_hash_table[26].next=NULL;//26
		head_of_long_hash_table[27].head_hash_table="TCG";head_of_long_hash_table[27].next=NULL;//27
		head_of_long_hash_table[28].head_hash_table="TGA";head_of_long_hash_table[28].next=NULL;//28
		head_of_long_hash_table[29].head_hash_table="TGT";head_of_long_hash_table[29].next=NULL;//29
		head_of_long_hash_table[30].head_hash_table="TGC";head_of_long_hash_table[30].next=NULL;//30
		head_of_long_hash_table[31].head_hash_table="TGG";head_of_long_hash_table[31].next=NULL;//31

		head_of_long_hash_table[32].head_hash_table="CAA";head_of_long_hash_table[32].next=NULL;//32
		head_of_long_hash_table[33].head_hash_table="CAT";head_of_long_hash_table[33].next=NULL;//33
		head_of_long_hash_table[34].head_hash_table="CAC";head_of_long_hash_table[34].next=NULL;//34
		head_of_long_hash_table[35].head_hash_table="CAG";head_of_long_hash_table[35].next=NULL;//35
		head_of_long_hash_table[36].head_hash_table="CTA";head_of_long_hash_table[36].next=NULL;//36
		head_of_long_hash_table[37].head_hash_table="CTT";head_of_long_hash_table[37].next=NULL;//37
		head_of_long_hash_table[38].head_hash_table="CTC";head_of_long_hash_table[38].next=NULL;//38
		head_of_long_hash_table[39].head_hash_table="CTG";head_of_long_hash_table[39].next=NULL;//39
		head_of_long_hash_table[40].head_hash_table="CCA";head_of_long_hash_table[40].next=NULL;//40
		head_of_long_hash_table[41].head_hash_table="CCT";head_of_long_hash_table[41].next=NULL;//41
		head_of_long_hash_table[42].head_hash_table="CCC";head_of_long_hash_table[42].next=NULL;//42
		head_of_long_hash_table[43].head_hash_table="CCG";head_of_long_hash_table[43].next=NULL;//43
		head_of_long_hash_table[44].head_hash_table="CGA";head_of_long_hash_table[44].next=NULL;//44
		head_of_long_hash_table[45].head_hash_table="CGT";head_of_long_hash_table[45].next=NULL;//45
		head_of_long_hash_table[46].head_hash_table="CGC";head_of_long_hash_table[46].next=NULL;//46
		head_of_long_hash_table[47].head_hash_table="CGG";head_of_long_hash_table[47].next=NULL;//47

		head_of_long_hash_table[48].head_hash_table="GAA";head_of_long_hash_table[48].next=NULL;//48
		head_of_long_hash_table[49].head_hash_table="GAT";head_of_long_hash_table[49].next=NULL;//49
		head_of_long_hash_table[50].head_hash_table="GAC";head_of_long_hash_table[50].next=NULL;//50
		head_of_long_hash_table[51].head_hash_table="GAG";head_of_long_hash_table[51].next=NULL;//51
		head_of_long_hash_table[52].head_hash_table="GTA";head_of_long_hash_table[52].next=NULL;//52
		head_of_long_hash_table[53].head_hash_table="GTT";head_of_long_hash_table[53].next=NULL;//53
		head_of_long_hash_table[54].head_hash_table="GTC";head_of_long_hash_table[54].next=NULL;//54
		head_of_long_hash_table[55].head_hash_table="GTG";head_of_long_hash_table[55].next=NULL;//55
		head_of_long_hash_table[56].head_hash_table="GCA";head_of_long_hash_table[56].next=NULL;//56
		head_of_long_hash_table[57].head_hash_table="GCT";head_of_long_hash_table[57].next=NULL;//57
		head_of_long_hash_table[58].head_hash_table="GCC";head_of_long_hash_table[58].next=NULL;//58
		head_of_long_hash_table[59].head_hash_table="GCG";head_of_long_hash_table[59].next=NULL;//59
		head_of_long_hash_table[60].head_hash_table="GGA";head_of_long_hash_table[60].next=NULL;//60
		head_of_long_hash_table[61].head_hash_table="GGT";head_of_long_hash_table[61].next=NULL;//61
		head_of_long_hash_table[62].head_hash_table="GGC";head_of_long_hash_table[62].next=NULL;//62
		head_of_long_hash_table[63].head_hash_table="GGG";head_of_long_hash_table[63].next=NULL;//63
	}
	if(len_index==4)
	{


	}

}

void init_medium_hash_table()
{
	int real_len_head_hash;
	real_len_head_hash=int(pow(4,len_index));
	if(len_index==1)
	{
		head_of_medium_hash_table[0].head_hash_table="A";head_of_medium_hash_table[0].next=NULL;
		head_of_medium_hash_table[1].head_hash_table="T";head_of_medium_hash_table[1].next=NULL;
		head_of_medium_hash_table[2].head_hash_table="C";head_of_medium_hash_table[2].next=NULL;
		head_of_medium_hash_table[3].head_hash_table="G";head_of_medium_hash_table[3].next=NULL;
	}
	if(len_index==2)
	{
		head_of_medium_hash_table[0].head_hash_table="AA";head_of_medium_hash_table[0].next=NULL;//0
		head_of_medium_hash_table[1].head_hash_table="AT";head_of_medium_hash_table[1].next=NULL;//1
		head_of_medium_hash_table[2].head_hash_table="AC";head_of_medium_hash_table[2].next=NULL;//2
		head_of_medium_hash_table[3].head_hash_table="AG";head_of_medium_hash_table[3].next=NULL;//3
		head_of_medium_hash_table[4].head_hash_table="TA";head_of_medium_hash_table[4].next=NULL;//4
		head_of_medium_hash_table[5].head_hash_table="TT";head_of_medium_hash_table[5].next=NULL;//5
		head_of_medium_hash_table[6].head_hash_table="TC";head_of_medium_hash_table[6].next=NULL;//6
		head_of_medium_hash_table[7].head_hash_table="TG";head_of_medium_hash_table[7].next=NULL;//7
		head_of_medium_hash_table[8].head_hash_table="CA";head_of_medium_hash_table[8].next=NULL;//8
		head_of_medium_hash_table[9].head_hash_table="CT";head_of_medium_hash_table[9].next=NULL;//9
		head_of_medium_hash_table[10].head_hash_table="CC";head_of_medium_hash_table[10].next=NULL;//10
		head_of_medium_hash_table[11].head_hash_table="CG";head_of_medium_hash_table[11].next=NULL;//11
		head_of_medium_hash_table[12].head_hash_table="GA";head_of_medium_hash_table[12].next=NULL;//12
		head_of_medium_hash_table[13].head_hash_table="GT";head_of_medium_hash_table[13].next=NULL;//13
		head_of_medium_hash_table[14].head_hash_table="GC";head_of_medium_hash_table[14].next=NULL;//14
		head_of_medium_hash_table[15].head_hash_table="GG";head_of_medium_hash_table[15].next=NULL;//15
	}
	if(len_index==3)
	{
		head_of_medium_hash_table[0].head_hash_table="AAA";head_of_medium_hash_table[0].next=NULL;//0
		head_of_medium_hash_table[1].head_hash_table="AAT";head_of_medium_hash_table[1].next=NULL;//1
		head_of_medium_hash_table[2].head_hash_table="AAC";head_of_medium_hash_table[2].next=NULL;//2
		head_of_medium_hash_table[3].head_hash_table="AAG";head_of_medium_hash_table[3].next=NULL;//3
		head_of_medium_hash_table[4].head_hash_table="ATA";head_of_medium_hash_table[4].next=NULL;//4
		head_of_medium_hash_table[5].head_hash_table="ATT";head_of_medium_hash_table[5].next=NULL;//5
		head_of_medium_hash_table[6].head_hash_table="ATC";head_of_medium_hash_table[6].next=NULL;//6
		head_of_medium_hash_table[7].head_hash_table="ATG";head_of_medium_hash_table[7].next=NULL;//7
		head_of_medium_hash_table[8].head_hash_table="ACA";head_of_medium_hash_table[8].next=NULL;//8
		head_of_medium_hash_table[9].head_hash_table="ACT";head_of_medium_hash_table[9].next=NULL;//9
		head_of_medium_hash_table[10].head_hash_table="ACC";head_of_medium_hash_table[10].next=NULL;//10
		head_of_medium_hash_table[11].head_hash_table="ACG";head_of_medium_hash_table[11].next=NULL;//11
		head_of_medium_hash_table[12].head_hash_table="AGA";head_of_medium_hash_table[12].next=NULL;//12
		head_of_medium_hash_table[13].head_hash_table="AGT";head_of_medium_hash_table[13].next=NULL;//13
		head_of_medium_hash_table[14].head_hash_table="AGC";head_of_medium_hash_table[14].next=NULL;//14
		head_of_medium_hash_table[15].head_hash_table="AGG";head_of_medium_hash_table[15].next=NULL;//15

		head_of_medium_hash_table[16].head_hash_table="TAA";head_of_medium_hash_table[16].next=NULL;//16
		head_of_medium_hash_table[17].head_hash_table="TAT";head_of_medium_hash_table[17].next=NULL;//17
		head_of_medium_hash_table[18].head_hash_table="TAC";head_of_medium_hash_table[18].next=NULL;//18
		head_of_medium_hash_table[19].head_hash_table="TAG";head_of_medium_hash_table[19].next=NULL;//19
		head_of_medium_hash_table[20].head_hash_table="TTA";head_of_medium_hash_table[20].next=NULL;//20
		head_of_medium_hash_table[21].head_hash_table="TTT";head_of_medium_hash_table[21].next=NULL;//21
		head_of_medium_hash_table[22].head_hash_table="TTC";head_of_medium_hash_table[22].next=NULL;//22
		head_of_medium_hash_table[23].head_hash_table="TTG";head_of_medium_hash_table[23].next=NULL;//23
		head_of_medium_hash_table[24].head_hash_table="TCA";head_of_medium_hash_table[24].next=NULL;//24
		head_of_medium_hash_table[25].head_hash_table="TCT";head_of_medium_hash_table[25].next=NULL;//25
		head_of_medium_hash_table[26].head_hash_table="TCC";head_of_medium_hash_table[26].next=NULL;//26
		head_of_medium_hash_table[27].head_hash_table="TCG";head_of_medium_hash_table[27].next=NULL;//27
		head_of_medium_hash_table[28].head_hash_table="TGA";head_of_medium_hash_table[28].next=NULL;//28
		head_of_medium_hash_table[29].head_hash_table="TGT";head_of_medium_hash_table[29].next=NULL;//29
		head_of_medium_hash_table[30].head_hash_table="TGC";head_of_medium_hash_table[30].next=NULL;//30
		head_of_medium_hash_table[31].head_hash_table="TGG";head_of_medium_hash_table[31].next=NULL;//31

		head_of_medium_hash_table[32].head_hash_table="CAA";head_of_medium_hash_table[32].next=NULL;//32
		head_of_medium_hash_table[33].head_hash_table="CAT";head_of_medium_hash_table[33].next=NULL;//33
		head_of_medium_hash_table[34].head_hash_table="CAC";head_of_medium_hash_table[34].next=NULL;//34
		head_of_medium_hash_table[35].head_hash_table="CAG";head_of_medium_hash_table[35].next=NULL;//35
		head_of_medium_hash_table[36].head_hash_table="CTA";head_of_medium_hash_table[36].next=NULL;//36
		head_of_medium_hash_table[37].head_hash_table="CTT";head_of_medium_hash_table[37].next=NULL;//37
		head_of_medium_hash_table[38].head_hash_table="CTC";head_of_medium_hash_table[38].next=NULL;//38
		head_of_medium_hash_table[39].head_hash_table="CTG";head_of_medium_hash_table[39].next=NULL;//39
		head_of_medium_hash_table[40].head_hash_table="CCA";head_of_medium_hash_table[40].next=NULL;//40
		head_of_medium_hash_table[41].head_hash_table="CCT";head_of_medium_hash_table[41].next=NULL;//41
		head_of_medium_hash_table[42].head_hash_table="CCC";head_of_medium_hash_table[42].next=NULL;//42
		head_of_medium_hash_table[43].head_hash_table="CCG";head_of_medium_hash_table[43].next=NULL;//43
		head_of_medium_hash_table[44].head_hash_table="CGA";head_of_medium_hash_table[44].next=NULL;//44
		head_of_medium_hash_table[45].head_hash_table="CGT";head_of_medium_hash_table[45].next=NULL;//45
		head_of_medium_hash_table[46].head_hash_table="CGC";head_of_medium_hash_table[46].next=NULL;//46
		head_of_medium_hash_table[47].head_hash_table="CGG";head_of_medium_hash_table[47].next=NULL;//47

		head_of_medium_hash_table[48].head_hash_table="GAA";head_of_medium_hash_table[48].next=NULL;//48
		head_of_medium_hash_table[49].head_hash_table="GAT";head_of_medium_hash_table[49].next=NULL;//49
		head_of_medium_hash_table[50].head_hash_table="GAC";head_of_medium_hash_table[50].next=NULL;//50
		head_of_medium_hash_table[51].head_hash_table="GAG";head_of_medium_hash_table[51].next=NULL;//51
		head_of_medium_hash_table[52].head_hash_table="GTA";head_of_medium_hash_table[52].next=NULL;//52
		head_of_medium_hash_table[53].head_hash_table="GTT";head_of_medium_hash_table[53].next=NULL;//53
		head_of_medium_hash_table[54].head_hash_table="GTC";head_of_medium_hash_table[54].next=NULL;//54
		head_of_medium_hash_table[55].head_hash_table="GTG";head_of_medium_hash_table[55].next=NULL;//55
		head_of_medium_hash_table[56].head_hash_table="GCA";head_of_medium_hash_table[56].next=NULL;//56
		head_of_medium_hash_table[57].head_hash_table="GCT";head_of_medium_hash_table[57].next=NULL;//57
		head_of_medium_hash_table[58].head_hash_table="GCC";head_of_medium_hash_table[58].next=NULL;//58
		head_of_medium_hash_table[59].head_hash_table="GCG";head_of_medium_hash_table[59].next=NULL;//59
		head_of_medium_hash_table[60].head_hash_table="GGA";head_of_medium_hash_table[60].next=NULL;//60
		head_of_medium_hash_table[61].head_hash_table="GGT";head_of_medium_hash_table[61].next=NULL;//61
		head_of_medium_hash_table[62].head_hash_table="GGC";head_of_medium_hash_table[62].next=NULL;//62
		head_of_medium_hash_table[63].head_hash_table="GGG";head_of_medium_hash_table[63].next=NULL;//63
	}
	if(len_index==4)
	{


	}

}

int get_index_num(char *ch)                                             //check the first len_index bases of kmer len_index to get its index_num��AA---0,AT---1,AC---2,AG---3...��
{
	if(len_index==1)
	{
		if(ch[0]=='A') index_num=0;
		if(ch[0]=='T') index_num=1;
		if(ch[0]=='C') index_num=2;
		if(ch[0]=='G') index_num=3;
	}
	if(len_index==2)
	{
		if(ch[0]=='A')
		{
			if(ch[1]=='A') index_num=0;
			if(ch[1]=='T') index_num=1;
			if(ch[1]=='C') index_num=2;
			if(ch[1]=='G') index_num=3;
		}
		if(ch[0]=='T')
		{
			if(ch[1]=='A') index_num=4;
			if(ch[1]=='T') index_num=5;
			if(ch[1]=='C') index_num=6;
			if(ch[1]=='G') index_num=7;
		}
		if(ch[0]=='C')
		{
			if(ch[1]=='A') index_num=8;
			if(ch[1]=='T') index_num=9;
			if(ch[1]=='C') index_num=10;
			if(ch[1]=='G') index_num=11;
		}
		if(ch[0]=='G')
		{
			if(ch[1]=='A') index_num=12;
			if(ch[1]=='T') index_num=13;
			if(ch[1]=='C') index_num=14;
			if(ch[1]=='G') index_num=15;
		}
	}
	if(len_index==3)
	{
		if(ch[0]=='A')
		{
			if(ch[1]=='A')
			{
				if(ch[2]=='A') index_num=0;
				if(ch[2]=='T') index_num=1;
				if(ch[2]=='C') index_num=2;
				if(ch[2]=='G') index_num=3;
			}
			if(ch[1]=='T')
			{
				if(ch[2]=='A') index_num=4;
				if(ch[2]=='T') index_num=5;
				if(ch[2]=='C') index_num=6;
				if(ch[2]=='G') index_num=7;
			}
			if(ch[1]=='C')
			{
				if(ch[2]=='A') index_num=8;
				if(ch[2]=='T') index_num=9;
				if(ch[2]=='C') index_num=10;
				if(ch[2]=='G') index_num=11;
			}
			if(ch[1]=='G')
			{
				if(ch[2]=='A') index_num=12;
				if(ch[2]=='T') index_num=13;
				if(ch[2]=='C') index_num=14;
				if(ch[2]=='G') index_num=15;
			}
		}
		if(ch[0]=='T')
		{
			if(ch[1]=='A')
			{
				if(ch[2]=='A') index_num=16;
				if(ch[2]=='T') index_num=17;
				if(ch[2]=='C') index_num=18;
				if(ch[2]=='G') index_num=19;
			}
			if(ch[1]=='T')
			{
				if(ch[2]=='A') index_num=20;
				if(ch[2]=='T') index_num=21;
				if(ch[2]=='C') index_num=22;
				if(ch[2]=='G') index_num=23;
			}
			if(ch[1]=='C')
			{
				if(ch[2]=='A') index_num=24;
				if(ch[2]=='T') index_num=25;
				if(ch[2]=='C') index_num=26;
				if(ch[2]=='G') index_num=27;
			}
			if(ch[1]=='G')
			{
				if(ch[2]=='A') index_num=28;
				if(ch[2]=='T') index_num=29;
				if(ch[2]=='C') index_num=30;
				if(ch[2]=='G') index_num=31;
			}
		}
		if(ch[0]=='C')
		{
			if(ch[1]=='A')
			{
				if(ch[2]=='A') index_num=32;
				if(ch[2]=='T') index_num=33;
				if(ch[2]=='C') index_num=34;
				if(ch[2]=='G') index_num=35;
			}
			if(ch[1]=='T')
			{
				if(ch[2]=='A') index_num=36;
				if(ch[2]=='T') index_num=37;
				if(ch[2]=='C') index_num=38;
				if(ch[2]=='G') index_num=39;
			}
			if(ch[1]=='C')
			{
				if(ch[2]=='A') index_num=40;
				if(ch[2]=='T') index_num=41;
				if(ch[2]=='C') index_num=42;
				if(ch[2]=='G') index_num=43;
			}
			if(ch[1]=='G')
			{
				if(ch[2]=='A') index_num=44;
				if(ch[2]=='T') index_num=45;
				if(ch[2]=='C') index_num=46;
				if(ch[2]=='G') index_num=47;
			}
		}
		if(ch[0]=='G')
		{
			if(ch[1]=='A')
			{
				if(ch[2]=='A') index_num=48;
				if(ch[2]=='T') index_num=49;
				if(ch[2]=='C') index_num=50;
				if(ch[2]=='G') index_num=51;
			}
			if(ch[1]=='T')
			{
				if(ch[2]=='A') index_num=52;
				if(ch[2]=='T') index_num=53;
				if(ch[2]=='C') index_num=54;
				if(ch[2]=='G') index_num=55;
			}
			if(ch[1]=='C')
			{
				if(ch[2]=='A') index_num=56;
				if(ch[2]=='T') index_num=57;
				if(ch[2]=='C') index_num=58;
				if(ch[2]=='G') index_num=59;
			}
			if(ch[1]=='G')
			{
				if(ch[2]=='A') index_num=60;
				if(ch[2]=='T') index_num=61;
				if(ch[2]=='C') index_num=62;
				if(ch[2]=='G') index_num=63;
			}
		}
	}



	return index_num;
}

void insert_to_short_hash_table()
{
	//num_short_k_mer=0;
	int index_num;
	char *ch=(char *)malloc(short_k*sizeof(char));                                                   //temp str for temp k-mer
	
	struct k_mer_node *new_node,*old_node;                                //new_node stands for the nodes haven't been setted into hash table
	                                                                      //old_node stands for the nodes have been setted into hash table
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-short_k+1;j++)
		{
			int k=0;
			for(k=0;k<short_k;k++)
			{
				ch[k]=sample[i][j+k];
			}
			ch[k]='\0';
			index_num=get_index_num(ch);
			int flag=NEW_KMER_NODE;
			                                                               //check if the node have been in hash table
			if (head_of_short_hash_table[index_num].next!=NULL)
			{
				old_node=head_of_short_hash_table[index_num].next;               //old_node points to the first node headed by index_num
				while(old_node!=NULL)
				{
					int m;
					for(m=0;m<short_k;m++)
					{
						if(ch[m]==sample[old_node->i_read][old_node->j_position+m])					
							continue;
						else
							break;
					}
					if(m==short_k)
					{
						flag=OLD_KMER_NODE;
					}
					if(flag==OLD_KMER_NODE)
					{
						break;                                              //if the node is found in hash tabld , break
					}
					else                                                    //if the node is not found in hash tabld , continue
					{
						old_node=old_node->next;
					}
				}
			}
			if(flag==NEW_KMER_NODE)                                         //if the node is decide to be new node, add it into hash table (to be the first node pointed by index_num),and count_of_node++;
			{
				new_node=(struct k_mer_node*)malloc(LEN);
				new_node->i_read=i;
				new_node->j_position=j;
				new_node->next=head_of_short_hash_table[index_num].next;
				head_of_short_hash_table[index_num].next=new_node;
				count_of_short_node[new_node->i_read][new_node->j_position]++;
				len_short_hash_table[index_num]++;
				every_short_kmer_times_position[i][j].i_position=i;
				every_short_kmer_times_position[i][j].j_position=j;
			}
			if(flag==OLD_KMER_NODE)                                         //if the node is decide to be old node,  count_of_node++, add the old position to temp kmer
			{
				count_of_short_node[old_node->i_read][old_node->j_position]++;
				every_short_kmer_times_position[i][j].i_position=old_node->i_read;
				every_short_kmer_times_position[i][j].j_position=old_node->j_position;
			}
			//num_short_k_mer++;
		}
	}
	free(ch);
}

void insert_to_long_hash_table()
{
	//num_short_k_mer=0;
	int index_num;
	char *ch=(char *)malloc(long_k*sizeof(char));                                                   //temp str for temp k-mer
	
	struct k_mer_node *new_node,*old_node;                                //new_node stands for the nodes haven't been setted into hash table
	                                                                      //old_node stands for the nodes have been setted into hash table
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-long_k+1;j++)
		{
			int k=0;
			for(k=0;k<long_k;k++)
			{
				ch[k]=sample[i][j+k];
			}
			ch[k]='\0';
			index_num=get_index_num(ch);
			int flag=NEW_KMER_NODE;
			                                                               //check if the node have been in hash table
			if (head_of_long_hash_table[index_num].next!=NULL)
			{
				old_node=head_of_long_hash_table[index_num].next;               //old_node points to the first node headed by index_num
				while(old_node!=NULL)
				{
					int m;
					for(m=0;m<long_k;m++)
					{
						if(ch[m]==sample[old_node->i_read][old_node->j_position+m])					
							continue;
						else
							break;
					}
					if(m==long_k)
					{
						flag=OLD_KMER_NODE;
					}
					if(flag==OLD_KMER_NODE)
					{
						break;                                              //if the node is found in hash tabld , break
					}
					else                                                    //if the node is not found in hash tabld , continue
					{
						old_node=old_node->next;
					}
				}
			}
			if(flag==NEW_KMER_NODE)                                         //if the node is decide to be new node, add it into hash table (to be the first node pointed by index_num),and count_of_node++;
			{
				new_node=(struct k_mer_node*)malloc(LEN);
				new_node->i_read=i;
				new_node->j_position=j;
				new_node->next=head_of_long_hash_table[index_num].next;
				head_of_long_hash_table[index_num].next=new_node;
				count_of_long_node[new_node->i_read][new_node->j_position]++;
				len_long_hash_table[index_num]++;
				every_long_kmer_times_position[i][j].i_position=i;
				every_long_kmer_times_position[i][j].j_position=j;
			}
			if(flag==OLD_KMER_NODE)                                         //if the node is decide to be old node,  count_of_node++, add the old position to temp kmer
			{
				count_of_long_node[old_node->i_read][old_node->j_position]++;
				every_long_kmer_times_position[i][j].i_position=old_node->i_read;
				every_long_kmer_times_position[i][j].j_position=old_node->j_position;
			}
			//num_short_k_mer++;
		}
	}
	free(ch);
}

void insert_to_medium_hash_table()
{
	//num_short_k_mer=0;
	int index_num;
	char *ch=(char *)malloc(medium_k*sizeof(char));                                                   //temp str for temp k-mer
	
	struct k_mer_node *new_node,*old_node;                                //new_node stands for the nodes haven't been setted into hash table
	                                                                      //old_node stands for the nodes have been setted into hash table
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-medium_k+1;j++)
		{
			int k=0;
			for(k=0;k<medium_k;k++)
			{
				ch[k]=sample[i][j+k];
			}
			ch[k]='\0';
			index_num=get_index_num(ch);
			int flag=NEW_KMER_NODE;
			                                                               //check if the node have been in hash table
			if (head_of_medium_hash_table[index_num].next!=NULL)
			{
				old_node=head_of_medium_hash_table[index_num].next;               //old_node points to the first node headed by index_num
				while(old_node!=NULL)
				{
					int m;
					for(m=0;m<long_k;m++)
					{
						if(ch[m]==sample[old_node->i_read][old_node->j_position+m])					
							continue;
						else
							break;
					}
					if(m==medium_k)
					{
						flag=OLD_KMER_NODE;
					}
					if(flag==OLD_KMER_NODE)
					{
						break;                                              //if the node is found in hash tabld , break
					}
					else                                                    //if the node is not found in hash tabld , continue
					{
						old_node=old_node->next;
					}
				}
			}
			if(flag==NEW_KMER_NODE)                                         //if the node is decide to be new node, add it into hash table (to be the first node pointed by index_num),and count_of_node++;
			{
				new_node=(struct k_mer_node*)malloc(LEN);
				new_node->i_read=i;
				new_node->j_position=j;
				new_node->next=head_of_medium_hash_table[index_num].next;
				head_of_medium_hash_table[index_num].next=new_node;
				count_of_medium_node[new_node->i_read][new_node->j_position]++;
				len_medium_hash_table[index_num]++;
				every_medium_kmer_times_position[i][j].i_position=i;
				every_medium_kmer_times_position[i][j].j_position=j;
			}
			if(flag==OLD_KMER_NODE)                                         //if the node is decide to be old node,  count_of_node++, add the old position to temp kmer
			{
				count_of_medium_node[old_node->i_read][old_node->j_position]++;
				every_medium_kmer_times_position[i][j].i_position=old_node->i_read;
				every_medium_kmer_times_position[i][j].j_position=old_node->j_position;
			}
			//num_short_k_mer++;
		}
	}
	free(ch);
}

void get_every_short_kmer_times()
{
	for(int i =0;i<num_read;i++)
	{
		for(int j=0 ;j<length[i]-short_k+1;j++)
		{
			count_of_short_node[i][j]=count_of_short_node[every_short_kmer_times_position[i][j].i_position][every_short_kmer_times_position[i][j].j_position];
		}
	}
}

void get_every_long_kmer_times()
{
	for(int i =0;i<num_read;i++)
	{
		for(int j=0 ;j<length[i]-long_k+1;j++)
		{
			count_of_long_node[i][j]=count_of_long_node[every_long_kmer_times_position[i][j].i_position][every_long_kmer_times_position[i][j].j_position];
		}
	}
}

void get_every_medium_kmer_times()
{
	for(int i =0;i<num_read;i++)
	{
		for(int j=0 ;j<length[i]-medium_k+1;j++)
		{
			count_of_medium_node[i][j]=count_of_medium_node[every_medium_kmer_times_position[i][j].i_position][every_medium_kmer_times_position[i][j].j_position];
		}
	}
}

void count_short_k_mer_times()
{
	int index_num;
	char *ch=(char *)malloc(short_k*sizeof(char));                                                     //temp str for temp k-mer
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-short_k+1;j++)
		{
			int k=0;
			for(k=0;k<short_k;k++)
			{
				ch[k]=sample[i][j+k];
			}
			ch[k]='\0';
			index_num=get_index_num(ch);
			struct k_mer_node *node;
			if (head_of_short_hash_table[index_num].next!=NULL)
			{
				node=head_of_short_hash_table[index_num].next;
				//int flag;
				
				while (node!=NULL)
				{
					int hamming=0;
					int m;
					for(m=len_index;m<short_k;m++)                         //get hamming distance between kmer in hash table and kmer in read
					{
						if(ch[m]!=sample[node->i_read][node->j_position+m])
						{
							hamming++;
						}
					}
					if(hamming<=len_k_mer_error)                           //the  kmer in hash table can match to  kmer in read
					{
						if(count_of_short_node[node->i_read][node->j_position]>=count_value)//if the count of node is a valid kmer, add it to times 
						{
							
							int n=0;
							for(n=0;n<short_k;n++)
							{
								short_k_mer_times[i][j+n].times++;
								short_k_mer_times[i][j+n].flag=1;
							}
								short_k_mer_times[i][j+n-1].flag=0;
							//k_mer_times[i][j+short_k-1].flag=1;//kmerβ����־
						}
						times[i][j]+=count_of_short_node[node->i_read][node->j_position];
					}
					node=node->next;
				}
			}
		}
	}
	free(ch);
}

void count_long_k_mer_times()
{
	int index_num;
	char *ch=(char *)malloc(long_k*sizeof(char));                                                     //temp str for temp k-mer
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-long_k+1;j++)
		{
			int k=0;
			for(k=0;k<long_k;k++)
			{
				ch[k]=sample[i][j+k];
			}
			ch[k]='\0';
			index_num=get_index_num(ch);
			struct k_mer_node *node;
			if (head_of_long_hash_table[index_num].next!=NULL)
			{
				node=head_of_long_hash_table[index_num].next;
				//int flag;
				
				while (node!=NULL)
				{
					int hamming=0;
					int m;
					for(m=len_index;m<long_k;m++)                         //get hamming distance between kmer in hash table and kmer in read
					{
						if(ch[m]!=sample[node->i_read][node->j_position+m])
						{
							hamming++;
						}
					}
					if(hamming<=len_k_mer_error)                           //the  kmer in hash table can match to  kmer in read
					{
						if(count_of_long_node[node->i_read][node->j_position]>=count_value)//if the count of node is a valid kmer, add it to times 
						{
							
							int n=0;
							for(n=0;n<long_k;n++)
							{
								long_k_mer_times[i][j+n].times++;
								long_k_mer_times[i][j+n].flag=1;
							}
								long_k_mer_times[i][j+n-1].flag=0;
							//k_mer_times[i][j+short_k-1].flag=1;//kmerβ����־
						}
						times[i][j]+=count_of_long_node[node->i_read][node->j_position];
					}
					node=node->next;
				}
			}
		}
	}
	free(ch);
}

void count_medium_k_mer_times()
{
	int index_num;
	char *ch=(char *)malloc(medium_k*sizeof(char));                                                     //temp str for temp k-mer
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-medium_k+1;j++)
		{
			int k=0;
			for(k=0;k<medium_k;k++)
			{
				ch[k]=sample[i][j+k];
			}
			ch[k]='\0';
			index_num=get_index_num(ch);
			struct k_mer_node *node;
			if (head_of_medium_hash_table[index_num].next!=NULL)
			{
				node=head_of_medium_hash_table[index_num].next;
				//int flag;
				
				while (node!=NULL)
				{
					int hamming=0;
					int m;
					for(m=len_index;m<medium_k;m++)                         //get hamming distance between kmer in hash table and kmer in read
					{
						if(ch[m]!=sample[node->i_read][node->j_position+m])
						{
							hamming++;
						}
					}
					if(hamming<=len_k_mer_error)                           //the  kmer in hash table can match to  kmer in read
					{
						if(count_of_medium_node[node->i_read][node->j_position]>=count_value)//if the count of node is a valid kmer, add it to times 
						{
							
							int n=0;
							for(n=0;n<medium_k;n++)
							{
								medium_k_mer_times[i][j+n].times++;
								medium_k_mer_times[i][j+n].flag=1;
							}
								medium_k_mer_times[i][j+n-1].flag=0;
							//k_mer_times[i][j+short_k-1].flag=1;//kmerβ����־
						}
						times[i][j]+=count_of_medium_node[node->i_read][node->j_position];
					}
					node=node->next;
				}
			}
		}
	}
	free(ch);
}

void get_short_base_score()
{
	for(int i = 0 ; i< num_read;i++)
	{
		for(int j =0 ;j<length[i];j++)
		{
			for(int k=0;k<short_k;k++)
			{
				short_base_score[i][j+k].times+=count_of_short_node[i][j];
			}
		}
	}
}

void get_long_base_score()
{
	for(int i = 0 ; i< num_read;i++)
	{
		for(int j =0 ;j<length[i];j++)
		{
			for(int k=0;k<long_k;k++)
			{
				long_base_score[i][j+k].times+=count_of_long_node[i][j];
			}
		}
	}
}

void get_medium_base_score()
{
	for(int i = 0 ; i< num_read;i++)
	{
		for(int j =0 ;j<length[i];j++)
		{
			for(int k=0;k<medium_k;k++)
			{
				medium_base_score[i][j+k].times+=count_of_medium_node[i][j];
			}
		}
	}
}

void save_as_short_base_score()
{
	ofstream fout1;
	fout1.open("E:\\short_base_score.txt",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout1<<"NO."<<i<<"\t";
		for(int j=0;j<length[i];j++)
		{
			fout1<<short_base_score[i][j].times<<"\t";
		}
		fout1<<endl;
	}

}

void save_as_long_base_score()
{
	ofstream fout1;
	fout1.open("E:\\long_base_score.txt",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout1<<"NO."<<i<<"\t";
		for(int j=0;j<length[i];j++)
		{
			fout1<<long_base_score[i][j].times<<"\t";
		}
		fout1<<endl;
	}

}

void save_as_medium_base_score()
{
	ofstream fout1;
	fout1.open("E:\\medium_base_score.txt",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout1<<"NO."<<i<<"\t";
		for(int j=0;j<length[i];j++)
		{
			fout1<<medium_base_score[i][j].times<<"\t";
		}
		fout1<<endl;
	}

}

void save_as_count_of_short_nodes()
{
	ofstream fout3;
	fout3.open("E:\\count_of_short_node.txt",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout3<<"NO."<<i<<"\t";
		for(int j=0;j<length[i];j++)
		{
			fout3<<"["<<j<<"]"<<count_of_short_node[i][j]<<"\t";
		}
		fout3<<endl;
	}

}

void save_as_count_of_long_nodes()
{
	ofstream fout3;
	fout3.open("E:\\count_of_long_node.txt",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout3<<"NO."<<i<<"\t";
		for(int j=0;j<length[i];j++)
		{
			fout3<<"["<<j<<"]"<<count_of_long_node[i][j]<<"\t";
		}
		fout3<<endl;
	}

}

void save_as_count_of_medium_nodes()
{
	ofstream fout3;
	fout3.open("E:\\count_of_medium_node.txt",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout3<<"NO."<<i<<"\t";
		for(int j=0;j<length[i];j++)
		{
			fout3<<"["<<j<<"]"<<count_of_medium_node[i][j]<<"\t";
		}
		fout3<<endl;
	}

}

void save_as_short_kmer_times()
{
	ofstream fout3;
	fout3.open("E:\\short_kmer_times.txt",ios::trunc);

	fout3<<">target"<<" "<<endl;
	for(int i=0;i<num_read;i++)
	{
		fout3<<"NO."<<i<<"   ";
		for(int j=0;j<length[i];j++)
		{
			fout3<<short_k_mer_times[i][j].times<<" ";
		}
		fout3<<endl;
		fout3<<endl;
	}
}

void save_as_long_kmer_times()
{
	ofstream fout3;
	fout3.open("E:\\long_kmer_times.txt",ios::trunc);

	fout3<<">target"<<" "<<endl;
	for(int i=0;i<num_read;i++)
	{
		fout3<<"NO."<<i<<"   ";
		for(int j=0;j<length[i];j++)
		{
			fout3<<long_k_mer_times[i][j].times<<" ";
		}
		fout3<<endl;
		fout3<<endl;
	}
}

void save_as_medium_kmer_times()
{
	ofstream fout3;
	fout3.open("E:\\medium_kmer_times.txt",ios::trunc);

	fout3<<">target"<<" "<<endl;
	for(int i=0;i<num_read;i++)
	{
		fout3<<"NO."<<i<<"   ";
		for(int j=0;j<length[i];j++)
		{
			fout3<<medium_k_mer_times[i][j].times<<" ";
		}
		fout3<<endl;
		fout3<<endl;
	}
}

extern void hash_table()
{
	init_short_hash_table();
	init_long_hash_table();
	init_medium_hash_table();

	init_k_mer_times();
	init_every_kmer_times_position();
	//short
	insert_to_short_hash_table();
	
	get_every_short_kmer_times();
	
	count_short_k_mer_times();
	
	get_short_base_score();
	
	save_as_short_base_score();
	
	save_as_short_kmer_times();
	
	save_as_count_of_short_nodes();
	
	//long
	insert_to_long_hash_table();
	get_every_long_kmer_times();
	count_long_k_mer_times();
	get_long_base_score();
	save_as_long_base_score();
	save_as_long_kmer_times();
	save_as_count_of_long_nodes();
	//medium
	insert_to_medium_hash_table();
	
	get_every_medium_kmer_times();
	
	count_medium_k_mer_times();
	
	get_medium_base_score();
	
	save_as_medium_base_score();
	
	save_as_medium_kmer_times();
	
	save_as_count_of_medium_nodes();
	
}