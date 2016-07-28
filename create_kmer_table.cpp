#include"main.h"
#include"create_kmer_table.h"

extern int read_count ;
extern int len_read ;         //read length 1000-10000

int index_num;
extern int *length;
extern char **sample;

extern node_position **every_short_kmer_times_position;          
extern node_position **every_medium_kmer_times_position;       

extern head_node_index head_of_short_hash_table;                  //hash table's head node
extern head_node_index head_of_medium_hash_table;                  //hash table's head node

/*
extern int len_short_hash_table[max_len_head_hash];                                  //number of hash table's head node
extern int len_long_hash_table[max_len_head_hash];                                  //number of hash table's head node
extern int len_medium_hash_table[max_len_head_hash];                                  //number of hash table's head node
*/
extern int **count_of_short_node;                           //array of counting every kmer (in read) times
extern int **count_of_medium_node;                           //array of counting every kmer (in read) times

int get_index_num(char *ch)                                             //check the first len_index bases of kmer len_index to get its index_num£¨AA---0,AT---1,AC---2,AG---3...£©
{
	int i;
	//0
	if(ch[0]=='A')
	{
		i=0;
	}
	else if(ch[0]=='T')
	{
		i=256;
	}
	else if(ch[0]=='C')
	{
		i=512;
	}
	else if(ch[0]=='G')
	{
		i=768;
	}
	//1
	if(ch[1]=='A')
	{
		i+=0;
	}
	else if(ch[1]=='T')
	{
		i+=64;
	}
	else if(ch[1]=='C')
	{
		i+=128;
	}
	else if(ch[1]=='G')
	{
		i+=192;
	}
	//2
	if(ch[2]=='A')
	{
		i+=0;
	}
	else if(ch[2]=='T')
	{
		i+=16;
	}
	else if(ch[2]=='C')
	{
		i+=32;
	}
	else if(ch[2]=='G')
	{
		i+=48;
	}
	//3
	if(ch[3]=='A')
	{
		i+=0;
	}
	else if(ch[3]=='T')
	{
		i+=4;
	}
	else if(ch[3]=='C')
	{
		i+=8;
	}
	else if(ch[3]=='G')
	{
		i+=12;
	}
	//4
	if(ch[4]=='A')
	{
		i+=0;
	}
	else if(ch[4]=='T')
	{
		i+=1;
	}
	else if(ch[4]=='C')
	{
		i+=2;
	}
	else if(ch[4]=='G')
	{
		i+=3;
	}
	return i;
}

void insert_to_short_hash_table()
{
	int index_num;
	char *ch=(char *)malloc((short_k+1)*sizeof(char));                                                   //temp str for temp k-mer
	struct k_mer_node *new_node,*old_node;                       
	//new_node stands for the nodes haven't been setted into hash table; old_node stands for the nodes have been setted into hash table
	for(int i=0;i<read_count;i++)
	{
		cout<<"searching No. "<<i<<" read..."<<endl;
		if(length[i]>short_k)
		{
			for(int j=0;j<length[i]-short_k+1;j++)
			{
				int k=0;
				for(k=0;k<len_index;k++)
				{
					ch[k]=sample[i][j+k];
				}
				//ch[k]='\0';
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
					new_node=(k_mer_node*)malloc(LEN);
					new_node->i_read=i;
					new_node->j_position=j;
					new_node->next=head_of_short_hash_table[index_num].next;
					head_of_short_hash_table[index_num].next=new_node;
					count_of_short_node[new_node->i_read][new_node->j_position]++;
					//len_short_hash_table[index_num]++;
					every_short_kmer_times_position[i][j].i_position=i;
					every_short_kmer_times_position[i][j].j_position=j;
				}
				if(flag==OLD_KMER_NODE)                                         //if the node is decide to be old node,  count_of_node++, add the old position to temp kmer
				{
					count_of_short_node[old_node->i_read][old_node->j_position]++;
					every_short_kmer_times_position[i][j].i_position=old_node->i_read;
					every_short_kmer_times_position[i][j].j_position=old_node->j_position;
				}
			}
		}
	}
	free(ch);
}

void insert_to_medium_hash_table()
{
	int index_num;
	char *ch=(char *)malloc((medium_k+1)*sizeof(char));                                                   //temp str for temp k-mer
	struct k_mer_node *new_node,*old_node;                                //new_node stands for the nodes haven't been setted into hash table
	                                                                      //old_node stands for the nodes have been setted into hash table
	for(int i=0;i<read_count;i++)
	{
				cout<<"searching No. "<<i<<" read..."<<endl;

		if(length[i]>medium_k)
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
					//len_medium_hash_table[index_num]++;
					every_medium_kmer_times_position[i][j].i_position=i;
					every_medium_kmer_times_position[i][j].j_position=j;
				}
				if(flag==OLD_KMER_NODE)                                         //if the node is decide to be old node,  count_of_node++, add the old position to temp kmer
				{
					count_of_medium_node[old_node->i_read][old_node->j_position]++;
					every_medium_kmer_times_position[i][j].i_position=old_node->i_read;
					every_medium_kmer_times_position[i][j].j_position=old_node->j_position;
				}
			}
		}
	}
	free(ch);
}

void get_every_short_kmer_times()
{
	for(int i =0;i<read_count;i++)
	{
		for(int j=0 ;j<length[i]-short_k+1;j++)
		{
			count_of_short_node[i][j]=count_of_short_node[every_short_kmer_times_position[i][j].i_position][every_short_kmer_times_position[i][j].j_position];
		}
	}
}

void get_every_medium_kmer_times()
{
	for(int i =0;i<read_count;i++)
	{
		for(int j=0 ;j<length[i]-medium_k+1;j++)
		{
			count_of_medium_node[i][j]=count_of_medium_node[every_medium_kmer_times_position[i][j].i_position][every_medium_kmer_times_position[i][j].j_position];
		}
	}
}

void save_as_count_of_short_nodes()
{
	ofstream fout3;
	fout3.open("./count_of_short_node.txt",ios::trunc);
	for(int i=0;i<read_count;i++)
	{
		fout3<<"NO."<<i<<"\t";
		for(int j=0;j<length[i];j++)
		{
			//fout3<<"["<<j<<"]"<<count_of_short_node[i][j]<<"\t";
			fout3<<count_of_short_node[i][j]<<" ";
		}
		fout3<<endl;
	}

}

void save_as_count_of_medium_nodes()
{
	ofstream fout3;
	fout3.open("./count_of_medium_node.txt",ios::trunc);
	for(int i=0;i<read_count;i++)
	{
		fout3<<"NO."<<i<<"\t";
		for(int j=0;j<length[i];j++)
		{
			//fout3<<"["<<j<<"]"<<count_of_medium_node[i][j]<<"\t";
			fout3<<count_of_medium_node[i][j]<<" ";

		}
		fout3<<endl;
	}
}

extern void hash_table()
{
	//short
	cout<<"short_hash_table"<<endl;
	insert_to_short_hash_table();
	get_every_short_kmer_times();
	save_as_count_of_short_nodes();
	
	//medium
	cout<<"medium_hash_table"<<endl;
	insert_to_medium_hash_table();
	get_every_medium_kmer_times();
	save_as_count_of_medium_nodes();
}