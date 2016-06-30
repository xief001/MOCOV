#include"dataset.h"
#include"hash_table.h"
#include"correct_errors.h"
extern struct times_and_flag;

extern times_and_flag short_k_mer_times[num_read][int(len_read*2)];
extern times_and_flag long_k_mer_times[num_read][int(len_read*2)];
extern times_and_flag medium_k_mer_times[num_read][int(len_read*2)];

extern int count_of_short_node[num_read][int(len_read*2)];
extern int count_of_long_node[num_read][int(len_read*2)];
extern int count_of_medium_node[num_read][int(len_read*2)];


extern int get_index_num(char *ch);
extern int num_of_insertion=0;
rawbase rawread[num_read][int(len_read*2)];

char corrected_read[num_read][int(len_read*2)]={0};
char deleted_read[num_read*100][int(len_read*2)]={0};
char deleted_read1[num_read*100][int(len_read*2)]={0};
int deleted_read_num=0;
int deleted_read_num1=0;
int deleted_ori[num_read*100]={0};
int deleted_len[num_read*100]={0};
int deleted_len1[num_read*100]={0};
int corrected_len[num_read]={0};

int insertion[num_read*100][int(len_read*2)]={0};
int deletion[num_read*100][int(len_read*2)]={0};
//com complex[num_read*100][int(len_read*2)]={0};

extern head_node_index head_of_short_hash_table[max_len_head_hash];
extern head_node_index head_of_long_hash_table[max_len_head_hash];
extern head_node_index head_of_medium_hash_table[max_len_head_hash];

extern int original_len[num_read];
extern char original[num_read][int(len_read*2)];
void init_rawread()
{
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i];j++)
		{
			rawread[i][j].base='0';
			rawread[i][j].new_stat=0;
		}
	}
}
void init_insertion()
{
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i];j++)
		{
			insertion[i][j]=0;
		}
	}
}
void init_deletion()
{
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i];j++)
		{
			deletion[i][j]=0;
		}
	}
}

int num_of_found_insertion=0;
int num_of_found_deletion=0;
int num_of_found_complex=0;

int kmer_value;                                                       //kmer appearance times , to decide if it is a valid match

int search_valid_short_kmer(char* kmer)
{
	int index_num;
	index_num=get_index_num(kmer);

	if (head_of_short_hash_table[index_num].next!=NULL)
	{
		struct k_mer_node *node;
		node=head_of_short_hash_table[index_num].next;        
		while(node!=NULL)
		{
			int m;
			for(m=0;m<short_k;m++)
			{
				if(kmer[m]==sample[node->i_read][node->j_position+m])					
					continue;
				else
					break;
			}
			if(m==short_k)
			{
				kmer_value=count_of_short_node[node->i_read][node->j_position];
				return kmer_value;
			}
			else
			{
				node=node->next;
			}
		}
	}
}

int search_valid_long_kmer(char* kmer)
{
	int index_num;
	index_num=get_index_num(kmer);

	if (head_of_long_hash_table[index_num].next!=NULL)
	{
		struct k_mer_node *node;
		node=head_of_long_hash_table[index_num].next;        
		while(node!=NULL)
		{
			int m;
			for(m=0;m<long_k;m++)
			{
				if(kmer[m]==sample[node->i_read][node->j_position+m])					
					continue;
				else
					break;
			}
			if(m==long_k)
			{
				kmer_value=count_of_long_node[node->i_read][node->j_position];
				return kmer_value;
			}
			else
			{
				node=node->next;
			}
		}
	}
}

int search_valid_medium_kmer(char* kmer)
{
	int index_num;
	index_num=get_index_num(kmer);

	if (head_of_medium_hash_table[index_num].next!=NULL)
	{
		struct k_mer_node *node;
		node=head_of_medium_hash_table[index_num].next;        
		while(node!=NULL)
		{
			int m;
			for(m=0;m<medium_k;m++)
			{
				if(kmer[m]==sample[node->i_read][node->j_position+m])					
					continue;
				else
					break;
			}
			if(m==medium_k)
			{
				kmer_value=count_of_medium_node[node->i_read][node->j_position];
				return kmer_value;
			}
			else
			{
				node=node->next;
			}
		}
	}
}

int distin_ins_del(int i ,int j )
{
	char ch1[short_k],ch2[short_k];                                  //ch1 stands for the first kmer that kmer_value=1,ch2 stands for the second kmer that kmer_value=1, if ch1 and ch2 both have high value in hash table, it is insertion error ,else it is deletion error
	int m;
	for( m = 0 ; m< short_k-1;m++)
	{
		ch1[m]=sample[i][j-short_k+1+m];
	}
	ch1[m]=sample[i][j+1];

	for( m = 0 ; m< short_k-2;m++)
	{
		ch2[m]=sample[i][j-short_k+2+m];
	}
	ch2[m]=sample[i][j+1];
	m++;
	ch2[m]=sample[i][j+2];

	int kmer_value1,kmer_value2;
	int flag=0;                                                       //flag=1 stands for insertion error ;flag=2 stands for deletion error
	kmer_value1=search_valid_short_kmer(ch1);
	kmer_value2=search_valid_short_kmer(ch2);
	if(kmer_value1>=count_value&&kmer_value2>=count_value)            //insertion error
	{
		flag=1;
		//insertion[i][j]=1;
	}
	else                                                              //deletion error
	{
		flag=2;
		//deletion[i][j]=1;                                           //move backwards subread
		//num_of_found_deletion++;
	}
	return flag;
}


/*
void get_error_position()
{
	int count=0;
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-short_k;j++)
		{
			if(short_k_mer_times[i][j].times<=valid_value)                       //use num of 1 to find error
			{
				count++;
				continue;
			}
			else
			{
				if(count==0)
				{
					continue;
				}
				if(count==short_k)                                     //num of 1 = short_k, insertion
				{
					num_of_found_insertion++;
					insertion[i][j-1]=1;
					count=0;
					continue;
				}
				if(count==short_k-2)
				{
					if(sample[i][j]!=sample[i][j+1])                   //deletion error
					{
						num_of_found_deletion++;
						deletion[i][j+1]=1;                            //add a random base to position j,and move bases backward 
						count=0;
						continue;
					}
					else                                               //need further judgement
					{
						num_of_found_insertion++;
						insertion[i][j]=1;
						count=0;
						continue;
					}
				}
				if(count==short_k-1)                                   //num of 1 = short_k-1, insertion/deletion
				{	
					if(sample[i][j-1]!=sample[i][j])                   //deletion error
					{
						num_of_found_deletion++;
						deletion[i][j]=1;                              //add a random base to position j,and move bases backward 
						count=0;
						continue;
					}
					else                                               //need further judgement
					{
						int flag=0;
						flag=distin_ins_del( i , j);
						if(flag==1)
						{
							insertion[i][j]=1;
						}
						if(flag==2)
						{
							deletion[i][j]=1;                          //add a random base to position j,and move backward bases
						}
					}
				}
				count=0;
			}
		}
	}
}

*/


void save_as_insertion()
{
	ofstream fout4;
	fout4.open("E:\\insertion.txt",ios::trunc);

	fout4<<">target"<<" "<<endl;
	for(int i=0;i<deleted_read_num1;i++)
	{
		fout4<<"NO."<<i<<"   ";
		for(int j=0;j<deleted_len1[i];j++)
		{
			if(insertion[i][j]!=0)
			{
				fout4<<i<<","<<j<<endl;
			}
		}
		fout4<<endl;
		fout4<<endl;
	}
}
void save_as_deletion()
{
	ofstream fout4;
	fout4.open("E:\\deletion.txt",ios::trunc);

	//fout4<<">target"<<" "<<endl;
	for(int i=0;i<deleted_read_num1;i++)
	{
		fout4<<"NO."<<i<<"   ";
		for(int j=0;j<deleted_len1[i];j++)
		{
			if(deletion[i][j]!=0)
			{
				fout4<<i<<","<<j<<" ";
			}
		}
		fout4<<endl;
		fout4<<endl;
	}
}
/*
void save_as_complex()
{
	ofstream fout4;
	fout4.open("E:\\complex.txt",ios::trunc);

	fout4<<">target"<<" "<<endl;
	for(int i=0;i<num_read;i++)
	{
		fout4<<"NO."<<i<<"   ";
		for(int j=0;j<length[i];j++)
		{
			if(complex[i][j].flag!=0)
			{
			fout4<<i<<","<<j<<","<<complex[i][j].len<<endl;
			}
		}
		fout4<<endl;
		fout4<<endl;
	}
}
*/


void subf_correct(int i)
{
	int pos=0;
		for(int j=0;j<deleted_len1[i];j++)
		{
			if(insertion[i][j]==1)
			{
				//sample[i][m]=sample[i][m+1];
				rawread[i][j].base=deleted_read1[i][j];
				rawread[i][j].new_stat=-1;
			}
			else if(deletion[i][j]==1)
			{
				char ch;                                   
				switch (rand()%4)
				{
				case 0:
					ch='A';
					break;
				case 1:
					ch='T';
					break;
				case 2:
					ch='C';
					break;
				case 3:
					ch='G';
					break;
				}
				rawread[i][j].base=deleted_read1[i][j];
				rawread[i][j].new_stat=-2;
				corrected_read[i][pos]=ch;
				pos++;
			}
			else
			{

				rawread[i][j].base=deleted_read1[i][j];
				rawread[i][j].new_stat=pos;
				corrected_read[i][pos]=deleted_read1[i][j];
				pos++;

			}
		}
		corrected_len[i]=pos;
}

void correct()
{
	for(int i=0;i<deleted_read_num1;i++)
	{
		subf_correct(i);
	}
}

void check_insertion(int i,int j,int count ,char* tmp)
{
	int flag=0;
	char ch1[long_k]={0},ch2[long_k]={0};
	int m;
	//==================ch1======================
	for( m = 0 ; m< long_k;m++)
	{
		if(m<count-1)
		{
		ch1[m]=sample[i][j+m];
		}
		else 
		{
			ch1[m]=sample[i][j+m+1];
		}
	}
	
	//==================ch2======================
	for( m = 0 ; m< long_k-1;m++)
	{
		ch2[m]=ch1[m+1];
	}
	ch2[m]=sample[i][j+long_k+1];
	//==================get k-mer value===================
	int kmer_value1,kmer_value2;
	kmer_value1=search_valid_long_kmer(ch1);
	kmer_value2=search_valid_long_kmer(ch2);
	if(kmer_value1>=valid_value&&kmer_value2>=valid_value)            //insertion error
	{
		flag=1;
		for(int n=0;n<count;n++)
		{
			corrected_read[i][corrected_len[i]]=ch1[n];
			corrected_len[i]++;
		}
	}
	else                                                             
	{
		for(int n=0;n<count+1;n++)
		{
			corrected_read[i][corrected_len[i]]=ch1[n];
			corrected_len[i]++;
		}
	}
}

void check_deletion(int i,int j,int count ,char* tmp)
{
	//================================================================get the corrected base at [i][j+count]==================================================================
	char ch0[short_k-1]={0};                          //first k-1 base
	int score_A=0,score_T=0,score_C=0,score_G=0;         //score of the deletion base
	char ch1[long_k]={0},ch2[long_k]={0},ch3[long_k]={0},ch4[long_k]={0};
	int flag=0;
	
	for(int n=0;n<long_k;n++)
	{
		if(n<count)
		{
			ch1[n]=tmp[n];
			ch2[n]=tmp[n];
			ch3[n]=tmp[n];
			ch4[n]=tmp[n];
		}
		else if(n==count)
		{
			ch1[n]='A';
			ch2[n]='T';
			ch3[n]='C';
			ch4[n]='G';
		}
		else
		{
			ch1[n]=sample[i][j+n-1];
			ch2[n]=sample[i][j+n-1];
			ch3[n]=sample[i][j+n-1];
			ch4[n]=sample[i][j+n-1];
		}
		
	}
	
	//==================get k-mer value===================
	score_A=search_valid_long_kmer(ch1);
	score_T=search_valid_long_kmer(ch2);
	score_C=search_valid_long_kmer(ch3);
	score_G=search_valid_long_kmer(ch4);
	//==================get corrected deleted base===========
	char deleted_base='0';
	if(score_A>=score_C  &&  score_A>=score_G  &&  score_A>=score_T  && score_A>=valid_value)
	{
		deleted_base='A';
		flag=1;
	}
	else if(score_T>=score_A  &&  score_T>=score_G  &&  score_T>=score_C  && score_T>=valid_value)
	{
		deleted_base='T';
		flag=2;
	}
	else if(score_C>=score_A  &&  score_C>=score_G  &&  score_C>=score_T  && score_C>=valid_value)
	{
		deleted_base='C';
		flag=3;
	}
	
	else if(score_G>=score_C  &&  score_G>=score_A  &&  score_G>=score_T  && score_G>=valid_value)
	{
		deleted_base='G';
		flag=4;
	}
	if(flag!=0)            //insertion error
	{
		for(int n=0;n<count+2;n++)
		{
			if(flag==1)
			{
				corrected_read[i][corrected_len[i]]=ch1[n];
			}
			if(flag==2)
			{
				corrected_read[i][corrected_len[i]]=ch2[n];
			}
			if(flag==3)
			{
				corrected_read[i][corrected_len[i]]=ch3[n];
			}
			if(flag==4)
			{
				corrected_read[i][corrected_len[i]]=ch4[n];
			}
			corrected_len[i]++;
		}
	}
	else                                                              //deletion error
	{
		for(int n=0;n<count+1;n++)
		{
			corrected_read[i][corrected_len[i]]=sample[i][j+n];
			corrected_len[i]++;
		}
	}
	flag=0;
}

char ch[short_k]={0};
int j_tail=0;

int get_next_short_k_mer(int i,int j_tail,char* complex)
{
	int n=0;
	for(n=0;n<short_k-1;n++)
	{
		ch[n]=complex[n+1];
	}
	ch[n]=sample[i][j_tail];
	return j_tail;
}

int find_left_complex(int i,int j ,char * complex ,int count)
{
	int tmp=count_value;
	int n=0;
	for(n=0;n<short_k-1;n++)
	{
		ch[n]=complex[medium_k-short_k+1+n];
	}
	ch[n]=sample[i][j+medium_k+1];
	tmp=search_valid_short_kmer(ch);
	j_tail=j+medium_k+2;
	while(tmp>=count_value)
	{
		j_tail=get_next_short_k_mer(i,j_tail,ch);
		tmp=search_valid_short_kmer(ch);
	}
	//base at position j is error. now distinguish insertion or deletion 
}

void check_complex(int i,int j,int count,char* tmp)
{
	int flag=0;
	char ch1[medium_k]={0},ch2[medium_k]={0};
	int m;
	//==================ch1======================
	for( m = 0 ; m< medium_k;m++)
	{
		if(m<short_k-1)
		{
			ch1[m]=sample[i][j+m];
		}
		else 
		{
			ch1[m]=sample[i][j+m+1];
		}
	}
	//==================ch2======================
	//for( m = 0 ; m< medium_k-1;m++)
	//{
	//	ch2[m]=ch1[m+1];
	//}
	//ch2[m]=sample[i][j+medium_k+1];
	//==================get k-mer value===================
	int kmer_value1,kmer_value2;
	kmer_value1=search_valid_medium_kmer(ch1);
	//kmer_value2=search_valid_medium_kmer(ch2);
	if(kmer_value1>=valid_value) 
	//if(kmer_value1>=valid_value&&kmer_value2>=valid_value)            //insertion error
	{
		flag=1;
		for(int n=0;n<short_k;n++)
		{
			corrected_read[i][corrected_len[i]]=ch1[n];
			corrected_len[i]++;
		}
		int left_count=0;
		left_count=find_left_complex(i,j,ch1,count);
	}
	else                                                             
	{
		for(int n=0;n<count+1;n++)
		{
			corrected_read[i][corrected_len[i]]=ch1[n];
			corrected_len[i]++;
		}
	}
}

void subf_position(int i)
{
	int len=0;
	int count=0;
	int flag=0;
	char tmp1[short_k*20]={0};
	for(int j =0 ;j<length[i];j++)
	{
		if(j<length[i]-short_k)
		{
			if(count_of_short_node[i][j]<count_value)
			{	
				tmp1[count]=sample[i][j];											//tmp1 used to store "count_of_short_node=1" region
				count++;
			}
			else//count_of_node[i][j]<count_value
			{
				if(count==0)
				{
					//deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
					corrected_read[i][corrected_len[i]]=sample[i][j];
					corrected_len[i]++;
					//len++;
					//deleted_len[deleted_read_num]++;
				}
				else
				{
					if(count>short_k)
					{
						//deleted_len[deleted_read_num]=len;
						//deleted_ori[deleted_read_num]=i;
						//deleted_read_num++;
						//deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
						//deleted_len[deleted_read_num]++;
						check_complex(i,j-count,count,tmp1);
						memset(tmp1,0,sizeof(char)*short_k*20);
						count=0;
						continue;
					}
					if(count==short_k)                                     //num of 1 = short_k, insertion
					{
						/*
						num_of_found_insertion++;

						for(int n=0;n<count;n++)
						{
							deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=tmp1[n];
							deleted_len[deleted_read_num]++;
						}
						insertion[deleted_read_num][deleted_len[deleted_read_num]-1]=1;
						//deleted_read_num++;
						deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
						deleted_len[deleted_read_num]++;
						*/
						check_insertion(i,j-count,count,tmp1);
						memset(tmp1,0,sizeof(char)*short_k*20);
						count=0;
						continue;
					}
					if(count==short_k-2)
					{
						if(sample[i][j]!=sample[i][j+1])                   //deletion error
						{
							/*
							num_of_found_deletion++;

							for(int n=0;n<count;n++)
							{
								deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=tmp1[n];
								deleted_len[deleted_read_num]++;
							}
							deletion[deleted_read_num][deleted_len[deleted_read_num]+1]=1;   							//add a random base to position j,and move bases backward 
							//deleted_read_num++;
							deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
							deleted_len[deleted_read_num]++;
							*/
							check_deletion(i,j-count,count,tmp1);
							memset(tmp1,0,sizeof(char)*short_k*20);

							//count=0;
							//continue;
						}
						else                                               //need further judgement
						{
							/*
							num_of_found_insertion++;

							for(int n=0;n<count;n++)
							{
								deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=tmp1[n];
								deleted_len[deleted_read_num]++;
							}
							insertion[deleted_read_num][deleted_len[deleted_read_num]]=1;
							//deleted_read_num++;
							deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
							deleted_len[deleted_read_num]++;*/
							check_insertion(i,j-count,count,tmp1);
							memset(tmp1,0,sizeof(char)*short_k*20);
							//count=0;
							//continue;
						}
					}
					if(count==short_k-1)                                   //num of 1 = short_k-1, insertion/deletion
					{	
						if(sample[i][j-1]!=sample[i][j])                   //deletion error
						{
							/*
							num_of_found_deletion++;

							for(int n=0;n<count;n++)
							{
								deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=tmp1[n];
								deleted_len[deleted_read_num]++;
							}
							deletion[deleted_read_num][deleted_len[deleted_read_num]]=1;      //add a random base to position j,and move bases backward 
							//deleted_read_num++;
							deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
							deleted_len[deleted_read_num]++;
							*/
							check_deletion(i,j-count,count,tmp1);

							memset(tmp1,0,sizeof(char)*short_k*20);
							//count=0;
							//continue;
						}
						else                                               //need further judgement
						{
							int flag=0;
							flag=distin_ins_del( i , j);
							if(flag==1)
							{
								/*
								for(int n=0;n<count;n++)
								{
									deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=tmp1[n];
									deleted_len[deleted_read_num]++;
								}
								insertion[deleted_read_num][deleted_len[deleted_read_num]]=1;
								//	deleted_read_num++;
								deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
								deleted_len[deleted_read_num]++;
								*/
								check_insertion(i,j-count,count,tmp1);
								memset(tmp1,0,sizeof(char)*short_k*20);
							}
							if(flag==2)
							{		
								/*
								for(int n=0;n<count;n++)
								{
									deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=tmp1[n];
									deleted_len[deleted_read_num]++;
								}
								deletion[deleted_read_num][deleted_len[deleted_read_num]]=1;                          //add a random base to position j,and move backward bases
								//	deleted_read_num++;
								deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
								deleted_len[deleted_read_num]++;
								*/
								check_deletion(i,j-count,count,tmp1);
								memset(tmp1,0,sizeof(char)*short_k*20);
							}
						}
					}
					else if(count<short_k-2)
					{
						for(int n=0;n<count;n++)
						{
							corrected_read[i][corrected_len[i]]=tmp1[n];
							corrected_len[i]++;
						}
						memset(tmp1,0,sizeof(char)*short_k*20);
					}
					count=0;			
				}		
			}			
		}
		//	=========================================================================需要处理count！=0情况
		else
		{
			if(count!=0)
			{
				for(int n=0;n<count;n++)
				{
					deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=tmp1[n];
					deleted_len[deleted_read_num]++;
				}
				count=0;
			}			
			deleted_read[deleted_read_num][deleted_len[deleted_read_num]]=sample[i][j];
			deleted_len[deleted_read_num]++;
		}
	}
	count=0;
	deleted_ori[deleted_read_num]=i;
	deleted_read_num++;
	len=0;	
}

void delete_invalid_region()
{
	for(int i =0 ;i<num_read;i++)
	{
		subf_position(i);
	}
}

void save_as_deleted_reads()
{
	ofstream fout1;
	fout1.open("E:\\deleted_reads_for_MOCOV.fasta",ios::trunc);
	int k=0;
	for(int i=0;i<deleted_read_num;i++)
	{
		if(deleted_len[i]>=len_read*0.2)
		{
			fout1<<">"<<k<<"_"<<deleted_ori[i]<<"/0_"<<deleted_len[i]<<endl;
			for(int j=0;j<deleted_len[i];j++)
			{
				deleted_len1[k]=deleted_len[i];
				deleted_read1[k][j]=deleted_read[i][j];
				
				fout1<<deleted_read[i][j];
				//fout1<<j<<corrected_read[i][j]<<" ";
			}
			fout1<<endl;
			k++;
		}
	}
	deleted_read_num1=k;
}

void save_as_reads()
{
	ofstream fout1;
	fout1.open("E:\\corrected_reads_for_MOCOV.fasta",ios::trunc);
	for(int i=0;i<deleted_read_num1;i++)
	{
		fout1<<">"<<i<<"/0_"<<corrected_len[i]<<endl;
		for(int j=0;j<corrected_len[i];j++)
		{
			fout1<<corrected_read[i][j];
		}
		fout1<<endl;
	}
}

void save_as_compare()
{
	ofstream fout1;
	fout1.open("E:\\compare_reads_for_MOCOV.fasta",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout1<<">"<<i<<"\t"<<original_len[i]<<"\t";
		for(int j=0;j<original_len[i];j++)
		{
			fout1<<original[i][j];
		}
		fout1<<endl;
		fout1<<">"<<i<<"\t"<<length[i]<<"\t";
		for(int j=0;j<length[i];j++)
		{
			fout1<<sample[i][j];
		}
		fout1<<endl;
		fout1<<">"<<i<<"\t"<<corrected_len[i]<<"\t";
		for(int j=0;j<corrected_len[i];j++)
		{
			fout1<<corrected_read[i][j];
		}
		fout1<<endl;
	}
}

extern void correct_errors()
{
	init_rawread();
	init_insertion();
	init_deletion();
	//init_complex();
	get_error_position();
	//delete_invalid_region();
	//save_as_deleted_reads();
	correct();
	//output_insertion();
	//output_deletion();
	//output_complex();
	save_as_insertion();
	save_as_deletion();
	//save_as_complex();
	
	save_as_reads();
	//save_as_compare();
	
	
}