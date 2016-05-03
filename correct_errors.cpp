#include"dataset.h"
#include"hash_table.h"
#include"correct_errors.h"
extern struct times_and_flag;

extern times_and_flag k_mer_times[num_read][int(len_read*2)];

int error_result[num_read][int(len_read*2)];
//int head;
//int tail;

typedef struct com
{
	int flag;
	int len;
}com;

int insertion[num_read][int(len_read*2)];
int deletion[num_read][int(len_read*2)];
com complex[num_read][int(len_read*2)];

extern head_node_index head_of_hash_table[max_len_head_hash];//哈希表头节点数组

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
void init_complex()
{
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i];j++)
		{
			complex[i][j].flag=0;
			complex[i][j].len=0;
		}
	}
}
/*
void find_head_and_tail(int pos_i,int pos_j,int len)
{
	if(k_mer_times[pos_i][pos_j]>invalid_value)
	{
		while((k_mer_times[pos_i][pos_j]>invalid_value)&&(pos_j<(len-short_k+1)))//如果开头可匹配
		{
			pos_j++;
			continue;
		}
		head=pos_j;
		while((k_mer_times[pos_i][pos_j]<=invalid_value)&&(pos_j<(len-short_k+1)))
		{
			pos_j++;
			continue;
		}
		tail=pos_j;
	}
	else
	{
		head=pos_j-1-short_k;
		while ((k_mer_times[pos_i][pos_j]<=invalid_value)&&(pos_j<(len-short_k+1)))
		{
			pos_j++;
			continue;
		}
		tail=pos_j;
	}
}

void get_error_position()
{
	//将匹配失败处的错误分类：insertion，deletion，complex
	int len_insertion;
	int begin_position;
	for(int i=0;i<num_read;i++)
	{
		int j=0;
		while(j<length[i])
		{
			find_head_and_tail(i ,j ,length[i]);
			len_insertion=tail-head-short_k;
			if(len_insertion<0)
			{
				continue;
			}
			else
			{
				if(len_insertion==0)//deletion
				{
					deletion[i][j]=1;
				}
				else
				{
					if(len_insertion==1)//insertion
					{
						insertion[i][j]=1;
					}
					else
					{
						complex[i][j]=len_insertion;
					}
				}
			}
			j=tail+1;
		}
	}
}
*/

int num_of_found_insertion=0;//检测到的insertion个数
int num_of_found_deletion=0;//检测到的deletion个数
int num_of_found_complex=0;//检测到的complex个数

/*
void get_error_position()
{
	int count=0;
	for(int i=0;i<num_read;i++)
	{
		for(int j=short_k;j<length[i]-1;j++)
		{
			if(k_mer_times[i][j].times==0)//0的个数判断插入/复杂错误
			{
				count++;
				continue;
			}
			else
			{
				if(count==1)//0的个数为1，是插入错误
				{
					num_of_found_insertion++;
					insertion[i][j-1]=1;
					count=0;
				}
				else
				{
					if(count>1)//0的个数>1，是复杂错误
					{	
						complex[i][j-count].flag=1;
						complex[i][j-count].len=count;
						num_of_found_complex++;
						count=0;
					}
					else//0的个数为0，判断是不是删除错误
					{
						if(k_mer_times[i][j].flag==0&&k_mer_times[i][j+1].flag==1)
						{
							deletion[i][j]=1;
							num_of_found_deletion++;
						}
						else
						{
							continue;
						}
					}
				}
			}
		}
	}
}
*/
extern int count_of_node[num_read][int(len_read*2)];//节点计数数组

extern int get_index_num(char *ch);

int kmer_value;//kmer出现的次数，用来区分是否有效匹配


int search_valid_kmer(char* kmer)
{
	int index_num;
	index_num=get_index_num(kmer);

	if (head_of_hash_table[index_num].next!=NULL)
	{
		struct k_mer_node *node;
		node=head_of_hash_table[index_num].next;//old_node指向index_num链后的第一个节点
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
				kmer_value=count_of_node[node->i_read][node->j_position];
				return kmer_value;
			}
			else//若没有扫描到相同k-mer，继续向后
			{
				node=node->next;
			}
		}
	}
}

int distin_ins_del(int i ,int j )
{
	char ch1[short_k],ch2[short_k];//ch1 和ch2 表示从第一个“1”开始和第二个“1”开始对应的kmer，如果都能在哈希表中对应到高频率，则是插入错误，否则是删除错误
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
	int flag=0;//flag=1表示插入，2表示删除
	kmer_value1=search_valid_kmer(ch1);
	kmer_value2=search_valid_kmer(ch2);
	if(kmer_value1>=valid_value&&kmer_value2>=valid_value)//此处是插入错误
	{
		flag=1;
		//insertion[i][j]=1;
	}
	else//是删除错误
	{
		flag=2;
		//deletion[i][j]=1;//在j位置插入，j及以后向后移动。
		//num_of_found_deletion++;
	}
	return flag;
}

void get_error_position()
{
	int count=0;
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-short_k;j++)
		{
			if(count_of_node[i][j]<=valid_value)//1的个数判断插入/复杂错误
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
				if(count==short_k)//1的个数为5，是插入错误
				{
					num_of_found_insertion++;
					insertion[i][j-1]=1;
					count=0;
					continue;
				}
				if(count==short_k-2)
				{
					if(sample[i][j]!=sample[i][j+1])//deletion error
					{
						num_of_found_deletion++;
						deletion[i][j+1]=1;//纠错时在j位置插入碱基，j+1及其以后的碱基向后移动。
						count=0;
						continue;
					}
					else//需要判断是插入还是删错
					{
						num_of_found_insertion++;
						insertion[i][j]=1;
						count=0;
						continue;
					}
				}
				if(count==short_k-1)//1的个数为4，是删除或插入错误
				{	
					if(sample[i][j-1]!=sample[i][j])//deletion error
					{
						num_of_found_deletion++;
						deletion[i][j]=1;//纠错时在j位置插入碱基，j及其以后的碱基向后移动。
						count=0;
						continue;
					}
					else//需要判断是插入还是删错
					{
						int flag=0;
						flag=distin_ins_del( i , j);
						if(flag==1)
						{
							insertion[i][j]=1;
						}
						if(flag==2)
						{
							deletion[i][j]=1;//在j位置插入，j及以后向后移动。
						}
					}
				}
				count=0;
			}
		}
	}
}



void save_as_insertion()
{
	ofstream fout4;
	fout4.open("E:\\insertion.txt",ios::trunc);

	fout4<<">target"<<" "<<endl;
	for(int i=0;i<num_read;i++)
	{
		fout4<<"NO."<<i<<"   ";
		for(int j=0;j<length[i];j++)
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
	for(int i=0;i<num_read;i++)
	{
		fout4<<"NO."<<i<<"   ";
		for(int j=0;j<length[i];j++)
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

extern int num_of_insertion=0;


void correct()
{
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i];j++)
		{
			if(insertion[i][j]==1)
			{
				for(int m=j;m<length[i];m++)//向插入错误位置删除一个碱基
				{
					sample[i][m]=sample[i][m+1];
				}
				for(int k=j+1;k<length[i];k++)//把该行后面的错误向前移
				{
					if(insertion[i][j]==1)
					{
						insertion[i][j-1]=1;
						insertion[i][j]=0;
					}
					if(deletion[i][j]==1)
					{
						deletion[i][j-1]=1;
						insertion[i][j]=0;
					}
				}

				num_of_insertion--;
				length[i]--;
			}
			if(deletion[i][j]==1)
			{
				char ch;//随机插入的碱基
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
				for(int m=length[i];m>j+1;m--)//将选定位置后面的碱基向后移动
				{
					sample[i][m]=sample[i][m-1];
				}
				sample[i][j+1]=ch;//插入一个碱基
				length[i]++;
			}
			if(complex[i][j].flag==1)
			{
				for(int m=j;m<length[i];m++)//向插入错误位置删除一个碱基
				{
					sample[i][m]=sample[i][m+1];
				}
				num_of_insertion--;
				length[i]--;
			}
		}
	}
}


void save_as_reads()
{
	ofstream fout1;
	fout1.open("E:\\reads_and_times.fasta",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout1<<">"<<i<<" "<<position[i]<<endl;
		for(int j=0;j<length[i];j++)
		{
			fout1<<j<<sample[i][j]<<" "<<count_of_node[i][j]<<" , ";
		}
		fout1<<endl;
	}
}


extern void correct_errors()
{
	init_insertion();
	init_deletion();
	init_complex();
	get_error_position();
	//output_insertion();
	//output_deletion();
	//output_complex();
	save_as_insertion();
	save_as_deletion();
	save_as_complex();
	save_as_reads();
}