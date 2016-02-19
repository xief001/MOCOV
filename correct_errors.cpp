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

void get_error_position()
{
	int count=0;

	for(int i=0;i<num_read;i++)
	{
		for(int j=short_k;j<length[i];j++)
		{
			count=0;
			if(k_mer_times[i][j].times==0)//0的个数判断插入/复杂错误
			{
				count++;
				continue;
			}
			else
			{
				if(count==1)//0的个数为1，是插入错误
				{
					insertion[i][j-1]=1;
				}
				else
				{
					if(count>1)//0的个数>1，是复杂错误
					{	
						complex[i][j-count].flag=1;
						complex[i][j-count].len=count;
					}
					else//0的个数为0，判断是不是删除错误
					{
						if(k_mer_times[i][j].flag==0&&k_mer_times[i][j+1].flag==1)
						{
							deletion[i][j]=1;
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
}