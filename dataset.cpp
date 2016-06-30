#include"main.h"
//#include"hash_table.h"
#include"dataset.h"

char genome[len_genome+1];                                 //genome
int position[num_read];                                    //array for long read starting positions
int length[num_read];                                      //array for every long read's length
char sample[num_read][int(len_read*2)];                    //array for store every long read with error
int total_num_of_bases=0 ;                                 //total num of bases in sample
//k_mer array_k_mer[len_array_k_mer];
int num_short_k_mer=0;                                     //number of k-mers for short match
char original[num_read][int(len_read*2)];  
int original_len[num_read]={0};

int read_count=0;
void initiate_sample()
{
	for(int i=0;i<num_read;i++)                           
	{
		for(int j=0;j<len_read;j++)
		{
			sample[i][j]='0';
		}
	}
}
/*
void initiate_array_k_mer()
{
	for(int i=0;i<len_array_k_mer;i++)
	{
		array_k_mer[i].ch.assign(" ");
		array_k_mer[i].count_ch=0;
	}
}
*/
unsigned long ulrand(void) {
    return (
     (((unsigned long)rand()<<24)&0xFF000000ul)
    |(((unsigned long)rand()<<12)&0x00FFF000ul)
    |(((unsigned long)rand()    )&0x00000FFFul));
}

char* rand_genome(char* str, const int len)
{
	int i;
	for(i=0;i<len;i++)
	{
		switch (rand()%4)
		{
		case 0:
			str[i]='A';
			break;
		case 1:
			str[i]='T';
			break;
		case 2:
			str[i]='C';
			break;
		case 3:
			str[i]='G';
			break;
		}
	}
	str[++i]='\0';
	
	return str;
}
void output_genome()
{
	cout<<"genome:"<<endl;
	for(int i=0;i<len_genome;i++)
	{
		cout<<i <<genome[i]<<' ';
	}
	cout<<endl;
	//cout<<str<<endl;
}

int* rand_position(int* str, const int len)//生成随机起始点
{
	int i;
	unsigned long ul;
	for(i=0;i<len;i++)
	{
		ul=ulrand();
		str[i]=ul%len_genome;
	}
	return str;
}
void output_position()
{
	cout<<"random position："<<endl;
	int i;
	for(i=0;i<num_read;i++)
	{
		cout<<position[i]<<' ';
	}
	cout<<endl;
}

void rand_length()
{
	int i;
	for(i=0;i<num_read;i++)
	{
		//length[i]=10000;
		original_len[i]=rand()%((int )(0.9*len_read))+((int )(0.1*len_read));
		//str[i]=rand()%(len_read-1000)+1000;//str[i]=rand()%9000+1000;此处为1~10
	}
}

void rand_read()
{
	int i,j;
	total_num_of_bases=0;
	for(i=0;i<num_read;i++)
	{
		if((position[i]+original_len[i])>len_genome)
		{
			original_len[i]=len_genome-position[i];
		}
		for(j=0;j<original_len[i];j++)
		{
			original[i][j]=genome[position[i]+j];
		}
		total_num_of_bases+=original_len[i];
		//rand_read(sample[i],genome,position[i],length[i]);
	}
}

void save_as_original()
{
	ofstream fout1;
	fout1.open("E:\\original_reads_for_MOCOV.fasta",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout1<<">No."<<i<<" position="<<position[i]<<" length="<<original_len[i]<<endl;
		for(int j=0;j<original_len[i];j++)
		{
			fout1<<original[i][j];
		}
		fout1<<endl;
	}
}

void copy_length()
{
	for(int i=0;i<num_read;i++)
	{
		length[i]=original_len[i];
	}
}
void copy_read()
{
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i];j++)
		{
			sample[i][j]=original[i][j];
		}
	}
}

void replace()
{
	int i=0,j=0;                                                        //the ith read ,the jth positiong,sample[i][j]is the base to be replaced
	int num_of_replacement=0;
	while(num_of_replacement<total_num_of_bases*0.01)
	{
		//cout<<"num_of_replace: "<<num_of_replace<<endl;
		i=rand()%num_read;
		j=rand()%length[i];
		char ch;                                                        //generate a random base to replace the original sample[i][j]
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
		if(sample[i][j]!=ch)
		{
			//cout<<"sample["<<i<<"]["<<j<<"] replace"<<sample[i][j]<<"->"<< ch<<endl;
			sample[i][j]=ch;
			num_of_replacement++;
		}
		else
		{
			continue;
		}
	}
}

void insert()
{
	ofstream fout4;
	fout4.open("E:\\insert.txt",ios::trunc);
	int num_of_insertion=0;

	int insert_i=0,insert_j=0;                                           //the ith read ,the jth positiong,sample[i][j]is the base to be inserted
	
	while(num_of_insertion<total_num_of_bases*0.03)
	{
		//cout<<"num_of_insertion: "<<num_of_insertion<<endl;
		insert_i=rand()%num_read;
		//insert_j=rand()%(length[insert_i]);
		if(length[insert_i]>(2*short_k))
		{
			insert_j=rand()%(length[insert_i]-short_k*2)+short_k;
		}
		else
			continue;
		char ch;                                                         //generate a random base to insert to the position original sample[i][j]
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
		fout4<<"sample["<<insert_i<<"]["<<insert_j<<"] insert"<<ch<<endl;
		int j=length[insert_i];
		for(;j>insert_j;j--)                                             //move the subread after insert-position backward
		{
			sample[insert_i][j]=sample[insert_i][j-1];
		}
		sample[insert_i][j]=ch;
		num_of_insertion++;
		length[insert_i]++;
	}
	fout4<<"num of insertion = "<<num_of_insertion<<endl;
}

void deletion()
{
	ofstream fout5;
	fout5.open("E:\\delete.txt",ios::trunc);

	int deletion_i=0,deletion_j=0;                                       //the ith read ,the jth positiong,sample[i][j]is the base to be inserted
	int num_of_deletion=0;
	//cout<<"total_num_of_bases: "<<total_num_of_bases<<endl;
	while(num_of_deletion<total_num_of_bases*0.03)
	{
		//cout<<"num_of_deletion: "<<num_of_deletion<<endl;
		deletion_i=rand()%num_read;
		if(length[deletion_i]>(2*short_k))
		{
			deletion_j=rand()%(length[deletion_i]-short_k*2)+short_k;
		}
		else 
			continue;
		fout5<<"sample["<<deletion_i<<"]["<<deletion_j<<"] delete"<<sample[deletion_i][deletion_j]<<endl;
		int j=deletion_j;
		for(;j<length[deletion_i]-1;j++)                                 //move the subread after delete-position forward
		{
			sample[deletion_i][j]=sample[deletion_i][j+1];
		}
		sample[deletion_i][j]='0';
		num_of_deletion++;
		length[deletion_i]--;
	}
	fout5<<"num of deletion = "<<num_of_deletion<<endl;
}



void save_as_query()
{
	ofstream fout1;
	fout1.open("E:\\query_read_for_MOCOV.fasta",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		//fout1<<">No."<<i<<" position="<<position[i]<<" length="<<length[i]<<endl;
		fout1<<">No."<<i<<"/0_"<<length[i]<<endl;

		for(int j=0;j<length[i];j++)
		{
			fout1<<sample[i][j];
		}
		fout1<<endl;
	}
}

void save_as_target()
{
	ofstream fout2;
	fout2.open("E:\\target.fasta",ios::trunc);

	fout2<<">target\n";
	for(int j=0;j<len_genome;j++)
	{
		fout2<<genome[j];
	}
	fout2<<endl;
}

void read_file()
{
	FILE * f_fp;
	char ch;
	int i=0,j=0;
	int flag=0;
	f_fp=fopen("E:\\query_read_for_MOCOV.fasta","r");
	while (!feof(f_fp))
	{		
		ch=fgetc(f_fp);
		if(ch=='\"')
		{
			continue;
		}
		if(ch=='\n')
		{
			if(flag==0)
			{
				flag=1;
				read_count++;
				continue;
			}
			else
			{
				length[read_count-1]=j;
				j=0;
				flag=0;
				continue;
			}
		}
		if(ch==' ')
		{
			continue;
		}
		if(ch>=48&&ch<=57)
		{
			continue;
		}
		if(ch=='A'||ch=='T'||ch=='C'||ch=='G')
		{
			sample[read_count-1][j]=ch;
			j++;
		}
	}
	if(ch==-1)
	{
		//length[read_count-1]=j;
	}
}



extern void dataset()
{
	srand((int)time(0));
	memset(position,0,sizeof(int)*num_read);
	initiate_sample();
	/*
	rand_genome(genome,len_genome);
	cout<<"rand_genome(genome,len_genome);"<<endl;
	
	rand_position(position, num_read);
	cout<<"rand_position(position, num_read);"<<endl;
	//output_position();
	rand_length();
	cout<<"rand_length(length,num_read);"<<endl;
	//output_read_length();	
	rand_read();
	cout<<"rand_read();"<<endl;
	save_as_original();
	copy_length();
	copy_read();
	//replace();
	insert();
	cout<<"insert();"<<endl;

	
	deletion();
	cout<<"deletion();"<<endl;

	
	save_as_query();
	save_as_target();
	//save_as_kmer_times();
	*/
	read_file();

	
}

