#include"main.h"

using namespace std;
string query_file_name="0";
string corrected_file_name="0";

int short_k=10;//-s
int len_index=5;//-i
int valid_value=5;//v
int len_k_mer_error=0;//e
int parseargs(int argc, char * argv[])
{
	if(argc<3)
	{
		cout<<"Invalid parameters!"<<endl;
		return -1;
	}
	else
	{
		query_file_name=argv[1];
		corrected_file_name=argv[2];
		for(int i=3;i<argc;i++)
		{
			if(argv[i][0]=='-')
			{
				if(argv[i][1]=='s')
				{
					short_k=atoi(argv[++i]);
					medium_k=short_k+4;
					continue;
				}
				if(argv[i][1]=='l')
				{
					long_k=atoi(argv[++i]);
					continue;
				}
				if(argv[i][1]=='v')
				{
					valid_value=atoi(argv[++i]);
					continue;
				}
				if(argv[i][1]=='e')
				{
					len_k_mer_error=atoi(argv[++i]);
					continue;
				}
			}
		}
	}
}

int read_count=0;
int len_read=0;

int *length=NULL;
int *name_length=NULL;
char **sample = NULL;
char **readname=NULL;

//create_kmer_table
node_position **every_short_kmer_times_position;          
node_position **every_medium_kmer_times_position;       

head_node_index head_of_short_hash_table[1024];                  
head_node_index head_of_medium_hash_table[1024];                  //hash table's head node

int **count_of_short_node;                           //array of counting every kmer (in read) times
int **count_of_medium_node;                           //array of counting every kmer (in read) times

char **corrected_read ;
int *corrected_len;//(int *)malloc(read_count*sizeof(int));

int **insertion;
int **deletion;

void init_every_kmer_times_position()
{
	every_short_kmer_times_position = new node_position*[read_count];
	if(every_short_kmer_times_position!=NULL)
	{
		for(int i =0; i<read_count; i++)
		{
			every_short_kmer_times_position[i] = new node_position[length[i]];
		}
	}
	every_medium_kmer_times_position = new node_position*[read_count];
	if(every_medium_kmer_times_position!=NULL)
	{
		for(int i =0; i<read_count; i++)
		{
			every_medium_kmer_times_position[i] = new node_position[length[i]];
		}
	}
	for(int i=0; i<read_count; i++)
	{
		for(int j=0; j<length[i]; j++)
		{
			every_short_kmer_times_position[i][j].i_position = 0;
			every_short_kmer_times_position[i][j].j_position=0;
			every_medium_kmer_times_position[i][j].i_position = 0;
			every_medium_kmer_times_position[i][j].j_position=0;
		}
	}
}

void init_count_of_node()
{
	count_of_short_node = new int*[read_count];//(int**)malloc(sizeof(int*)*read_count);
	for(int i =0; i<read_count; i++)
	{
		count_of_short_node[i] = new int[length[i]];//(int*)malloc(sizeof(int)*length[i]);
	}
	count_of_medium_node = new int*[read_count];//(int**)malloc(sizeof(int*)*read_count);
	for(int i =0; i<read_count; i++)
	{
		count_of_medium_node[i] = new int[length[i]];//(int*)malloc(sizeof(int)*length[i]);
	}
    for(int i=0; i<read_count; i++)
    {
		for(int j=0; j<length[i]; j++)
        {
			count_of_short_node[i][j] = 0;
			count_of_medium_node[i][j] = 0;
        }
    }
}

void init_head_of_hash_table()
{	
	for(int i =0;i<1024;i++)
	{
		//0
		if(i/256==0)
		{
			head_of_short_hash_table[i].head_hash_table[0]='A';
			head_of_medium_hash_table[i].head_hash_table[0]='A';
		}
		else if(i/256==1)
		{
			head_of_short_hash_table[i].head_hash_table[0]='T';
			head_of_medium_hash_table[i].head_hash_table[0]='T';
		}
		else if(i/256==2)
		{
			head_of_short_hash_table[i].head_hash_table[0]='C';	
			head_of_medium_hash_table[i].head_hash_table[0]='C';
		}
		else if(i/256==3)
		{
			head_of_short_hash_table[i].head_hash_table[0]='G';		
			head_of_medium_hash_table[i].head_hash_table[0]='G';
		}
		//1
		if((i%256)/64==0)
		{
			head_of_short_hash_table[i].head_hash_table[0]='A';
			head_of_medium_hash_table[i].head_hash_table[0]='A';
		}
		else if((i%256)/64==1)
		{
			head_of_short_hash_table[i].head_hash_table[0]='T';
			head_of_medium_hash_table[i].head_hash_table[0]='T';	
		}
		else if((i%256)/64==2)
		{
			head_of_short_hash_table[i].head_hash_table[0]='C';
			head_of_medium_hash_table[i].head_hash_table[0]='C';	
		}
		else if((i%256)/64==3)
		{
			head_of_short_hash_table[i].head_hash_table[0]='G';
			head_of_medium_hash_table[i].head_hash_table[0]='G';	
		}
		//2
		if((i%64)/16==0)
		{
			head_of_short_hash_table[i].head_hash_table[0]='A';
			head_of_medium_hash_table[i].head_hash_table[0]='A';	
		}
		else if((i%64)/16==1)
		{
			head_of_short_hash_table[i].head_hash_table[0]='T';
			head_of_medium_hash_table[i].head_hash_table[0]='T';	
		}
		else if((i%64)/16==2)
		{
			head_of_short_hash_table[i].head_hash_table[0]='C';
			head_of_medium_hash_table[i].head_hash_table[0]='C';	
		}
		else if((i%64)/16==3)
		{
			head_of_short_hash_table[i].head_hash_table[0]='G';
			head_of_medium_hash_table[i].head_hash_table[0]='G';	
		}
		//3
		if((i%16)/4==0)
		{
			head_of_short_hash_table[i].head_hash_table[0]='A';
			head_of_medium_hash_table[i].head_hash_table[0]='A';	 
		}
		else if((i%16)/4==1)
		{
			head_of_short_hash_table[i].head_hash_table[0]='T';
			head_of_medium_hash_table[i].head_hash_table[0]='T';	
		}
		else if((i%16)/4==2)
		{
			head_of_short_hash_table[i].head_hash_table[0]='C';
			head_of_medium_hash_table[i].head_hash_table[0]='C';	
		}
		else if((i%16)/4==3)
		{
			head_of_short_hash_table[i].head_hash_table[0]='G';
			head_of_medium_hash_table[i].head_hash_table[0]='G';	
		}
		//4
		if((i%4)==0)
		{
			head_of_short_hash_table[i].head_hash_table[0]='A';
			head_of_medium_hash_table[i].head_hash_table[0]='A';	
		}
		else if((i%4)==1)
		{
			head_of_short_hash_table[i].head_hash_table[0]='T';
			head_of_medium_hash_table[i].head_hash_table[0]='T';	
		}
		else if((i%4)==2)
		{
			head_of_short_hash_table[i].head_hash_table[0]='C';
			head_of_medium_hash_table[i].head_hash_table[0]='C';	
		}
		else if((i%4)==3)
		{
			head_of_short_hash_table[i].head_hash_table[0]='G';
			head_of_medium_hash_table[i].head_hash_table[0]='G';	
		}
	}
}

void init_corrected_read()
{
	corrected_read = new char*[read_count];            // (char**)malloc(sizeof(char*)*read_count);
	for(int i =0; i<read_count; i++)
	{
		int size=length[i]*2;
		corrected_read[i] = new char[size];        //(char*)malloc(sizeof(char)*(length[i]*2));
	}
	for(int i=0; i<read_count; i++)
	{
		for(int j=0; j<(length[i]*2); j++)
		{
			corrected_read[i][j] ='0';
		}
	}
}

void init_corrected_len()
{
	corrected_len=new int[read_count];
	for(int i=0;i<read_count;i++)
	{
		corrected_len[i]=0;
	}
}

void create_sample(FILE *f_fp)
{
	char read[100000]={'0'};
	char name[200]={'0'};
	char ch;
	int pos=0;
	
	while(!feof(f_fp))
	{
		ch=fgetc(f_fp);
		if(ch==-1)
		{
			break;
		}
		if(ch=='>')
		{
			while(!feof(f_fp))
			{
				ch=fgetc(f_fp);
				if(ch!='\n')
				{
					name[pos]=ch;
					pos++;
				}
				else
				{
					name_length[read_count]=pos;
					readname[read_count]=new char[pos];//(char*)malloc(sizeof(char)*(pos+1));
					for(int j =0 ;j<pos;j++)
					{
						readname[read_count][j]=name[j];	
					}
					read_count++;
					pos=0;
					memset(name,'0',sizeof(name));
					break;
				}
			}
		}
		while(!feof(f_fp))
		{
			ch=fgetc(f_fp);
			if(ch!='\n')
			{
				read[pos]=ch;
				pos++;
			}
			else
			{
				length[read_count-1]=pos;
				sample[read_count-1]= new char[int(pos*1.1)];//(char *)malloc(sizeof(char)*int(pos*1.1));
				for(int j =0 ;j<int(pos*1.1);j++)
				{
					if(j<pos)
					{
						sample[read_count-1][j]=read[j];
					}
					else
					{
						sample[read_count-1][j]='0';
					}
				}
				memset(read,'0',10000);
				pos=0;
				break;
			}
		}
	}
}

void read_file()
{
	sample = new char*[20000];
	//sample = (char**)malloc(sizeof(char*)*20000);
	readname = new char*[20000];
	//readname = (char**)malloc(sizeof(char*)*20000);
	name_length=new int[20000];
	length=new int[20000];
	FILE * f_fp;
	char ch;
	int i=0,j=0;
	int flag=0;
	const char* filename;
	filename=query_file_name.c_str();
	f_fp=fopen(filename,"r");
	if(f_fp==NULL)
	{
		perror("open");
	}
	create_sample(f_fp);
	fclose(f_fp);
}

void save_as_corrected_read()
{
	const char* filename;
	filename=corrected_file_name.c_str();
	ofstream fout1;
	fout1.open(filename,ios::trunc);
	for(int i=0;i<read_count;i++)
	{
		//fout1<<">"<<i<<"/0_"<<corrected_len[i]<<endl;
		for(int j=0;j<name_length[i];j++)
		{
			fout1<<readname[i][j];
		}
		fout1<<endl;
		for(int j=0;j<corrected_len[i];j++)
		{
			fout1<<corrected_read[i][j];
		}
		fout1<<endl;
	}
}

void free_sample()
{
    int i,j;
	for(i=0; i<read_count; i++)
	{
        delete []sample[i];
		sample[i]=NULL;
	}
    delete []sample;
	sample=NULL;
}

void free_readname()
{
    int i,j;
	for(i=0; i<read_count; i++)
	{
        delete []readname[i];
		readname[i]=NULL;
	}
    delete []readname;
	readname=NULL;
}

void free_Array_node_position()
{
    int i,j;
	for(i=0; i<read_count; i++)
	{
		delete []every_short_kmer_times_position[i];
		every_short_kmer_times_position[i]=NULL;
	}
	delete []every_short_kmer_times_position;
	every_short_kmer_times_position=NULL;

	for(i=0; i<read_count; i++)
	{
		delete []every_medium_kmer_times_position[i];
		every_medium_kmer_times_position[i]=NULL;
	}
	delete []every_medium_kmer_times_position;
	every_medium_kmer_times_position=NULL;
}

void free_count_of_node()
{
    int i,j;
	for(i=0; i<read_count; i++)
	{
		delete []count_of_short_node[i];
		count_of_short_node[i]=NULL;
	}
   delete []count_of_short_node;
   count_of_short_node=NULL;

   for(i=0; i<read_count; i++)
	{
		delete []count_of_medium_node[i];
		count_of_medium_node[i]=NULL;
	}
   delete []count_of_medium_node;
   count_of_medium_node=NULL;
}

void free_corrected_read()
{
    int i,j;
	for(i=0; i<read_count; i++)
	{
        delete []corrected_read[i];
		corrected_read[i]=NULL;
	}
    delete []corrected_read;
	corrected_read=NULL;
}

void init_dataset()
{
	init_every_kmer_times_position();
	init_count_of_node();
	init_head_of_hash_table();
	init_corrected_read();
	init_corrected_len();
}

void free_arrays()
{
	free_sample();
	free_readname();
	delete []length;
	length=NULL;
	delete []name_length;    
	name_length=NULL;
	free_Array_node_position();
	free_count_of_node();
	free_corrected_read();
	delete []corrected_len;
	corrected_len=NULL;
}

int main(int argc,char * argv[])
{
	parseargs(argc,argv);
	read_file();
	cout<<"finish reading "<<read_count<<" reads."<<endl;
	init_dataset();
	cout<<"finish initiating dataset."<<endl;
	hash_table();
	//correct_errors();
	//save_as_corrected_read();
	free_arrays();
	//cout<<"finished..."<<endl;
	return 0;
}

