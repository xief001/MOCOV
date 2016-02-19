#include"main.h"
//#include"hash_table.h"
#include"dataset.h"

char genome[len_genome+1];//基因组字符串
int position[num_read];//长read起始位置数组
int length[num_read];//1000条read的随机长度。1K-10K
char sample[num_read][int(len_read*2)];//存放长read，长度扩大为1.1倍，用来存放碱基插入错误
int total_num_of_bases=0 ;//生成的长read中所有碱基总数
k_mer array_k_mer[len_array_k_mer];
int num_short_k_mer=0; //标记短k-mer组中出现的k-mer的数量

void initiate_sample()
{
	for(int i=0;i<num_read;i++) //sample初始化为0
	{
		for(int j=0;j<len_read;j++)
		{
			sample[i][j]='0';
		}
	}
}

void initiate_array_k_mer()
{
	for(int i=0;i<len_array_k_mer;i++)
	{
		array_k_mer[i].ch.assign(" ");
		array_k_mer[i].count_ch=0;
	}
}

unsigned long ulrand(void) {
    return (
     (((unsigned long)rand()<<24)&0xFF000000ul)
    |(((unsigned long)rand()<<12)&0x00FFF000ul)
    |(((unsigned long)rand()    )&0x00000FFFul));
}

char* rand_genome(char* str, const int len)//生成随机基因组
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
	cout<<"随机起点："<<endl;
	int i;
	for(i=0;i<num_read;i++)
	{
		cout<<position[i]<<' ';
	}
	cout<<endl;
}
/*
int* rand_length(int* str, const int num)//生成随机read长度
{
	total_num_of_bases=0;
	int i;
	for(i=0;i<num;i++)
	{
		str[i]=rand()%(len_read-1)+10;
		//str[i]=rand()%(len_read-1000)+1000;//str[i]=rand()%9000+1000;此处为1~100
		total_num_of_bases+=str[i];
	}
	return str;
}
*/
void rand_length()//生成随机read长度
{
	int i;
	for(i=0;i<num_read;i++)
	{
		length[i]=500;
		//length[i]=rand()%(len_read-10)+10;
		//str[i]=rand()%(len_read-1000)+1000;//str[i]=rand()%9000+1000;此处为1~10
	}
}
void output_read_length()
{
	cout<<"每条read长度："<<endl;
	int i;
	for(i=0;i<num_read;i++)
	{
		cout<<length[i]<<' ';
	}
	cout<<endl;
}

void rand_read()//生成第i个长read
{
	int i,j;
	total_num_of_bases=0;
	for(i=0;i<num_read;i++)
	{
		if((position[i]+length[i])>len_genome)
		{
			length[i]=len_genome-position[i];
		}
		for(j=0;j<length[i];j++)
		{
			sample[i][j]=genome[position[i]+j];
		}
		total_num_of_bases+=length[i];
		//rand_read(sample[i],genome,position[i],length[i]);
	}
}

void replace()//添加3%的碱基替换错误
{
	int i=0,j=0;//i表示sample数组的行，j表示sample数组的列，sample[i][j]表示进行替换的碱基
	int num_of_replacement=0;
	while(num_of_replacement<total_num_of_bases*0.01)
	{
		//cout<<"num_of_replace: "<<num_of_replace<<endl;
		i=rand()%num_read;
		j=rand()%length[i];
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

	int insert_i=0,insert_j=0;//i表示sample数组的行，j表示sample数组的列，sample[i][j]表示进行插入的位置
	int num_of_insertion=0;
	while(num_of_insertion<total_num_of_bases*0.01)
	{
		//cout<<"num_of_insertion: "<<num_of_insertion<<endl;
		insert_i=rand()%num_read;
		insert_j=rand()%length[insert_i];
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
		fout4<<"sample["<<insert_i<<"]["<<insert_j<<"] insert"<<ch<<endl;
		int j=length[insert_i];
		for(;j>insert_j;j--)//将选定位置后面的碱基向后移动
		{
			sample[insert_i][j]=sample[insert_i][j-1];
		}
		sample[insert_i][j]=ch;//插入一个碱基
		num_of_insertion++;
		length[insert_i]++;
	}
}

void deletion()
{
	ofstream fout5;
	fout5.open("E:\\delete.txt",ios::trunc);

	int deletion_i=0,deletion_j=0;//i表示sample数组的行，j表示sample数组的列，sample[i][j]表示进行删除的位置
	int num_of_deletion=0;
	//cout<<"total_num_of_bases: "<<total_num_of_bases<<endl;
	while(num_of_deletion<total_num_of_bases*0.01)
	{
		//cout<<"num_of_deletion: "<<num_of_deletion<<endl;
		deletion_i=rand()%num_read;
		if(length[deletion_i]!=0)
		{
			deletion_j=rand()%length[deletion_i];
		}
		else 
			continue;
		
		fout5<<"sample["<<deletion_i<<"]["<<deletion_j<<"] delete"<<sample[deletion_i][deletion_j]<<endl;
		int j=deletion_j;
		for(;j<length[deletion_i]-1;j++)//将选定位置后面的碱基向后移动
		{
			sample[deletion_i][j]=sample[deletion_i][j+1];
		}
		sample[deletion_i][j]='0';//插入一个碱基
		num_of_deletion++;
		length[deletion_i]--;
	}
}

void output_sample()
{
	for(int i=0;i<num_read;i++)//输出添加碱基替换错误后的所有长read
	{
		cout<<"第"<<i<<"条长read：";
		for(int j=0;j<length[i];j++)
		{
			cout<<j<<sample[i][j]<<" ";
		}
		cout<<endl;
	}
}

void get_short_k_mer()
{
	num_short_k_mer=0;
	initiate_array_k_mer();
	char ch[short_k+1];            //用来暂时存储当前k-mer
	int flag=0;            //用来标记当前k-mer是否存入array_k_mer[]，1为已存在，0为不存在，需添加
	for(int i=0;i<num_read;i++)
	{
		for(int j=0;j<length[i]-short_k+1;j++)
		{
			for(int k=0;k<short_k;k++)
			{
				ch[k]=sample[i][j+k];
				ch[k+1]='\0';
			}

			//此处应添加判断array_k_mer中不含有ch
			for(int m=0;m<num_short_k_mer;m++)
			{
				if(array_k_mer[m].ch==ch)
				{
					array_k_mer[m].count_ch++;
					flag=1;         //表示已有k-mer无需添加
					break;
				}
				else
				{
					continue;
				}
			}
			if(flag==0)
			{
				array_k_mer[num_short_k_mer].ch=ch;
				array_k_mer[num_short_k_mer].count_ch++;
				num_short_k_mer++;
			}
			flag=0;
		}
	}
	
}

void output_short_k_mer()
{
	cout<<num_short_k_mer<<endl;
	for(int i=0;i<num_short_k_mer;i++)
	{
		cout<<array_k_mer[i].ch<<" "<<array_k_mer[i].count_ch<<endl;
	}
}

void save_as_query()
{
	ofstream fout1;
	fout1.open("E:\\original.fasta",ios::trunc);
	for(int i=0;i<num_read;i++)
	{
		fout1<<">"<<i<<" "<<endl;
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
	
		fout2<<">target"<<" "<<endl;
		for(int j=0;j<len_genome;j++)
		{
			fout2<<genome[j];
		}
		fout2<<endl;
	
}



extern void dataset()
{

	//static char genome[len_genome+1];//基因组字符串
	srand((int)time(0));
	memset(position,0,sizeof(int)*num_read);//position 初始化为0
	initiate_sample();
	initiate_array_k_mer();

	rand_genome(genome,len_genome);//生成随机基因组
	cout<<"rand_genome(genome,len_genome);"<<endl;
	//output_genome();
	rand_position(position, num_read);//产生1000个随机起点
	cout<<"rand_position(position, num_read);//产生1000个随机起点"<<endl;
	//output_position();
	rand_length();//生成随机长度
	cout<<"rand_length(length,num_read);//生成随机长度"<<endl;
	//output_read_length();	
	rand_read();
	cout<<"rand_read();"<<endl;
	//cout<<"生成的初始read："<<endl;
	//output_sample();
	//cout<<"引入碱基替换错误："<<endl;
	//replace();//引入碱基替换错误
	//cout<<"replace();//引入碱基替换错误"<<endl;

	//output_read_length();
	//output_sample();
	//cout<<"引入碱基插入错误："<<endl;
	insert();//引入插入错误
	cout<<"insert();//引入插入错误"<<endl;

	//output_read_length();
	//output_sample();
	//cout<<"引入碱基删除错误："<<endl;
	deletion();//引入删除错误
	cout<<"deletion();//引入删除错误"<<endl;

	//output_read_length();
	//output_sample();
	//get_short_k_mer();
	//cout<<"用数组存k-mer："<<endl;
	//output_short_k_mer();
	//save_as_query();
	//save_as_target();
	//save_as_kmer_times();

	
}

