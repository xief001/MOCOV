#include"main.h"

extern int  read_count ;
extern int  len_read ;         //read length 1000-10000
extern int *length;

extern int **count_of_short_node;                           //array of counting every kmer (in read) times
extern int **count_of_medium_node;                           //array of counting every kmer (in read) times

extern int get_index_num(char *ch);

extern char **corrected_read ;
extern int *corrected_len;//(int *)malloc(read_count*sizeof(int));

char *ch_add1[4]={"0"};
char *ch_add2[16]={"0"};
//char *ch_add3[64]={"0"};

int circle_times=0;//if circle_times>=3 give up this region
int circle_value1=0,circle_value2=0;

extern int **insertion;
extern int **deletion;

extern head_node_index head_of_short_hash_table[1024];
extern head_node_index head_of_long_hash_table[1024];
extern head_node_index head_of_medium_hash_table[1024];

int kmer_value=0;                                                       //kmer appearance times , to decide if it is a valid match

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
	return 0;
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
	return 0;
}

int distin_ins_del(int i ,int j )
{
	char *ch1=(char *)malloc((short_k)*sizeof(char));                                            
	char *ch2=(char *)malloc((short_k)*sizeof(char));
	//char ch1[short_k],ch2[short_k];                                  //ch1 stands for the first kmer that kmer_value=1,ch2 stands for the second kmer that kmer_value=1, if ch1 and ch2 both have high value in hash table, it is insertion error ,else it is deletion error
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
	if(kmer_value1>=valid_value&&kmer_value2>=valid_value)            //insertion error
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
	free(ch1);
	free(ch2);
	return flag;
}

void check_insertion(int i,int j,int count ,char* tmp)
{
	int flag=0;
	char *ch0_1=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_2=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_3=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_4=(char *)malloc((medium_k)*sizeof(char));
	char *ch1=(char *)malloc((medium_k)*sizeof(char));
	char *ch2=(char *)malloc((medium_k)*sizeof(char));
	int m;
	int score_A=0,score_T=0,score_C=0,score_G=0;         //score of the deletion base

	//========================ch0===============================

	for(int n=0;n<medium_k;n++)
	{
		if(n<count-1)
		{
			ch0_1[n]=tmp[n];
			ch0_2[n]=tmp[n];
			ch0_3[n]=tmp[n];
			ch0_4[n]=tmp[n];
		}
		else if(n==count-1)
		{
			ch0_1[n]='A';
			ch0_2[n]='T';
			ch0_3[n]='C';
			ch0_4[n]='G';
		}
		else
		{
			ch0_1[n]=sample[i][j+n];
			ch0_2[n]=sample[i][j+n];
			ch0_3[n]=sample[i][j+n];
			ch0_4[n]=sample[i][j+n];
		}
	}
	//==================get k-mer value===================
	score_A=search_valid_medium_kmer(ch0_1);
	score_T=search_valid_medium_kmer(ch0_2);
	score_C=search_valid_medium_kmer(ch0_3);
	score_G=search_valid_medium_kmer(ch0_4);
	//==================get corrected deleted base===========
	if(score_A>=score_C  &&  score_A>=score_G  &&  score_A>=score_T  && score_A>=valid_value)
	{
		flag=1;
	}
	else if(score_T>=score_A  &&  score_T>=score_G  &&  score_T>=score_C  && score_T>=valid_value)
	{
		flag=2;
	}
	else if(score_C>=score_A  &&  score_C>=score_G  &&  score_C>=score_T  && score_C>=valid_value)
	{
		flag=3;
	}
	else if(score_G>=score_C  &&  score_G>=score_A  &&  score_G>=score_T  && score_G>=valid_value)
	{
		flag=4;
	}
	if(flag!=0)            //insertion error
	{
		for(int n=0;n<count+1;n++)
		{
			if(flag==1)
			{
				corrected_read[i][corrected_len[i]]=ch0_1[n];
			}
			if(flag==2)
			{
				corrected_read[i][corrected_len[i]]=ch0_2[n];
			}
			if(flag==3)
			{
				corrected_read[i][corrected_len[i]]=ch0_3[n];
			}
			if(flag==4)
			{
				corrected_read[i][corrected_len[i]]=ch0_4[n];
			}
			corrected_len[i]++;
		}
	}
	else
	{
		//==================ch1======================
		for( m = 0 ; m< medium_k;m++)
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
		//==================get k-mer value===================
		int kmer_value1;
		kmer_value1=search_valid_medium_kmer(ch1);
		if(kmer_value1>=valid_value)            //insertion error
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
			for( m = 0 ; m< medium_k;m++)
			{
				if(m<count-1)
				{
					ch2[m]=sample[i][j+m];
				}
				else 
				{
					ch2[m]=sample[i][j+m+2];
				}
			}		
			int kmer_value2;
			kmer_value2=search_valid_medium_kmer(ch2);
			if(kmer_value2>=valid_value)            //insertion error
			{
				flag=1;
				for(int n=0;n<count;n++)
				{
					corrected_read[i][corrected_len[i]]=ch2[n];
					corrected_len[i]++;
				}
			}
			else                                                             
			{
				for(int n=0;n<count+1;n++)
				{
					corrected_read[i][corrected_len[i]]=sample[i][j+n];
					corrected_len[i]++;
				}
			}
		}
	}
	free(ch0_1);
	free(ch0_2);
	free(ch0_3);
	free(ch0_4);
	free(ch1);
	free(ch2);
}

void check_deletion(int i,int j,int count ,char* tmp)
{
	char *ch0_1=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_2=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_3=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_4=(char *)malloc((medium_k)*sizeof(char));
	int flag=0;
	int m;
	int score_A=0,score_T=0,score_C=0,score_G=0;         //score of the deletion base

	//========================ch0===============================

	for(int n=0;n<medium_k;n++)
	{
		if(n<count)
		{
			ch0_1[n]=tmp[n];
			ch0_2[n]=tmp[n];
			ch0_3[n]=tmp[n];
			ch0_4[n]=tmp[n];
		}
		else if(n==count)
		{
			ch0_1[n]='A';
			ch0_2[n]='T';
			ch0_3[n]='C';
			ch0_4[n]='G';
		}
		else
		{
			ch0_1[n]=sample[i][j+n];
			ch0_2[n]=sample[i][j+n];
			ch0_3[n]=sample[i][j+n];
			ch0_4[n]=sample[i][j+n];
		}
	}
	//==================get k-mer value===================
	score_A=search_valid_medium_kmer(ch0_1);
	score_T=search_valid_medium_kmer(ch0_2);
	score_C=search_valid_medium_kmer(ch0_3);
	score_G=search_valid_medium_kmer(ch0_4);
	//==================get corrected deleted base===========
	if(score_A>=score_C  &&  score_A>=score_G  &&  score_A>=score_T  && score_A>=valid_value)
	{
		flag=1;
	}
	else if(score_T>=score_A  &&  score_T>=score_G  &&  score_T>=score_C  && score_T>=valid_value)
	{
		flag=2;
	}
	else if(score_C>=score_A  &&  score_C>=score_G  &&  score_C>=score_T  && score_C>=valid_value)
	{
		flag=3;
	}
	else if(score_G>=score_C  &&  score_G>=score_A  &&  score_G>=score_T  && score_G>=valid_value)
	{
		flag=4;
	}
	if(flag!=0)            
	{
		for(int n=0;n<count+1;n++)
		{
			if(flag==1)
			{
				corrected_read[i][corrected_len[i]]=ch0_1[n];
			}
			if(flag==2)
			{
				corrected_read[i][corrected_len[i]]=ch0_2[n];
			}
			if(flag==3)
			{
				corrected_read[i][corrected_len[i]]=ch0_3[n];
			}
			if(flag==4)
			{
				corrected_read[i][corrected_len[i]]=ch0_4[n];
			}
			corrected_len[i]++;
		}
	}
	else
	{
		//================================================================get the corrected base at [i][j+count]==================================================================
		int score_A=0,score_T=0,score_C=0,score_G=0;         //score of the deletion base
		char *ch1_1=(char *)malloc((medium_k)*sizeof(char));
		char *ch1_2=(char *)malloc((medium_k)*sizeof(char));
		char *ch1_3=(char *)malloc((medium_k)*sizeof(char));
		char *ch1_4=(char *)malloc((medium_k)*sizeof(char));
		for(int n=0;n<medium_k;n++)
		{
			if(n<count)
			{
				ch1_1[n]=tmp[n];
				ch1_2[n]=tmp[n];
				ch1_3[n]=tmp[n];
				ch1_4[n]=tmp[n];
			}
			else if(n==count)
			{
				ch1_1[n]='A';
				ch1_2[n]='T';
				ch1_3[n]='C';
				ch1_4[n]='G';
			}
			else
			{
				ch1_1[n]=sample[i][j+n-1];
				ch1_2[n]=sample[i][j+n-1];
				ch1_3[n]=sample[i][j+n-1];
				ch1_4[n]=sample[i][j+n-1];
			}
		}
		//==================get k-mer value===================
		score_A=search_valid_medium_kmer(ch1_1);
		score_T=search_valid_medium_kmer(ch1_2);
		score_C=search_valid_medium_kmer(ch1_3);
		score_G=search_valid_medium_kmer(ch1_4);
		//==================get corrected deleted base===========
		if(score_A>=score_C  &&  score_A>=score_G  &&  score_A>=score_T  && score_A>=valid_value)
		{
			flag=1;
		}
		else if(score_T>=score_A  &&  score_T>=score_G  &&  score_T>=score_C  && score_T>=valid_value)
		{
			flag=2;
		}
		else if(score_C>=score_A  &&  score_C>=score_G  &&  score_C>=score_T  && score_C>=valid_value)
		{
			flag=3;
		}
		else if(score_G>=score_C  &&  score_G>=score_A  &&  score_G>=score_T  && score_G>=valid_value)
		{
			flag=4;
		}
		if(flag!=0)            //insertion error
		{
			for(int n=0;n<count+2;n++)
			{
				if(flag==1)
				{
					corrected_read[i][corrected_len[i]]=ch1_1[n];
				}
				if(flag==2)
				{
					corrected_read[i][corrected_len[i]]=ch1_2[n];
				}
				if(flag==3)
				{
					corrected_read[i][corrected_len[i]]=ch1_3[n];
				}
				if(flag==4)
				{
					corrected_read[i][corrected_len[i]]=ch1_4[n];
				}
				corrected_len[i]++;
			}
		}
		else                                                              //deletion error
		{
			int score=0;         //score of the deletion base
			char *ch2=(char *)malloc((medium_k+1)*sizeof(char));
			for(int m=0;m<16;m++)
			{
				for(int n=0;n<medium_k;n++)
				{
					if(n<count)
					{
						ch2[n]=tmp[n];
					}
					else if(n==count)
					{
						ch2[n]=ch_add2[m][0];
					}
					else if(n==count+1)
					{
						ch2[n]=ch_add2[m][1];
					}
					else
					{
						ch2[n]=sample[i][j+n-1];
					}
				}
				score=search_valid_medium_kmer(ch2);
				if(score>=valid_value)            //insertion error
				{
					flag=m+1;
					break;
				}
			}
			if(flag!=0)            //insertion error
			{
				for(int n=0;n<count+3;n++)
				{
					corrected_read[i][corrected_len[i]]=ch2[n];
					corrected_len[i]++;
				}
			}
			else
			{
				for(int n=0;n<count+1;n++)
				{
					corrected_read[i][corrected_len[i]]=sample[i][j+n];
					corrected_len[i]++;
				}
			}
			free(ch2);
		}
		free(ch1_1);
		free(ch1_2);
		free(ch1_3);
		free(ch1_4);
	}
	flag=0;
	free(ch0_1);
	free(ch0_2);
	free(ch0_3);
	free(ch0_4);
}

int check_next_complex(int i,int j,int count);
int find_left_complex(int i,int j ,int count)
{
	int return_value=0;
	if(count<=0)
	{
		return j;
	}
	char *ch=(char *)malloc((medium_k)*sizeof(char));
	//char ch[medium_k]={0};
	int tmp=valid_value;
	int n=0;
	for(n=0;n<medium_k;n++)
	{
		if(n<medium_k-1)
		{
			ch[n]=corrected_read[i][corrected_len[i]-(short_k)+n];
		}
		else
		{
			ch[n]=sample[i][j];
		}
	}
	tmp=search_valid_medium_kmer(ch);

	if(tmp>=valid_value)
	{
		corrected_read[i][corrected_len[i]]=ch[medium_k-1];
		corrected_len[i]++;
		count--;
		free(ch);
		return_value=find_left_complex(i,j+1,count);
		if(circle_times>3)
		{
			return j;
		}
	}
	else
	{
		free(ch);
		circle_value2=j-short_k+1;
		if(circle_value2==circle_value1)
		{
			circle_times++;
		}
		circle_value1=circle_value2;
		if(circle_times>3)
		{
			return j;
		}
		return_value=check_next_complex(i,j-short_k+1,count+short_k-1);
	}
	return return_value;
}


void init_ch_add()
{
	ch_add1[0]="A";ch_add1[1]="T";ch_add1[2]="C";ch_add1[3]="G";
	
	ch_add2[0]="AA";ch_add2[1]="AT";ch_add2[2]="AC";ch_add2[3]="AG";
	ch_add2[4]="TA";ch_add2[5]="TT";ch_add2[6]="TC";ch_add2[7]="TG";
	ch_add2[8]="CA";ch_add2[9]="CT";ch_add2[10]="CC";ch_add2[11]="CG";
	ch_add2[12]="GA";ch_add2[13]="GT";ch_add2[14]="GC";ch_add2[15]="GG";

}

int complex_error(int i,int j,char* ch)
{
	int flag=0;
	int kmer_value[3][64]={0};
	for(int m=0;m<1;m++) 
	{
		if(flag!=0)
		{
			break;
		}
		char *ch1=(char *)malloc((short_k)*sizeof(char));
		char *ch2=(char *)malloc((medium_k)*sizeof(char));
		char *ch3=(char *)malloc((medium_k+1)*sizeof(char));
		if(m==0)
		{
			for(int x=0;x<short_k-1;x++)
			{
				ch1[x]=ch[x];
			}
			for(int n=0;n<4;n++)
			{
				ch1[short_k-1]=ch_add1[n][0];
				kmer_value[m][n]=search_valid_short_kmer(ch1);
				if(kmer_value[m][n]>=valid_value)
				{
					corrected_read[i][corrected_len[i]]=ch_add1[n][0];
					corrected_len[i]++;
					flag=1;
					break;
				}
			}
		}
		free(ch1);
		free(ch2);
		free(ch3);
	}
	return flag;
}

int check_complex(int i,int j,int count)
{
	int return_value=0;
	int flag=0;
	char *ch0_1=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_2=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_3=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_4=(char *)malloc((medium_k)*sizeof(char));
	int m;
	int score_A=0,score_T=0,score_C=0,score_G=0;         //score of the deletion base

	//========================ch0===============================
	
	for(int n=0;n<medium_k;n++)
	{
		if(n<short_k-1)
		{
			ch0_1[n]=sample[i][j+n];
			ch0_2[n]=sample[i][j+n];
			ch0_3[n]=sample[i][j+n];
			ch0_4[n]=sample[i][j+n];
		}
		else if(n==short_k-1)
		{
			ch0_1[n]='A';
			ch0_2[n]='T';
			ch0_3[n]='C';
			ch0_4[n]='G';
		}
		else
		{
			ch0_1[n]=sample[i][j+n];
			ch0_2[n]=sample[i][j+n];
			ch0_3[n]=sample[i][j+n];
			ch0_4[n]=sample[i][j+n];
		}
	}
	//==================get k-mer value===================
	score_A=search_valid_medium_kmer(ch0_1);
	score_T=search_valid_medium_kmer(ch0_2);
	score_C=search_valid_medium_kmer(ch0_3);
	score_G=search_valid_medium_kmer(ch0_4);
	//==================get corrected replace base===========
	if(score_A>=score_C  &&  score_A>=score_G  &&  score_A>=score_T  && score_A>=valid_value)
	{
		flag=1;
	}
	else if(score_T>=score_A  &&  score_T>=score_G  &&  score_T>=score_C  && score_T>=valid_value)
	{
		flag=2;
	}
	else if(score_C>=score_A  &&  score_C>=score_G  &&  score_C>=score_T  && score_C>=valid_value)
	{
		flag=3;
	}
	else if(score_G>=score_C  &&  score_G>=score_A  &&  score_G>=score_T  && score_G>=valid_value)
	{
		flag=4;
	}
	if(flag!=0)            
	{
		for(int n=0;n<short_k;n++)
		{
			if(flag==1)
			{
				corrected_read[i][corrected_len[i]]=ch0_1[n];
			}
			if(flag==2)
			{
				corrected_read[i][corrected_len[i]]=ch0_2[n];
			}
			if(flag==3)
			{
				corrected_read[i][corrected_len[i]]=ch0_3[n];
			}
			if(flag==4)
			{
				corrected_read[i][corrected_len[i]]=ch0_4[n];
			}
			corrected_len[i]++;
		}
		return_value=find_left_complex(i,j+short_k,count-short_k);
	}
	else
	{
		char *ch1=(char *)malloc((medium_k)*sizeof(char));
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
		//==================insertion1/2===================
		int kmer_value1;
		kmer_value1=search_valid_medium_kmer(ch1);
		if(kmer_value1>=valid_value)    //insertion error
		{
			flag=1;
			for(int n=0;n<short_k;n++)
			{
				corrected_read[i][corrected_len[i]]=ch1[n];
				corrected_len[i]++;
			}
			return_value=find_left_complex(i,j+medium_k,count-(medium_k));
			if(circle_times>3)
			{
				flag=-1;
			}
		}
		else
		{
			int m=0;
			char *ch2=(char *)malloc((medium_k)*sizeof(char));
			for( m = 0 ; m< medium_k;m++)
			{
				if(m<short_k-1)
				{
					ch2[m]=sample[i][j+m];
				}
				else 
				{
					ch2[m]=sample[i][j+m+2];
				}
			}		
			int kmer_value2;
			kmer_value2=search_valid_medium_kmer(ch2);
			if(kmer_value2>=valid_value)            //insertion error
			{
				flag=1;
				for(int n=0;n<medium_k;n++)
				{
					corrected_read[i][corrected_len[i]]=ch2[n];
					corrected_len[i]++;
				}
				return_value=find_left_complex(i,j+medium_k+1,count-(medium_k+1));
				if(circle_times>3)
				{
					flag=-1;
				}
			}
			free(ch2);
		}
		if(flag==0)
		{
			//==================ch2======================
			int score_A=0,score_T=0,score_C=0,score_G=0,flag1=0;         //score of the deletion base
			char *ch2_1=(char *)malloc((medium_k)*sizeof(char));
			char *ch2_2=(char *)malloc((medium_k)*sizeof(char));
			char *ch2_3=(char *)malloc((medium_k)*sizeof(char));
			char *ch2_4=(char *)malloc((medium_k)*sizeof(char));
			//char ch2_1[medium_k]={0},ch2_2[medium_k]={0},ch2_3[medium_k]={0},ch2_4[medium_k]={0};
			for( m = 0 ; m< medium_k;m++)
			{
				if(m<short_k-1)
				{
					//ch2[m]=sample[i][j+m];
					ch2_1[m]=sample[i][j+m];
					ch2_2[m]=sample[i][j+m];
					ch2_3[m]=sample[i][j+m];
					ch2_4[m]=sample[i][j+m];
				}
				else if(m==short_k-1)
				{
					ch2_1[m]='A';
					ch2_2[m]='T';
					ch2_3[m]='C';
					ch2_4[m]='G';
				}
				else
				{
					ch2_1[m]=sample[i][j+m-1];
					ch2_2[m]=sample[i][j+m-1];
					ch2_3[m]=sample[i][j+m-1];
					ch2_4[m]=sample[i][j+m-1];
				}
			}

			//==================get k-mer value2===================
			score_A=search_valid_medium_kmer(ch2_1);
			score_T=search_valid_medium_kmer(ch2_2);
			score_C=search_valid_medium_kmer(ch2_3);
			score_G=search_valid_medium_kmer(ch2_4);

			//==================get corrected deleted base===========
			if(score_A>=score_C  &&  score_A>=score_G  &&  score_A>=score_T  && score_A>=valid_value)
			{
				flag1=1;
			}
			else if(score_T>=score_A  &&  score_T>=score_G  &&  score_T>=score_C  && score_T>=valid_value)
			{
				flag1=2;
			}
			else if(score_C>=score_A  &&  score_C>=score_G  &&  score_C>=score_T  && score_C>=valid_value)
			{
				flag1=3;
			}
			else if(score_G>=score_C  &&  score_G>=score_A  &&  score_G>=score_T  && score_G>=valid_value)
			{
				flag1=4;
			}
			if(flag1!=0)            //insertion error
			{
				flag=2;
				int left_count=0;
				if(flag1==1)
				{
					for(int n=0;n<short_k;n++)
					{
						corrected_read[i][corrected_len[i]]=ch2_1[n];
						corrected_len[i]++;
					}
				}

				else if(flag1==2)
				{
					for(int n=0;n<short_k;n++)
					{
						corrected_read[i][corrected_len[i]]=ch2_2[n];
						corrected_len[i]++;
					}
				}

				else if(flag1==3)
				{
					for(int n=0;n<short_k;n++)
					{
						corrected_read[i][corrected_len[i]]=ch2_3[n];
						corrected_len[i]++;
					}
				}

				else if(flag1==4)
				{
					for(int n=0;n<short_k;n++)
					{
						corrected_read[i][corrected_len[i]]=ch2_4[n];
						corrected_len[i]++;
					}
				}
				return_value=find_left_complex(i,j+short_k-1,count-short_k+1);
				if(circle_times>3)
				{
					flag=-1;
				}
			}
			if(flag==0)
			{
				int score=0;         //score of the deletion base
				char *ch2=(char *)malloc((medium_k+1)*sizeof(char));
				for(int m=0;m<16;m++)
				{
					for(int n=0;n<medium_k;n++)
					{
						if(n<short_k-1)
						{
							ch2[n]=sample[i][j+n];
						}
						else if(n==short_k-1)
						{
							ch2[n]=ch_add2[m][0];
						}
						else if(n==short_k)
						{
							ch2[n++]=ch_add2[m][1];
						}
						else
						{
							ch2[n]=sample[i][j+n-1];
						}
					}
					score=search_valid_medium_kmer(ch2);
					if(score>=valid_value)            //insertion error
					{
						flag=m+1;
						break;
					}
				}
				if(flag!=0)            //insertion error
				{
					for(int n=0;n<medium_k;n++)
					{
						corrected_read[i][corrected_len[i]]=ch2[n];
						corrected_len[i]++;
					}
					return_value=find_left_complex(i,j+short_k-1,count-short_k+1);
					if(circle_times>3)
					{
						flag=-1;
					}
				}
				free(ch2);
			}
			free(ch2_1);
			free(ch2_2);
			free(ch2_3);
			free(ch2_4);
		}
		//complex error
		if(flag==0)
		{
			int left_count;
			char *ch=(char *)malloc((short_k-1)*sizeof(char));
			for(int n=0;n<short_k-1;n++)
			{
				corrected_read[i][corrected_len[i]]=ch1[n];
				corrected_len[i]++;
				ch[n]=ch1[n];
			}
			flag=complex_error(i,j,ch);
			if(flag==0)
			{
				for(int n=0;n<short_k;n++)
				{
					if(n<short_k-1)
					{
						corrected_read[i][corrected_len[i]]=ch[n];
						corrected_len[i]++;
					}
					corrected_read[i][corrected_len[i]]=sample[i][j+short_k];
					corrected_len[i]++;
				}
			}
			else if(flag==1&&count==1)
			{
				for(int n=0;n<short_k-1;n++)
				{
					corrected_read[i][corrected_len[i]]=ch[n];
					corrected_len[i]++;
				}
				free(ch);
				return_value=j+short_k;
			}
			else if(flag==2&&count==2)
			{
				for(int n=0;n<short_k-1;n++)
				{
					corrected_read[i][corrected_len[i]]=ch[n];
					corrected_len[i]++;
				}
				free(ch);
				return_value=j+short_k+1;
			}
			else if(flag==3&&count==3)
			{
				for(int n=0;n<short_k-1;n++)
				{
					corrected_read[i][corrected_len[i]]=ch[n];
					corrected_len[i]++;
				}
				free(ch);
				return_value=j+short_k+2;
			}
			else
			{
				return_value=find_left_complex(i,j+short_k,count-short_k);
				if(circle_times>3)
				{
					flag=-1;
				}
			}
			free(ch);
		}
		free(ch1);
	}
	free(ch0_1);
	free(ch0_2);
	free(ch0_3);
	free(ch0_4);
	if(circle_times>3)
	{
		for(int n=0;n<j+count-return_value+1;n++)
		{
			corrected_read[i][corrected_len[i]]=sample[i][return_value+n];
			corrected_len[i]++;
		}
		circle_times=0;
		return j+count+1;
	}
	return return_value;
}

int check_next_complex(int i,int j,int count)
{
	int return_value=0;
	for(int n=0;n<short_k-1;n++)
	{
		sample[i][j+n]=corrected_read[i][corrected_len[i]-short_k+1+n];
	}
	int flag=0;

	char *ch0_1=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_2=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_3=(char *)malloc((medium_k)*sizeof(char));
	char *ch0_4=(char *)malloc((medium_k)*sizeof(char));
	int m;
	int score_A=0,score_T=0,score_C=0,score_G=0;         //score of the deletion base

	//========================ch0===============================
	
	for(int n=0;n<medium_k;n++)
	{
		if(n<short_k-1)
		{
			ch0_1[n]=sample[i][j+n];
			ch0_2[n]=sample[i][j+n];
			ch0_3[n]=sample[i][j+n];
			ch0_4[n]=sample[i][j+n];
		}
		else if(n==short_k-1)
		{
			ch0_1[n]='A';
			ch0_2[n]='T';
			ch0_3[n]='C';
			ch0_4[n]='G';
		}
		else
		{
			ch0_1[n]=sample[i][j+n];
			ch0_2[n]=sample[i][j+n];
			ch0_3[n]=sample[i][j+n];
			ch0_4[n]=sample[i][j+n];
		}
	}
	//==================get k-mer value===================
	score_A=search_valid_medium_kmer(ch0_1);
	score_T=search_valid_medium_kmer(ch0_2);
	score_C=search_valid_medium_kmer(ch0_3);
	score_G=search_valid_medium_kmer(ch0_4);
	//==================get corrected replace base===========
	if(score_A>=score_C  &&  score_A>=score_G  &&  score_A>=score_T  && score_A>=valid_value)
	{
		flag=1;
	}
	else if(score_T>=score_A  &&  score_T>=score_G  &&  score_T>=score_C  && score_T>=valid_value)
	{
		flag=2;
	}
	else if(score_C>=score_A  &&  score_C>=score_G  &&  score_C>=score_T  && score_C>=valid_value)
	{
		flag=3;
	}
	else if(score_G>=score_C  &&  score_G>=score_A  &&  score_G>=score_T  && score_G>=valid_value)
	{
		flag=4;
	}
	if(flag!=0)            //insertion error
	{
		if(flag==1)
		{
			corrected_read[i][corrected_len[i]]=ch0_1[short_k-1];
			corrected_len[i]++;
		}
		if(flag==2)
		{
			corrected_read[i][corrected_len[i]]=ch0_2[short_k-1];
			corrected_len[i]++;
		}
		if(flag==3)
		{
			corrected_read[i][corrected_len[i]]=ch0_3[short_k-1];
			corrected_len[i]++;
		}
		if(flag==4)
		{
			corrected_read[i][corrected_len[i]]=ch0_4[short_k-1];
			corrected_len[i]++;		
		}		
		return_value=find_left_complex(i,j+short_k,count-(short_k));
		if(circle_times>3)
		{
			flag=-1;
		}
	}
	else
	{
		char *ch1=(char *)malloc((medium_k)*sizeof(char));
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
		//==================get k-mer value1=============================================================================
		int kmer_value1;
		kmer_value1=search_valid_medium_kmer(ch1);
		if(kmer_value1>=valid_value)    //insertion error
		{
			flag=1;
			corrected_read[i][corrected_len[i]]=ch1[short_k-1];
			corrected_len[i]++;
			return_value=find_left_complex(i,j+medium_k,count-(medium_k));
		}
		else
		{
			int m=0;
			char *ch2=(char *)malloc((medium_k)*sizeof(char));
			for( m = 0 ; m< medium_k;m++)
			{
				if(m<short_k-1)
				{
					ch2[m]=sample[i][j+m];
				}
				else 
				{
					ch2[m]=sample[i][j+m+2];
				}
			}		
			int kmer_value2;
			kmer_value2=search_valid_medium_kmer(ch2);
			if(kmer_value2>=valid_value)            //insertion error
			{
				flag=1;
				corrected_read[i][corrected_len[i]]=ch1[short_k-1];
				corrected_len[i]++;
				return_value=find_left_complex(i,j+medium_k+1,count-(medium_k+1));
				if(circle_times>3)
				{
					flag=-1;
				}
			}
			free(ch2);
		}
		if(flag==0)
		{
			//==================ch2========================================================================================================
			int score_A=0,score_T=0,score_C=0,score_G=0,flag1=0;         //score of the deletion base
			char *ch2_1=(char *)malloc((medium_k)*sizeof(char));
			char *ch2_2=(char *)malloc((medium_k)*sizeof(char));
			char *ch2_3=(char *)malloc((medium_k)*sizeof(char));
			char *ch2_4=(char *)malloc((medium_k)*sizeof(char));
			//char ch2_1[medium_k]={0},ch2_2[medium_k]={0},ch2_3[medium_k]={0},ch2_4[medium_k]={0};
			for( m = 0 ; m< medium_k;m++)
			{
				if(m<short_k-1)
				{
					//ch2[m]=sample[i][j+m];
					ch2_1[m]=sample[i][j+m];
					ch2_2[m]=sample[i][j+m];
					ch2_3[m]=sample[i][j+m];
					ch2_4[m]=sample[i][j+m];
				}
				else if(m==short_k-1)
				{
					ch2_1[m]='A';
					ch2_2[m]='T';
					ch2_3[m]='C';
					ch2_4[m]='G';
				}
				else
				{
					ch2_1[m]=sample[i][j+m-1];
					ch2_2[m]=sample[i][j+m-1];
					ch2_3[m]=sample[i][j+m-1];
					ch2_4[m]=sample[i][j+m-1];
				}
			}

			//==================get k-mer value2===================
			score_A=search_valid_medium_kmer(ch2_1);
			score_T=search_valid_medium_kmer(ch2_2);
			score_C=search_valid_medium_kmer(ch2_3);
			score_G=search_valid_medium_kmer(ch2_4);

			//==================get corrected deleted base===========
			if(score_A>=score_C  &&  score_A>=score_G  &&  score_A>=score_T  && score_A>=valid_value)
			{
				flag1=1;
			}
			else if(score_T>=score_A  &&  score_T>=score_G  &&  score_T>=score_C  && score_T>=valid_value)
			{
				flag1=2;
			}
			else if(score_C>=score_A  &&  score_C>=score_G  &&  score_C>=score_T  && score_C>=valid_value)
			{
				flag1=3;
			}
			else if(score_G>=score_C  &&  score_G>=score_A  &&  score_G>=score_T  && score_G>=valid_value)
			{
				flag1=4;
			}
			if(flag1!=0)            //insertion error
			{
				flag=2;
				if(flag1==1)
				{
					corrected_read[i][corrected_len[i]]=ch2_1[short_k-1];
					corrected_len[i]++;
				}
				if(flag1==2)
				{
					corrected_read[i][corrected_len[i]]=ch2_2[short_k-1];
					corrected_len[i]++;
				}
				if(flag1==3)
				{
					corrected_read[i][corrected_len[i]]=ch2_3[short_k-1];
					corrected_len[i]++;
				}
				if(flag1==4)
				{
					corrected_read[i][corrected_len[i]]=ch2_4[short_k-1];
					corrected_len[i]++;
				}		
				return_value=find_left_complex(i,j+short_k-1,count-short_k+1);
				if(circle_times>3)
				{
					flag=-1;
				}
			}
			if(flag==0)
			{
				int score=0;         //score of the deletion base
				char *ch2=(char *)malloc((medium_k+1)*sizeof(char));
				for(int m=0;m<16;m++)
				{
					for(int n=0;n<medium_k;n++)
					{
						if(n<short_k-1)
						{
							ch2[n]=sample[i][j+n];
						}
						else if(n==short_k-1)
						{
							ch2[n]=ch_add2[m][0];
						}
						else if(n==short_k)
						{
							ch2[n++]=ch_add2[m][1];
						}
						else
						{
							ch2[n]=sample[i][j+n-1];
						}
					}
					score=search_valid_medium_kmer(ch2);
					if(score>=valid_value)            //insertion error
					{
						flag=m+1;
						break;
					}
				}
				if(flag!=0)            //insertion error
				{
					corrected_read[i][corrected_len[i]]=ch2[short_k-1];
					corrected_len[i]++;
					corrected_read[i][corrected_len[i]]=ch2[short_k];
					corrected_len[i]++;
					return_value=find_left_complex(i,j+short_k-1,count-short_k+1);
					if(circle_times>3)
					{
						flag=-1;
					}
				}
				free(ch2);
			}
			free(ch2_1);
			free(ch2_2);
			free(ch2_3);
			free(ch2_4);
		}
		//complex error
		if(flag==0)
		{
			corrected_read[i][corrected_len[i]]=sample[i][j+short_k];
			corrected_len[i]++;
			/*
			int left_count;
			char *ch=(char *)malloc((short_k-1)*sizeof(char));
			for(int n=0;n<short_k-1;n++)
			{
			ch[n]=ch1[n];
			}
			flag=complex_error(i,j,ch);

			if(flag==0)
			{

			}
			else if(flag==1&&count==(short_k))
			{
			free(ch);
			return_value=j+short_k;
			}
			else if(flag==2&&count==(short_k+1))
			{
			free(ch);
			return_value=j+short_k+1;
			}
			else if(flag==3&&count==(short_k+2))
			{
			free(ch);
			return_value=j+short_k+2;
			}
			*/

			return_value=find_left_complex(i,j+short_k,count-short_k);
			if(circle_times>3)
			{
				flag=-1;

			}
		}
		free(ch1);
	}
	free(ch0_1);
	free(ch0_2);
	free(ch0_3);
	free(ch0_4);
	if(circle_times>3)
	{
		return j;
	}
	return return_value;
}

void subf_correct(int i)
{
	int len=0;
	int count=0;
	int flag=0;
	char *tmp1=(char *)malloc((short_k*20)*sizeof(char));
	//char tmp1[short_k*20]={0};
	for(int j =0 ;j<length[i];j++)
	{
		if(j<length[i]-short_k)
		{
			if(count_of_short_node[i][j]<valid_value)
			{	
				tmp1[count]=sample[i][j];											//tmp1 used to store "count_of_short_node=1" region
				count++;
			}
			else//count_of_node[i][j]<count_value
			{
				if(count==0)
				{
					corrected_read[i][corrected_len[i]]=sample[i][j];
					corrected_len[i]++;
				}
				else
				{
					if(count>short_k)
					{
						j=check_complex(i,j-count,count)-1;
						memset(tmp1,0,sizeof(char)*short_k*20);
						count=0;
						continue;
					}
					if(count==short_k)                                     //num of 1 = short_k, insertion
					{

						check_insertion(i,j-count,count,tmp1);
						memset(tmp1,0,sizeof(char)*short_k*20);
						count=0;
						continue;
					}
					if(count==short_k-2)
					{
						if(sample[i][j]!=sample[i][j+1])                   //deletion error
						{

							check_deletion(i,j-count,count,tmp1);
							memset(tmp1,0,sizeof(char)*short_k*20);
						}
						else                                               //need further judgement
						{
							check_insertion(i,j-count,count,tmp1);
							memset(tmp1,0,sizeof(char)*short_k*20);
						}
					}
					if(count==short_k-1)                                   //num of 1 = short_k-1, insertion/deletion
					{	
						if(sample[i][j-1]!=sample[i][j])                   //deletion error
						{
							check_deletion(i,j-count,count,tmp1);
							memset(tmp1,0,sizeof(char)*short_k*20);
						}
						else                                               //need further judgement
						{
							int flag=0;
							flag=distin_ins_del( i , j);
							if(flag==1)
							{
								check_insertion(i,j-count,count,tmp1);
								memset(tmp1,0,sizeof(char)*short_k*20);
							}
							if(flag==2)
							{
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
						corrected_read[i][corrected_len[i]]=sample[i][j];
						corrected_len[i]++;
						memset(tmp1,0,sizeof(char)*short_k*20);
					}
					count=0;			
				}		
			}			
		}
		//	=========================================================================需要处理count！=0情况
		else
		{
			corrected_read[i][corrected_len[i]]=sample[i][j];
			corrected_len[i]++;
		}
	}

	count=0;

	len=0;	
}

void correct()
{
	for(int i =0 ;i<read_count;i++)
	{
		if(length[i]>short_k)
		{
			cout<<"correcting No. "<<i<<" read..."<<endl;
			subf_correct(i);
		}
	}
}

extern void correct_errors()
{
	init_ch_add();
	cout<<"correct..."<<endl;
	correct();
	save_as_insertion();
	save_as_deletion();
}