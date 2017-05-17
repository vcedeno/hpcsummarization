//
//  main.c
//  summarization
//
//  Created by Vanessa Cedeño on 5/27/16.
//  Copyright © 2016 Vanessa Cedeño. All rights reserved.
//
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int m,n,mm,nn;
//MaxInt
int MaxInt=99999;
int MinInt=-99999;

static int mat[1500][1500];

void insertionSort(int matrix[mm][3])
{
    int i,j;
    int temp,temp2,temp3;
    for(i=2;i<=m;i++)
    {
        temp=matrix[i][1];
        temp2=matrix[i][0];
        temp3=matrix[i][2];
        j=i-1;
        while ((temp > matrix[j][1]) && (j >= 1))
        {
            matrix[j+1][0]=matrix[j][0];
            matrix[j+1][1]=matrix[j][1];
            matrix[j+1][2]=matrix[j][2];
            j=j-1;
        }
        matrix[j+1][0]=temp2;
        matrix[j+1][1]=temp;
        matrix[j+1][2]=temp3;
    }
}
//file into the matrix
void generateMatrix(FILE *fstream, int mat[1500][1500])
{
    char buffer[16384] ;
    //char buffer[4096] ;
    char *record,*line;
    int i=0,j=0;
    
    while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
    {
        record = strtok(line,",");
        while(record != NULL)
        {
            mat[i][j++] = atof(record) ;
            record = strtok(NULL,",");
        }
        ++i ;
        n=j;
        nn=n+1;
        j=0;
    }
    m=i;
    mm=m+1;
    printf("m:%d n:%d\n",m,n);
}
//repetition matrix ne based on the event matrix
void repetitionMatrix(int mat[][1500],int ne[][nn])
{
    int i,j;
    
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=m; j++)
        {
            ne[j][i]=ne[j][i-1]+mat[j-1][i-1];
        }
    }
}
float maxValue(int mat[1500][1500], int c, int b, int e)
{
    int i;
    int max=0;
    for (i = b; i <= e ; i++)
    {
        if (mat[c][i] > max)
        {
            max = mat[c][i];
        }
    }
    
    return max;
}
float localmodelcost(int summ[mm][3],int i, int j,float rr[mm],int ss[mm],float sum)
{
    int a,b,c,d,e,counting,uniquev=1,occurrences,valor;
    float totalValue,localdata,p1,probability,pd,temp,pdi,days,dayso;
    int interval[n];
    float dd[m];
    memset(dd, 0,(m)*sizeof(float));
    for(a=1;a<=m;a++)
    {
        totalValue=MaxInt;
        dayso=MaxInt;
        localdata=0;
        days=0;
        for(b=1;b<=a;b++)
        {
            occurrences=0;
            valor=summ[a-b+1][1];
            counting=0;
            pd=0;
            memset(interval, 0,(n)*sizeof(int));
            for(c=(a-b+1);c<=a;c++)
            {
                occurrences=occurrences+(summ[c][1]);
                e=j;
                for (d=0;d<(i-j);d++)
                {
                    interval[d]=(interval[d]+mat[m-(int)summ[c][0]][e]);
                    e++;
                }

                if(summ[c][1]==valor)
                {
                    counting=counting+1;
                }
                if(summ[c][1]!=valor )
                {
                    p1=0;
                    probability=counting/(float)b;
                    if(probability!=0 && probability != 1)
                    {
                        p1=log2(probability);
                    }
                    pd=pd-probability*p1;
                    valor=summ[c][1];
                    counting=1;
                    uniquev=uniquev+1;
                }
                if(c==a)
                {
                    p1=0;
                    probability=counting/(float)b;
                    if(probability!=0 && probability != 1)
                    {
                        p1=log2(probability);
                    }
                    pd=pd-probability*p1;
                }
                
            }
            
            localdata=b*pd;
            for(c=1;c<i;c++)
            {
                temp=interval[c];
                d=c-1;
                
                while ((temp > interval[d]) && (d >= 0))
                {
                    interval[d+1]=interval[d];
                    d=d-1;
                }
                interval[d+1]=temp;
                
            }
            valor=interval[0];
            counting=0;
            pdi=0;
            for(d=0;d<(i-j);d++)
            {
                if(interval[d]==valor)
                    counting=counting+1;
                
                if(interval[d]!=valor)
                {
                    p1=0;
                    probability=(counting/(float)(i-j));
                    if(probability!=0 && probability != 1)
                    {
                        p1=log2(probability);
                    }
                    pdi=pdi-probability*p1;
                    valor=interval[d];
                    counting=1;
                    
                }
                if(d==(i-j-1))
                {
                    p1=0;
                    probability=(counting/(float)(i-j));
                    
                    if(probability!=0 && probability != 1)
                    {
                        p1=log2(probability);
                    }
                    pdi=pdi-probability*p1;
                    
                }
            }
            
            days=(i-j)*pdi;
            days=days+dd[a-b];
            localdata=localdata+log2(m);
            localdata=localdata+rr[a-b];//+days;
        
            if(localdata<totalValue)
            {
                totalValue=localdata;
                ss[a]=b;
                dayso=days;
            }
            
        }
        dd[a]=dayso;
        rr[a]=totalValue;
        
    }
    totalValue=rr[m]+m*log2(m)+log2(n)+dd[m];
    return(totalValue);
}
//DP segmenting
void globalSegment(int ne[mm][nn],float r[nn],int s[nn]){
    
    int i,j,c,e;
    double tvalue,tcost;
    float sum,q;
    int summ[mm][3];
    memset( summ, 0, (mm)*(3)*sizeof(int));
    
    //optimal values for cuts of size 0..n
    float rr[mm];
    memset(rr, 0,(mm)*sizeof(float));
    //optimal first cut for cuts of length 0..n
    int ss[mm];
    memset(ss, 0,(mm)*sizeof(int));
    
    float interval[n];
    memset(interval, 0,(n)*sizeof(float));
    for(j=1;j<=n;j++)
    {
        printf("%d\n",j);
        q=MaxInt;
        for (i=1;i<=j;i++)
        {
            for (c=1;c<=m;c++)
            {
                summ[c][0]=m-c+1;
                summ[c][1]=(ne[c][j]-ne[c][j-i]);
                
                e=0;
            }
            
            insertionSort(summ);
            tcost=localmodelcost(summ,j,j-i,rr,ss,sum);
            
            tvalue= tcost+r[j-i];
            if (tvalue<q)
            {
                q=tvalue;
                s[j]=i;
            }
        }
        r[j]=q;
    }
    //show the results for the global model cuts

}
void printsolution_xml(int ne[mm][nn])
{
    FILE *xmlFile = fopen(sys.argv[2],"w");
    
    if(xmlFile == NULL)
    {
        printf("\n xml file opening failed ");
    }
    float interval[n];
    memset( interval, 0, (n)*sizeof(float));
    
    FILE *csvFile = fopen("visualization.csv","w");
    if(csvFile == NULL)
    {
        printf("\n csv file opening failed ");
    }
    
    char buff[100];
    int i,j,k,c,d,e, elimit, blimit,elocal,blocal,counting;
    float probability,probability1,sum=0,totalcost=0,localcost,temp,valor,pdi,p1,occurrences,pd,localdatavalue,days;
    //optimal values for cuts of size 0..n
    fprintf(xmlFile,"<summarization>\n");
    float r[nn];
    memset(r, 0,(nn)*sizeof(float));
    //optimal first cut for cuts of legnth 0..n
    int s[nn];
    memset(s, 0,(nn)*sizeof(int));
    float rrr[mm];
    memset(rrr, 0,(mm)*sizeof(float));
    //optimal first cut for cuts of legnth 0..n
    int sss[mm];
    memset(sss, 0,(mm)*sizeof(int));
    
    
    globalSegment(ne,r,s);
    
    fprintf(csvFile,"\"date\",\"bucket\",\"count\"\n");
    elimit=n;
    
    time_t now = time (0);
    struct tm* tm = localtime(&now);
    
    while(elimit>0)
    {
        blimit=elimit-s[elimit]+1;
        fprintf(xmlFile,"<time_segment>\n<start_time>%d</start_time>\n<end_time>%d</end_time>\n", blimit, elimit);
        //generate sums
        int summ[mm][3];
        memset( summ, 0, (mm)*(3)*sizeof(int));
        for (c=1;c<=m;c++)
        {
            summ[c][0]=m-c+1;
            summ[c][1]=(ne[c][elimit]-ne[c][blimit-1]);
        }
        
        insertionSort(summ);
        
       
        localcost=localmodelcost(summ,elimit,blimit-1,rrr,sss,sum);
        localcost=localcost;//+log2(n);
        fprintf(xmlFile,"<cost>%f</cost>\n", localcost);
        //printf("Segment cost: %f\n",localcost);
        totalcost=totalcost+localcost;
        elocal=m;
        while(elocal>0)
        {
            blocal=elocal-sss[elocal]+1;
            probability=0;
            memset( interval, 0, (n)*sizeof(float));
            fprintf(xmlFile,"<local_segment>\n");
            valor=summ[blocal][1];
            counting=0;
            pd=0;
            localdatavalue=0;
            occurrences=0;
            days=0;
            for(k=blocal; k<=elocal; k++)
            {
                probability=probability+summ[k][1];
                int event_id = m-summ[k][0]+1;
                fprintf(xmlFile,"<event_id>%d</event_id>\n", event_id);
                occurrences=occurrences+(summ[k][1]);
                e=blimit-1;
                for (d=0;d<(elimit-blimit+1);d++)
                {
                    interval[d]=(interval[d]+mat[m-(int)summ[k][0]][e]);
                    e++;
                }
                if(summ[k][1]==valor)
                {
                    counting=counting+1;
                }
                if(summ[k][1]!=valor )
                {
                    p1=0;
                    probability1=counting/(float)(elocal-blocal+1);
                    if(probability1!=0 && probability1 != 1)
                    {
                        p1=log2(probability1);
                    }
                    pd=pd-probability1*p1;
                    valor=summ[k][1];
                    counting=1;
                }
                if(k==elocal)
                {
                    p1=0;
                    probability1=counting/(float)(elocal-blocal+1);
                    if(probability1!=0 && probability1 != 1)
                    {
                        p1=log2(probability1);
                    }
                    pd=pd-probability1*p1;
                }
                
            }
            localdatavalue=(elocal-blocal+1)*pd;
            for(c=1;c<elimit;c++)
            {
                temp=interval[c];
                d=c-1;
                
                while ((temp > interval[d]) && (d >= 0))
                {
                    interval[d+1]=interval[d];
                    d=d-1;
                }
                interval[d+1]=temp;
                
            }
            valor=interval[0];
            counting=0;
            pdi=0;
            for(d=0;d<(elimit-blimit+1);d++)
            {
                if(interval[d]==valor)
                    counting=counting+1;
                
                if(interval[d]!=valor)
                {
                    p1=0;
                    probability=(counting/(float)(elimit-blimit+1));
                    if(probability!=0 && probability != 1)
                    {
                        p1=log2(probability);
                    }
                    pdi=pdi-probability*p1;
                    valor=interval[d];
                    counting=1;
                    
                }
                if(d==(elimit-blimit))
                {
                    p1=0;
                    probability=(counting/(float)(elimit-blimit+1));
                    
                    if(probability!=0 && probability != 1)
                    {
                        p1=log2(probability);
                    }
                    pdi=pdi-probability*p1;
                    
                }
            }
            days=(elimit-blimit+1)*pdi;
            
            localdatavalue=localdatavalue+log2(m)+days;
            fprintf(xmlFile,"<cost>%f</cost>\n",localdatavalue);
            fprintf(xmlFile,"</local_segment>\n");
            probability=probability/(elocal-blocal+1);
            
            for(i=blocal; i<=elocal; i++)
            {
                for(j=elimit; j>=blimit; j--)
                {
                    tm->tm_mday -=1;
                    time_t next = mktime(tm);
                    strftime (buff, sizeof(buff), "%Y-%m-%d", localtime (&next));
                    fprintf(csvFile,"%s,%d,%d\n",buff,summ[i][0]*100,(int)(probability*10));
                }
                
                tm->tm_mday +=elimit-blimit+1;
            }
            
            elocal=blocal-1;
            //printf("\n");
        }
        fprintf(xmlFile,"</time_segment>\n");
        tm->tm_mday -=elimit-blimit+1;
        elimit=blimit-1;
        
    }
    fprintf(xmlFile,"<total_cost>%f</total_cost>\n",totalcost);
    fprintf(xmlFile,"</summarization>\n");
}
int main(int argc, char** argv)
{
    
    FILE *fstream = fopen(sys.argv[1],"r");
    
    if(fstream == NULL)
    {
        printf("\n file opening failed ");
        return -1 ;
    }
    generateMatrix(fstream,mat);
    int ne[mm][nn];
    memset(ne, 0, (mm)*(nn)*sizeof(int));
    
    repetitionMatrix(mat,ne);
    printsolution_xml(ne);
    return 0;
}
