//
//  main.c
//  eventSummarization
//
//  Created by Vanessa Cedeño on 10/5/16.
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

void insertionSort(int matrix[mm][4], int interval)
{
    int i,j,temp,temp2,temp3,temp4;
    float tempf,tempf2;
    for(i=2;i<=m;i++)
    {
        temp=matrix[i][1];
        temp2=matrix[i][0];
        temp3=matrix[i][2];
        temp4=matrix[i][3];
        if(matrix[i][3]==0)
            tempf=0;
        else
            tempf=((matrix[i][1]/(float)interval)*(matrix[i][2]/(float)(matrix[i][3]*interval)));
       
        j=i-1;
        if(matrix[j][3]==0)
            tempf2=0;
        else
            tempf2=((matrix[j][1]/(float)interval)*(matrix[j][2]/(float)(matrix[j][3]*interval)));
        while ((tempf > tempf2) && (j >= 1))
        {
            matrix[j+1][0]=matrix[j][0];
            matrix[j+1][1]=matrix[j][1];
            matrix[j+1][2]=matrix[j][2];
            matrix[j+1][3]=matrix[j][3];
            j=j-1;
        }
        matrix[j+1][0]=temp2;
        matrix[j+1][1]=temp;
        matrix[j+1][2]=temp3;
        matrix[j+1][3]=temp4;
    }
}


//file into the matrix
void generateMatrix(FILE *fstream, int mat[][2500])
{
    char buffer[8192] ;
    char *record,*line;
    int i=0,j=0;
    
    while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
    {
        record = strtok(line,",");
        while(record != NULL)
        {
            mat[i][j++] = atoi(record) ;
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
void repetitionMatrix(int mat[][2500],int ne[mm][nn])
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
//returns optimal localmodelcost
float localmodelcost(int summ[mm][4],int i, int j,float rr[mm],int ss[mm])
{
    int a,b,c;
    double totalValue,localdata,occurrences,occurrences2,p1,p2,probability,expected,maxv,probabilityevents;
    for(a=1;a<=m;a++)
    {
        totalValue=MaxInt;
        localdata=0;
        //printf("a:%d \n", a);
        for(b=1;b<=a;b++)
        {
            occurrences=0;
            occurrences2=0;
            maxv=0;
            for(c=(a-b+1);c<=a;c++)
            {
                occurrences=occurrences+summ[c][1];
                occurrences2=occurrences2+summ[c][2];
                maxv=maxv+summ[c][3];
                
            }
            expected=(i-j)*(b);

            if((maxv/(b))>1)
                probabilityevents=occurrences2/(maxv*(i-j));
            else
                probabilityevents=1;
            
            probability=(occurrences/expected)*probabilityevents;
            p1=0;
            p2=0;
            
            if(probability!=0 && probability != 1)
            {
                p1=log2(probability);
                p2=log2(1-probability);
            }
            localdata=0;
            for(c=(a-b+1);c<=a;c++)
            {
                localdata=localdata-(summ[c][1]*p1)-((i-j-summ[c][1])*p2);
            }
            localdata=localdata+log2(m)+log2(m);
            localdata=localdata+rr[a-b];
            
            if(localdata<totalValue)
            {
                totalValue=localdata;
                ss[a]=b;
            }
            
        }
        rr[a]=totalValue;
        
    }
    totalValue=rr[m]+m*log2(m);
    return(totalValue);
}


//DP global segmenting
void globalSegment(int ne[mm][nn],float r[nn],int s[nn])
{
    
    int i,j,c,d,uvalue,maxv;//,q
    double q,tvalue,tcost;
    int summ[mm][4];

    float rr[mm];
    rr[0]=0;
    memset(rr, 0,(mm)*sizeof(float));
    //optimal first cut for cuts of length 0..n
    int ss[mm];
    memset(ss, 0,(mm)*sizeof(int));
    
    for(j=1;j<=n;j++)
    {
        q=MaxInt;
        for (i=1;i<=j;i++)
        {
            for (c=1;c<=m;c++)
            {
                summ[c][0]=m-c+1;
                summ[c][2]=ne[c][j]-ne[c][j-i];
                uvalue=ne[c][j-i];
                maxv=0;
                summ[c][1]=0;
                for(d=(j-i+1);d<=j;d++)
                {
                    if(ne[c][d]-ne[c][d-1]> maxv)
                    {
                        maxv=ne[c][d]-ne[c][d-1];
                    }
                    if(ne[c][d]!=uvalue)
                    {
                        summ[c][1]=summ[c][1]+1;
                        uvalue=ne[c][d];
                    }
                }
                summ[c][3]=maxv;
            }
            insertionSort(summ,i);
            tcost=localmodelcost(summ,j,j-i,rr,ss);
            //tcost=localmodelcost(summ,j,j-i+1,rr,ss);
            tvalue= tcost+r[j-i]+(1*log2(n));
            if (tvalue<q)
            {
                q=tvalue;
                s[j]=i;
                //printf("j: %d, i:%d, cost: %f\n",j,i, q);
            }
        }
        r[j]=q;
    }

    
}
float localmodelvalue(int frequency1, int frequency2, int interval, int sumevents, int events)
{
    float probability=0,probabilityevents=1,p1,p2,localdata;
    //printf("%f\n",(sumevents/events));
    if((sumevents/events)>1)
    {
        probabilityevents=frequency2/(float)(sumevents*interval);
    }
    probability=(frequency1/(float)(interval*events))*probabilityevents;
    p1=0;
    p2=0;
    if(probability!=0 && probability != 1)
    {
        p1=log2(probability);
        p2=log2(1-probability);
    }
    localdata=-(frequency1*p1)-(((interval*events)-frequency1)*p2)+log2(m)+log2(m);
    return(localdata);
    
}

void PrintSolution(int ne[mm][nn])
{
    FILE *xmlFile = fopen(sys.argv[2],"w");
    
    if(xmlFile == NULL)
    {
        printf("\n xml file opening failed ");
    }
    
    char buff[100];
    int i,j,k,c,d,uvalue, elimit, blimit,elocal,blocal,probability,maxv,occurrences1,occurrences2,sumevents;
    float localcost,localvalue;
    //optimal values for cuts of size 0..n
    fprintf(xmlFile,"<summarization>\n");
    float r[nn];
    memset(r, 0,(nn)*sizeof(float));
    //optimal first cut for cuts of legnth 0..n
    int s[nn];
    memset(s, 0,(nn)*sizeof(int));
    globalSegment(ne,r,s);
    
    float rrr[mm];
    memset(rrr, 0,(mm)*sizeof(float));
    //optimal first cut for cuts of legnth 0..n
    int sss[mm];
    memset(sss, 0,(mm)*sizeof(int));
    
    FILE *csvFile = fopen("visualization.csv","w");
    if(csvFile == NULL)
    {
        printf("\n csv file opening failed ");
    }
    
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
        int summ[mm][4];
        //memset( summ, 0, (mm)*(4)*sizeof(int));
        for (c=1;c<=m;c++)
        {
            summ[c][0]=m-c+1;
            summ[c][2]=ne[c][elimit]-ne[c][blimit-1];
            uvalue=ne[c][blimit-1];
            maxv=0;
            summ[c][1]=0;
            for(d=(blimit);d<=elimit;d++)
            {
                if(ne[c][d]-ne[c][d-1]> maxv)
                {
                    maxv=ne[c][d]-ne[c][d-1];
                }
                if(ne[c][d]!=uvalue)
                {
                    summ[c][1]=summ[c][1]+1;
                    uvalue=ne[c][d];
                }
            }
            summ[c][3]=maxv;
        }

        insertionSort(summ,elimit-blimit+1);

        fprintf(xmlFile,"<cost>%f</cost>\n", localmodelcost(summ,elimit,blimit-1,rrr,sss));
        localcost=0;
        elocal=m;
        while(elocal>0)
        {
            blocal=elocal-sss[elocal]+1;
            probability=0;
            fprintf(xmlFile,"<local_segment>\n");
            occurrences1=0;
            occurrences2=0;
            sumevents=0;
            for(k=blocal; k<=elocal; k++)
            {
                occurrences1=occurrences1+summ[k][1];
                occurrences2=occurrences2+summ[k][2];
                sumevents=sumevents+summ[k][3];
                fprintf(xmlFile,"<event_id>%d</event_id>\n", m-summ[k][0]+1);
            }
            
            probability=occurrences1/(elocal-blocal+1);
            localvalue=localmodelvalue(occurrences1, occurrences2, elimit-blimit+1, sumevents, elocal-blocal+1);
            fprintf(xmlFile,"<cost>%f</cost>\n",localvalue);
            fprintf(xmlFile,"</local_segment>\n");
            localcost=localcost+localvalue;
            for(i=blocal; i<=elocal; i++)
            {
                for(j=elimit; j>=blimit; j--)
                {
                    tm->tm_mday -=1;// elimit-j;
                    time_t next = mktime(tm);
                    strftime (buff, sizeof(buff), "%Y-%m-%d", localtime (&next));
                    fprintf(csvFile,"%s,%d,%d\n",buff,summ[i][0]*100,probability*10);
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
    fprintf(xmlFile,"</summarization>\n");
}

int main()
{
    
    static int mat[2500][2500];
    
    //reading the datafile
    FILE *fstream = fopen(sys.argv[1],"r");
    if(fstream == NULL)
    {
        printf("\n file opening failed ");
        return -1 ;
    }
    generateMatrix(fstream,mat);
    
    int ne[mm][nn];
    memset( ne, 0, (mm)*(nn)*sizeof(int));
    
    repetitionMatrix(mat,ne);
    PrintSolution(ne);
    return 0;
}
