/*
  Created by Le Li, May 11, 2023
  Complement to inversion caller of COMSV pipeline
*/
#include "Header.h"

FILE* inputBed;
FILE* inputFragileSite;
FILE* outputF;


struct stdSV
{
	char chr[10];
	LL start;
	LL end;
	char type[20];
	double size = 0;
	char zyg[20];
};

void readSVBed(char* inputSV)
{
        int numberOfDL = 0;
        char temp[10000];
        while(fgetc(inputBed)=='#')
        {
                numberOfDL++;
                fgets(temp,10000,inputBed);
        }
        fclose(inputBed);
        if ((inputBed = fopen(inputSV,"r"))==NULL) perror("Wrong in reading standard SVs!");
        for (int i = 0; i < numberOfDL; i++)
        {
                fgets(temp,10000,inputBed);
        }
        stdSV tempSV,preSV;
	char temp1[100];
	char temp2[100];
	char temp3[1000];
        int cnt =0;
        while(fscanf(inputBed,"%s\t%lld\t%lld\t%s\t%s\t%s\t%s\tsizeChange=%lf",tempSV.chr, &tempSV.start, &tempSV.end, tempSV.type, temp1, temp2, temp3, &tempSV.size)==8)
        {
		if (cnt !=0){
                    if (strcmp(tempSV.chr,preSV.chr)==0&&preSV.end+20000>=tempSV.start&&tempSV.type[2]+preSV.type[2]=='s'+'l'&&abs(abs(preSV.size)-abs(tempSV.size))<0.05*min(abs(preSV.size),abs(tempSV.size))){
                        fprintf(outputF,"%s\t%lld\t%lld\tInversion\t-\t%s:%lld-%lld:%lf\t%s:%lld-%lld:%lf\tsizeChange=%lf\tNone\tInferred-Inversion\t-\t-\t-\n",tempSV.chr,preSV.start,tempSV.end,preSV.chr,preSV.start,preSV.end,preSV.size,tempSV.chr,tempSV.start,tempSV.end,tempSV.size,1.0*tempSV.end-preSV.start);
                    }
                }
                strcpy(preSV.chr,tempSV.chr);
                strcpy(preSV.type,tempSV.type);
                strcpy(preSV.zyg,tempSV.zyg);
                preSV.start = tempSV.start;
                preSV.end = tempSV.end;
                preSV.size = tempSV.size;
		fgets(temp,10000,inputBed);
                cnt++;
        }
        cout << "Read " << cnt << " SVs\n";
        fclose(inputBed);
	fclose(outputF);
}



int main(int argc, char* argv[])
{
        char inputSV[500];
        char outFile[500];
        memset(inputSV,0,sizeof(inputSV));
        memset(outFile,0,sizeof(outFile));
                if (argc==1)
                {
                        printf("\nThe valid parameters are described as follows:\n");
                        printf("\t-inputSVFile: \n\t\t The input file name of the detected SVs.\n");
                        printf("\t-outputFile: \n\t\t The output file name.\n");
                        return -1;
                }
        //}
        for (int i = 1; i < argc; i=i+2)
        {
                string temp(argv[i]);
                if (temp.compare("-inputSVFile")==0)
                        strcpy(inputSV,argv[i+1]);
                else if (temp.compare("-outputFile")==0)
                        strcpy(outFile,argv[i+1]);
                else
                {
                        printf("No such parameter or wrong : %s\n",argv[i]);
                        printf("\nThe format to arrange parameters is:\n\t./convertFromCaoToStandard -param_1 val_1 -param_2 val_2 ... -param_n val_n\n");
                        return -1;
                }
        }
        int ccnt = 0;
        int dotPos = -1;
        while(outFile[ccnt]!='\0')
        {
                if (outFile[ccnt]=='.') dotPos = ccnt;
                ccnt++;
        }
	char outFiles[500];
        memset(outFiles,0,sizeof(outFiles));
        if (dotPos == -1) dotPos = ccnt;
        for (int i = 0; i < dotPos; i++)
        {
                outFiles[i] = outFile[i];
        }
        strcat(outFiles,".bed");
        if (((inputBed = fopen(inputSV,"r"))==NULL)||((outputF = fopen(outFiles,"w+"))==NULL))
                perror("Failed to open files!");
	fprintf(outputF,"#");
	for (int par = 0; par < argc; par++)
		fprintf(outputF,"%s  ",argv[par]);
	fprintf(outputF,"\n");
	readSVBed(inputSV);
        return 0;
}

