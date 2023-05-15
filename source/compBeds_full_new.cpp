/*
  Created by Le Li, May 11, 2023
  tool to compare two lists of SVs other than translocations
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
vector<stdSV> SVbed;
vector<stdSV> fragileSite;

struct outRow
{
	char chr[10];
	LL startSV;
	LL endSV;
	char typeSV[20];
	double sizeSV = 0;
	char zygSV[20];
	LL startFS = -1;
	LL endFS = -1;
	char typeFS[20] = "-";
	double sizeFS = 0;
	LL score = 0;
	char zygFS[20] = "-";
};
vector<outRow> resList;

void swapIfRev(LL &num1, LL &num2)
{
	LL temp;
	if (num1 > num2)
	{
		temp = num2;
		num2 = num1;
		num1 = temp;
	}
}

int relation(LL start1, LL stop1, LL start2, LL stop2, int shift, char *type1, char *type2, char *zyg1, char *zyg2)
{//1 is the SV bed, 2 is the FS bed
	swapIfRev(start1,stop1);
	swapIfRev(start2,stop2);
	int typSame = 0;
	if ((type1[0]==type2[0]&&type1[1]==type2[1]&&type1[2]==type2[2])||(type1[0]=='B'&&type1[1]=='o')||(type2[0]=='B'&&type2[1]=='o')) typSame += 20;
	if (strcmp(zyg1,zyg2)==0) typSame += 10;
	if (start2 > stop1+shift || start1 > stop2+shift) return 0; // no relation
	shift = 0;
	if (abs(start1-start2)<=shift && abs(stop1-stop2)<=shift) return 4+typSame; // complete
	if ((start2 <= start1-shift && stop2 > stop1+shift)||(start2 < start1-shift && stop2 >= stop1+shift)) return 3+typSame; // covered
	if ((start1 <= start2-shift && stop1 > stop2+shift)||(start1 < start2-shift && stop1 >= stop2+shift)) return 2+typSame; // covering
	return 1+typSame; // overlap
}

void readFS(char* inputFS, int mode=0)
{
	fragileSite.clear();
	int numberOfDL = 0;
	char temp[10000];
	while(fgetc(inputFragileSite)=='#') 
	{
		numberOfDL++;
		fgets(temp,10000,inputFragileSite);
	}
	fclose(inputFragileSite);
	if ((inputFragileSite = fopen(inputFS,"r"))==NULL) perror("Wrong in reading standard SVs!");
	for (int i = 0; i < numberOfDL; i++)
	{
		fgets(temp,10000,inputFragileSite);
	}
	stdSV tempSV;
	char temp1[100];
	char temp2[100];
	char temp3[1000];
        LL templl;
        if (mode == 0){
	while(fscanf(inputFragileSite,"%s\t%lld\t%lld\t%s\t%s\t%s\t%s\tsizeChange=%lf\t%s",tempSV.chr, &tempSV.start, &tempSV.end, tempSV.type, temp1, temp2, temp3, &tempSV.size, tempSV.zyg)==9)
	{
		fragileSite.push_back(tempSV);
		fgets(temp,10000,inputFragileSite);
	}
        }else{
            
	while(fscanf(inputFragileSite,"%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%s\t%s\t%s\tsizeChange=%lf\t%s",tempSV.chr, &tempSV.start, &tempSV.end, &templl,&templl,&templl,&templl,tempSV.type, temp1, temp2, temp3, &tempSV.size, tempSV.zyg)==13)
	{
		fragileSite.push_back(tempSV);
		fgets(temp,10000,inputFragileSite);
	}
        }
	fclose(inputFragileSite);
}

void readSVBed(char* inputSV, int mode=0)
{
        SVbed.clear();
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
        stdSV tempSV;
	char temp1[100];
	char temp2[100];
	char temp3[1000];
        LL templl;
        if (mode ==0){
        while(fscanf(inputBed,"%s\t%lld\t%lld\t%s\t%s\t%s\t%s\tsizeChange=%lf\t%s",tempSV.chr, &tempSV.start, &tempSV.end, tempSV.type, temp1, temp2, temp3, &tempSV.size, tempSV.zyg)==9)
        {
                SVbed.push_back(tempSV);
		fgets(temp,10000,inputBed);
        }
        }else{
        while(fscanf(inputBed,"%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%s\t%s\t%s\tsizeChange=%lf\t%s",tempSV.chr, &tempSV.start, &tempSV.end, &templl,&templl,&templl,&templl,tempSV.type, temp1, temp2, temp3, &tempSV.size, tempSV.zyg)==13)
        {
                SVbed.push_back(tempSV);
		fgets(temp,10000,inputBed);
        }
        }
        fclose(inputBed);
}

void crossValidate(int shift)
{
	resList.clear();
	for (LL i = 0; i < (LL)SVbed.size(); i++)
	{
		outRow tempOR;
		strcpy(tempOR.chr,SVbed[i].chr);
		tempOR.startSV = SVbed[i].start;
		tempOR.endSV = SVbed[i].end;
		tempOR.sizeSV = SVbed[i].size;
		strcpy(tempOR.zygSV,SVbed[i].zyg);
		strcpy(tempOR.typeSV,SVbed[i].type);
		int oldScore = 0;
		for (LL j = 0; j < (LL)fragileSite.size(); j++)
		{
			if (strcmp(SVbed[i].chr,fragileSite[j].chr)!=0)continue;
//			if (fragileSite[j].end + shift < SVbed[i].start) continue;
//			if (fragileSite[j].start > SVbed[i].end + shift) break;
			int score = relation(SVbed[i].start, SVbed[i].end, fragileSite[j].start, fragileSite[j].end, shift, SVbed[i].type, fragileSite[j].type,SVbed[i].zyg,fragileSite[j].zyg);
			if (score>oldScore) {oldScore = score; tempOR.startFS = fragileSite[j].start; tempOR.endFS = fragileSite[j].end; strcpy(tempOR.typeFS,fragileSite[j].type); tempOR.score = score; tempOR.sizeFS = fragileSite[j].size;strcpy(tempOR.zygFS,fragileSite[j].zyg);}
		}
		resList.push_back(tempOR);
	}
}

void printFile()
{
	fprintf(outputF,"#chr\tstart\tend\ttype\tsize\tstartFS\tendFS\ttypeFS\tsizeFS\tValidation_Result\n");
	for (LL i = 0; i < (LL)resList.size(); i++)
	{
		fprintf(outputF,"%s\t%lld\t%lld\t%s\t%lf\t%s\t%lld\t%lld\t%s\t%lf\t%s\t%d\n",resList[i].chr, resList[i].startSV, resList[i].endSV, resList[i].typeSV,resList[i].sizeSV,resList[i].zygSV,resList[i].startFS,resList[i].endFS, resList[i].typeFS,resList[i].sizeFS,resList[i].zygFS,resList[i].score);
	}
}



int main(int argc, char* argv[])
{
        char inputSV[500];
	char inputFS[500];
        char outFile[500];
	int shift = 0;
        int mode_sv=0, mode_fs=0;
        memset(inputSV,0,sizeof(inputSV));
	memset(inputFS,0,sizeof(inputFS));
        memset(outFile,0,sizeof(outFile));
                if (argc==1)//||temp.compare("--help")==0||temp.compare("--Help")==0||temp.compare("--HELP")==0)
                {
                        printf("\nThe valid parameters are described as follows:\n");
                        printf("\t-inputSVFile: \n\t\t The input file name of the detected SVs.\n");
                        printf("\t-modeSV: \n\t\t The mode of input SV file.\n");
			printf("\t-inputFSFile: \n\t\t The input file name of the fragile site list.\n");
                        printf("\t-modeFS: \n\t\t The mode of input FS file.\n");
                        printf("\t-outputFile: \n\t\t The output file name.\n");
                        printf("\t-shift: \n\t\t The tolerated position error of SVs. The default value is: %lf.\n",shift);
                        return -1;
                }
        //}
        for (int i = 1; i < argc; i=i+2)
        {
                string temp(argv[i]);
                if (temp.compare("-inputSVFile")==0)
                        strcpy(inputSV,argv[i+1]);
		else if (temp.compare("-inputFSFile")==0)
			strcpy(inputFS,argv[i+1]);
                else if (temp.compare("-outputFile")==0)
                        strcpy(outFile,argv[i+1]);
                else if (temp.compare("-shift")==0)
                        shift = atoi(argv[i+1]);
                else if (temp.compare("-modeSV")==0)
                        mode_sv = atoi(argv[i+1]);
                else if (temp.compare("-modeFS")==0)
                        mode_fs = atoi(argv[i+1]);
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
        if (((inputBed = fopen(inputSV,"r"))==NULL)||((inputFragileSite = fopen(inputFS,"r"))==NULL)||((outputF = fopen(outFiles,"w+"))==NULL))
                perror("Failed to open files!");
	fprintf(outputF,"#");
	for (int par = 0; par < argc; par++)
		fprintf(outputF,"%s  ",argv[par]);
	fprintf(outputF,"\n");
	readSVBed(inputSV,mode_sv);
	printf("Done 1.\n");
	readFS(inputFS,mode_fs);
	printf("Done 2.\n");
	crossValidate(shift);
	printf("Done 3.\n");
	printFile();
        return 0;
}

