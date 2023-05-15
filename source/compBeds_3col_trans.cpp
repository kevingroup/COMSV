#include "Header.h"

FILE* inputBed;
FILE* inputFragileSite;
FILE* outputF;


struct stdSV
{
	char chr1[10];
	char chr2[10];
	LL start1;
	LL end1;
	LL start2;
	LL end2;
	char type[20]="-";
	double size = 0;
	char zyg[20]="-";
};
vector<stdSV> SVbed;
vector<stdSV> fragileSite;

struct outRow
{
	char chr1[10];
	char chr2[10];
	LL startSV1;
	LL endSV1;
	LL startSV2;
	LL endSV2;
	char typeSV[20];
	double sizeSV = 0;
	char zygSV[20] = "-";
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

int relation(char* chr1, char* chr2, LL start1, LL stop1, LL start2, LL stop2, int shift)
{//1 is the SV bed, 2 is the FS bed
	swapIfRev(start1,stop1);
	swapIfRev(start2,stop2);
        if (strcmp(chr1,chr2)!=0)return 0;
	int typSame = 0;
	if (start2 > stop1+shift || start1 > stop2+shift) return 0; // no relation
	shift = 0;
	if (abs(start1-start2)<=shift && abs(stop1-stop2)<=shift) return 4+typSame; // complete
	if ((start2 <= start1-shift && stop2 > stop1+shift)||(start2 < start1-shift && stop2 >= stop1+shift)) return 3+typSame; // covered
	if ((start1 <= start2-shift && stop1 > stop2+shift)||(start1 < start2-shift && stop1 >= stop2+shift)) return 2+typSame; // covering
	return 1+typSame; // overlap
}

void readFS(char* inputFS)
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
	while(fscanf(inputFragileSite,"%s\t%lld\t%lld",tempSV.chr1, &tempSV.start1, &tempSV.end1)==3)
	{
		strcpy(tempSV.chr2,tempSV.chr1);
                tempSV.start2 = tempSV.start1;
                tempSV.end2 = tempSV.end1;
		fragileSite.push_back(tempSV);
		fgets(temp,10000,inputFragileSite);
	}
	fclose(inputFragileSite);
}

void readSVBed(char* inputSV)
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
        while(fscanf(inputBed,"%s\t%s\t%s\t%s", tempSV.chr1, temp1, tempSV.chr2, temp2)==4)
        {
                char* pos1 = strtok(temp1,".");
                int cnt=0;
                while(pos1){
                    cnt++;
                    if (cnt==1)tempSV.start1 = atol(pos1);
                    else tempSV.end1 = atol(pos1);
                    pos1 = strtok(NULL,".");
                }
                if (cnt==1)tempSV.end1 = tempSV.start1;

                char* pos2 = strtok(temp2,"-");
                cnt = 0;
                while(pos2){
                    cnt++;
                    if (cnt==1)tempSV.start2 = atol(pos2);
                    else tempSV.end2 = atol(pos2);
                    pos2 = strtok(NULL,"-");
                }
                if (cnt==1)tempSV.end2 = tempSV.start2;

                SVbed.push_back(tempSV);
		fgets(temp,10000,inputBed);
        }
        fclose(inputBed);
}

void crossValidate(int shift)
{
	resList.clear();
	for (LL i = 0; i < (LL)SVbed.size(); i++)
	{
		outRow tempOR;
		strcpy(tempOR.chr1,SVbed[i].chr1);
		strcpy(tempOR.chr2,SVbed[i].chr2);
		tempOR.startSV1 = SVbed[i].start1;
		tempOR.endSV1 = SVbed[i].end1;
		tempOR.startSV2 = SVbed[i].start2;
		tempOR.endSV2 = SVbed[i].end2;
		tempOR.sizeSV = SVbed[i].size;
		strcpy(tempOR.zygSV,SVbed[i].zyg);
		strcpy(tempOR.typeSV,SVbed[i].type);
		int oldScore = 0;
                tempOR.score = 0;
		for (LL j = 0; j < (LL)fragileSite.size(); j++)
		{
//			if (strcmp(SVbed[i].chr1,fragileSite[j].chr1)!=0 && strcmp(SVbed[i].chr2,fragileSite[j].chr1)!=0)continue;
                        tempOR.score += relation(SVbed[i].chr1,fragileSite[j].chr1,SVbed[i].start1, SVbed[i].end1, fragileSite[j].start1, fragileSite[j].end1, shift);
                        tempOR.score += relation(SVbed[i].chr2,fragileSite[j].chr1,SVbed[i].start2, SVbed[i].end2, fragileSite[j].start1, fragileSite[j].end1, shift);
                        tempOR.score += relation(SVbed[i].chr1,fragileSite[j].chr2,SVbed[i].start1, SVbed[i].end1, fragileSite[j].start2, fragileSite[j].end2, shift);
                        tempOR.score += relation(SVbed[i].chr2,fragileSite[j].chr2,SVbed[i].start2, SVbed[i].end2, fragileSite[j].start2, fragileSite[j].end2, shift);
                        if (tempOR.score>0)break;
		}
		resList.push_back(tempOR);
	}
}

void printFile()
{
	fprintf(outputF,"#chr1\tpoint1\tchr2\tpoint2\tValidation_Result\n");
	for (LL i = 0; i < (LL)resList.size(); i++)
	{
                if (resList[i].startSV1!=resList[i].endSV1 || resList[i].startSV2!=resList[i].endSV2)
			fprintf(outputF,"%s\t%lld.%lld\t%s\t%lld-%lld\t%d\n",resList[i].chr1, resList[i].startSV1, resList[i].endSV1, resList[i].chr2, resList[i].startSV2, resList[i].endSV2,resList[i].score);
                else
			fprintf(outputF,"%s\t%lld\t%s\t%lld\t%d\n",resList[i].chr1, resList[i].startSV1, resList[i].chr2, resList[i].endSV2,resList[i].score);
	}
}



int main(int argc, char* argv[])
{
        char inputSV[500];
	char inputFS[500];
        char outFile[500];
	int shift = 0;
        memset(inputSV,0,sizeof(inputSV));
	memset(inputFS,0,sizeof(inputFS));
        memset(outFile,0,sizeof(outFile));
                if (argc==1)//||temp.compare("--help")==0||temp.compare("--Help")==0||temp.compare("--HELP")==0)
                {
                        printf("\nThe valid parameters are described as follows:\n");
                        printf("\t-inputSVFile: \n\t\t The input file name of the detected SVs.\n");
			printf("\t-inputFSFile: \n\t\t The input file name of the fragile site list.\n");
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
	readSVBed(inputSV);
	printf("Done 1.\n");
//        cout << SVbed.size() << endl;
	readFS(inputFS);
	printf("Done 2.\n");
//        cout << fragileSite.size() << endl;
	crossValidate(shift);
	printf("Done 3.\n");
//        cout << resList.size() << endl;
	printFile();
        return 0;
}
