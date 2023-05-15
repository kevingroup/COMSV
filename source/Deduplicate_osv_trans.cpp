/*
  Created by Le Li, May 11, 2023
  tool to deduplicate a list of translocations
*/
#include "Header.h"
#include <float.h>
#include <unistd.h>
#include <numeric>   
#include <algorithm> 

struct OMSV{
    LL chr1;
    LL start1;
    LL end1;
    LL chr2;
    LL start2;
    LL end2;
    int id;
    string type;
    string method;
    string part1;
    string part2;
    string molecules;
    string tp1;
    string tp2;
    string tp3;
    char extra[10000]="";
    LL coverage=1;
    bool operator<(const OMSV &x) const{
        return (chr1 < x.chr1 || (chr1 == x.chr1 && chr2 < x.chr2) || (chr1 == x.chr1 && chr2 == x.chr2 && min(start1,end1) < min(x.start1,x.end1)) || (chr1 == x.chr1 && chr2 == x.chr2 && min(start1,end1) == min(x.start1,x.end1) && min(start2,end2) < min(x.start2,x.end2)));
    }
};
vector<OMSV> allSVs;


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

int relation(LL chr1, LL chr2, LL start1, LL stop1, LL start2, LL stop2, LL shift)
{//1 is the SV bed, 2 is the FS bed
        swapIfRev(start1,stop1);
        swapIfRev(start2,stop2);
        if (chr1 != chr2)return 0;
        int typSame = 0;
        if (start2 > stop1+shift || start1 > stop2+shift) return 0; // no relation
        shift = 0;
        if (abs(start1-start2)<=shift && abs(stop1-stop2)<=shift) return 4+typSame; // complete
        if ((start2 <= start1-shift && stop2 > stop1+shift)||(start2 < start1-shift && stop2 >= stop1+shift)) return 3+typSame; // covered
        if ((start1 <= start2-shift && stop1 > stop2+shift)||(start1 < start2-shift && stop1 >= stop2+shift)) return 2+typSame; // covering
        return 1+typSame; // overlap
}

int call_relation(OMSV sv1, OMSV sv2, LL shift){
	//here we only consider the same chromosome order
	return relation(sv1.chr1,sv2.chr1,sv1.start1,sv1.end1,sv2.start1,sv2.end1,shift)*relation(sv1.chr2,sv2.chr2,sv1.start2,sv1.end2,sv2.start2,sv2.end2,shift);
}


LL max(LL x, LL y){
    return x > y ? x : y;
}

LL min(LL x, LL y){
    return x > y ? y : x;
}

struct triEle{
    string zyg;
    string type;
    int flagIdx;
    int selIdx;
    double coverage;
    int num;
};


void readSVBed(char* inputSV)
{
        allSVs.clear();
        FILE* inputBed = fopen(inputSV,"r");
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
        OMSV tempSV;
        char temp1[100];
        char temp2[100];
        char temp3[1000];
        char tstart[100];
        char tend[100];
        char tpart1[100];
        char tpart2[100];
        char tmed[100];
        char tmol[100000];
        char typ[1000];
        char zyg[1000];
        memset(typ,0,sizeof(typ));
        memset(zyg,0,sizeof(zyg));
        while(fscanf(inputBed,"%lld\t%s\t%lld\t%s\t%s\t%s\t%s\t%s\t%s\t%lf\t%lld\t%s\t%s\t%s\t%s\t%s",&tempSV.chr1, temp1, &tempSV.chr2, temp2, typ, tpart1, tpart2, temp3, temp3, temp3, &tempSV.coverage, temp3, temp3, temp3, tmed, tmol)==16)
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

                tempSV.type = typ;
                tempSV.part1 = tpart1;
                tempSV.part2 = tpart2;
                tempSV.method = tmed;
                tempSV.molecules = tmol;
                memset(typ,0,sizeof(typ));
                memset(tpart1,0,sizeof(tpart1));
                memset(tpart2,0,sizeof(tpart2));
                memset(temp1,0,sizeof(temp1));
                memset(temp2,0,sizeof(temp2));
                memset(tmed,0,sizeof(tmed));
                memset(tmol,0,sizeof(tmol));

                fgets(tempSV.extra,10000,inputBed);
                allSVs.push_back(tempSV);
        }

        sort(allSVs.begin(),allSVs.end());
        printf("There are %lld SVs read!\n",(LL)allSVs.size());
        fclose(inputBed);
}



void correctOverlapVariant(LL shift){
    sort(allSVs.begin(),allSVs.end());


    vector<OMSV> uniqSVs;
    uniqSVs.clear();
    vector<OMSV> grpSVs;
    LL tempRealCnt = (LL)allSVs.size();
    if (tempRealCnt < 2) return;
    vector<bool> ocpy(tempRealCnt,false);
    vector<triEle> triVec;
    for (int i = 0; i < tempRealCnt; i++){
        if (ocpy[i]) continue;
        grpSVs.clear();
        grpSVs.push_back(allSVs[i]);
        ocpy[i] = true;
        OMSV lstSV = allSVs[i];
        for (int j = i+1; j < tempRealCnt; j++){
            if (ocpy[j])continue;
            if (allSVs[j].type == lstSV.type && call_relation(allSVs[j],lstSV,shift)>0){
                        grpSVs.push_back(allSVs[j]);
                        ocpy[j] = true;
                        lstSV = allSVs[j];
            }
        }
        if (grpSVs.size() == 1){
            uniqSVs.push_back(grpSVs[0]);
        }
        else{
            triEle tpEle;
            tpEle.coverage = grpSVs[0].coverage;
            tpEle.flagIdx = 0;
            tpEle.selIdx = 0;
            tpEle.num = tpEle.coverage;
            
            for (int s = 0; s < grpSVs.size(); s++){
                    if (tpEle.coverage < grpSVs[s].coverage){
                        tpEle.coverage = grpSVs[s].coverage;
                        tpEle.flagIdx = s;
                        tpEle.selIdx = s;
                        tpEle.num += grpSVs[s].coverage;
                    }
            }

            grpSVs[tpEle.selIdx].coverage = tpEle.num;
            uniqSVs.push_back(grpSVs[tpEle.selIdx]);
        }
    }
    allSVs = uniqSVs;
    printf("After removing the redundant SVs, there are %lld left\n",(LL)allSVs.size());
}


void printFile(char* outF,int minDep)
{
        FILE* outputF = fopen(outF,"w");
        char tp[1000];
        fprintf(outputF,"#chr1\tregion1\tchr2\tregion2\ttype\tsplit-alignment1\tsplit-alignment2\tcoverage\tmethods\tmolecules\n");
        for (LL i = 0; i < (LL)allSVs.size(); i++)
        {
		if (allSVs[i].coverage<minDep && (allSVs[i].type=="Inter-TL-BP" || allSVs[i].type=="Intra-TL-BP"|| allSVs[i].type=="Intra-TL-BP/Deletion"))continue;
                if (allSVs[i].start1!=allSVs[i].end1 || allSVs[i].start2!=allSVs[i].end2)
	                fprintf(outputF,"%lld\t%lld.%lld\t%lld\t%lld-%lld\t%s\t%s\t%s\t%lld\t%s\t%s\n",allSVs[i].chr1, allSVs[i].start1, allSVs[i].end1, allSVs[i].chr2, allSVs[i].start2, allSVs[i].end2, allSVs[i].type.c_str(),allSVs[i].part1.c_str(),allSVs[i].part2.c_str(),allSVs[i].coverage,allSVs[i].method.c_str(),allSVs[i].molecules.c_str());
		else
	                fprintf(outputF,"%lld\t%lld\t%lld\t%lld\t%s\t%s\t%s\t%lld\t%s\t%s\n",allSVs[i].chr1, allSVs[i].start1, allSVs[i].chr2, allSVs[i].start2, allSVs[i].type.c_str(),allSVs[i].part1.c_str(),allSVs[i].part2.c_str(),allSVs[i].coverage,allSVs[i].method.c_str(),allSVs[i].molecules.c_str());
        }
        fclose(outputF);
}

int main(int argv, char *argc[]){
    char outputFile[1000];
    char inputSV[1000];
    LL shift = 500000;
    LL minDep = 5;
    memset(outputFile,0,sizeof(outputFile));
    memset(inputSV,0,sizeof(inputSV));
    if (argv == 1)
    {
        printf("\nThe valid parameters are described as follows:\n");
        printf("\t-inputSV: \n\t\t The input file of SVs!\n");
        printf("\t-outputFile: \n\t\t Default value: %s. The file name of SVs (.osv).\n",outputFile);
        printf("\t-shift: \n\t\t Default value: %lld. The mode of checking overlap.\n",shift);
        printf("\t-minDepth: \n\t\t Default value: %d. The mode of checking overlap.\n",minDep);
        return -1;
    }
    bool paraWrongFlag = false;
    for (int i = 1; i < argv; i=i+2)
    {
        string temp(argc[i]);
        if (temp.compare("-inputSV")==0)
        {
            memset(inputSV,0,sizeof(inputSV));
            strcpy(inputSV, argc[i+1]);
        }
        else if (temp.compare("-outputFile")==0)
        {
            memset(outputFile,0,sizeof(outputFile));
            strcpy(outputFile,argc[i+1]);
        }
        else if (temp.compare("-shift")==0)
        {
            shift = atol(argc[i+1]);
        }
        else if (temp.compare("-minDepth")==0)
        {
            minDep = atoi(argc[i+1]);
        }
        else
        {
            printf("No such parameter or wrong : %s\n",argc[i]);
            paraWrongFlag = true;
        }
    }
    if (paraWrongFlag)
    {
        printf("\nThe format to arrange parameters is:\n\t./makeRefine_en -param_1 val_1 -param_2 val_2 ... -param_n val_n\n");
        printf("\nPlease check your input parameters!\n\n");
        return -1;
    }
    allSVs.clear();
    readSVBed(inputSV);
    correctOverlapVariant(shift);
    printFile(outputFile,minDep);
    return 0;
}

