/*
  Created by Le Li, May 11, 2023
  tools to deduplicate a list of inversions and/or duplications
*/
#include "Header.h"
#include <float.h>
#include <unistd.h>
#include <numeric>  
#include <algorithm>

struct OMSV{
    LL chr;
    LL start;
    LL end;
    int id;
    string type;
    string zyg;
    string tp1;
    string tp2;
    string tp3;
    char extra[10000]="";
    double support=1;
    double coverage=1;
    double score=0;
    double difSiz1=0;
    double difSiz2=0;
    bool operator<(const OMSV &x) const{
        return (chr < x.chr || (chr == x.chr && start < x.start) || (chr == x.chr && start == x.start && end < x.end));
    }
};
vector<OMSV> allSVs;

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
        char typ[1000];
        char zyg[1000];
        memset(typ,0,sizeof(typ));
        memset(zyg,0,sizeof(zyg));
        while(fscanf(inputBed,"%lld\t%lld\t%lld\t%s\t%s\t%s\t%s\tsizeChange=%lf\t%s",&tempSV.chr, &tempSV.start, &tempSV.end, typ, temp1, temp2, temp3, &tempSV.difSiz1, zyg)==9)
        {
                tempSV.type = typ;
                tempSV.zyg = zyg;
                tempSV.tp1 = temp1;
                tempSV.tp2 = temp2;
                tempSV.tp3 = temp3;
                memset(typ,0,sizeof(typ));
                memset(zyg,0,sizeof(zyg));
                tempSV.difSiz2 = tempSV.difSiz1;
                fgets(tempSV.extra,10000,inputBed);
                allSVs.push_back(tempSV);
        }

        sort(allSVs.begin(),allSVs.end());
        printf("There are %lld SVs read!\n",(LL)allSVs.size());
        fclose(inputBed);
}



void correctOverlapVariant(int mode=1){
    vector<OMSV> trSVs;
    trSVs.clear();
    for (auto ele:allSVs){
        if (ele.type != "Reference")trSVs.push_back(ele);
    }
    allSVs = trSVs;
    vector<OMSV>().swap(trSVs);
    sort(allSVs.begin(),allSVs.end());
    printf("After removing non-SV cases, there are %lld left\n",(LL)allSVs.size());


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
            if (lstSV.chr != allSVs[j].chr || lstSV.end <= allSVs[j].start)break;
            LL opL = lstSV.end - allSVs[j].start;
            bool ovlp_cond;
            if (mode == 1){ovlp_cond = (lstSV.start == allSVs[j].start)&&(lstSV.end == allSVs[j].end);}
            // 50% intra-overlapping
            else if (mode == 2) {ovlp_cond = (opL > min(lstSV.end - lstSV.start, allSVs[j].end - allSVs[j].start)*0.5 && allSVs[i].end >= allSVs[j].start);}
            else if (mode == 3) {ovlp_cond = (opL > min(lstSV.end - lstSV.start, allSVs[j].end - allSVs[j].start)*0.5 && lstSV.type == allSVs[j].type && allSVs[i].end >= allSVs[j].start);}

            if (ovlp_cond){
                        grpSVs.push_back(allSVs[j]);
                        ocpy[j] = true;
                        lstSV = allSVs[j];
            }
        }
        if (grpSVs.size() == 1){
            uniqSVs.push_back(grpSVs[0]);
        }
        else{
            triVec.clear();
            triEle tpEle;
            tpEle.coverage = grpSVs[0].coverage;
            tpEle.flagIdx = 0;
            tpEle.selIdx = 0;
            tpEle.num = 1;
            tpEle.zyg = grpSVs[0].zyg;
            tpEle.type = grpSVs[0].type;
            triVec.push_back(tpEle);
            
            bool rep = false;
            for (int s = 1; s < grpSVs.size(); s++){
                rep = false;
                for (int t = 0; t < triVec.size(); t++){
                    if (triVec[t].zyg == grpSVs[s].zyg && triVec[t].type == grpSVs[s].type){
                        rep = true;
                        if (triVec[t].coverage < grpSVs[s].coverage){
                            triVec[t].coverage = grpSVs[s].coverage;
                            triVec[t].flagIdx = s;
                            triVec[t].selIdx = s;
                        }
                        triVec[t].num++;
                        break;
                    }
                }
                if (!rep){
                    tpEle.zyg = grpSVs[s].zyg;
                    tpEle.type = grpSVs[s].type;
                    tpEle.flagIdx = s;
                    tpEle.selIdx = s;
                    tpEle.num = 1;
                    tpEle.coverage = grpSVs[s].coverage;
                    triVec.push_back(tpEle);
                }
            }

            for (int s = 0; s < grpSVs.size(); s++){
                for (int t = 0; t < triVec.size(); t++){
                    if (triVec[t].zyg == grpSVs[s].zyg && triVec[t].type == grpSVs[s].type){
                        int fidx = triVec[t].flagIdx;
                        int sidx = triVec[t].selIdx;
                        if ( min(max(abs(grpSVs[s].difSiz1),abs(grpSVs[s].difSiz2)),max(abs(grpSVs[fidx].difSiz1),abs(grpSVs[fidx].difSiz2))) >= min(max(abs(grpSVs[s].difSiz1),abs(grpSVs[s].difSiz2)),max(abs(grpSVs[fidx].difSiz1),abs(grpSVs[fidx].difSiz2)))*0.8 ){
                            if (grpSVs[s].end - grpSVs[s].start < grpSVs[sidx].end - grpSVs[sidx].start){
                                triVec[t].selIdx = s;
                            }
                        }
                        break;
                    }
                }
            }

            int maxN = 0;
            for (int s = 0; s < triVec.size(); s++){
                maxN = max(maxN,triVec[s].num);
            }
            double maxCov = 0;
            int flIdx = -1;
            for (int s = 0; s < triVec.size(); s++){
                if (triVec[s].num == maxN){
                    if (flIdx == -1 || maxCov < triVec[s].coverage){
                        maxCov = triVec[s].coverage;
                        flIdx = triVec[s].selIdx;
                    }
                }
            }
            uniqSVs.push_back(grpSVs[flIdx]);
        }
    }
    allSVs = uniqSVs;
    printf("After removing the redundant SVs, there are %lld left\n",(LL)allSVs.size());
}


void printFile(char* outF)
{
        FILE* outputF = fopen(outF,"w");
        char tp[1000];
        fprintf(outputF,"#chr\tstart\tend\ttype\tsplit-alignment1\tsplit-alignment2\tDup-BP\tsize\tzygosity\tsupp-rate1\tsupp-rate2\tcov\tsupport\n");
        for (LL i = 0; i < (LL)allSVs.size(); i++)
        {
                fprintf(outputF,"%lld\t%lld\t%lld\t%s\t%s\t%s\t%s\tsizeChange=%lf\t%s\t%s",allSVs[i].chr, allSVs[i].start, allSVs[i].end, allSVs[i].type.c_str(),allSVs[i].tp1.c_str(),allSVs[i].tp2.c_str(),allSVs[i].tp3.c_str(),allSVs[i].difSiz1,allSVs[i].zyg.c_str(),allSVs[i].extra);
        }
        fclose(outputF);
}

int main(int argv, char *argc[]){
    char outputFile[1000];
    char inputSV[1000];
    int md = 1;
    memset(outputFile,0,sizeof(outputFile));
    memset(inputSV,0,sizeof(inputSV));
    if (argv == 1)
    {
        printf("\nThe valid parameters are described as follows:\n");
        printf("\t-inputSV: \n\t\t The input file of SVs!\n");
        printf("\t-outputFile: \n\t\t Default value: %s. The file name of SVs (.osv).\n",outputFile);
        printf("\t-mode: \n\t\t Default value: %d. The mode of checking overlap.\n",md);
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
        else if (temp.compare("-mode")==0)
        {
            md = atoi(argc[i+1]);
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
    correctOverlapVariant(md);
    printFile(outputFile);
    return 0;
}

