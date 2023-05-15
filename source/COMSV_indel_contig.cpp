/*
  Created by Le Li, May 11, 2023
  Indel caller based on contig alignment of COMSV pipeline
*/

#include "Header.h"
#include <float.h>
#include <unistd.h>
#include <numeric>      
#include <algorithm>    

FILE* inputChrConsen; 
FILE* inputOptAlign; 
FILE* outputVariantResultFile; 
char inputAlignmentFileName[1000]; 
char inputRepeatFileName[1000]; 
char outputFileLocation[500]; 
FILE *inputAlignmentFile; 
FILE *outputFileList[500]; 
LL numberOfOpticalMap; 
char tempString[10000]; 
LL tempLongLong; 
double tempDouble; 
bool canOpenFile = true;
double minIndelSize; 
double minSVDiffSize; 
double minClsSizRat;
int minClsSI;
int minClsSD;
int minCluster;
int minVar;
double minSVSuppRat;
double minIndelRatio; 
LL inputFlag; 
double confidenceLimit; 
LL numberOfSupportIndelMolecule; 
LL distancePairCount = 0; 
LL numberOfSupportedSV = 0;
LL curId; 
int chrId;
vector <LL> listOfChromosome;

struct OMSV{
    LL chr;
    LL start;
    LL end;
    int id;
    string type;
    string zyg;
    double support;
    double coverage;
    double score;
    double difSiz1=0;
    double difSiz2=0;
    bool operator<(const OMSV &x) const{
        return (chr < x.chr || (chr == x.chr && start < x.start) || (chr == x.chr && start == x.start && end < x.end));
    }
};
vector<OMSV> allSVs;

struct opticalMapType{
    string mapId;
    LL refStartIndex;
    LL refEndIndex;
    LL refStart;
    LL refEnd;
    LL optStart;
    LL optEnd;
    LL chrId;
    bool orientation;
    double score;
    double confidence;
    string hitEnum;
    double fpr=0;////
    double fnr=0;////
    double alignRate=0;////
    vector<int> position;
    void print(){
        printf("%s %lf %lf %lf %lf %lf %lld %lld %lld %lld %lld %s %lld ", mapId.c_str(), score, confidence, fpr,fnr,alignRate, chrId, optStart, optEnd, refStart, refEnd, hitEnum.c_str(), (LL)position.size());
        for (LL i=0; i<(LL)position.size(); i++)
            printf("%d ", position[i]);
        printf("\n");
    }
    void print(FILE* targetFile){
        sprintf(tempString, "%lld", inputFlag);
        fprintf(targetFile, "%s %s %lf %lf %lf %lf %lf %lld %lld %lld %lld %lld %s %lld ", tempString, mapId.c_str(), score, confidence, fpr, fnr, alignRate, chrId, optStart, optEnd, refStart, refEnd, hitEnum.c_str(), (LL)position.size());
        for (LL i=0; i<(LL)position.size(); i++)
            fprintf(targetFile, "%d ", position[i]);
        fprintf(targetFile, "\n");
    }
};

bool ss(opticalMapType q, opticalMapType w){
    return ((q.mapId < w.mapId) || (q.mapId==w.mapId && q.chrId < w.chrId) || (q.mapId == w.mapId && q.chrId == w.chrId && min(q.optStart,q.optEnd) < min(w.optStart,w.optEnd)));
}

vector<opticalMapType> opticalMap1;

struct bedReg{
    LL chr;
    LL start;
    LL stop;
};
vector<bedReg> highDens;


void readHighDensity(char* inputAlignmentFileName){
    cout << "Start to read high density regions!: " << inputAlignmentFileName << endl;
    if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL)return;
    int numberOfDL = 0;
    char ttt;
    char ttts[100000];
    while(fgetc(inputAlignmentFile) == '#'){
        numberOfDL++;
        fgets(ttts,100000,inputAlignmentFile);
    }
    fclose(inputAlignmentFile);
    if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
    for (LL i=0; i<numberOfDL; i++)
        fgets(ttts, 5000, inputAlignmentFile);
    bedReg tpBed;
    while(fscanf(inputAlignmentFile, "%lld\t%lld\t%lld", &tpBed.chr,&tpBed.start,&tpBed.stop) == 3){
        cout << tpBed.chr << "   " << tpBed.start << "   " << tpBed.stop << endl;
        highDens.push_back(tpBed);
        fgets(ttts,5000,inputAlignmentFile);
    }
    cout << "Finish reading high density regions!\n";
    fclose(inputAlignmentFile);
}


void readSourceFile(){
    if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
    int numberOfDL = 0;
    char ttt;
    char ttts[100000];
    while(fgetc(inputAlignmentFile) == '#')
    {
        numberOfDL++;
        fgets(ttts,100000,inputAlignmentFile);
    }
    LL totSize = 0;
    while(fscanf(inputAlignmentFile,"%s",ttts)==1){
        totSize++;
        fgets(ttts,100000,inputAlignmentFile);
    }
    fclose(inputAlignmentFile);
    if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
    opticalMap1.clear();
    opticalMap1.resize(totSize+10000);
    for (LL i=0; i<numberOfDL; i++)
        fgets(tempString, 10000, inputAlignmentFile);
    numberOfOpticalMap = 0;
    char mmid[10000];
    memset(mmid,0,sizeof(mmid));
    LL rem_cnt = 0;
    while (fscanf(inputAlignmentFile, "%s", mmid) == 1){
        opticalMap1[numberOfOpticalMap].mapId = mmid;
        memset(mmid,0,sizeof(mmid));
        LL tempNumberOfSites, tempLL;
        fscanf(inputAlignmentFile, "%lld", &tempNumberOfSites);
        LL tempIndex = 0;
        LL totLen = 0;
        for (LL i=0; i<tempNumberOfSites; i++){
            if (i != tempNumberOfSites - 1){
                fscanf(inputAlignmentFile, "%lld;", &tempLL);
                if (i != 0) totLen += tempLL;
            }
            else fscanf(inputAlignmentFile, "%lld", &tempLL);
            tempIndex += tempLL;
            opticalMap1[numberOfOpticalMap].position.push_back((int)tempIndex);
        }
        fscanf(inputAlignmentFile, "%s", tempString);
        if (tempString[0] == 'c'){
            if (tempString[3] == 'X')
                opticalMap1[numberOfOpticalMap].chrId = 23;
            else if (tempString[3] == 'Y')
                opticalMap1[numberOfOpticalMap].chrId = 24;
            else if (tempString[3] == 'M')
                opticalMap1[numberOfOpticalMap].chrId = 25;
            else
                sscanf(tempString, "chr%lld", &opticalMap1[numberOfOpticalMap].chrId);
        }
        else{
            if (tempString[0] == 'X')
                opticalMap1[numberOfOpticalMap].chrId = 23;
            else if (tempString[0] == 'Y')
                opticalMap1[numberOfOpticalMap].chrId = 24;
            else if (tempString[0] == 'M')
                opticalMap1[numberOfOpticalMap].chrId = 25;
            else
                sscanf(tempString, "%lld", &opticalMap1[numberOfOpticalMap].chrId);
        }
        fscanf(inputAlignmentFile, "%s", tempString);
        if (tempString[0] == 'r' || tempString[0] == '-') opticalMap1[numberOfOpticalMap].orientation = false; 
        else opticalMap1[numberOfOpticalMap].orientation = true;
        fscanf(inputAlignmentFile, "%lf", &opticalMap1[numberOfOpticalMap].score);
        fscanf(inputAlignmentFile, "%lf", &opticalMap1[numberOfOpticalMap].confidence);
        fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStartIndex, &opticalMap1[numberOfOpticalMap].refEndIndex);
        fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].optStart, &opticalMap1[numberOfOpticalMap].optEnd);
        fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStart, &opticalMap1[numberOfOpticalMap].refEnd);
        char hitE[500000];
        memset(hitE,0,sizeof(hitE));
        fscanf(inputAlignmentFile, "%s", hitE);



        opticalMap1[numberOfOpticalMap].hitEnum = hitE;

        numberOfOpticalMap++;
        if (numberOfOpticalMap==totSize){
            totSize+=1000000;
            opticalMap1.resize(totSize);
        }
    }
    opticalMap1.resize(numberOfOpticalMap);
    sort(opticalMap1.begin(), opticalMap1.end(), ss);
    printf("Number of dense regions: %lld\n", (LL)highDens.size());
    printf("Number of optical map: %lld (removed %lld molecules overlapping bad regions)\n", numberOfOpticalMap,rem_cnt);
    LL old_i = 0, new_i=0;////
    for (LL i=0; i<numberOfOpticalMap; i++){
        if (opticalMap1[i].orientation){
            opticalMap1[i].optStart = opticalMap1[i].position[opticalMap1[i].optStart - 1];
            opticalMap1[i].optEnd = opticalMap1[i].position[opticalMap1[i].optEnd];
        } else {
            opticalMap1[i].optStart = opticalMap1[i].position[opticalMap1[i].optStart];
            opticalMap1[i].optEnd = opticalMap1[i].position[opticalMap1[i].optEnd - 1];
        }
        if (!opticalMap1[i].orientation){
            LL tempMaxDis = opticalMap1[i].position[(int)opticalMap1[i].position.size() - 1];
            for (LL j=0; j<(LL)opticalMap1[i].position.size(); j++)
                opticalMap1[i].position[j] = (int)tempMaxDis - opticalMap1[i].position[j];
            sort(opticalMap1[i].position.begin(), opticalMap1[i].position.end());
            opticalMap1[i].optStart = tempMaxDis - opticalMap1[i].optStart;
            opticalMap1[i].optEnd = tempMaxDis - opticalMap1[i].optEnd;
            opticalMap1[i].orientation = true;
        }
    }
}

void initOutput(int chr, const char *type){
    char buffer[200];
    char nameOfFile[1000];
    strcpy(nameOfFile, outputFileLocation);
    sprintf(buffer, "%lld_%d",inputFlag,chr);
    strcat(nameOfFile, buffer);
    strcat(nameOfFile, ".bmap");
    if ((outputFileList[chr] = fopen(nameOfFile, type)) == NULL)
    {
        perror("error opening file(): nameOfFile");
    }
}

void addSplitedMap(){
    LL tempCC = numberOfOpticalMap;
    sort(opticalMap1.begin(), opticalMap1.end(), ss);
    LL totSize = numberOfOpticalMap+50000;
    opticalMap1.resize(totSize);
    for (LL i=0; i<numberOfOpticalMap-1; i++){
        if (opticalMap1[i].mapId==opticalMap1[i+1].mapId && opticalMap1[i].chrId == opticalMap1[i+1].chrId && opticalMap1[i+1].refStart - opticalMap1[i].refEnd > 0 && opticalMap1[i+1].refStart - opticalMap1[i].refEnd < 100000 && (opticalMap1[i].position[0] == opticalMap1[i+1].position[0]) && (opticalMap1[i].position[1] == opticalMap1[i+1].position[1]) && (opticalMap1[i].position[2] == opticalMap1[i+1].position[2])){
            opticalMap1[tempCC] = opticalMap1[i];
            opticalMap1[tempCC].optStart = opticalMap1[i].optEnd;
            opticalMap1[tempCC].optEnd = opticalMap1[i+1].optStart;
            opticalMap1[tempCC].refStart = opticalMap1[i].refEnd;
            opticalMap1[tempCC].refEnd = opticalMap1[i+1].refStart;
            opticalMap1[tempCC].score = 10.0;
            opticalMap1[tempCC].confidence = 0.1;
            opticalMap1[tempCC].orientation = true;
            opticalMap1[tempCC].fpr = opticalMap1[i].fpr + opticalMap1[i+1].fpr;////
            opticalMap1[tempCC].fnr = opticalMap1[i].fnr + opticalMap1[i+1].fnr;////
            opticalMap1[tempCC].position.clear();
            opticalMap1[tempCC].hitEnum = "FFF";
            if (opticalMap1[tempCC].optStart > opticalMap1[tempCC].optEnd)
                tempCC--;
            tempCC++;
            if (tempCC==totSize){
                totSize+=50000;
                opticalMap1.resize(totSize);
            }
        }
    }
    numberOfOpticalMap = tempCC;
    opticalMap1.resize(numberOfOpticalMap);
}

void outputToBillMapDestinationFile(){
    for (LL i=0; i<numberOfOpticalMap; i++){
        if (outputFileList[opticalMap1[i].chrId] == NULL)
            initOutput(opticalMap1[i].chrId,"w+");
        else
            initOutput(opticalMap1[i].chrId,"a+");
        opticalMap1[i].print(outputFileList[opticalMap1[i].chrId]);
        fclose(outputFileList[opticalMap1[i].chrId]);
    }
    opticalMap1.clear();
    vector<opticalMapType>().swap(opticalMap1);
}


void setDefault() {
    confidenceLimit = 9;
    numberOfSupportIndelMolecule = 1; //coverage of the indel
    minIndelSize = 2000;
    minSVDiffSize = 5000;
    minClsSI = 1; // 0.2 for cancer samples, 0.3 for normal samples
    minClsSD = 1; // 0.2 for cancer samples, 0.3 for normal samples
    minVar = 1; // minimum number of molecules having variation distance to be considered
    minClsSizRat = 0.01; // 0.2 for cancer samples, 0.3 for normal samples
    minSVSuppRat = 0.01;
    minIndelRatio = 0.05;
    inputFlag = 0507;
}




bool feq(double x, double y){
    return fabs(x - y) <= eps*(min(fabs(x),fabs(y)));//should be normalized by the values of x and y
}
bool feq1(double x, double y){
    return fabs(x - y) <= eps;
}


LL max(LL x, LL y){
    return x > y ? x : y;
}

LL min(LL x, LL y){
    return x > y ? y : x;
}

double average(vector <double> a){
    double ans = 0;
    for (LL i=0; i<(LL)a.size(); i++)
        ans += a[i];
    ans /= a.size();
    return ans;
}


struct variantType{
    LL people;
    LL chr;
    LL start;
    LL end;
    LL size;
    LL sizeExt;
    LL support;
    LL oppoSupp;////
    double ratio;
    double ratioExt;
    double likelihood;
    bool isSignal;
    bool isDel;
    bool isHomo;
    bool isSupported;
    void update(LL inputPeople, LL inputChr, LL inputStart, LL inputEnd, LL inputSize, LL inputSizeExt, LL inputSignal, LL inputDel, LL inputHomo, LL inputSupport, double inputRatio, double inputRatioExt, double likely, LL supp, LL IoppoSupp)
    {//version for SV
        oppoSupp = IoppoSupp;////
        people = inputPeople;
        chr = inputChr;
        start = inputStart;
        end = inputEnd;
        size = inputSize;
        sizeExt = inputSizeExt;
        ratio = inputRatio;
        ratioExt = inputRatioExt;
        likelihood = likely;
        support = supp;
        isSignal = (inputSignal == 1);
        isDel = (inputDel == 1);
        isHomo = (inputHomo == 1);
        isSupported = (inputSupport == 1);
    };
    void update(LL inputPeople, LL inputChr, LL inputStart, LL inputEnd, LL inputSize, LL inputSignal, LL inputDel, LL inputHomo, LL inputSupport, double inputRatio, double likely, LL supp,  LL IoppoSupp)
    {//version for signal changes
        oppoSupp = IoppoSupp;////
        people = inputPeople;
        chr = inputChr;
        start = inputStart;
        end = inputEnd;
        size = inputSize;
        sizeExt = 0;
        ratio = inputRatio;
        ratioExt = 0;
        likelihood = likely;
        support = supp;
        isSignal = (inputSignal == 1);
        isDel = (inputDel == 1);
        isHomo = (inputHomo == 1);
        isSupported = (inputSupport == 1);
    };
    bool operator<(const variantType &x) const{
        return (chr < x.chr || (chr == x.chr && start < x.start) || (chr == x.chr && start == x.start && end < x.end) || (chr == x.chr && start == x.start && end == x.end && isSupported));
    }

};

struct chrType{
    LL length;
    LL numberOfSites;

    LL* position = new LL[10];
    LL* distance = new LL[10];
    int* coverage = new int[10];
    int* occurrence = new int[10];
    int* gapCount = new int[10];
    int* gapCoverage = new int[10];
    bool* gapSigni = new bool[10];
    bool* signi = new bool[10];
};
chrType chromosome;

struct optAlignType{
    LL belongs;
    string mapId;
    LL optStart;
    LL optEnd;
    LL refStart;
    LL refEnd;
    LL numberOfSites;
    double score;
    double confidence;
    string hitEnum;

    vector<int> position;
    vector<int> oldPosition;
    double fpr;////
    double fnr;////
    double alignRate;////
    bool operator<(const optAlignType &x) const{
        return (refStart < x.refStart || (refStart == x.refStart && refEnd < x.refEnd));
    }

};

struct distanceType{
    string mapId;
    double uniqId;
    LL start;
    LL end;
    bool mP;
    double distance;
    LL cnt;
    bool operator<(const distanceType &x) const{
        return (start < x.start || (start == x.start && end < x.end));
    }
    void print(){
        printf("chrId:%d, id:%s start:%lld end:%lld oldDis:%lld, molDis:%lf cnt:%lld\n", chrId, mapId.c_str(), chromosome.position[start], chromosome.position[end], chromosome.position[end]-chromosome.position[start], distance, cnt);
    }
    void printToFile(FILE* outputFile){
        fprintf(outputFile, "chrId:%d, id:%s start:%lld end:%lld oldDis:%lld, molDis:%lf cnt:%lld\n", chrId, mapId.c_str(), chromosome.position[start], chromosome.position[end], chromosome.position[end] - chromosome.position[start], distance, cnt);
    }
};

vector<optAlignType> opticalMap;
vector<distanceType> distancePair;
vector<variantType> variant;


void init(){
    variant.clear();
    variant.resize(100000);
}

char outputFolder[1000];

double statInfo(vector<double> vec, int mod);
bool iterCompletePair(vector<distanceType> & moreDP, const vector<LL> &headVec, int curStart, int curEnd);
bool advanceLikelihoodDistanceCalculation2(LL start, LL end, LL &cnt, double &svRat);
bool postProcess(vector<int>& finalIdx, vector<LL>& finalMap, vector<double> finalDat, LL refDist, OMSV & resSV);
int advanceLikelihoodDistanceCalculation3(LL start, LL end, double cnt, vector<int> &idx, vector<LL> &maps, vector<double> &distanceRatio, int maxClsNum, int ft_cov);


void getChromosomeList(char* chrMapFile){
    listOfChromosome.clear();
    char nameOfFile[1000];
    strcpy(nameOfFile, chrMapFile);
    if ((inputChrConsen = fopen(nameOfFile, "r")) == NULL){
        perror("error opening file(): ChromosomeInfo");
    }
    fscanf(inputChrConsen,"%s",tempString);
    while(tempString[0]=='#'){
        fgets(tempString, 10000, inputChrConsen);
        fscanf(inputChrConsen,"%s",tempString);
    }
    LL tempChrId, tempNumberOfSites;
    double tempDouble1, tempDouble2;
    tempChrId = atoll(tempString);
    if (fscanf(inputChrConsen, "%lf %lld %lld %lld %lf %lf %lf %lf", &tempDouble1, &tempNumberOfSites, &tempLongLong, &tempLongLong, &tempDouble2, &tempDouble, &tempDouble, &tempDouble) != 8) perror("Wrong in read chromosome maps!");
    listOfChromosome.push_back(tempChrId);
        fgets(tempString, 10000, inputChrConsen);
    while (fscanf(inputChrConsen, "%lld %lf %lld %lld %lld %lf %lf %lf %lf", &tempChrId, &tempDouble1, &tempNumberOfSites, &tempLongLong, &tempLongLong, &tempDouble2, &tempDouble, &tempDouble, &tempDouble) == 9){
        listOfChromosome.push_back(tempChrId);
        fgets(tempString, 10000, inputChrConsen);
    }
    fclose(inputChrConsen);
    sort(listOfChromosome.begin(), listOfChromosome.end());
    vector<LL>::iterator tempIt;
    tempIt = unique(listOfChromosome.begin(), listOfChromosome.end());
    listOfChromosome.resize(distance(listOfChromosome.begin(), tempIt));
    printf("Number of Chromosomes: %lld\n", (LL)listOfChromosome.size());
}

void readChromosomeInfo(int chr,char* chrMapFile){
    char nameOfFile[1000];
    strcpy(nameOfFile, chrMapFile);
    if ((inputChrConsen = fopen(nameOfFile, "r")) == NULL){
        perror("error opening file(): ChromosomeInfo");
    }
    fscanf(inputChrConsen,"%s",tempString);
    while(tempString[0]=='#'){
        fgets(tempString, 10000, inputChrConsen);
        fscanf(inputChrConsen,"%s",tempString);
    }
    LL cc = 0, tempChrId, tempNumberOfSites;
    double tempDouble1, tempDouble2;
    tempChrId = atoll(tempString);
    bool done= false;
    delete[] chromosome.position;
    delete[] chromosome.distance;
    delete[] chromosome.coverage;
    delete[] chromosome.occurrence;
    delete[] chromosome.gapCount;
    delete[] chromosome.gapCoverage;
    delete[] chromosome.gapSigni;
    delete[] chromosome.signi;
    if (fscanf(inputChrConsen, "%lf %lld %lld %lld %lf %lf %lf %lf", &tempDouble1, &tempNumberOfSites, &tempLongLong, &tempLongLong, &tempDouble2, &tempDouble, &tempDouble, &tempDouble) != 8) perror("Wrong in read chromosome maps!");
    if (tempChrId == chr){
        chromosome.numberOfSites = tempNumberOfSites;
        chromosome.position = new LL[tempNumberOfSites+10];
        chromosome.distance = new LL[tempNumberOfSites+10];
        chromosome.coverage = new int[tempNumberOfSites+10];
        chromosome.occurrence = new int[tempNumberOfSites+10];
        chromosome.gapCount = new int[tempNumberOfSites+10];
        chromosome.gapCoverage = new int[tempNumberOfSites+10];
        chromosome.gapSigni = new bool[tempNumberOfSites+10];
        chromosome.signi = new bool[tempNumberOfSites+10];
        done = true;
        chromosome.position[cc] = (LL) tempDouble2;
        chromosome.length = (LL) tempDouble1;
        chromosome.gapCount[cc] = 0;
        chromosome.coverage[cc] = 0;
        chromosome.gapCoverage[cc] = 0;        chromosome.occurrence[cc] = 0;
        cc++;
    }
    fgets(tempString, 10000, inputChrConsen);
    while (fscanf(inputChrConsen, "%lld %lf %lld %lld %lld %lf %lf %lf %lf", &tempChrId, &tempDouble1, &tempNumberOfSites, &tempLongLong, &tempLongLong, &tempDouble2, &tempDouble, &tempDouble, &tempDouble) == 9){
        if (tempChrId == chr){
            chromosome.numberOfSites = tempNumberOfSites;
            if (!done){
                chromosome.position = new LL[tempNumberOfSites+10];
                chromosome.distance = new LL[tempNumberOfSites+10];
                chromosome.coverage = new int[tempNumberOfSites+10];
                chromosome.occurrence = new int[tempNumberOfSites+10];
                chromosome.gapCount = new int[tempNumberOfSites+10];
                chromosome.gapCoverage = new int[tempNumberOfSites+10];
                chromosome.gapSigni = new bool[tempNumberOfSites+10];
                chromosome.signi = new bool[tempNumberOfSites+10];
                done = true;
            }
            chromosome.length = (LL) tempDouble1;
            chromosome.position[cc] = (LL) tempDouble2;
            chromosome.gapCount[cc] = 0;
            chromosome.coverage[cc] = 0;
            chromosome.gapCoverage[cc] = 0;
            chromosome.occurrence[cc] = 0;
            cc++;
        }
        fgets(tempString, 10000, inputChrConsen);
    }
    printf("Finish assign\n");
    chromosome.numberOfSites++;
    for (LL i=0; i<chromosome.numberOfSites-1; i++){
        chromosome.distance[i] = chromosome.position[i+1] - chromosome.position[i];
    }
    printf("chromosome's Sites: %lld\n", chromosome.numberOfSites);
    fclose(inputChrConsen);
}
void readOpticalAlign(int chr, char* outputFileLocation){
    char buffer[200];
    char nameOfFile[1000];
    strcpy(nameOfFile, outputFileLocation);
    sprintf(buffer, "%lld_%d",inputFlag,chr);
    strcat(nameOfFile, buffer);
    strcat(nameOfFile, ".bmap");
    if ((inputOptAlign = fopen(nameOfFile, "r")) == NULL){
        canOpenFile = false;
        return;
    }
    char temp[10000],temp2[100000];
    LL totSiz = 0;
    while (fscanf(inputOptAlign,"%s",temp)==1){
        totSiz++;
        fgets(temp2,100000,inputOptAlign);
    }
    fclose(inputOptAlign);
    opticalMap.clear();
    opticalMap.resize(totSiz+1);
    inputOptAlign = fopen(nameOfFile, "r");
    LL cc = 0;
    LL cr = 0;
    char hitE[500000];
    memset(hitE,0,sizeof(hitE));
    char mmid[100000];
    memset(mmid,0,sizeof(mmid));
    while (fscanf(inputOptAlign, "%lld %s %lf %lf %lf %lf %lf %lld %lld %lld %lld %lld %s %lld", &opticalMap[cc].belongs, mmid, &opticalMap[cc].score, &opticalMap[cc].confidence, &opticalMap[cc].fpr, &opticalMap[cc].fnr, &opticalMap[cc].alignRate, &tempLongLong, &opticalMap[cc].optStart, &opticalMap[cc].optEnd, &opticalMap[cc].refStart, &opticalMap[cc].refEnd, hitE, &opticalMap[cc].numberOfSites) == 14){
        opticalMap[cc].hitEnum = hitE;
        memset(hitE,0,sizeof(hitE));
        opticalMap[cc].mapId = mmid;
        memset(mmid,0,sizeof(mmid));
        vector<int> tempPosition;
        tempPosition.clear();
        tempPosition.resize(opticalMap[cc].numberOfSites);
        opticalMap[cc].position.clear();
        opticalMap[cc].position.resize(opticalMap[cc].numberOfSites);
        opticalMap[cc].oldPosition.clear();
        opticalMap[cc].oldPosition.resize(opticalMap[cc].numberOfSites);
        for (LL i=0; i<opticalMap[cc].numberOfSites; i++){
            fscanf(inputOptAlign, "%d", &opticalMap[cc].position[i]);
            opticalMap[cc].oldPosition[i] = opticalMap[cc].position[i];
            tempPosition[i] = opticalMap[cc].position[i];
        }
        for (LL i=0; i<opticalMap[cc].numberOfSites; i++){
            if (i == 0)
                opticalMap[cc].position[i] = 0;
            else
                opticalMap[cc].position[i] = tempPosition[i] - tempPosition[i-1];
        }
        bool isIn = false;
        for (int k = 0; k < highDens.size(); k++){
            if (highDens[k].chr == chr && ((opticalMap[cc].refStart - highDens[k].stop)*(opticalMap[cc].refEnd - highDens[k].start))<=0){
                isIn = true;
                break;
            }
        }
        if (isIn)cr++;
        else cc++;
    }
    opticalMap[cc].refStart = -1;
    opticalMap[cc].refEnd = -1;
    //Change in April 13
    sort(opticalMap.begin(),opticalMap.begin()+cc);
    printf("Number of kept optical map(cc) of chr %d: %lld; removed: %lld\n", chr, cc, cr);
    fclose(inputOptAlign);
}



map<int,vector<int>> aggl_cluster(vector<double> dat, int & maxNum, int mode = 1){
// Here is the agglomerative hierarchical clustering for 1D dataset
// mode denotes the linkage methods: 1)average linkage (average distance, default), 2)single linkage (minimum distance), 3)complete linkage (maximum distance)
    // Initialization
    if (dat.size() < 1)perror("The number of samples in the agglomerative clustering is less than 1!");
    map<int,vector<int>> aggl_idx;
    aggl_idx.clear();
    int n = dat.size();
    vector<double> curDat = dat;    
    vector<int> curIdx;
    curIdx.clear();
    vector<int> idx;
    idx.clear();
    vector<vector<double>> distMat;
    distMat.clear();
    for (int i = 0; i < dat.size(); i++){
        vector<double> distVec(dat.size(),0.0);
        distMat.push_back(distVec);
    }
    int r_tag, c_tag;
    double maxL = 10000000;
    double minEle;
    int minI, maxI, cur_s, pre_s;
    map<int,double> pairs;
    map<int,int> cnts;
    int cls_inc = 1;
    minEle = maxL;
    for (int i = 0; i < dat.size(); i++){
        curIdx.push_back(i); // Only here use this command
        idx.push_back(i+1);
        for (int j = i+1; j < dat.size(); j++){
            distMat[i][j] = abs(dat[i] - dat[j]); // Only here use this command
            if (distMat[i][j] < minEle){
                r_tag = i; 
                c_tag = j; 
                minEle = distMat[r_tag][c_tag];
            }
        }
    }
    vector<vector<double>> distMat_bak = distMat;
    
    //Here to update final index, current index, and distance matrix
    minI = min(idx[curIdx[r_tag]],idx[curIdx[c_tag]]);
    maxI = max(idx[curIdx[r_tag]],idx[curIdx[c_tag]]);
    for (int i = idx.size()-1; i >= 0; i--){
        if (idx[i] == maxI)idx[i] = minI;
    }
    pairs.clear();
    cnts.clear();
    for (int i = 0; i < idx.size(); i++){
        int s = idx[i];
        if (i > 0 && s == idx[i-1])continue;
        if (mode == 1) pairs[s] = 0;
        else if (mode == 2) pairs[s] = 10000000;
        else if (mode == 3) pairs[s] = -10000000;
        cnts[s] = 0;
    }
    for (int i = 0; i < idx.size()-1; i++){
        for (int j = i+1; j < idx.size(); j++){
            if ((idx[j] != minI && idx[i] != minI) || (idx[i] == idx[j]))continue;
            int othI;
            if (idx[i] == minI) othI = idx[j];
            else othI = idx[i];
            double d_ij = distMat_bak[i][j];
            if (mode == 1){// average linkage, d(X,Y)/|X|/|Y|
                pairs[othI] += d_ij;
                cnts[othI]++;
            }else if (mode == 2){// minimum linkage, min(d(X,Y))
                pairs[othI] = min(pairs[othI], d_ij);
            }else if (mode == 3){// maximum linkage, max(d(X,Y))
                pairs[othI] = max(pairs[othI], d_ij);
            }
        }
    }
    if (mode == 1){
        cur_s = -1; 
        pre_s = -1;
        for (int s = 0; s < idx.size(); s++){
            if (idx[s] == minI)continue;
            cur_s = s;
            if (pre_s == -1 || idx[cur_s] != idx[pre_s]){
                pairs[idx[s]] /= cnts[idx[s]];
                cnts[idx[s]] = 1;
            }
            pre_s = cur_s;
        }
    }
    for (int i = 0; i < curIdx.size(); i++){
        if (curIdx[i] < curIdx[c_tag])
            distMat[curIdx[i]][curIdx[c_tag]] = -1;
        else
            distMat[curIdx[c_tag]][curIdx[i]] = -1;
        if (i == r_tag || i == c_tag)continue;
        if (curIdx[i] < curIdx[r_tag])
            distMat[curIdx[i]][curIdx[r_tag]] = pairs[idx[curIdx[i]]];
        else 
            distMat[curIdx[r_tag]][curIdx[i]] = pairs[idx[curIdx[i]]];
    }
    curIdx.erase(curIdx.begin()+c_tag);

      maxNum = min(10,((int)(log2((double)idx.size()))));

    int steps = 1;
    while(curIdx.size() > 1){

        minEle = 10000000;
        for (int i = 0; i < curIdx.size()-1; i++){
            for (int j = i+1; j < curIdx.size(); j++){
                if (distMat[curIdx[i]][curIdx[j]] < minEle){
                    r_tag = i; 
                    c_tag = j; 
                    minEle = distMat[curIdx[r_tag]][curIdx[c_tag]];
                }
            }
        }
        //Here to update final index, current index, and distance matrix
        minI = min(idx[curIdx[r_tag]],idx[curIdx[c_tag]]);
        maxI = max(idx[curIdx[r_tag]],idx[curIdx[c_tag]]);
        for (int i = idx.size()-1; i >= 0; i--){
            if (idx[i] == maxI)idx[i] = minI;
        }
        pairs.clear();
        cnts.clear();
        for (int i = 0; i < idx.size(); i++){
            int s = idx[i];
            if (i > 0 && s == idx[i-1])continue;
            if (mode == 1) pairs[s] = 0;
            else if (mode == 2) pairs[s] = 10000000;
            else if (mode == 3) pairs[s] = -10000000;
            cnts[s] = 0;
        }
        for (int i = 0; i < idx.size()-1; i++){
            for (int j = i+1; j < idx.size(); j++){
                if ((idx[j] != minI && idx[i] != minI) || (idx[i] == idx[j]))continue;
                            int othI;
                            if (idx[i] == minI) othI = idx[j];
                            else othI = idx[i];

                double d_ij = distMat_bak[i][j];
                if (mode == 1){// average linkage, d(X,Y)/|X|/|Y|
                    pairs[othI] += d_ij;
                    cnts[othI]++;
                }else if (mode == 2){// minimum linkage, min(d(X,Y))
                    pairs[othI] = min(pairs[othI], d_ij);
                }else if (mode == 3){// maximum linkage, max(d(X,Y))
                    pairs[othI] = max(pairs[othI], d_ij);
                }
            }
        }
        if (mode == 1){
            cur_s = -1;
            pre_s = -1;
            for (int s = 0; s < idx.size(); s++){
                if (idx[s] == minI)continue;
                cur_s = s;
                if (pre_s == -1 || idx[cur_s] != idx[pre_s]){
                    pairs[idx[s]] /= cnts[idx[s]];
                    cnts[idx[s]] = 1;
                }
                pre_s = cur_s;
            }
        }

        for (int i = 0; i < curIdx.size(); i++){
            if (curIdx[i] < curIdx[c_tag])
                            distMat[curIdx[i]][curIdx[c_tag]] = -1;
                    else
                            distMat[curIdx[c_tag]][curIdx[i]] = -1;
            if (i == r_tag || i == c_tag)continue;
            if (curIdx[i] < curIdx[r_tag])
                distMat[curIdx[i]][curIdx[r_tag]] = pairs[idx[curIdx[i]]];
            else 
                distMat[curIdx[r_tag]][curIdx[i]] = pairs[idx[curIdx[i]]];
        }

        curIdx.erase(curIdx.begin()+c_tag);
        
        if (curIdx.size() <= maxNum){ 
            // here to adjust index and store it, adjust index to be from 1 to K
            vector<int> idx_store = idx;
            vector<int> idx_save = idx;
            sort(idx_store.begin(),idx_store.end());
            int cls_num = 1;
            map<int,int> cls_map;
            cls_map[idx_store[0]] = cls_num++;

            for (int i = 1; i < idx_store.size(); i++){
                if (idx_store[i] == idx_store[i-1])continue;
                cls_map[idx_store[i]] = cls_num++;
            }
            for (auto &ele:idx_save){
                ele = cls_map[ele];
            }
            aggl_idx[curIdx.size()] = idx_save;
        }
    }
    return aggl_idx;
}

double eucDist(const vector<double> & v1, const vector<double> & v2){
    if (v1.size()!=v2.size()){
        perror("The dimensions of samples are not consistent!");
        exit(0);
    }
    double dist = 0;
    double dif = 0;
    for (int i = 0; i < v1.size(); i++){
        dif = (v1[0] - v2[0]);
        dist += (dif*dif);
    }
    return dist;
}

map<int,vector<int>> aggl_cluster(vector<vector<double>> dat, int & maxNum, int mode = 1){
// Here is the agglomerative hierarchical clustering for 1D dataset
// mode denotes the linkage methods: 1)average linkage (average distance, default), 2)single linkage (minimum distance), 3)complete linkage (maximum distance)
    // Initialization
    map<int,vector<int>> aggl_idx;
    aggl_idx.clear();
    int n = dat.size();
    vector<int> curIdx;
    curIdx.clear();
    vector<int> idx;
    idx.clear();
    vector<vector<double>> distMat;
    distMat.clear();
    for (int i = 0; i < dat.size(); i++){
        vector<double> distVec(dat.size(),0.0);
        distMat.push_back(distVec);
    }
    int r_tag, c_tag;
    double maxL = 10000000;
    double minEle;
    int minI, maxI, cur_s, pre_s;
    map<int,double> pairs;
    map<int,int> cnts;
    int cls_inc = 1;
    minEle = maxL;
    for (int i = 0; i < dat.size(); i++){
        curIdx.push_back(i); // Only here use this command
        idx.push_back(i+1);
        for (int j = i+1; j < dat.size(); j++){
            distMat[i][j] = eucDist(dat[i],dat[j]);
            if (distMat[i][j] < minEle){
                r_tag = i; 
                c_tag = j; 
                minEle = distMat[r_tag][c_tag];
            }
        }
    }
    vector<vector<double>> distMat_bak = distMat;
    
    //Here to update final index, current index, and distance matrix
    minI = min(idx[curIdx[r_tag]],idx[curIdx[c_tag]]);
    maxI = max(idx[curIdx[r_tag]],idx[curIdx[c_tag]]);
    for (int i = idx.size()-1; i >= 0; i--){
        if (idx[i] == maxI)idx[i] = minI;
    }
    pairs.clear();
    cnts.clear();
    for (int i = 0; i < idx.size(); i++){
        int s = idx[i];
        if (i > 0 && s == idx[i-1])continue;
        if (mode == 1) pairs[s] = 0;
        else if (mode == 2) pairs[s] = 10000000;
        else if (mode == 3) pairs[s] = -10000000;
        cnts[s] = 0;
    }
    for (int i = 0; i < idx.size()-1; i++){
        for (int j = i+1; j < idx.size(); j++){
            if ((idx[j] != minI && idx[i] != minI) || (idx[i] == idx[j]))continue;
            int othI;
            if (idx[i] == minI) othI = idx[j];
            else othI = idx[i];
            double d_ij = distMat_bak[i][j];
            if (mode == 1){// average linkage, d(X,Y)/|X|/|Y|
                pairs[othI] += d_ij;
                cnts[othI]++;
            }else if (mode == 2){// minimum linkage, min(d(X,Y))
                pairs[othI] = min(pairs[othI], d_ij);
            }else if (mode == 3){// maximum linkage, max(d(X,Y))
                pairs[othI] = max(pairs[othI], d_ij);
            }
        }
    }
    if (mode == 1){
        cur_s = -1; 
        pre_s = -1;
        for (int s = 0; s < idx.size(); s++){
            if (idx[s] == minI)continue;
            cur_s = s;
            if (pre_s == -1 || idx[cur_s] != idx[pre_s]){
                pairs[idx[s]] /= cnts[idx[s]];
                cnts[idx[s]] = 1;
            }
            pre_s = cur_s;
        }
    }
    for (int i = 0; i < curIdx.size(); i++){
        if (curIdx[i] < curIdx[c_tag])
            distMat[curIdx[i]][curIdx[c_tag]] = -1;
        else
            distMat[curIdx[c_tag]][curIdx[i]] = -1;
        if (i == r_tag || i == c_tag)continue;
        if (curIdx[i] < curIdx[r_tag])
            distMat[curIdx[i]][curIdx[r_tag]] = pairs[idx[curIdx[i]]];
        else 
            distMat[curIdx[r_tag]][curIdx[i]] = pairs[idx[curIdx[i]]];
    }
    curIdx.erase(curIdx.begin()+c_tag);

    int maxNumX = (int)(idx.size()/3);//added feb 25, 2019
    maxNum = -1;
    int steps = 1;
    while(curIdx.size() > 1){

        minEle = 10000000;
        for (int i = 0; i < curIdx.size()-1; i++){
            for (int j = i+1; j < curIdx.size(); j++){
                if (distMat[curIdx[i]][curIdx[j]] < minEle){
                    r_tag = i; 
                    c_tag = j; 
                    minEle = distMat[curIdx[r_tag]][curIdx[c_tag]];
                }
            }
        }
        minI = min(idx[curIdx[r_tag]],idx[curIdx[c_tag]]);
        maxI = max(idx[curIdx[r_tag]],idx[curIdx[c_tag]]);
        for (int i = idx.size()-1; i >= 0; i--){
            if (idx[i] == maxI)idx[i] = minI;
        }
        pairs.clear();
        cnts.clear();
        for (int i = 0; i < idx.size(); i++){
            int s = idx[i];
            if (i > 0 && s == idx[i-1])continue;
            if (mode == 1) pairs[s] = 0;
            else if (mode == 2) pairs[s] = 10000000;
            else if (mode == 3) pairs[s] = -10000000;
            cnts[s] = 0;
        }
        for (int i = 0; i < idx.size()-1; i++){
            for (int j = i+1; j < idx.size(); j++){
                if ((idx[j] != minI && idx[i] != minI) || (idx[i] == idx[j]))continue;
                            int othI;
                            if (idx[i] == minI) othI = idx[j];
                            else othI = idx[i];

                double d_ij = distMat_bak[i][j];
                if (mode == 1){// average linkage, d(X,Y)/|X|/|Y|
                    pairs[othI] += d_ij;
                    cnts[othI]++;
                }else if (mode == 2){// minimum linkage, min(d(X,Y))
                    pairs[othI] = min(pairs[othI], d_ij);
                }else if (mode == 3){// maximum linkage, max(d(X,Y))
                    pairs[othI] = max(pairs[othI], d_ij);
                }
            }
        }
        if (mode == 1){
            cur_s = -1;
            pre_s = -1;
            for (int s = 0; s < idx.size(); s++){
                if (idx[s] == minI)continue;
                cur_s = s;
                if (pre_s == -1 || idx[cur_s] != idx[pre_s]){
                    pairs[idx[s]] /= cnts[idx[s]];
                    cnts[idx[s]] = 1;
                }
                pre_s = cur_s;
            }
        }

        for (int i = 0; i < curIdx.size(); i++){
            if (curIdx[i] < curIdx[c_tag])
                            distMat[curIdx[i]][curIdx[c_tag]] = -1;
                    else
                            distMat[curIdx[c_tag]][curIdx[i]] = -1;
            if (i == r_tag || i == c_tag)continue;
            if (curIdx[i] < curIdx[r_tag])
                distMat[curIdx[i]][curIdx[r_tag]] = pairs[idx[curIdx[i]]];
            else 
                distMat[curIdx[r_tag]][curIdx[i]] = pairs[idx[curIdx[i]]];
        }

        curIdx.erase(curIdx.begin()+c_tag);
        
        vector<int> cekIdx = idx;//add Feb 22, 2019
        sort(cekIdx.begin(),cekIdx.end());
        int curClsSiz = 1;
        bool stp = true;
        for (int xt = 1; xt < cekIdx.size(); xt++){
            if (cekIdx[xt] == cekIdx[xt-1])curClsSiz++;
            else{
                if (curClsSiz*10 < idx.size()) {
                    stp = false;
                    break;
                }
                curClsSiz = 1;
            }
        }// add until here
        if (curClsSiz*10 < idx.size())stp = false;
        if (stp && (maxNum == -1)) maxNum = curIdx.size(); // Add Feb 22, 2019


        if (curIdx.size() <= maxNumX){ // here to adjust index and store it, adjust index to be from 1 to K
            vector<int> idx_store = idx;
            vector<int> idx_save = idx;
            sort(idx_store.begin(),idx_store.end());
            int cls_num = 1;
            map<int,int> cls_map;
            cls_map[idx_store[0]] = cls_num++;

            for (int i = 1; i < idx_store.size(); i++){
                if (idx_store[i] == idx_store[i-1])continue;
                cls_map[idx_store[i]] = cls_num++;
            }
            for (auto &ele:idx_save){
                ele = cls_map[ele];
            }
            aggl_idx[curIdx.size()] = idx_save;
        }
    }
    return aggl_idx;
}


vector<distanceType> completeGivenPairs(LL ost, LL oed){
// Function to complete any given distance pair based on the global distance pair vector distancePair
    vector<distanceType> extraDP;
    extraDP.clear();
    distanceType tempPair;
    if (oed-ost == 1)return extraDP;
    LL i;
    for (LL j = 0; j < distancePair.size(); j++)//this search the possible sub-regions
    {
        if (distancePair[j].start < ost) continue;
        else if (distancePair[j].start > ost) break;
        else
        {
            i = j - 1;
            tempPair = distancePair[j];
            bool rep = true;
            for (LL k = j+1; distancePair[k].start<oed && k<distancePair.size(); k++)
            {
                if (distancePair[k].mapId== distancePair[j].mapId && distancePair[k].uniqId== distancePair[j].uniqId){///Here modified in Dec 3, 2019
                    if (distancePair[k].start == tempPair.end){
                        tempPair.end = distancePair[k].end;
                        tempPair.distance += distancePair[k].distance;
                    }else if (distancePair[k].start > tempPair.end)
                        break;
                }
                if (tempPair.end >= oed || distancePair[k].start > tempPair.end) break;
            }
            if (tempPair.end == oed) //no need
            {
                rep = false;
                for (LL l = max((LL)0,i); l < distancePairCount; l++)
                {
                    if (distancePair[l].start > ost || (distancePair[l].start == ost && distancePair[l].end > oed))break;
                    else if (distancePair[l].mapId==tempPair.mapId && distancePair[l].uniqId==tempPair.uniqId && tempPair.start == distancePair[l].start && tempPair.end == distancePair[l].end) {rep = true; break;}
                }
            }
            if (!rep) extraDP.push_back(tempPair);
        }
    }
    return extraDP;
}


bool iterCompletePair(vector<distanceType> & moreDP, const vector<LL> &headVec, int curStart, int curEnd){
    vector<distanceType> extraDP;
    extraDP.clear();
    bool ext = false;
    int oldLast, newLast;
    oldLast = moreDP.size();
    int cnt = 0;
    /**consider to add all possible internal regions**/
    LL ost = distancePair[headVec[curStart]].start, oed = distancePair[headVec[curStart]].end;
    for (int i = curStart + 1; i <= curEnd; i++){
        ost = min(ost,distancePair[headVec[i]].start);
        oed = max(oed,distancePair[headVec[i]].end);
    }
    int maxRg = 5; // here use 
    for (LL i = ost; i <= oed; i++){
        for (LL j = i + 1; j <= min(oed,i+maxRg); j++){
            for (int s = max((LL)0,curStart-1000); s <= min(curEnd+1000,(LL)headVec.size()-1); s++){
                if ((distancePair[headVec[s]].start == i) && (distancePair[headVec[s]].end == j)){
                                            // if the continuous region already exists in current set, then skip
                    ext = true;
                    break;
                }
            }
            if (!ext){
                extraDP = completeGivenPairs(i,j);
                if (extraDP.size() > 0)cnt++;
                for (auto ele:extraDP)moreDP.push_back(ele);
                extraDP.clear();
            }
            ext = false;
        }
    }

    return true;
}

void completePairs(int ft_cov){
    sort(distancePair.begin(), distancePair.end());
    vector<distanceType> extraDP;
    extraDP.clear();
    extraDP.reserve(distancePairCount*5);
    LL ost = -1, oed = -1;
    LL max_range = 0, max_depth=0;
    for (LL i = 0; i < distancePairCount; i++)
    {
        if (distancePair[i].start!= ost || distancePair[i].end!=oed)//this locate the non-redundant reference regions (NRRR)
        {
            distanceType tempPair;
            ost = distancePair[i].start;
            oed = distancePair[i].end;
            if (oed-ost <= 1)continue;
            if ((oed-ost > 10) || (chromosome.position[oed]-chromosome.position[ost]>100000)){
                max_range++;
                continue;//here remove the large regions
            }
            int cur_cov = 0;
            
            for (LL j = 0; j <= i-1; j++)//this search the possible sub-regions
            {
                if (distancePair[j].start < ost) continue;
                else if (distancePair[j].start > ost) break;
                else cur_cov++;
            }
            if (cur_cov>ft_cov){
                max_depth++;
                continue; //do not complete pairs in the very high coverage regions (e.g. 4 times of median coverage)
            }

            for (LL j = 0; j <= i-1; j++)//this search the possible sub-regions
            {
                if (distancePair[j].start < ost) continue;
                else if (distancePair[j].start > ost) break;
                else//for each sub-region holding the same start point with NRRR, we explore the consecutive regions on the same molecule
                {
                    tempPair = distancePair[j];
                    bool rep = true;
                    for (LL k = j+1; distancePair[k].start<oed&&k<distancePairCount; k++)
                    {
                        if (distancePair[k].mapId==distancePair[j].mapId && distancePair[k].uniqId==distancePair[j].uniqId){ // Dec 3, 2019
                            if (distancePair[k].start == tempPair.end)
                            {
                                tempPair.end = distancePair[k].end;
                                tempPair.distance += distancePair[k].distance;
                            }
                            else if (distancePair[k].start > tempPair.end)
                                                        break;
                        }
                        if (tempPair.end >= oed || distancePair[k].start > tempPair.end) break;
                    }
                    if (tempPair.end == oed) //no need
                    {
                        rep = false;
                        if (tempPair.start != ost)printf("Wrong completion!\n");
                        for (LL l = i; l < distancePairCount; l++)
                        {
                            if (distancePair[i] < distancePair[l])break;
                            else if (distancePair[l].mapId==tempPair.mapId && distancePair[l].uniqId==tempPair.uniqId) {rep = true; break;} // here consider to remove the uniqId identity condition
                        }
                    }
                    if (!rep) {
                        extraDP.push_back(tempPair);
                    }
                }
            }
        }
        else
            continue;
    }
    printf("There are %lld extra distance pairs are added, %lld are removed due to maximum range, %lld removed due to high depth.\n",(LL)extraDP.size(), max_range, max_depth);
    // also need to include the new distance pairs into the original set
    distancePair.resize(distancePairCount+extraDP.size());
    for (LL i = 0; i < (LL)extraDP.size(); i++){
        distancePair[distancePairCount++]=extraDP[i];
    }
}


struct reg{
    LL start;
    LL end;
};

bool advanceLikelihoodDistanceCalculation1(LL start, LL end, int fltLv){
    double referenceDistance = chromosome.position[distancePair[start].end] - chromosome.position[distancePair[start].start];

    vector<double> distanceRatio;
    distanceRatio.clear();

    LL suppCnt = 0;

    vector<LL> maps;
    vector<double> dists;
    maps.clear();
    dists.clear();
    for (LL i=start; i<=end; i++)
    {
        bool fal = true;
        for (LL j = 0; j < (LL)maps.size(); j++)
        {
            try{
                if (maps[j] == stol(distancePair[i].mapId) && abs(dists[j] - distancePair[i].distance) < 1){fal = false;break;}
            }catch(std::invalid_argument&) {
                cout << 1111111 << "###" << distancePair[i].mapId << "###" << endl;
            }
        }
        if (fal)
        {
            try
            {
                maps.push_back(stol(distancePair[i].mapId));
            }catch(std::invalid_argument&){
                cout << 1111112 << "###" << distancePair[i].mapId << "###" << endl;
            }
            dists.push_back(distancePair[i].distance);
            distanceRatio.push_back(distancePair[i].distance/referenceDistance);
            if (abs(distancePair[i].distance - referenceDistance) >= minIndelSize && abs(distancePair[i].distance/referenceDistance-1)>=minIndelRatio)    suppCnt++;////
        }
    }
    LL newId = -1;
    if (suppCnt<=fltLv){
        for (LL i=start; i<=end; i++){
            distanceRatio.push_back(distancePair[i].distance/referenceDistance);
            if (abs(distancePair[i].distance - referenceDistance) >= minIndelSize && abs(distancePair[i].distance/referenceDistance-1)>=minIndelRatio){

                //new in April 21                
                distancePair[i].distance = referenceDistance;
                string tpCA = "-";
                tpCA += distancePair[i].mapId;
                distancePair[i].mapId = tpCA; 
            }
        }
    }
    if (suppCnt>0)
        return true;
    else return false;
}


void detectByDistancePair(int maxClsNum, int fltLv, bool cancer){

    //change in april 19

    sort(distancePair.begin(), distancePair.end());//    sort(distancePair, distancePair + distancePairCount);
    LL tempHead = 0, tempTail = 0, totalCov = 0, totalNum = 0, tmpCov;
    vector<LL> covVec;
    covVec.clear();
    LL nbNum = 0; // Key parameter, set to 0 means no ensemble, otherwise (>0) means ensembles
    vector<int> suspVec; // 2 is the center regions, 1 is the neighboring regions, 0 is the ref regions
    vector<LL> headVec; // 2 is the center regions, 1 is the neighboring regions, 0 is the ref regions
    vector<LL> tailVec; // 2 is the center regions, 1 is the neighboring regions, 0 is the ref regions
    suspVec.clear();
    headVec.clear();
    tailVec.clear();
    bool tmpSusp;
    int susNum = 0;
    LL lastHitEnd, lastHitStart;
    vector<double> svRatVec;
    svRatVec.clear();
    double svRat;


    
    // changed in April 13-14
    distancePair.resize(distancePairCount+1);// printDistancePair();
    distancePair[distancePairCount].start = -1;
    distancePair[distancePairCount].end = -1;
    distancePair[distancePairCount].mapId = "-1";
    distancePairCount++;
    vector<int> cov_list;
    // here add a vector to monitor the coverage distribution
    cov_list.clear();
    printf("For the chromosome %lld, we applied level-%d filter on the distance ratios\n",chrId,fltLv);
    for (LL i=0; i<distancePairCount; i++){
        if (i!=distancePairCount-2 && distancePair[i].start == distancePair[i+1].start && distancePair[i].end == distancePair[i+1].end){
        } else {
            tempTail = i;
            if (tempTail-tempHead+1 >= numberOfSupportIndelMolecule)
                cov_list.push_back(tempTail-tempHead+1);
            bool notPure = advanceLikelihoodDistanceCalculation1(tempHead, tempTail, fltLv);
            if (notPure){
            }
            tempHead = tempTail + 1;
        }
    }
    sort(cov_list.begin(),cov_list.end());
    cout << cov_list[0] << "  " << cov_list[cov_list.size()/4] << "  " << cov_list[cov_list.size()/2] << "  " << cov_list[cov_list.size()-1] << endl;
    int median_cov = cov_list[cov_list.size()/2];
    printf("Start to erase with median_cov %d\n",median_cov);
    sort(distancePair.begin(),distancePair.end());
    LL st;
    for (LL i = 0; i < distancePair.size(); i++){
        st = i;
        if (distancePair[i].start != -1)break;
    }
    distancePair.erase(distancePair.begin(),distancePair.begin()+st);

    distancePairCount = distancePair.size();
    printf("After noise removal(distance ratio level), there are %lld distance pairs!\n",distancePairCount);

    completePairs(median_cov*4);

    sort(distancePair.begin(),distancePair.end());


    int btClsSize = 1;

    distancePair.resize(distancePairCount+1);// printDistancePair();
    distancePair[distancePairCount].start = -1;
    distancePair[distancePairCount].end = -1;
    distancePair[distancePairCount].mapId = "-1";
    distancePairCount++;
    tempHead = 0;
    tempTail = 0;
    for (LL i=0; i<distancePairCount; i++){
        if (i!=distancePairCount-2 && distancePair[i].start == distancePair[i+1].start && distancePair[i].end == distancePair[i+1].end){
        } else {
            tempTail = i;
            tmpSusp = advanceLikelihoodDistanceCalculation2(tempHead, tempTail, tmpCov, svRat);
            if (tmpCov>=btClsSize)svRatVec.push_back(svRat);
            if (tmpSusp){
                suspVec.push_back(2);
                susNum++;
            }
            else{
                suspVec.push_back(0);
            }
            headVec.push_back(tempHead);
            tailVec.push_back(tempTail);
            totalCov += tmpCov;
            covVec.push_back(tmpCov);
            totalNum++;
            tempHead = tempTail + 1;
        }
    }
    distancePair.resize(--distancePairCount);
    // here we extract the median cov as the background coverage
    vector<double> ratVec;
    ratVec.clear();
    LL insCnt = 0, delCnt = 0, refCnt = 0;
    for (LL i=0; i<distancePairCount; i++){
        double referenceDistance = chromosome.position[distancePair[i].end] - chromosome.position[distancePair[i].start];
        double rati = distancePair[i].distance / referenceDistance;
        if (distancePair[i].distance-referenceDistance > 2000 && rati > 1.05)insCnt++;
        else if (distancePair[i].distance-referenceDistance < -2000 && rati < 0.95)delCnt++;
        else if (rati >= 0.95 && rati <= 1.05) refCnt++;
    }
    printf("Statistics \t%lld\t%lld\t%lld\t%lld\n",chrId,delCnt,insCnt,refCnt);


    sort(svRatVec.begin(),svRatVec.end());
    statInfo(svRatVec, 3);
    printf("First pick the suspecting cases: %d cases!\n",susNum);
    double meanCov = totalCov*1.0/totalNum;

    int curSsStart = -1, curSsEnd = -1;
    LL ost = -1, oed = -1;
    LL maxIdx, maxRange = 0;
    vector<distanceType> extraDP;
    vector<distanceType> moreDP;
    moreDP.clear();
    /***************Open part 1: here consider first locate the continuous SV candidate regions, and then exhaustly extract all inside distance regions****************/
    /***************Now, only the original continuous regions as well as the whole region are considered****************/

    
    for (LL i = 0; i < suspVec.size(); i++){
        if (suspVec[i] == 2 && (oed == -1 || oed>= distancePair[headVec[i]].start)){
            if (curSsStart == -1){
                curSsStart = i;
            }
            curSsEnd = i;
            oed = distancePair[headVec[curSsEnd]].end;
        }else if (curSsEnd != -1 && oed < distancePair[headVec[i]].start){
            if (curSsStart<=curSsEnd){
                iterCompletePair(moreDP, headVec,curSsStart,curSsEnd);
            }
            curSsStart = i;
            curSsEnd = i;
            oed = distancePair[headVec[curSsEnd]].end;
        }else{//suspVec[i] != 2
        }
    }
    if (curSsStart<=curSsEnd){
        iterCompletePair(moreDP, headVec,curSsStart,curSsEnd);
    }
    
    sort(moreDP.begin(),moreDP.end());

    printf("After initial pair complete, there are %lld pairs!\n",(LL)moreDP.size());
    for (LL i = 0; i < headVec.size(); i++){
        if (suspVec[i] > 0){
            //here to remove redundant distance pairs
            bool exst = false;
            auto dp_ele = distancePair[headVec[i]];
            for (LL t = 0; t < moreDP.size(); t++){
                if (moreDP[t].start < dp_ele.start || (moreDP[t].start == dp_ele.start && moreDP[t].end < dp_ele.end))continue;
                else if (moreDP[t].start > dp_ele.start || (moreDP[t].start == dp_ele.start && moreDP[t].end > dp_ele.end))break;
                exst = true;
            }
            if (exst)continue;
            for (LL j = headVec[i]; j <= tailVec[i]; j++){
                moreDP.push_back(distancePair[j]);
            }
        }
    }
    vector<distanceType>().swap(distancePair);
    distancePair = moreDP;
    vector<distanceType>().swap(moreDP);
    vector<distanceType>().swap(extraDP);
    distancePairCount = (LL)distancePair.size();

    printf("Added internal continuous regions: %lld cases!\n",distancePairCount);

    tempHead = 0;
    tempTail = 0;
    headVec.clear();
    tailVec.clear();
    vector<bool> centerVec;
    centerVec.clear();
    
    int curNum = 0;
    int center = 0;
    vector<int> clsIdx;
    vector<LL> clsMaps;
    vector<vector<int>> allClsIdx;
    vector<vector<LL>> allClsMaps;
    vector<vector<double>> allClsDat;
    allClsIdx.clear();
    allClsMaps.clear();
    allClsDat.clear();
    vector<double> clsDat;
    vector<int> clsNumVec;
    clsNumVec.clear();
    for (LL i=0; i<distancePairCount; i++){
        if (i!=distancePairCount-1&&distancePair[i].start == distancePair[i+1].start && distancePair[i].end == distancePair[i+1].end){
        } else {
            tempTail = i;
            int rep_time = 0;
            int limit;
            do{
                // now center stores the size of unique maps
                center = advanceLikelihoodDistanceCalculation3(tempHead, tempTail, meanCov, clsIdx, clsMaps, clsDat, maxClsNum, median_cov*4);
                if (center>0)clsNumVec.push_back((int)*max_element(clsIdx.begin(),clsIdx.end()));
                headVec.push_back(tempHead);
                tailVec.push_back(tempTail);
                centerVec.push_back(center);
                allClsIdx.push_back(clsIdx);
                allClsMaps.push_back(clsMaps);
                allClsDat.push_back(clsDat);
                rep_time++;
                limit = int(ceil(log2(center/median_cov/4.0)));
//                limit = 0;
            }while(rep_time < limit);
            tempHead = tempTail + 1;
            if (rep_time>1)cout << "Repeat times is " << rep_time << endl;
        }
    }
    if (clsNumVec.size() == 0)return;

    printf("After get the suspecting pairs only, there are %lld pairs\n",(LL)headVec.size());

    // At this step, all centers and their neighbors have been clustered, and the maps of distances and index of clusters are stored, as well as the vector indicating cneters or not,
    // and the heads and tails of their distance pairs in the whole set of distance pairs
    // Now need to break all results into groups based on each center, and then extract the clustering results of the molecules in the centers, and make an ensemble


    LL lastEnd, firstStart;
    OMSV resSV;
    LL refDistV;
    vector<vector<int>> keptClsIdx;
    vector<vector<LL>> keptClsMaps;
    keptClsIdx.clear();
    keptClsMaps.clear();
    vector<LL> startVec;
    vector<LL> endVec;
    startVec.clear();
    endVec.clear();
    for (int i = 0; i < centerVec.size(); i++){
        cout << i << " out of " << centerVec.size() << " with start " << headVec[i] << " and end " << tailVec[i] << endl;
        if (centerVec[i]){
            // changed in May 30, 2021
            refDistV = chromosome.position[distancePair[headVec[i]].end] - chromosome.position[distancePair[headVec[i]].start];
            if (allClsIdx[i].size() > 1000) {
                printf("Region (%lld:%lld-%lld) has extremely high coverage(>1000)!\n",chrId,headVec[i],tailVec[i]);
                continue;
            }
            resSV.chr = chrId;
            resSV.id = i;
            resSV.start = chromosome.position[distancePair[headVec[resSV.id]].start];
            resSV.end = chromosome.position[distancePair[headVec[resSV.id]].end];
            cout << "There are " << allClsIdx[i].size() << " items and " << allClsIdx[i][allClsIdx[i].size()-1] << " class labels (median coverage is " << median_cov << ")\n";
            bool nonR = postProcess(allClsIdx[i], allClsMaps[i], allClsDat[i], refDistV, resSV);

            if (nonR)allSVs.push_back(resSV);
            cout << "Here the finalIdx need to be changed if really want to use it!" << endl;
            
            keptClsMaps.push_back(allClsMaps[i]);
            startVec.push_back(resSV.start);
            endVec.push_back(resSV.end);
            cout << "Done all" << endl;
        }else if (allClsIdx[i].size()>=15){
            refDistV = chromosome.position[distancePair[headVec[i]].end] - chromosome.position[distancePair[headVec[i]].start];
        }
    }
}


struct OMSV_cls{
    int siz;
    int clsLabel;
    double botmVal;
    double topVal;
    double meanVal;
    double midVal;
    vector<double> dat;
    bool operator<(const OMSV_cls &x) const{
        return (meanVal < x.meanVal);
    }
};


bool postProcess(vector<int>& finalIdx, vector<LL>& finalMap, vector<double> finalDat, LL refDist, OMSV & resSV){
    /********************Here insert the part to identify the SV types according to the clustering index*****************************/
    /*        Post-process: drop off the outlier clusters, and merge the clusters with very close centers            */
    //First merge the close clusters

    double minRatGap = 0;
    double minRatGapStp = 0;
    double glbMinRatio = max(minIndelRatio, minIndelSize*1.0/refDist);
    double glbMinRatioSVs = max(minIndelRatio, minSVDiffSize*1.0/refDist);
    int selI, selJ, k_Num;
    int mods = 3;

    int sNum = finalIdx.size();

    printf("\nSV region (%lld-%lld), the Ref dist: %lld\t the minimum ratio: %g\n",resSV.start, resSV.end,refDist,glbMinRatio);
    k_Num = *max_element(finalIdx.begin(),finalIdx.end());
    if (k_Num>=2){
        do{
            k_Num = *max_element(finalIdx.begin(),finalIdx.end());
            vector<vector<double>> maxInterDist(k_Num, vector<double>(k_Num,0.0));
            vector<vector<double>> cntInterDist(k_Num, vector<double>(k_Num,0.0));
            minRatGap = glbMinRatio;
            minRatGapStp = glbMinRatio;
            selI = -1;
            selJ = -1;
        /**Open part x: here can replace the maximum distance by average distance between groups**/


            if (mods == 3){
                for (int i = 0; i < finalIdx.size(); i++){
                    for (int j = 0; j < finalIdx.size(); j++){
                        maxInterDist[finalIdx[i]-1][finalIdx[j]-1] = max(abs(finalDat[i]-finalDat[j]),maxInterDist[finalIdx[i]-1][finalIdx[j]-1]);
                    }
                }
            }else if (mods == 2){
                for (int i = 0; i < finalIdx.size(); i++){
                    for (int j = 0; j < finalIdx.size(); j++){
                        maxInterDist[finalIdx[i]-1][finalIdx[j]-1] = min(abs(finalDat[i]-finalDat[j]),maxInterDist[finalIdx[i]-1][finalIdx[j]-1]);
                    }
                }
            }else if (mods == 1){
                for (int i = 0; i < finalIdx.size(); i++){
                    for (int j = 0; j < finalIdx.size(); j++){
                        maxInterDist[finalIdx[i]-1][finalIdx[j]-1] += abs(finalDat[i]-finalDat[j]);
                        cntInterDist[finalIdx[i]-1][finalIdx[j]-1]++;
                    }
                }
                for (int i = 1; i <= k_Num; i++){
                    for (int j = 1; j <= k_Num; j++){
                        if (cntInterDist[i-1][j-1]>0){
                            maxInterDist[i-1][j-1] /= cntInterDist[i-1][j-1];
                        }
                    }
                }
            }


            minRatGap = maxInterDist[0][1];
            selI = 0;
            selJ = 1;


            for (int i = 0; i < k_Num; i++){
                for (int j = i+1; j<k_Num; j++){
                    if (maxInterDist[i][j] < minRatGap){
                        minRatGap = maxInterDist[i][j];
                        selI = i;
                        selJ = j;
                    }
                }
            }

            if (minRatGap >= glbMinRatio)break;
            for (int i = finalIdx.size() - 1; i >=0; i--){
                if(finalIdx[i] == selJ+1){
                    finalIdx[i] = selI + 1;
                }else if (finalIdx[i] > selJ +1)
                    finalIdx[i]--;
            }
            if (k_Num == 2){
                k_Num--;
                break;
            }
        }while(1); // Each time pick the most adjacent (with minimum max-inter cluster distance) cluster pair
    }
    printf("Merged adjacent clusters\n");
    // then remove the clusters whose sizes are less than minClsSize
    vector<int> clsSize(k_Num,0);
    for (int i = 0; i < finalIdx.size(); i++){
        clsSize[finalIdx[i]-1]++;
    }
    int minClsSize = 1;



    glbMinRatio = max(minIndelRatio, minIndelSize*1.0/refDist);
    glbMinRatioSVs = max(minIndelRatio, minSVDiffSize*1.0/refDist);
    for (int i = finalIdx.size()-1; i >= 0; i--){
        if (clsSize[finalIdx[i]-1]< minClsSize){
            finalMap.erase(finalMap.begin()+i);
            finalIdx.erase(finalIdx.begin()+i);
            finalDat.erase(finalDat.begin()+i);
        }
    }
    for (int i = k_Num-1; i >= 0; i--){
        if (clsSize[i] < minClsSize){
            for (int j = 0; j < finalIdx.size(); j++){
                if (finalIdx[j]>i+1)finalIdx[j]--;
            }
        }
    }


    if (finalIdx.size() == 0){
        resSV.type = "Reference";
        return false;
    }

    cout << "The finalidx size is " << finalIdx.size() << endl;

    int sel_m = 1; // here 1 denote median linkage, otherwise average linkage
    k_Num = *max_element(finalIdx.begin(), finalIdx.end());
    while (k_Num > 1){
        vector<vector<double>> grpDist(k_Num, vector<double>(k_Num,0.0));
        vector<double> grpCent(k_Num,0.0);
        vector<double> clsCps;
        for (int j = 0; j < k_Num; j++){
            clsCps.clear();
            for (int i = 0; i < finalDat.size(); i++){
                if (finalIdx[i] == j+1){
                    clsCps.push_back(finalDat[i]);
                }
            }
            sort(clsCps.begin(),clsCps.end());
            if (sel_m == 1){
                grpCent[j] = clsCps[(int)clsCps.size()/2];
            }else {
                grpCent[j] = accumulate(clsCps.begin(),clsCps.end(),0.0) / clsCps.size();
            }
        }
        for (int i = 0; i < k_Num; i++)
            for (int j = 0; j < k_Num; j++)
                grpDist[i][j] = abs(grpCent[i] - grpCent[j]);
        
        minRatGap = grpDist[0][1];
        selI = 0;
        selJ = 1;
        

        for (int i = 0; i < k_Num; i++){
            for (int j = i+1; j < k_Num; j++){
                if (grpDist[i][j] < minRatGap){
                    minRatGap = grpDist[i][j];
                    selI = i;
                    selJ = j;
                }
            }
        }
        /*If the selected two groups are all SV-supporting groups, we need to use SV size gap but not SV size thre to decide if they should be merged*/
        cout << "CHECKCHECK:  " << grpCent[selI] << "\t" << grpCent[selJ] << "\t" << glbMinRatio << "\t" << minRatGap << "\t" << glbMinRatioSVs << endl;
        if (abs(grpCent[selI]-1)>=glbMinRatio && abs(grpCent[selJ]-1)>=glbMinRatio){
            if (minRatGap >= glbMinRatioSVs)break;
        }else{
            if (minRatGap >= glbMinRatio)break;
        }
        for (int i = finalIdx.size() - 1; i >=0; i--){
            if(finalIdx[i] == selJ+1){
                finalIdx[i] = selI + 1;
            }else if (finalIdx[i] > selJ +1)
                finalIdx[i]--;
        }
        cout << 4 << endl;
        

        k_Num = *max_element(finalIdx.begin(), finalIdx.end());
    }



    printf("Then the cluster information is: \n");

    vector<OMSV_cls> allClsInfo;
    allClsInfo.clear();
    int red_Num = 0;
    for (int i = 1; i <= k_Num; i++){
        OMSV_cls tmpCls;
        tmpCls.siz = 0;
        tmpCls.botmVal = 100000;
        tmpCls.topVal = 0;
        tmpCls.meanVal = 0;
        tmpCls.clsLabel = i;
        tmpCls.dat.clear();
        for (int j = 0; j < finalIdx.size(); j++){
            if (finalIdx[j] == i){
                tmpCls.siz++;
                tmpCls.botmVal = min(finalDat[j],tmpCls.botmVal);
                tmpCls.topVal = min(finalDat[j],tmpCls.topVal);
                tmpCls.meanVal += finalDat[j];
                tmpCls.dat.push_back(finalDat[j]);
            }
        }
        sort(tmpCls.dat.begin(),tmpCls.dat.end());
        tmpCls.meanVal /= tmpCls.siz;
        tmpCls.midVal = tmpCls.dat[round((int)tmpCls.dat.size()/2)];
        if (sel_m == 1)tmpCls.meanVal = tmpCls.midVal;
        if (((1 - tmpCls.meanVal)>glbMinRatio) && tmpCls.siz<minClsSD)
            red_Num++;
        else if (((tmpCls.meanVal - 1) > glbMinRatio) && tmpCls.siz<minClsSI)
            red_Num++;
        else if (((1 - tmpCls.meanVal)<=glbMinRatio)&&((tmpCls.meanVal - 1) <= glbMinRatio)&& tmpCls.siz<max(minClsSD,minClsSI))
            red_Num++;
        else{
            allClsInfo.push_back(tmpCls);
            printf("kept!\n");
        }
    }
    k_Num -= red_Num;
    cout << allClsInfo.size() << "  " << k_Num << endl;
    
    if (allClsInfo.size() == 0) return false;
    sort(allClsInfo.begin(),allClsInfo.end());
    /*    SV classification: identify the types of each clusters, and then decide the SV type, zygosity, and/or copy number    */
    resSV.support = 0;
    resSV.coverage = 0;
    if (k_Num > 2){ // complex SVs, e.g. mixed indels, or CNVs
        char sizT[10];
        resSV.zyg = " ";
        resSV.difSiz1 = 1000000000;
        resSV.difSiz2 = -1000000000;
        bool insHas=false;
        bool delHas=false;
        for (int i = 0; i < k_Num; i++){
            resSV.support += allClsInfo[i].siz;
            resSV.coverage += allClsInfo[i].siz;
            sprintf(sizT,"%.1fkb ",(allClsInfo[i].midVal - 1)*refDist/1000);
            resSV.zyg += sizT;
            resSV.difSiz1 = min((allClsInfo[i].meanVal - 1)*refDist,resSV.difSiz1);
            resSV.difSiz2 = max((allClsInfo[i].meanVal - 1)*refDist,resSV.difSiz2);
        }
        if (resSV.difSiz1>(0-minIndelSize) && resSV.difSiz2>minIndelSize)resSV.type = "Insertion-ND";
        else if (resSV.difSiz1<(0-minIndelSize) && resSV.difSiz2<minIndelSize)resSV.type = "Deletion-ND";
        else resSV.type = "Both-ND";
    }else{ // when there is only 1 or 2 clusters, then call homozygous or heterozygous or reference type
        if (k_Num == 2){
            resSV.support = allClsInfo[0].siz + allClsInfo[1].siz;
            resSV.coverage = allClsInfo[0].siz + allClsInfo[1].siz;
            if (allClsInfo[0].meanVal - 1>glbMinRatio){
                resSV.type = "Insertions";
                resSV.zyg = "HI1I2";
                for (int i = 0 ; i < allClsInfo[0].dat.size(); i++){
                    if (allClsInfo[0].dat[i] - 1 <= glbMinRatio)resSV.support--;
                }
            }else if (1 - allClsInfo[1].meanVal>glbMinRatio){
                resSV.type = "Deletions";
                resSV.zyg = "HD1D2";
                for (int i = allClsInfo[1].dat.size() - 1 ; i >= 0; i--){
                    if (1- allClsInfo[1].dat[i] <= glbMinRatio)resSV.support--;
                }
                
            }else if (1 - allClsInfo[0].meanVal>glbMinRatio && allClsInfo[1].meanVal - 1>glbMinRatio){
                resSV.type = "Both";
                resSV.zyg = "HDI";
                for (int i = 0 ; i < allClsInfo[0].dat.size(); i++){
                    if (1 - allClsInfo[0].dat[i] <= glbMinRatio)resSV.support--;
                }
                for (int i = allClsInfo[1].dat.size() - 1 ; i >= 0; i--){
                    if (allClsInfo[1].dat[i] - 1 <= glbMinRatio)resSV.support--;
                }
            }else if (allClsInfo[1].meanVal - 1>glbMinRatio){
                resSV.type = "Insertion";
                resSV.zyg = "Heterozygous";
                for (int i = 0 ; i < allClsInfo[0].dat.size(); i++){
                    if (allClsInfo[0].dat[i] - 1 <= glbMinRatio)resSV.support--;
                }
                for (int i = allClsInfo[1].dat.size() - 1 ; i >= 0; i--){
                     if (allClsInfo[1].dat[i] - 1 <= glbMinRatio)resSV.support--;
                }
            }else if (1 - allClsInfo[0].meanVal>glbMinRatio){
                resSV.type = "Deletion";
                resSV.zyg = "Heterozygous";
                for (int i = 0 ; i < allClsInfo[0].dat.size(); i++){
                    if (1 - allClsInfo[0].dat[i] <= glbMinRatio)resSV.support--;
                }
                for (int i = allClsInfo[1].dat.size() - 1 ; i >= 0; i--){
                    if (1 - allClsInfo[1].dat[i] <= glbMinRatio)resSV.support--;
                }
            }else {// reference type into 2 groups
                resSV.type = "Reference";
                resSV.zyg = "Heterozygous";
                return false;
                for (int i = 0 ; i < allClsInfo[0].dat.size(); i++){
                    if (1 - allClsInfo[0].dat[i] > glbMinRatio)resSV.support--;
                }
                for (int i = allClsInfo[1].dat.size() - 1 ; i >= 0; i--){
                    if (allClsInfo[1].dat[i] - 1 > glbMinRatio)resSV.support--;
                }
            }
            resSV.difSiz1 = (allClsInfo[0].meanVal - 1)*refDist;
            resSV.difSiz2 = (allClsInfo[1].meanVal - 1)*refDist;                
        }else if (k_Num == 1){
            resSV.support = allClsInfo[0].siz;
            resSV.coverage = allClsInfo[0].siz;
            resSV.zyg = "Homozygous";
            if (allClsInfo[0].meanVal - 1>glbMinRatio){
                resSV.type = "Insertion";
                for (int i = 0 ; i < allClsInfo[0].dat.size(); i++){
                    if (allClsInfo[0].dat[i] - 1 <= glbMinRatio)resSV.support--;
                }
            }else if (1 - allClsInfo[0].meanVal > glbMinRatio){
                resSV.type = "Deletion";
                if (allClsInfo[0].siz < minClsSD)return false;
                for (int i = allClsInfo[0].dat.size() - 1 ; i >= 0; i--){
                    if (1 - allClsInfo[0].dat[i] <= glbMinRatio)resSV.support--;
                }
            }else{
                resSV.type = "Reference";
                return false;
                for (int i = 0 ; i < allClsInfo[0].dat.size(); i++){
                    if (1 - allClsInfo[0].dat[i] > glbMinRatio || allClsInfo[0].meanVal - 1>glbMinRatio)resSV.support--;
                }
            }
            resSV.difSiz1 = (allClsInfo[0].meanVal - 1)*refDist;
            resSV.difSiz2 = (allClsInfo[0].meanVal - 1)*refDist;
        }
    }
    resSV.score = resSV.support / resSV.coverage;
    printf("Done post-process!\n");
    return true;
}

bool advanceLikelihoodDistanceCalculation2(LL start, LL end, LL &cnt, double &svRat){
    //change in april 24
    double minimumRatio = minSVSuppRat; // here fix the minimum ratio as 0.2 for the 1st-run check

    double referenceDistance = chromosome.position[distancePair[start].end] - chromosome.position[distancePair[start].start];

    cnt = 0;
    vector<double> distanceRatio;
    distanceRatio.clear();


    LL suppCnt = 0;

    vector<LL> maps;
    vector<double> dists;
    maps.clear();
    dists.clear();
    for (LL i=start; i<=end; i++)
    {
        bool fal = true;
        for (LL j = 0; j < (LL)maps.size(); j++)
        {
            try{
                if (maps[j] == stol(distancePair[i].mapId) && abs(dists[j] - distancePair[i].distance) < 1){fal = false;break;}
            }catch(std::invalid_argument&){
                cout << 1111113 << "###" << distancePair[i].mapId << "###" << endl;
            }
        }
        if (fal)
        {
            try{
                maps.push_back(stol(distancePair[i].mapId));
            }catch(std::invalid_argument&){
                cout << 1111114 << "###" << distancePair[i].mapId << "###" << endl;
            }
            dists.push_back(distancePair[i].distance);
            distanceRatio.push_back(distancePair[i].distance/referenceDistance);
            cnt++;
            if (abs(distancePair[i].distance - referenceDistance) >= minIndelSize && abs(distancePair[i].distance/referenceDistance-1)>=minIndelRatio) suppCnt++;////
        }
    }
    int outNum = 0;
    svRat = suppCnt*1.0/cnt;
    LL minGap = referenceDistance - *min_element(dists.begin(),dists.end());
    LL maxGap = -(referenceDistance - *max_element(dists.begin(),dists.end()));
    double suppRate = suppCnt*1.0/cnt; // here removed the leading and tailing one distance pair
    if (suppRate < minimumRatio)
        return false;
    else
        return true;
}


template <typename T>
vector<int> sort_indexes(const vector<T> &v) {
  // initialize original index locations
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}


int advanceLikelihoodDistanceCalculation3(LL start, LL end, double cnt, vector<int> &idx, vector<LL> &maps, vector<double> &distanceRatio, int maxClsNum, int ft_cov=200){
    //change in april 24 3PM Indel100
    double minimumRatio = minSVSuppRat;

    double referenceDistance = chromosome.position[distancePair[start].end] - chromosome.position[distancePair[start].start];
    cnt = 0;
    distanceRatio.clear();
    LL suppCnt = 0;
    vector<double> dists;
    maps.clear();
    dists.clear();
    for (LL i=start; i<=end; i++)
    {
        bool fal = true;
        for (LL j = 0; j < (LL)maps.size(); j++)
        {
            try{
                if (maps[j] == stol(distancePair[i].mapId) && abs(dists[j] - distancePair[i].distance) < 1){fal = false;break;}
            }catch(std::invalid_argument&){
                cout << 1111115 << "###" << distancePair[i].mapId << "###" << endl;
            }
        }
        if (fal)
        {
            try{
                maps.push_back(stol(distancePair[i].mapId));
            }catch(std::invalid_argument&){
                cout << 1111116 << "###" << distancePair[i].mapId << "###" << endl;
            }
            dists.push_back(distancePair[i].distance);
            cnt++;
            if (abs(distancePair[i].distance - referenceDistance) >= minIndelSize && abs(distancePair[i].distance/referenceDistance-1)>=minIndelRatio) suppCnt++;////
        }
    }
    int real_size = maps.size();



    vector<LL> maps_sf;
    vector<double> dists_sf;

    maps_sf = maps;
    dists_sf = dists;


// here we sort the distances
    vector<LL> s_maps;
    vector<double> s_dists;
    s_maps.clear();
    s_dists.clear();
    for (auto i: sort_indexes(dists_sf)){
        s_maps.push_back(maps_sf[i]);
        s_dists.push_back(dists_sf[i]);
        distanceRatio.push_back(dists_sf[i]/referenceDistance);
    }
    maps = s_maps;
    dists = s_dists;
    int outNum = 0;
    LL minGap = referenceDistance - dists[0];
    LL maxGap = -(referenceDistance - dists[dists.size()-1]);
    if (minGap > max(minIndelSize,(LL)(referenceDistance*minIndelRatio)))
        outNum++;
    if (maxGap > max(minIndelSize,(LL)(referenceDistance*minIndelRatio)))
        outNum++;    
    double suppRate = suppCnt*1.0/cnt; // here removed the leading and tailing one distance pair

    int maxNum;
    int selK1, selK2;
    map<int,vector<int>> idxs;
    idx.clear();
    if (cnt < numberOfSupportIndelMolecule || (suppRate < minimumRatio)){
        int stIdx=1;
        if (cnt>=numberOfSupportIndelMolecule)
            for (int i = 0; i < distanceRatio.size(); i++){
                    idx.push_back(stIdx++);
            }
        return 0;
    }else{
        int stIdx = 1;
        for (int i = 0; i < distanceRatio.size(); i++){
                idx.push_back(stIdx++);
        }
        return real_size;
    }
}


double statInfo(vector<double> vec, int mod){
    double minL, maxL, stp, aver;
    if (mod == 1){
        minL = 0.9;
        maxL = 1.1;
        stp = (maxL - minL)/100;
        aver = accumulate(vec.begin(),vec.end(),0.0)/vec.size();
        printf("The info for mean ratio vector (molecule level): (global: %g)\n",aver);
    }
    else if (mod == 2){
        minL = 0;
        maxL = 1;
        stp = 0.01;
        aver = accumulate(vec.begin(),vec.end(),0.0)/vec.size();
        printf("The info for SV ratio vector (molecule level): (global: %g)\n",aver);
    }
    else {
        minL = 0;
        maxL = 1;
        stp = 0.01;
        aver = accumulate(vec.begin(),vec.end(),0.0)/vec.size();
        printf("The info for distance ratio vector (SV level): (global: %g)\n",aver);
    }
    int cnt = 1;
    vector<int> stCnt(102,0);
    
    for (int j = 0; j < vec.size(); j++){
        if (vec[j] <= minL){
            stCnt[cnt-1]++;
            continue;
        }
        if (vec[j] > maxL){
            stCnt[101]++;
            continue;
        }
        if (vec[j] > minL + stp*cnt)
            cnt++;
        stCnt[cnt]++;
    }
    if (mod == 4){
        printf("<=%g\t%d\n",minL,stCnt[0]);

        for (int i = 1; i <= 100; i++){
            printf("%g-%g\t%d\n",minL+(i-1)*stp,minL+i*stp,stCnt[i]);
        }
        printf(">%g\t%d\n",maxL,stCnt[101]);
    }
    return aver;
}

struct tpMol{
    LL idx;
    LL start;
    LL end;
    LL depth;
    vector<distanceType> distP;
    bool operator<(const tpMol &x) const{
        return (start < x.start || (start == x.start && end < x.end));
    }
};

/*Here insert a struct to store the sequence of markers for molecules*/
/*Use HitEnum to get the signal, two bits for each site: M->1(bit 1), D->0(bit 1), and I->1(bit 2), no I->0(bit 2)*/
/*Can also consider the distance between sites, one way is also encoding the distance into the above sequence, the other is using ensemble way to combine the clustering results based on both bit sequences and distances*/
/*First extract this sequence vector, and then select the molecules covering the candidate region and extract a feature matrix for clustering*/
struct markVec{
    LL mapId;
    LL start;
    LL end;
    LL startIdx;
    LL endIdx;
    vector<double> feaVec;
    bool operator<(const markVec &x) const{
        return (start < x.start || (start == x.start && end < x.end));
    }
};
vector<markVec> markSeq;
/*Added in Feb 4, 2019*/

double parseHitEnum(){
    LL tempCount = 0;
    distancePair.clear();
    distancePair.resize(5000000);
    LL molStart, molEnd;
    vector<double> ratioVec;    
    vector<double> ratioVecG5;
    vector<LL> refVec;
    vector<double> ratioVecSort;
    vector<double> ratioMeanVec;
    ratioMeanVec.clear();
    vector<double> SvRatioVec;
    SvRatioVec.clear();
    LL countF = 0;
    LL countF2 = 0;
    LL countF3 = 0;
    double varRate = 0.75;


    vector<distanceType> dP;
    vector<tpMol> allRmMol;
    allRmMol.clear();
    tpMol rmMol;
    LL minFNS = 5;
    markSeq.clear();
    
    cout << "There are " << opticalMap.size() << " molecules extracted!\n";

    uniform_real_distribution<double> unif(0.0,1.0);
    default_random_engine re;
    double a_random_double;

    for (LL i=0; opticalMap[i].refStart != -1 || opticalMap[i].refEnd != -1; i++){
        a_random_double = unif(re);
        if (opticalMap[i].score < confidenceLimit) continue; 
        LL RefIndex = upper_bound(chromosome.position, chromosome.position + chromosome.numberOfSites, opticalMap[i].refStart) - chromosome.position - 1;
        LL OpticalIndex = upper_bound(opticalMap[i].oldPosition.begin(), opticalMap[i].oldPosition.end(), opticalMap[i].optStart) - opticalMap[i].oldPosition.begin() - 1;
        LL BeginRefIndex = RefIndex;
        LL BeginOpticalIndex = OpticalIndex;
        LL EndRefIndex = lower_bound(chromosome.position, chromosome.position + chromosome.numberOfSites, opticalMap[i].refEnd - 1) - chromosome.position;
        LL countM = 0;
        LL countD = 0;
        LL countI = 0;
        char previousChar = 'M';
        double tempDistance = 0;
        LL PreviousIndex = 0;
        if (opticalMap[i].hitEnum[0] == 'F'){
            distancePair[distancePairCount].start = BeginRefIndex;
            distancePair[distancePairCount].end = EndRefIndex;
            distancePair[distancePairCount].mP = false;
            distancePair[distancePairCount].distance = opticalMap[i].optEnd - opticalMap[i].optStart;
            distancePair[distancePairCount].mapId = opticalMap[i].mapId;
            distancePair[distancePairCount].uniqId = a_random_double;
            distancePairCount++;
            continue;
        }
        /*added Feb4,2019*/
        markVec tpMV;
        vector<double>().swap(tpMV.feaVec);
        tpMV.end = chromosome.position[EndRefIndex];
        tpMV.start = chromosome.position[BeginRefIndex];
        tpMV.startIdx = BeginRefIndex;
        tpMV.endIdx = EndRefIndex;
        try{
            tpMV.mapId = stol(opticalMap[i].mapId);
        }catch(std::invalid_argument&){
            cout << 1111117 << "###" << opticalMap[i].mapId << "###" << endl;
        }
        LL queryLen = 0;
        LL endTrim = 1, endCnt = 0;
        molStart = distancePairCount;
        for (LL j=0; j < opticalMap[i].hitEnum.size(); j++){
            if (opticalMap[i].hitEnum[j] >= '0' && opticalMap[i].hitEnum[j] <= '9'){
                tempCount = tempCount * 10 + (opticalMap[i].hitEnum[j]-'0');
            }
            else {
                if (opticalMap[i].hitEnum[j] == 'M'){
                    for (LL k=0; k<tempCount; k++){ 
                        if (tpMV.feaVec.size()>0&&previousChar == 'M')tpMV.feaVec.push_back(0);
                        else if (tpMV.feaVec.size()>0&&previousChar == 'I')tpMV.feaVec.push_back(1);//Here do not count the number of insert sites
                        else if (tpMV.feaVec.size()>0&&previousChar == 'D')tpMV.feaVec.push_back(0);
                        tpMV.feaVec.push_back(1);
                        if (OpticalIndex != BeginOpticalIndex){
                            tempDistance += opticalMap[i].position[OpticalIndex];
                            distancePair[distancePairCount].start = PreviousIndex;
                            distancePair[distancePairCount].end = RefIndex;
                            if (previousChar == 'M')
                                distancePair[distancePairCount].mP = true;
                            else
                                distancePair[distancePairCount].mP = false;
                            distancePair[distancePairCount].distance = tempDistance;
                            distancePair[distancePairCount].mapId =  opticalMap[i].mapId;
                            distancePair[distancePairCount].uniqId =  a_random_double;
                            endCnt++;
                            if (endCnt>endTrim)/*here use endTrim to trim the start end sites*/
                                distancePairCount++;
                            tempDistance = 0;
                        }
                        if (OpticalIndex != BeginOpticalIndex)chromosome.coverage[RefIndex]++;
                        chromosome.occurrence[RefIndex]++;
                        if (OpticalIndex > BeginOpticalIndex+1)chromosome.gapCoverage[RefIndex]++;
                        previousChar = 'M';
                        OpticalIndex++;

                        PreviousIndex = RefIndex;
                        countM += tempCount;
                        RefIndex++;
                    }
                } else if (opticalMap[i].hitEnum[j] == 'D'){
                    for (LL k=0; k<tempCount; k++){
                        /*Added in Feb4, 2019*/
                        if (previousChar == 'M')tpMV.feaVec.push_back(0);
                        else if (previousChar == 'I')tpMV.feaVec.push_back(1);//Here do not count the number of insert sites
                        else if (previousChar == 'D')tpMV.feaVec.push_back(0);
                        tpMV.feaVec.push_back(0);
                        if (OpticalIndex != BeginOpticalIndex)chromosome.coverage[RefIndex+k]++;
                        if (OpticalIndex > BeginOpticalIndex+1)chromosome.gapCoverage[RefIndex+k]++;
                        previousChar = 'D';
                    }
                    RefIndex += tempCount;
                    countD += tempCount;
                } else if (opticalMap[i].hitEnum[j] == 'I'){
                    double tempDouble = 0;
                    for (LL k=0; k<tempCount; k++){
                        if (previousChar != 'D'){
                            tempDouble += opticalMap[i].position[OpticalIndex];
                        }            //is there any problem? The length of this Inserted segment will be padded if followed by 'M'
                        tempDistance += opticalMap[i].position[OpticalIndex];
                        previousChar = 'I';
                        OpticalIndex++;
                    }
                    countI += tempCount;
                }
                tempCount = 0;
            }
        }
        markSeq.push_back(tpMV);
        if (RefIndex>=2){
            chromosome.gapCoverage[RefIndex-1]--;
            chromosome.gapCoverage[RefIndex-2]--;
            chromosome.coverage[RefIndex-1]--;
        }
        if (distancePairCount>distancePair.size()-10000)distancePair.resize(distancePair.size()+100000);


        /*Here use endTrim to remove the stop ends*/ 
        for (LL k = 0; k < min(endTrim,endCnt-endTrim); k++)distancePairCount--;
        if (endCnt <= 2*endTrim) continue;
    ////////////here is the new component to normalize distances per molecule, use median distance ratio on the molecule to do the normalization

        molEnd = distancePairCount;
        refVec.clear();
        ratioVec.clear();
        ratioVecG5.clear();
        for (LL k = molStart; k < molEnd; k++){
            queryLen += distancePair[k].distance;
                LL refDist = chromosome.position[distancePair[k].end] - chromosome.position[distancePair[k].start];
                ratioVec.push_back(distancePair[k].distance*1.0/refDist);
                refVec.push_back(refDist);
                if (refDist >= 5000 || distancePair[k].distance > 5000){
                    ratioVecG5.push_back(distancePair[k].distance*1.0/refDist);
                }
        }
        vector<double>().swap(ratioVecSort);
        ratioVecSort = ratioVec;
        sort(ratioVecSort.begin(),ratioVecSort.end());
        sort(ratioVecG5.begin(),ratioVecG5.end());
        double divRatio;
        if (ratioVecG5.size()<5)
            divRatio = 1;
        else
            divRatio = ratioVecG5[(int)ratioVecG5.size()/2]; //here use the round of integer to locate the index
        
        double meanRatio = accumulate(ratioVecG5.begin(), ratioVecG5.end(), 0.0)/ratioVecSort.size()/divRatio; //here use the round of integer to locate the index
        ratioMeanVec.push_back(meanRatio);
    ///////////If over half of the distances are supporting distance changes, then treat this molecule as badly aligned, and remove it
    /*******************Here add the component to make the alignment quality filter*****************/


        int cntRef = ratioVec.size();
        for (int t = 0; t < ratioVec.size(); t++){
            auto ele = ratioVec[t];
            if (abs((ele / divRatio)-1.0) > max(minIndelRatio,minIndelSize/refVec[t]))cntRef--;
        }
        /**Open part y: here can replace 0.5 by a larger number, and consider to remove the molecules based on the pair number**/
        //Here also add one other filter to remove the molecule being aligned with extremely low score 
        
        bool qltFlt = countM<2*(countI+countD);
        if ( qltFlt ){
            countF++;
        }
        if (cntRef*1.0 < ratioVec.size()*varRate)countF2++;
        for (LL k = molStart; k < molEnd; k++){
            distancePair[k].distance /= divRatio;
        } 
        if ( cntRef*1.0 < ratioVec.size()*varRate || qltFlt ){
            countF3++;
            distancePairCount = molStart;

            rmMol.idx = i;            
            dP.clear();
            for (LL k = molStart; k < molEnd; k++)dP.push_back(distancePair[k]);
            rmMol.distP = dP;
            allRmMol.push_back(rmMol);
        }else{
            SvRatioVec.push_back((double)ratioVec.size()-cntRef);
        }    

    }

    sort(markSeq.begin(),markSeq.end());

    LL countF4 = 0;
    distancePair.resize(distancePairCount);

    for (LL k = 0; k < allRmMol.size(); k++){
        for (LL l = 0; l < allRmMol[k].distP.size(); l++){
            LL refDist = chromosome.position[allRmMol[k].distP[l].end] - chromosome.position[allRmMol[k].distP[l].start];
            if (abs(allRmMol[k].distP[l].distance*1.0/refDist-1) > max(minIndelRatio,minIndelSize/refDist)){
                allRmMol[k].start = chromosome.position[allRmMol[k].distP[l].start];
                break;
            }
        }
        for (LL l = allRmMol[k].distP.size()-1; l >=0;  l--){
            LL refDist = chromosome.position[allRmMol[k].distP[l].end] - chromosome.position[allRmMol[k].distP[l].start];
            if (abs(allRmMol[k].distP[l].distance*1.0/refDist-1) > max(minIndelRatio,minIndelSize/refDist)){
                allRmMol[k].end = chromosome.position[allRmMol[k].distP[l].end];
                break;
            }
        }
        allRmMol[k].depth = 0;
        for (LL l = max(allRmMol[k].idx-10000,(LL)0); l < min(allRmMol[k].idx+10000,(LL)opticalMap.size()); l++){
            if (opticalMap[l].refStart < allRmMol[k].end && opticalMap[l].refEnd > allRmMol[k].end)allRmMol[k].depth++;
            if (opticalMap[l].refStart > allRmMol[k].end)break;
        }
    }
    sort(allRmMol.begin(),allRmMol.end());
    for (LL k = 0; k < allRmMol.size();){
        LL minFN = max(ceil(allRmMol[k].depth*minClsSizRat/2),minFNS);
        LL cont = 1;
        LL grpEnd = allRmMol[k].end;
        for (LL l = k+1; l < allRmMol.size(); l++){
            if (allRmMol[l].start < grpEnd){
                cont++;
                grpEnd = min(grpEnd,allRmMol[l].end);
            }
            else 
                break;
        }
        if (cont>=minFN){
            for (LL m = k; m < k + cont; m++){
                countF4++;
                for (LL n = 0; n < allRmMol[m].distP.size(); n++){
                    distancePair.push_back(allRmMol[m].distP[n]);
                }
            }
            k += cont;
        }else{
            k++;
        }
    }
    sort(ratioMeanVec.begin(),ratioMeanVec.end());
    sort(SvRatioVec.begin(),SvRatioVec.end());
    if (SvRatioVec.size() == 0)return 0;
    double mxNum = *max_element(SvRatioVec.begin(),SvRatioVec.end());
    vector<int> distSVCov(mxNum+1,0);
    
    for (LL i = 0; i < SvRatioVec.size(); i++){
        distSVCov[(int)SvRatioVec[i]]++;
    }
    
    int totN = 0;
    int fltLv = 0;
    double averSVRat = statInfo(SvRatioVec, 2);
    
    printf("Before Resume, there are %lld distance pairs, After resume, there are %lld distance pairs\n",distancePairCount,(LL)distancePair.size());
    
    distancePairCount = distancePair.size();
    printf("There are %lld molecules being filtered due to low alignment accuracy, %lld molecules being filtered due to much wrong distances, and totally %lld being filtered!\n",countF, countF2, countF3);
    printf("Then %lld have been resumed!\n",countF4);
    opticalMap.clear();
    vector<optAlignType>().swap(opticalMap);

    return fltLv;
}



void initData(){
    memset(chromosome.position, 0, sizeof(chromosome.position));
    memset(chromosome.coverage, 0, sizeof(chromosome.coverage));
    memset(chromosome.occurrence, 0, sizeof(chromosome.occurrence));
    memset(chromosome.gapCoverage, 0, sizeof(chromosome.gapCoverage));
    memset(chromosome.gapCount, 0, sizeof(chromosome.gapCount));
    chromosome.length = 0;
    chromosome.numberOfSites = 0;
    distancePairCount = 0;
}

struct triEle{
    string zyg;
    string type;
    int flagIdx;
    int selIdx; // store the index of the SV with the narrowest regions and similar size with the one having highest coverage
    double coverage;
    int num;
};

void correctOverlapVariant(){
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
            if (opL > min(lstSV.end - lstSV.start, allSVs[j].end - allSVs[j].start)*0.5 && allSVs[i].end >= allSVs[j].start){
            // 50% intra-overlapping 
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
                        if ( min(max(abs(grpSVs[s].difSiz1),abs(grpSVs[s].difSiz2)),max(abs(grpSVs[fidx].difSiz1),abs(grpSVs[fidx].difSiz2))) > min(max(abs(grpSVs[s].difSiz1),abs(grpSVs[s].difSiz2)),max(abs(grpSVs[fidx].difSiz1),abs(grpSVs[fidx].difSiz2)))*0.8 ){
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
                        //here use the one with highest coverage to represent
                        flIdx = triVec[s].flagIdx;
                    }
                }
            }
            uniqSVs.push_back(grpSVs[flIdx]);
        }
    }
    allSVs = uniqSVs;
    printf("After removing the redundant SVs, there are %lld left\n",(LL)allSVs.size());
}


void outputVariantResult(int argv, char* argc[], char* outputFolder, int chrID){
    char fnameTP[100];
    if (chrID != 0)
        sprintf(fnameTP,"%s/Label%lld_Chr%lld_%lld.osv",outputFolder,inputFlag,chrID,(LL)minIndelSize);
    else
        sprintf(fnameTP,"%s/Label%lld_All_%lld.osv",outputFolder,inputFlag,chrID,(LL)minIndelSize);
    outputVariantResultFile = fopen(fnameTP,"w");
    fprintf(outputVariantResultFile, "#");
    for (int i = 0; i < argv; i++)
        fprintf(outputVariantResultFile,"%s ",argc[i]);
    fprintf(outputVariantResultFile, "\n#Chr\tStart\tStop\tType\tBlank\tSize1(for diploid)\tID\tSize2(for diploid)\tZygosity(Sizes for >2 haplotypes)\tConfidence\tScore\tCoverage\tSupport\n");
    LL variantCallId = 1;
    numberOfSupportedSV = 0;
    for (LL i=0; i<allSVs.size(); i++){
        if (allSVs[i].type!="Reference"){
            if (allSVs[i].type=="Deletion"){
                double tpSiz = allSVs[i].difSiz2;
                allSVs[i].difSiz2 = allSVs[i].difSiz1;
                allSVs[i].difSiz1 = tpSiz;
            }
            fprintf(outputVariantResultFile, "%lld\t%lld\t%lld\t%s\t%d\t%.1f\t%d\tsizeChange=%.1f\t%s\t%g\t%g\t%g\t%g\n", allSVs[i].chr, allSVs[i].start, allSVs[i].end, allSVs[i].type.c_str(), 0,  allSVs[i].difSiz1, 0, allSVs[i].difSiz2, allSVs[i].zyg.c_str(), allSVs[i].score, allSVs[i].score,allSVs[i].coverage, allSVs[i].support);
            numberOfSupportedSV++;
        }
    }
    fclose(outputVariantResultFile);
}

int main(int argv, char *argc[]){
    setDefault();
    memset(outputFolder,0,sizeof(outputFolder));
    char chrMapFile[1000];
    int chrID = 0;
    int maxClsNum = -1;
    int cMode = 0;
    bool compul = false;
    memset(chrMapFile,0,sizeof(chrMapFile));
    if (argv == 1)
    {
        printf("\nThe valid parameters are described as follows:\n");
        
        printf("\t-inputLabel: \n\t\t Default value: %lld. The integer flag label of the run.\n",inputFlag);
        printf("\t-outputFolder: \n\t\t The path of the folder to store the output fils. This folder must be exist!\n");
        printf("\t-cancerSample: \n\t\t Default value: %d. 1 for cancer; otherwise non cancer.\n",cMode);
        printf("\t-chrMapFile: \n\t\t The file name of the reference map(.cmap).\n");
        printf("\t-optAlignFile: \n\t\t The file name of the alignment map file(.oma).\n");
        printf("\t-inputRepeatFileName: \n\t\t The file containing highly repeated regions.\n");
        printf("\t-optTempFolder: \n\t\t The folder to store the processed alignment maps by chromosomes.\n");
        printf("\t-minClusterSizeRatio: \n\t\t Default value: %g. The minimum ratio of molecules forming a cluster.\n",minClsSizRat);
        printf("\t-minInsClusterSize: \n\t\t Default value: %d. The minimum number of molecules forming a insertion cluster.\n",minClsSI);
        printf("\t-minDelClusterSize: \n\t\t Default value: %d. The minimum number of molecules forming a deletion cluster.\n",minClsSD);
        printf("\t-minSVSupportRatio: \n\t\t Default value: %g. The minimum ratio of molecules supporting SVs.\n",minSVSuppRat);
        printf("\t-numberOfSupportIndelMolecule: \n\t\t Default value: %lld. The minimum coverage of a segment being called SVs. The default value changes along with the experiment data.\n",numberOfSupportIndelMolecule);
        printf("\t-minIndelSize: \n\t\t Default value (b): %g. The minimum length of a segment to call SVs.\n",minIndelSize);
        printf("\t-minSVSizeDiff: \n\t\t Default value (b): %g. The minimum difference of sizes of two SVs to be considered separately.\n",minSVDiffSize);
        printf("\t-minIndelRatio: \n\t\t Default value: %g. The length proportion of a minimum SV could be detected on a segment. E.g. segment = 10000b, then the length of the minimum SV should be larger than 10000*0.05=500b.\n",minIndelRatio);
        printf("\t-confidenceLimit: \n\t\t Default value: %g. The lowest alignment confidence for molecules (optical maps) to call SVs or signal variations.\n",confidenceLimit);
        printf("\t-chromosomeID: \n\t\t Default value: %d (-1 indicates an initialization of all temp files). The first n chromosomes to detect SVs.\n",chrID);
        printf("\t-maxClusterNumber: \n\t\t Default value: %d. The maximum number of clusters.\n",maxClsNum);
        return -1;
    }
    bool paraWrongFlag = false;
    for (int i = 1; i < argv; i=i+2)
    {
        string temp(argc[i]);
        //printf("What is the parameter?? %s \n",temp.c_str());
        if (temp.compare("-inputLabel")==0)
            inputFlag = atol(argc[i+1]);
        else if (temp.compare("-chromosomeID")==0){
            chrID = atoi(argc[i+1]);
        }
        else if (temp.compare("-maxClusterNumber")==0)
            maxClsNum = atoi(argc[i+1]);
        else if (temp.compare("-cancerSample")==0)
            cMode = atoi(argc[i+1]);
        else if (temp.compare("-minClusterSizeRatio")==0)
            minClsSizRat = atof(argc[i+1]);
        else if (temp.compare("-minInsClusterSize")==0)
            minClsSI = atoi(argc[i+1]);
        else if (temp.compare("-minDelClusterSize")==0)
            minClsSD = atoi(argc[i+1]);
        else if (temp.compare("-minSVSupportRatio")==0)
            minSVSuppRat = atof(argc[i+1]);
        else if (temp.compare("-optAlignFile")==0)
        {
            memset(inputAlignmentFileName,0,sizeof(inputAlignmentFileName));
            strcpy(inputAlignmentFileName, argc[i+1]);
        }
        else if (temp.compare("-inputRepeatFileName")==0)
        {
            memset(inputRepeatFileName,0,sizeof(inputRepeatFileName));
            strcpy(inputRepeatFileName, argc[i+1]);
        }
        else if (temp.compare("-optTempFolder")==0)
        {
            memset(outputFileLocation,0,sizeof(outputFileLocation));
            strcpy(outputFileLocation, argc[i+1]);
        }
        else if (temp.compare("-chrMapFile")==0)
        {
            memset(chrMapFile,0,sizeof(chrMapFile));
            strcpy(chrMapFile,argc[i+1]);
        }
        else if (temp.compare("-outputFolder")==0)
            strcpy(outputFolder,argc[i+1]);
        else if (temp.compare("-numberOfSupportIndelMolecule")==0)
        {
            numberOfSupportIndelMolecule = atoi(argc[i+1]);
            if (numberOfSupportIndelMolecule < 0)
            {
                printf("Error! An negative numberOfSupportIndelMolecule is inputted!\n");
                paraWrongFlag = true;
            }
            if (numberOfSupportIndelMolecule >= 1000)
            {
                printf("Warning! An improper large (>=1000) numberOfSupportIndelMolecule is inputted!\n");
            }
        }
        else if (temp.compare("-minIndelSize")==0)
        {
            minIndelSize = atof(argc[i+1]);
            if (minIndelSize < 1)
            {
                printf("Error! Please make sure the minIndelSize is a positive integer!\n");
                paraWrongFlag = true;
            }
            if (minIndelSize >= 100000)
            {
                printf("Warning! An improper large number (>=100000b) Of minIndelSize is inputted!\n");
            }
        }
        else if (temp.compare("-minSVSizeDiff")==0)
        {
            minSVDiffSize = atof(argc[i+1]);
            if (minSVDiffSize < 1)
            {
                printf("Error! Please make sure the minSVSizeDiff is a positive integer!\n");
                paraWrongFlag = true;
            }
            if (minSVDiffSize >= 100000)
            {
                printf("Warning! An improper large number (>=100000b) Of minIndelSize is inputted!\n");
            }
        }
        else if (temp.compare("-minIndelRatio")==0)
        {
            minIndelRatio = atof(argc[i+1]);
            if (minIndelRatio <= 0||minIndelRatio >= 1)
            {
                printf("Error! An improper number Of minIndelRatio is inputted, should be in between 0 and 1!\n");
                paraWrongFlag = true;
            }
        }
        else if (temp.compare("-confidenceLimit")==0)
        {
            confidenceLimit = atof(argc[i+1]);
            if (confidenceLimit < 0)
            {
                printf("Error! An negative confidenceLimit is inputted!\n");
                paraWrongFlag = true;
            }
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
    
    init();
    minCluster=min(minClsSI,minClsSD);
    highDens.clear();
    readHighDensity(inputRepeatFileName);
    getChromosomeList(chrMapFile);
    bool hitNot = false;
    for (LL tempChr = 0; tempChr < (LL)listOfChromosome.size(); tempChr++){
        LL curChr = listOfChromosome[tempChr];
        if (curChr == chrID) hitNot = true;
        char buff[200];
        char nameOF[1000];
        strcpy(nameOF, outputFileLocation);
        sprintf(buff, "%lld_%d",inputFlag,curChr);
        strcat(nameOF, buff);
        strcat(nameOF, ".bmap");
        if ((inputOptAlign = fopen(nameOF, "r")) == NULL)
        {
            continue;
        }
    }
    if ((!hitNot) && (chrID > 0)){
        printf("The input chromosomeID is wrong or not listed in the reference file, please check it again!\n");
        return -1;
    }
    if (chrID == -1){
        printf("Start the initialization for all chromosome!\n");
        readSourceFile();
        puts("DONE readSourceFile");
        addSplitedMap();
        puts("DONE addSplitedMap");
        outputToBillMapDestinationFile();
        printf("Initialization finished! Please change the input chromosomeID for normal processing!\n");
        return 0;
    }
    allSVs.clear();
    bool cM = (cMode==1);
    for (LL tempChrId = 0; tempChrId < (LL)listOfChromosome.size(); tempChrId++){
        chrId = listOfChromosome[tempChrId];
        canOpenFile = true;
        if (chrID > 0 && chrID != chrId)continue;
        printf("Start the chr%lld:\n",tempChrId+1);
        initData();
        printf("Doing %d, chr file: \n\t%s\n", chrId, chrMapFile);
        readChromosomeInfo(chrId,chrMapFile);
        printf("Done Read Chromosome\n");
        readOpticalAlign(chrId,outputFileLocation);
        if (!canOpenFile)
        {
            printf("Open File fails?\n");
            continue;
        }
        printf("Done Read Optical\n");
        printf("Distance Pair Count-1: %lld\n", distancePairCount);
        parseHitEnum();
        printf("Distance Pair Count-2: %lld\n", distancePairCount);
        detectByDistancePair(maxClsNum,0,cM);
        printf("Distance Pair Count-3: %lld\n", distancePairCount);
        printf("Now get %d SVs\n",(int)allSVs.size());
    }
    correctOverlapVariant();
    outputVariantResult(argv,argc,outputFolder,chrID);
    printf("Number of support SV: %lld\n", numberOfSupportedSV);
    printf("\n");
    return 0;
}

