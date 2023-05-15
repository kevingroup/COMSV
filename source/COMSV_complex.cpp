/*
  Created by Le Li, May 11, 2023
  Complex SV caller of COMSV pipeline
*/
#include "Header.h"

FILE *inputAlignmentFile;  
LL tempLongLong;
double tempDouble;
double scoreLimit = -1; 
vector <LL> listOfChromosome;
char tempString[10000];

struct	opticalMapType{
        char mapId[1000];
        LL refStartIndex;
        LL refEndIndex;
        LL optStartIndex;
        LL optEndIndex;
        LL refStart;
        LL refEnd;
        LL optStart;
        LL optEnd;
        LL chrId;
        bool orientation;
	double score;
        double confidence;
	string hitEnum;
	double fpr;////
        double fnr;////
        double alignRate;////
        vector<LL> position;
	bool operator<(const opticalMapType &w) const{
		LL e = strcmp(mapId, w.mapId);
		return ((e<0?true:false) || (e == 0 && min(optStart,optEnd) < min(w.optStart,w.optEnd)));
	}
        void print(){
        printf("%s %lf %lf %lf %lf %lf %lld %lld %lld %lld %lld %s %lld ", mapId, score, confidence,fpr,fnr,alignRate, chrId, optStart, optEnd, refStart, refEnd, hitEnum.c_str(), (LL)position.size());
        for (LL i=0; i<(LL)position.size(); i++)
            printf("%d ", position[i]);
        printf("\n");
    }

};

struct bedReg{
	LL chr;
	LL start;
	LL stop;
};
vector<bedReg> highDens;


void readHighDensity(char* inputAlignmentFileName){
        highDens.clear();
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
	while(fscanf(inputAlignmentFile, "%lld\t%lld\t%lld\n", &tpBed.chr,&tpBed.start,&tpBed.stop) == 3){
		highDens.push_back(tpBed);
	}
	fclose(inputAlignmentFile);
}


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
vector<chrType> genome;
vector<LL> chrList;


void getChromosomeList(char* chrMapFile){
    FILE* inputChrConsen;
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

chrType returnEmptyChr(){
    chrType chrom;
    delete[] chrom.position;
    delete[] chrom.distance;
    delete[] chrom.coverage;
    delete[] chrom.occurrence;
    delete[] chrom.gapCount;
    delete[] chrom.gapCoverage;
    delete[] chrom.gapSigni;
    delete[] chrom.signi;
    return chrom;
}

chrType readChromosomeInfo(int chr,char* chrMapFile){
    chrType chromosome;
    char nameOfFile[1000];
    FILE* inputChrConsen;
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
        chromosome.gapCoverage[cc] = 0;
        chromosome.occurrence[cc] = 0;
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
    return chromosome;
}




double normalizeMol(opticalMapType OM, chrType chromosome, int prt = 0){
    if (OM.position.size()<5)return 1;
    LL RefIndex = OM.refStartIndex;
    LL OpticalIndex = OM.optStartIndex;
    LL BeginRefIndex = RefIndex;
    LL BeginOpticalIndex = OpticalIndex;
    LL EndRefIndex = OM.refEndIndex;

    char previousChar = 'M';
    double tempDistance = 0;
    LL PreviousIndex = 0;
    LL PreviousOptIndex = 0;
    LL tempCount = 0;
    vector<double> ratioVecG1;
    ratioVecG1.clear();
    double optDist, refDist;
    LL inc = 0;
    if (OM.orientation) inc = 1;
    else inc = -1;
    for (LL j=0; j < OM.hitEnum.size(); j++){
        if (OM.hitEnum[j] >= '0' && OM.hitEnum[j] <= '9'){
            tempCount = tempCount * 10 + (OM.hitEnum[j]-'0');
        }
        else {
            if (OM.hitEnum[j] == 'M'){
                for (LL k=0; k<tempCount; k++){ //#if defined(Real_Assembly_Hg38) // fprintf(outputAssembly, "%s$
                    if (OpticalIndex != BeginOpticalIndex){
                        optDist = abs(OM.position[OpticalIndex] - OM.position[PreviousOptIndex]);
                        refDist = chromosome.position[RefIndex] - chromosome.position[PreviousIndex];
                        if (optDist>1000 && refDist>1000)ratioVecG1.push_back(optDist/refDist);
                        if (prt>0)printf("%g:%g__%g\t",optDist/refDist,optDist,refDist);
                    }
                    PreviousOptIndex = OpticalIndex;
                    OpticalIndex+=inc;
                    PreviousIndex = RefIndex;
                    RefIndex++;
                }
            } else if (OM.hitEnum[j] == 'D'){
                RefIndex += tempCount;
            } else if (OM.hitEnum[j] == 'I'){
                OpticalIndex += (tempCount*inc);
            }
            tempCount = 0;
        }
    }
    if (prt>0)printf("\n");
    sort(ratioVecG1.begin(),ratioVecG1.end());
    if (ratioVecG1.size()<5)return 1;
    else return ratioVecG1[(int)ratioVecG1.size()/2];
}


int countMs(string hitEnum){
	int Nd=0, Ni=0, Nm=0;
        int tempCount = 0;
	for (LL j = 0; j < hitEnum.size(); j++){
		if (hitEnum[j] >= '0' && hitEnum[j] <= '9')
        		tempCount = tempCount * 10 + (hitEnum[j]-'0');
		else{
        		if (hitEnum[j] == 'M'){
				Nm+=tempCount;
			}
			else if (hitEnum[j] == 'D')
				Nd+=tempCount;
			else
				Ni+=tempCount;
			tempCount = 0;
		}
	}
        return Nm;
}


vector<opticalMapType> readSourceFile(char* inputAlignmentFileName, double minSegLen, int minSites){
	vector<opticalMapType> opticalMap1;
        if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
        int numberOfDL = 0;
        char ttt;
        char ttts[100000];
        while(fgetc(inputAlignmentFile) == '#'){
                numberOfDL++;
                fgets(ttts,100000,inputAlignmentFile);
        }
        fclose(inputAlignmentFile);
        if ((inputAlignmentFile = fopen(inputAlignmentFileName, "r")) == NULL) puts("ERROR IN READ SOURCE");
        opticalMap1.clear();
	opticalMap1.resize(5000000);
        for (LL i=0; i<numberOfDL; i++)
                fgets(ttts, 50000, inputAlignmentFile);
        LL numberOfOpticalMap = 0;
        LL cntRm = 0;
        while (fscanf(inputAlignmentFile, "%s", opticalMap1[numberOfOpticalMap].mapId) == 1){
                if (numberOfOpticalMap % 100000 == 0) printf("%lld %s\n", numberOfOpticalMap, opticalMap1[numberOfOpticalMap].mapId);
                LL tempNumberOfSites, tempLL;
                fscanf(inputAlignmentFile, "%lld", &tempNumberOfSites);
                LL tempIndex = 0;
		opticalMap1[numberOfOpticalMap].position.clear();
                for (LL i=0; i<tempNumberOfSites; i++){
                        if (i != tempNumberOfSites - 1)
                                fscanf(inputAlignmentFile, "%lld;", &tempLL);
                        else fscanf(inputAlignmentFile, "%lld", &tempLL);
                        tempIndex += tempLL;
                        opticalMap1[numberOfOpticalMap].position.push_back(tempIndex);
                }
                fscanf(inputAlignmentFile, "%s", ttts);
                if (ttts[0] == 'c'){
                        if (ttts[3] == 'X')
                                opticalMap1[numberOfOpticalMap].chrId = 23;
                        else if (ttts[3] == 'Y')
                                opticalMap1[numberOfOpticalMap].chrId = 24;
                        else if (ttts[3] == 'M')
                                opticalMap1[numberOfOpticalMap].chrId = 25;
                        else
                                sscanf(ttts, "chr%lld", &opticalMap1[numberOfOpticalMap].chrId);
                }
                else{
                        if (ttts[0] == 'X')
                                opticalMap1[numberOfOpticalMap].chrId = 23;
                        else if (ttts[0] == 'Y')
                                opticalMap1[numberOfOpticalMap].chrId = 24;
                        else if (ttts[0] == 'M')
                                opticalMap1[numberOfOpticalMap].chrId = 25;
			else
                        	sscanf(ttts, "%lld", &opticalMap1[numberOfOpticalMap].chrId);
                }
                fscanf(inputAlignmentFile, "%s", ttts);
                if (ttts[0] == 'r' || ttts[0] == '-') opticalMap1[numberOfOpticalMap].orientation = false; else opticalMap1[numberOfOpticalMap].orientation = true;
                fscanf(inputAlignmentFile, "%lf", &opticalMap1[numberOfOpticalMap].score);
                fscanf(inputAlignmentFile, "%lf", &opticalMap1[numberOfOpticalMap].confidence);
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStartIndex, &opticalMap1[numberOfOpticalMap].refEndIndex);
                opticalMap1[numberOfOpticalMap].refStartIndex--;
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].optStart, &opticalMap1[numberOfOpticalMap].optEnd);
                fscanf(inputAlignmentFile, "%lld%lld", &opticalMap1[numberOfOpticalMap].refStart, &opticalMap1[numberOfOpticalMap].refEnd);
		char hitEnum[500000];
		memset(hitEnum,0,sizeof(hitEnum));
                fscanf(inputAlignmentFile, "%s", hitEnum);
		opticalMap1[numberOfOpticalMap].hitEnum = hitEnum;
                LL lenAlign;
                if (opticalMap1[numberOfOpticalMap].orientation)lenAlign = abs(opticalMap1[numberOfOpticalMap].position[opticalMap1[numberOfOpticalMap].optStart-1] - opticalMap1[numberOfOpticalMap].position[opticalMap1[numberOfOpticalMap].optEnd]);
                else lenAlign = abs(opticalMap1[numberOfOpticalMap].position[opticalMap1[numberOfOpticalMap].optStart] - opticalMap1[numberOfOpticalMap].position[opticalMap1[numberOfOpticalMap].optEnd-1]);
                if (opticalMap1[numberOfOpticalMap].score >= scoreLimit && lenAlign>minSegLen && countMs(opticalMap1[numberOfOpticalMap].hitEnum)>=minSites)
                    numberOfOpticalMap++;
                else
                    cntRm++;
        }
	opticalMap1.resize(numberOfOpticalMap);
	opticalMap1.shrink_to_fit();
        printf("Number of optical map: %lld, and filtered %lld with low score\n", (LL)opticalMap1.size(),cntRm);
        for (LL i=0; i<(LL)opticalMap1.size(); i++){
                if (opticalMap1[i].orientation){
			opticalMap1[i].optStartIndex = opticalMap1[i].optStart - 1;
			opticalMap1[i].optEndIndex = opticalMap1[i].optEnd;
                        opticalMap1[i].optStart = opticalMap1[i].position[opticalMap1[i].optStart - 1];
                        opticalMap1[i].optEnd = opticalMap1[i].position[opticalMap1[i].optEnd];
                }else {
			opticalMap1[i].optStartIndex = opticalMap1[i].optStart;
			opticalMap1[i].optEndIndex = opticalMap1[i].optEnd - 1;
                        opticalMap1[i].optStart = opticalMap1[i].position[opticalMap1[i].optStart];
                        opticalMap1[i].optEnd = opticalMap1[i].position[opticalMap1[i].optEnd - 1];
                }
        }
        double medianRatio;
        double ratMod = 2;
        for (LL i = 0; i < numberOfOpticalMap; i++){
            if (i==0 || i == numberOfOpticalMap - 1){
                medianRatio = normalizeMol(opticalMap1[i], genome[opticalMap1[i].chrId-1],1);
            }
            else medianRatio = normalizeMol(opticalMap1[i], genome[opticalMap1[i].chrId-1],0);
            if (medianRatio >= ratMod || medianRatio <= 1/ratMod){
                medianRatio = 1;
            }
            for (LL j = 0; j < opticalMap1[i].position.size(); j++){
                opticalMap1[i].position[j] /= medianRatio;
            }
        }
        printf("Done normalization!\n");
	return opticalMap1;
}

struct breakpoint{//actually is break point pair
	LL chr1; // if chr1 != chr2, let chr1 < chr2
	LL point1; // if chr1 == chr2, let point1 < point2
	LL point1_rag;
	LL point1_bd; // boundary of alignment of segment 1
	LL pointIdx1;
	LL chr2;
	LL point2;
	LL point2_rag;
	LL point2_bd; // boundary of alignment of segment 2
	LL pointIdx2;
	LL dupStart = 0;
	LL dupEnd = 0;
	LL dupStart_rag = 0;
	LL dupEnd_rag = 0;
	LL invStart;
	LL invEnd;
	LL gapSiz = 0;
	LL cnvSiz = 0;
	LL cnvNum = 0;
	char cnvPat[10000];
	double size = 0;
	char type[100];
	char extInfo[10000];
	char orient1;
	char orient2;
	vector<LL> supportMapId;
	vector<double> supportScores;
	vector<LL> nonSupportMapId;
	vector<LL> nonSupportMapId1;
	vector<LL> nonSupportMapId2;
	vector<double> nonSuppScores;
	vector<double> nonSuppScores1;	
	vector<double> nonSuppScores2;	
	bool operator<(const breakpoint &x) const{
		return (chr1 < x.chr1 || (chr1 == x.chr1 && point1 < x.point1) || (chr1 == x.chr1 && point1 == x.point1 && chr2 < x.chr2) || (chr1 == x.chr1 && point1 == x.point1 && chr2 == x.chr2 && point2 < x.point2));
	}
	double score = 0;
};

bool ss(breakpoint y, breakpoint x){
	return (y.chr1 < x.chr1 || (y.chr1 == x.chr1 && y.chr2 < x.chr2) || (y.chr1 == x.chr1 && y.chr2 == x.chr2 && (y.point1+y.point2) < (x.point2+x.point1)));
}
bool st(breakpoint y, breakpoint x){
	return (y.chr1 < x.chr1 || (y.chr1 == x.chr1 && y.chr2 < x.chr2) || (y.chr1 == x.chr1 && y.chr2 == x.chr2 && strcmp(y.type,x.type) < 0) || (y.chr1 == x.chr1 && y.chr2 == x.chr2 && strcmp(y.type,x.type) == 0  && (y.point1+y.point2) < (x.point2+x.point1)));
}
bool sm(breakpoint y, breakpoint x){
	LL mapComp = y.supportMapId[0] - x.supportMapId[0];
	return (mapComp<0 || (mapComp == 0 && y.chr1 < x.chr1) || (mapComp == 0 && y.chr1 == x.chr1 && y.chr2 < x.chr2) || (mapComp == 0 && y.chr1 == x.chr1 && y.chr2 == x.chr2 && strcmp(y.type,x.type) < 0) || (mapComp == 0 && y.chr1 == x.chr1 && y.chr2 == x.chr2 && strcmp(y.type,x.type) == 0  && (y.point1+y.point2) < (x.point2+x.point1)));
}


double calc_score(string hitEnum, int miniM){
	int tempCount = 0;
	double N=0, n=0;
	for (LL j = 0; j < hitEnum.size(); j++)
	{
		if (hitEnum[j] >= '0' && hitEnum[j] <= '9')
	        	tempCount = tempCount * 10 + (hitEnum[j]-'0');
	    	else{
        		if (hitEnum[j] == 'M')
				n+=tempCount;
			N+=tempCount;
			tempCount = 0;
		}
	}
	if (n<miniM) return 0;
	return log(n)*n/N;
}

void swapNums(LL &x, LL &y)
{
	auto temp = x;
	x = y;
	y = temp;
}
void swapNums(char &x, char &y)
{
	auto temp = x;
	x = y;
	y = temp;
}

bool checkConsist(vector<LL> vec, LL idx1, LL idx2, LL kmer){
	if (idx1 < 0 || idx2 < 0 || idx1 >= vec.size()- (kmer - 1) || idx2 >= vec.size()-(kmer-1)){
		printf("The index for forward is wrong, exit!\n");
		exit(0);
	}
	bool pos = true;
	for (LL i = 0; i < kmer; i++)
		pos = pos && (abs((vec[idx1+i]-vec[idx1+i-1])-(vec[idx2+i]-vec[idx2+i-1])) < min((vec[idx1+i]-vec[idx1+i-1]),(vec[idx2+i]-vec[idx2+i-1]))*0.05);
	return pos;
}

bool checkConsistR(vector<LL> vec, LL idx1, LL idx2, LL kmer){
	if (idx1 < 0 || idx2 < kmer || idx1 >= vec.size()-(kmer-1) || idx2 >= vec.size()){
		printf("The index for reverse is wrong, exit!\n");
		exit(0);
	}
	bool pos = true;
	for (LL i = 0; i < kmer; i++)
		pos = pos && (abs((vec[idx1+i]-vec[idx1+i-1])-(vec[idx2-i]-vec[idx2-i-1])) < min((vec[idx1+i]-vec[idx1+i-1]),(vec[idx2-i]-vec[idx2-i-1]))*0.05);
	return pos;
}
bool checkConsistRc(vector<LL> vec, LL idx1, LL idx2, LL kmer){
	if (idx1 < 0 || idx2 < kmer || idx1 >= vec.size()-(kmer-1) || idx2 >= vec.size()){
		printf("The index for reverse is wrong, exit!\n");
		exit(0);
	}
	bool pos = true;
	for (LL i = 0; i < kmer; i++)
		pos = pos && (abs((vec[idx1+i]-vec[idx1+i-1])-(vec[idx2-i]-vec[idx2-i-1])) < min((vec[idx1+i]-vec[idx1+i-1]),(vec[idx2-i]-vec[idx2-i-1]))*0.05);
	return pos;
}

bool checkConsistInv(vector<LL> vec, LL chrId, LL idx1, LL idx2, LL kmer){
        if (idx1 < 0 || idx2 < 0 || idx1 >= vec.size()- (kmer - 1) || idx2 >= genome[chrId-1].numberOfSites-kmer){
                printf("The inversion index is wrong, exit!\n");
                exit(0);
        }
        bool pos = true;
        for (LL i = 0; i < kmer; i++)
                pos = pos && (abs((vec[idx1+i]-vec[idx1+i-1])-(genome[chrId-1].position[idx2+i]-genome[chrId-1].position[idx2+i-1])) < min((vec[idx1+i]-vec[idx1+i-1]),(genome[chrId-1].position[idx2+i]-genome[chrId-1].position[idx2+i-1]))*0.1);
        return pos;
}

bool checkConsistInvR(vector<LL> vec1, LL chrId, LL idx1, LL idx2, LL kmer){
        if (idx1 < 0 || idx2 < kmer || idx1 >= vec1.size()-(kmer-1) || idx2 >= genome[chrId-1].numberOfSites){
                printf("The inversion reserve index is wrong, exit!\n");
                exit(0);
        }
        bool pos = true;
        for (LL i = 0; i < kmer; i++)
                pos = pos && (abs((vec1[idx1+i]-vec1[idx1+i-1])-(genome[chrId-1].position[idx2-i]-genome[chrId-1].position[idx2-i-1])) < min((vec1[idx1+i]-vec1[idx1+i-1]),(genome[chrId-1].position[idx2-i]-genome[chrId-1].position[idx2-i-1]))*0.1);
        return pos;
}
bool checkConsistInvRTest(vector<LL> vec1, LL chrId, LL idx1, LL idx2, LL kmer){
        if (idx1 < 0 || idx2 < kmer || idx1 >= vec1.size()-(kmer-1) || idx2 >= genome[chrId-1].numberOfSites){
                printf("The inversion reserve index is wrong, exit!\n");
                exit(0);
        }
        bool pos = true;
        for (LL i = 0; i < kmer; i++){
                pos = pos && (abs((vec1[idx1+i]-vec1[idx1+i-1])-(genome[chrId-1].position[idx2-i]-genome[chrId-1].position[idx2-i-1])) < min((vec1[idx1+i]-vec1[idx1+i-1]),(genome[chrId-1].position[idx2-i]-genome[chrId-1].position[idx2-i-1]))*0.1);
        }
        return pos;
}



void clarify_insert(opticalMapType OM1, opticalMapType OM2, bool& dup_flag, bool& CNV_flag, bool& inv_flag, bool& BP_flag, breakpoint& tempBP, LL stpOptIdx, int s_mode = 1){
	LL kmer=3;
	LL opt1St, opt1End, opt2St, opt2End;
	LL ref1St, ref1End, ref2St, ref2End, refGapSt, refGapEnd;
	LL chr1, chr2;
	tempBP.cnvSiz = 0;
	tempBP.cnvNum = 0;
	memset(tempBP.cnvPat,0,sizeof(tempBP.cnvPat));
	bool orient_sm1, orient_sm2, orient_gap;
	if (s_mode == 1){
		opt1St = min(OM1.optStartIndex,OM1.optEndIndex);
		opt1End = max(OM1.optStartIndex,OM1.optEndIndex);
		opt2St = min(OM2.optStartIndex,OM2.optEndIndex);
		opt2End = max(OM2.optStartIndex,OM2.optEndIndex);
		ref2St = OM2.refStart;
		ref2End = OM2.refEnd;
                orient_gap = OM1.orientation;
                if (orient_gap) {
                    refGapSt = OM1.refEndIndex+1;
                    refGapEnd = OM2.refStartIndex;
                }else{
                    refGapSt = OM2.refEndIndex+1;
                    refGapEnd = OM1.refStartIndex;
                }
		orient_sm1 = OM1.orientation;
		orient_sm2 = OM1.orientation;
		chr1 = OM1.chrId;
		chr2 = OM2.chrId;
	}else if (s_mode == 2){
		//the left split map is not aligned, only check the closing regions with <=50 site and <=500Kb
		opt2St = min(OM2.optStartIndex,OM2.optEndIndex);
		opt2End = max(OM2.optStartIndex,OM2.optEndIndex);
		for (opt1End = opt2St; opt1End>=stpOptIdx; opt1End--){
			if (OM1.position[opt2St] - OM1.position[opt1End] >= 500000 || opt2St - opt1End >= 50)break;
		}
		opt1St = opt1End-2;
		if (opt2St < kmer)return;
                orient_gap = OM2.orientation;
                if (orient_gap){
                    refGapEnd = OM2.refStartIndex;
                    refGapSt = max((LL)2,OM2.refStartIndex-opt2St-20);
                }else{
                    refGapSt = OM2.refEndIndex + 1;
                    refGapEnd = min(genome[OM2.chrId-1].numberOfSites-2, (LL)(opt2St + OM2.refEndIndex + 20));
                }
		ref2St = OM2.refStart;
		ref2End = OM2.refEnd;
		orient_sm1 = OM1.orientation;
		orient_sm2 = OM2.orientation;
		chr1 = OM1.chrId;
		chr2 = OM2.chrId;
	}else if (s_mode == 3){
		//the right split map is not aligned
		opt1St = min(OM1.optStartIndex,OM1.optEndIndex);
		opt1End = max(OM1.optStartIndex,OM1.optEndIndex);

		for (opt2St = opt1End; opt2St <= stpOptIdx; opt2St++){
			if (OM1.position[opt2St] - OM1.position[opt1End] >= 500000 || opt2St - opt1End >= 50)break;
		}
		opt2End = opt2St+2;
		if (OM1.position.size() - opt1End < kmer)return;

                orient_gap = OM1.orientation;
                if (orient_gap){
                    refGapSt = OM1.refEndIndex + 1;
                    refGapEnd = min(genome[OM1.chrId-1].numberOfSites-2, (LL)(OM1.position.size() - opt1End + OM1.refEndIndex + 20));
                }else{
                    refGapEnd = OM1.refStartIndex;
                    refGapSt = max((LL)2, (LL)(OM1.refStartIndex - (OM1.position.size() - opt1End) - 20));
                }
  
		ref2St = OM1.refStart;
		ref2End = OM1.refEnd;
		orient_sm1 = OM1.orientation;
		orient_sm2 = OM2.orientation;
		chr1 = OM1.chrId;
		chr2 = OM2.chrId;
	}
	LL gapNum = opt2St - opt1End;
	LL gapNum_ref = refGapEnd - refGapSt;
	LL CNV_st, CNV_ed;
	if (gapNum < kmer){
		CNV_flag = false;
		inv_flag = false;
	}else {
		//*Here to detect CNVs (self-repeat pattern) from the insert segment*//
		//consider the continuous cases
		vector<int> CNV_cnt(gapNum,0);

		for (LL l = 1; l < gapNum/2; l++){
			for (LL s = opt1End+1; s <= opt2St-l+1; s++){
				for (LL t = opt1End+1; t <= opt2St; t++) CNV_cnt[t-opt1End-1] = 0;
				for (LL t = opt1End+1; t <= opt2St-l+1;){
					if (s>t-l&&s<t+l){
						t++;
						continue;
					}
					bool ident_cond = true;
					for (LL lx = 0; lx < l; lx++)
						ident_cond = ident_cond && (abs((OM1.position[s+lx] - OM1.position[s+lx-1]) - (OM1.position[t+lx] - OM1.position[t+lx-1])) < min(min((OM1.position[s+lx] - OM1.position[s+lx-1]), (OM1.position[t+lx] - OM1.position[t+lx-1]))*0.05,2000.0) );
					if (ident_cond){
						for (LL lx = 0; lx < l; lx++)CNV_cnt[t-opt1End-1+lx] = 1;
						t+=l;
					}else{
						t++;
					}
				}
				if ((accumulate(CNV_cnt.begin(),CNV_cnt.end(),0)+l)>CNV_cnt.size()*0.8){
					CNV_flag = true;
					LL CNV_num = accumulate(CNV_cnt.begin(),CNV_cnt.end(),0)/l+1;
						for (LL t = 0; t < l; t++){
							sprintf(tempBP.cnvPat,"%lld;",OM1.position[t+s]-OM1.position[t+s-1]);
							tempBP.cnvSiz += (OM1.position[t+s]-OM1.position[t+s-1]);
						}
						tempBP.cnvNum = CNV_num;
					break;
				}
			}
			if (CNV_flag)break;
		}
		LL gapNum2 = gapNum;
		if (s_mode > 1){
			gapNum2 = max(gapNum, (LL)50); // scan the adjacent 50 sites for the unaligned maps
		}
		LL cnt = 0, numP=0, numN = 0;
		vector<int> occ(gapNum,0);
		vector<int> occR(gapNum,0);
		vector<int> occInv(gapNum,0);
		vector<int> occInvR(gapNum,0);
		vector<int> occ_Ref1(opt1End-max(opt1St+1,opt1End-gapNum2),0);
		vector<int> occ_Ref2(min(opt2St+gapNum2,opt2End)-opt2St-1,0);
                vector<int> occR_Ref1(opt1End-max(opt1St+1,opt1End-gapNum2),0); 
                vector<int> occR_Ref2(min(opt2St+gapNum2,opt2End)-opt2St-1,0);

		int lb = 0;
		for (LL s = opt1End+1; s <= opt2St-(kmer-1); s++){
			bool ocpP = false;
			cnt++;
			for (LL t = opt1End-(kmer-1); t >= max(opt1St+1,opt1End-gapNum2)+1; t--){
                                bool tp_bl = checkConsist(OM1.position,s,t,kmer);
				if (tp_bl){
					ocpP = true;
					for (LL x = 0; x < kmer; x++)occ[s-opt1End-1+x] = 1;
					for (LL x = 0; x < kmer; x++)occ_Ref1[t-max(opt1St+1,opt1End-gapNum2)-1+x] = 1;
					break;
				}
			}
			if (!ocpP){
				for (LL t = opt2St; t < min(opt2St+gapNum2,opt2End)-(kmer-1)-1; t++){
                                        bool tp_bl = checkConsist(OM1.position,s,t,kmer);
		        		if (tp_bl){
						for (LL x = 0; x < kmer; x++)occ[s-opt1End-1+x] = 1;
						for (LL x = 0; x < kmer; x++)occ_Ref2[t-opt2St+x] = 1;
						ocpP = true;
						break;
					}							
				}
			}
			if (ocpP){
				numP++;
			}
                        ocpP = false;
                        for (LL t = opt1End; t >= max(opt1St+1,opt1End-gapNum2)+(kmer-1)+1; t--){
                                if (checkConsistRc(OM1.position,s,t,kmer)){
                                        for (LL x = 0; x < kmer; x++)occR[s-opt1End-1+x] = 1;
                                        for (LL x = 0; x > -kmer; x--)occR_Ref1[t-max(opt1St+1,opt1End-gapNum2)-1+x] = 1;
                                        ocpP = true;
                                        break;
                                }
                        }
                        if (!ocpP){
                                for (LL t = opt2St+(kmer-1); t < min(opt2St+gapNum2,opt2End)-1; t++){
                                        if (checkConsistRc(OM1.position,s,t,kmer)){
                                                for (LL x = 0; x < kmer; x++)occR[s-opt1End-1+x] = 1;
                                                for (LL x = 0; x > -kmer; x--)occR_Ref2[t-opt2St+x] = 1;
                                                ocpP = true;
                                                break;
                                        }
                                }
                        }
                        if (ocpP){
                                numN++;
                        }
                        lb++;

                }
                LL refStSeg = 3000000000, refEdSeg = -1;
                if (gapNum_ref >= kmer){
			for (LL s = opt1End+1; s <= opt2St-(kmer-1); s++){
        	                if (orient_gap){
					for (LL t = refGapSt+kmer-1; t <= refGapEnd; t++){
		                                if (strcmp(OM1.mapId,"99098")==0)checkConsistInvRTest(OM1.position,OM1.chrId,s,t,kmer);
						if (checkConsistInvR(OM1.position,OM1.chrId,s,t,kmer)){
							for (LL x = 0; x < kmer; x++)occInv[s-opt1End-1+x] = 1;
	                                                refStSeg = min(refStSeg,genome[OM1.chrId-1].position[t]);
        	                                        refEdSeg = max(refEdSeg,genome[OM1.chrId-1].position[t+kmer-1]);
							break;
    						}
    					}
					for (LL t = refGapSt; t <= refGapEnd-(kmer-1); t++){
						if (checkConsistInv(OM1.position,OM1.chrId,s,t,kmer)){
							for (LL x = 0; x < kmer; x++)occInvR[s-opt1End-1+x] = 1;
							break;
    						}
    					}
	                        }else{
					for (LL t = refGapSt; t <= refGapEnd-(kmer-1); t++){
						if (checkConsistInv(OM1.position,OM1.chrId,s,t,kmer)){
							for (LL x = 0; x < kmer; x++)occInv[s-opt1End-1+x] = 1;
                                                	refStSeg = min(refStSeg,genome[OM1.chrId-1].position[t]);
	                                                refEdSeg = max(refEdSeg,genome[OM1.chrId-1].position[t+kmer-1]);
							break;
    						}
    					}
					for (LL t = refGapSt+kmer-1; t <= refGapEnd; t++){
		                                if (strcmp(OM1.mapId,"99098")==0)checkConsistInvRTest(OM1.position,OM1.chrId,s,t,kmer);
						if (checkConsistInvR(OM1.position,OM1.chrId,s,t,kmer)){
							for (LL x = 0; x < kmer; x++)occInvR[s-opt1End-1+x] = 1;
							break;
    						}
    					}
	                        }
			}
                }
                //here remove the maps overlap with those aligned to the reference (i.e. occInvR) in the rescue
                vector<LL> invRseg;
                invRseg.clear();
                for (LL s = 0; s < occInvR.size(); s++){
                    if (occInvR[s]>0 && (s==0 || occInvR[s-1]==0)){
                        invRseg.push_back(s);
                    }
                    if (occInvR[s]>0 && (s==occInvR.size()-1 || occInvR[s+1]==0)){
                        invRseg.push_back(s);
                    }
                }
                for (LL s = 0; s < invRseg.size(); s+=2){
                    bool allSum = false;
                    bool estAll = false;
                    if (((invRseg[s]>0)&&(occ[invRseg[s]-1]>0))||((invRseg[s+1]<occ.size()-1) && (occ[invRseg[s+1]+1]>0))) estAll=true;
                    for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                        if (occ[t]==0){estAll=false;break;}
                    if (!estAll){
                        for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                            occ[t]=0;
                    }
                    allSum |= estAll;
                    estAll = false;
                    if (((invRseg[s]>0)&&(occR[invRseg[s]-1]>0))||((invRseg[s+1]<occR.size()-1) && (occR[invRseg[s+1]+1]>0))) estAll=true;
                    for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                        if (occR[t]==0){estAll=false;break;}
                    if (!estAll){
                        for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                            occR[t]=0;
                    }
                    allSum |= estAll;
                    estAll = false;
                    if (((invRseg[s]>0)&&(occInv[invRseg[s]-1]>0))||((invRseg[s+1]<occInv.size()-1) && (occInv[invRseg[s+1]+1]>0))) estAll=true;
                    for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                        if (occInv[t]==0){estAll=false;break;}
                    if (!estAll){
                        for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                            occInv[t]=0;
                    }
                    allSum |= estAll;
                    if (allSum){
                    //if the segment of reference rescue is fully covered by another type of alignment, then we need to remove it
                        for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                            occInvR[t]=0;
                    }
                }

                //here to remove the small fragment of alignments (<kmer) after removing the parts overlapped with occInvR
                invRseg.clear();
                for (LL s = 0; s < occ.size(); s++){
                    if (occ[s]>0 && (s==0 || occ[s-1]==0)){
                        invRseg.push_back(s);
                    }
                    if (occ[s]>0 && (s==occ.size()-1 || occ[s+1]==0)){
                        invRseg.push_back(s);
                    }
                }
		for (LL s = 0; s < invRseg.size(); s+=2){
		    if (invRseg[s+1]-invRseg[s]<kmer-1){
                        for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                            occ[t]=0;
                    }
                }
                invRseg.clear();
                for (LL s = 0; s < occR.size(); s++){
                    if (occR[s]>0 && (s==0 || occR[s-1]==0)){
                        invRseg.push_back(s);
                    }
                    if (occR[s]>0 && (s==occR.size()-1 || occR[s+1]==0)){
                        invRseg.push_back(s);
                    }
                }
		for (LL s = 0; s < invRseg.size(); s+=2){
		    if (invRseg[s+1]-invRseg[s]<kmer-1){
                        for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                            occR[t]=0;
                    }
                }
                invRseg.clear();
                for (LL s = 0; s < occInv.size(); s++){
                    if (occInv[s]>0 && (s==0 || occInv[s-1]==0)){
                        invRseg.push_back(s);
                    }
                    if (occInv[s]>0 && (s==occInv.size()-1 || occInv[s+1]==0)){
                        invRseg.push_back(s);
                    }
                }
		for (LL s = 0; s < invRseg.size(); s+=2){
		    if (invRseg[s+1]-invRseg[s]<kmer-1){
                        for (LL t = invRseg[s]; t<= invRseg[s+1]; t++)
                            occInv[t]=0;
                    }
                }
                

		LL occ_dis = 0, occR_dis = 0, occInv_dis = 0, occInvR_dis = 0;
		for (LL s = 0; s < occ.size(); s++){
			if (occ[s]>0)
				occ_dis += (OM1.position[s+opt1End+1] - OM1.position[s+opt1End]);
			if (occR[s]>0)
				occR_dis += (OM1.position[s+opt1End+1] - OM1.position[s+opt1End]);
			if (occInv[s]>0)
				occInv_dis += (OM1.position[s+opt1End+1] - OM1.position[s+opt1End]);
			if (occInvR[s]>0)
				occInvR_dis += (OM1.position[s+opt1End+1] - OM1.position[s+opt1End]);
		}
		LL tp1_chr = chr1;
		LL tp1_dupS = 0;
		LL tp1_dupE = 0;
		LL tp2_chr = chr2;
		LL tp2_dupS = 0;
		LL tp2_dupE = 0;

		//here to fill small gaps and remove small noises between duplicated regions
		vector<int> oneMark;
		bool updat;
		int sel_Reg, max_Reg;

		for (int typ = 0; typ < 4; typ++){
                        vector<int> occRef;
			if (typ == 0)occRef = occ_Ref1;
			else if (typ == 1)occRef = occR_Ref1;
			else if (typ == 2)occRef = occ_Ref2;
			else if (typ == 3)occRef = occR_Ref2;
			updat = true;
			while(updat){
				updat = false;
				oneMark.clear();
				for (int t = 0; t < occRef.size(); t++){
					if (occRef[t] > 0){
						if (t == 0 || occRef[t-1] == 0){
							oneMark.push_back(t);
							oneMark.push_back(t);
						}else
							oneMark[oneMark.size()-1] = t;
					}
				}
				int st1, ed1, st2, ed2;
				if (oneMark.size()<4)continue;
				st1 = oneMark[0];
				ed1 = oneMark[1];
				for (int u = 2; u < oneMark.size(); u+=2){
					//if the gap size is less than the minimal flanking dups, then fill the gap, 
					st2 = oneMark[u];
					ed2 = oneMark[u+1];
					if (min(ed1-st1+1, ed2-st2+1) >= (st2-ed1-1)){
						for (int v = ed1+1; v < st2; v++)occRef[v] = 1;
						ed1 = ed2;
						updat = true;
					}
				}
			}
			max_Reg = 0;
			if (typ < 2){
				for (int t = oneMark.size()-2; t >=0; t-=2){
					if (oneMark[t+1] - oneMark[t] + 1 > max_Reg){
						sel_Reg = t;
						max_Reg = oneMark[t+1] - oneMark[t] + 1;
					}
				}
			}else{
				for (int t = 0; t < oneMark.size(); t+=2){
					if (oneMark[t+1] - oneMark[t] + 1 > max_Reg){
						sel_Reg = t;
						max_Reg = oneMark[t+1] - oneMark[t] + 1;
					}
				}
			}
			if (max_Reg > 0){
				for (int t = 0; t < oneMark[sel_Reg]; t++)occRef[t]=0;
				for (int t = oneMark[sel_Reg+1]+1; t < occRef.size(); t++)occRef[t]=0;
				// if the duplicated patterns is too far away from the split point, we remove them (range < distance)
				if (typ < 2){
					if ((oneMark[sel_Reg+1] - oneMark[sel_Reg] + 1)*1.5 < (occRef.size() - oneMark[sel_Reg+1])){
						for (int t = oneMark[sel_Reg]; t <= oneMark[sel_Reg+1]; t++)occRef[t]=0;
					}
				}else {
					if ((oneMark[sel_Reg+1] - oneMark[sel_Reg] + 1)*1.5 < (oneMark[sel_Reg]) && (oneMark[sel_Reg+1] - oneMark[sel_Reg]<=5)){
						for (int t = oneMark[sel_Reg]; t <= oneMark[sel_Reg+1]; t++)occRef[t]=0;
					}
				}
			}
			if (typ == 0)occ_Ref1 = occRef;
			else if (typ == 1)occR_Ref1 = occRef;
			else if (typ == 2)occ_Ref2 = occRef;
			else if (typ == 3)occR_Ref2 = occRef;
		}

		vector<int> cnvCnt(11,0);
                if (OM1.position[opt2St] - OM1.position[opt1End] - occInvR_dis==0)cout << OM1.mapId << "'s unmapped part can be fully aligned to the reference!\n";
                else if (OM1.position[opt2St] - OM1.position[opt1End] - occInvR_dis<0)cout << OM1.mapId << " has a wrong rescue result!\n";
		if (occ_dis >= occR_dis*0.8 && occ_dis > (OM1.position[opt2St] - OM1.position[opt1End] - occInvR_dis)*0.3 && occ_dis > 5000){
			for (LL ix = 0; ix <= occ_Ref1.size()-1; ix++){
				if (occ_Ref1[ix]>0)cnvCnt[occ_Ref1[ix]]++;
				if (orient_sm1){
					if (occ_Ref1[ix]>0){
						if (tp1_dupS == 0)
							tp1_dupS = OM1.position[max(opt1St+1,opt1End-gapNum2)+1+ix-1] - OM1.position[opt1End] + OM1.refEnd;
						tp1_dupE = OM1.position[max(opt1St+1,opt1End-gapNum2)+1+ix] - OM1.position[opt1End] + OM1.refEnd;
					}
				}else{
					if (occ_Ref1[ix]>0){
						if (tp1_dupE == 0)
							tp1_dupE = -OM1.position[max(opt1St+1,opt1End-gapNum2)+1+ix-1] + OM1.position[opt1End] + OM1.refStart;
						tp1_dupS = -OM1.position[max(opt1St+1,opt1End-gapNum2)+1+ix] + OM1.position[opt1End] + OM1.refStart;
					}
				}
			}
			for (LL ix = 0; ix <= occ_Ref2.size()-1; ix++){
				if (occ_Ref2[ix]>0)cnvCnt[occ_Ref2[ix]]++;
				if (orient_sm2){
					if (occ_Ref2[ix]>0){								
						if (tp2_dupS == 0)
							tp2_dupS = OM1.position[opt2St+ix-1] - OM1.position[opt2St] + ref2St;
						tp2_dupE = OM1.position[opt2St+ix] - OM1.position[opt2St] + ref2St;
					}
				}else{
					if (occ_Ref2[ix]>0){
						if (tp2_dupE == 0)
							tp2_dupE = -OM1.position[opt2St+ix-1] + OM1.position[opt2St] + ref2End;
						tp2_dupS = -OM1.position[opt2St+ix] + OM1.position[opt2St] + ref2End;
					}
				}
			}
			if (tempBP.dupStart == 0){
				if (tp1_dupS == 0){
					tempBP.dupStart = tp2_dupS;
					tempBP.dupEnd = tp2_dupE;
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}else if (tp2_dupS == 0){
					tempBP.dupStart = tp1_dupS;
					tempBP.dupEnd = tp1_dupE;
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}else if ((max(tp1_dupS,tp2_dupS) - min(tp1_dupE,tp2_dupE)<10000)&&tp2_chr==tp1_chr){
					tempBP.dupStart = min(tp1_dupS,tp2_dupS);
					tempBP.dupEnd = max(tp1_dupE,tp2_dupE);
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}else if ((tp1_dupE - tp1_dupS)>(tp2_dupE-tp2_dupS)){
					tempBP.dupStart = tp1_dupS;
					tempBP.dupEnd = tp1_dupE;
					sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tp2_chr,tp2_dupS,tp2_dupE);
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}else{
					tempBP.dupStart = tp2_dupS;
					tempBP.dupEnd = tp2_dupE;
					sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tp1_chr,tp1_dupS,tp1_dupE);
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}
				dup_flag = true;
			}else{
				//merge the duplication segments if their distance < 10Kb
				if ((( min(tempBP.dupStart,tempBP.dupStart_rag) - tp1_dupE - 10000 ) * (tp1_dupS - max(tempBP.dupStart,tempBP.dupStart_rag) - 10000) >= 0)&&tempBP.chr1 == tp1_chr){
					if (tp1_dupS <= tempBP.dupStart){
						tempBP.dupStart = tp1_dupS;
						tempBP.dupStart_rag = tp1_dupS;
					}
					if (tp1_dupE >= tempBP.dupEnd){
						tempBP.dupEnd = tp1_dupE;
						tempBP.dupEnd_rag = tp1_dupE;
					}
					if (( min(tempBP.dupStart,tempBP.dupStart_rag) - tp2_dupE - 10000 ) * (tp2_dupS - max(tempBP.dupStart,tempBP.dupStart_rag) - 10000) >= 0){
						if (tp2_dupS <= tempBP.dupStart){
							tempBP.dupStart = tp2_dupS;
							tempBP.dupStart_rag = tp2_dupS;
						}
						if (tp2_dupE >= tempBP.dupEnd){
							tempBP.dupEnd = tp2_dupE;
							tempBP.dupEnd_rag = tp2_dupE;
						}
					}else if (tp2_dupS > 0){
						sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tp2_chr,tp2_dupS,tp2_dupE);
					}
				}else if ((( min(tempBP.dupStart,tempBP.dupStart_rag) - tp2_dupE - 10000 ) * (tp2_dupS - max(tempBP.dupStart,tempBP.dupStart_rag) - 10000) >= 0)&&tempBP.chr1 == tp2_chr){
					if (tp2_dupS <= tempBP.dupStart){
						tempBP.dupStart = tp2_dupS;
						tempBP.dupStart_rag = tp2_dupS;
					}
					if (tp2_dupE >= tempBP.dupEnd){
						tempBP.dupEnd = tp2_dupE;
						tempBP.dupEnd_rag = tp2_dupE;
					}
					if (( min(tempBP.dupStart,tempBP.dupStart_rag) - tp1_dupE - 10000 ) * (tp1_dupS - max(tempBP.dupStart,tempBP.dupStart_rag) - 10000) >= 0){
						if (tp1_dupS <= tempBP.dupStart){
							tempBP.dupStart = tp1_dupS;
							tempBP.dupStart_rag = tp1_dupS;
						}
						if (tp1_dupE >= tempBP.dupEnd){
							tempBP.dupEnd = tp1_dupE;
							tempBP.dupEnd_rag = tp1_dupE;
						}
					}else if (tp1_dupS > 0){
						sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tp1_chr,tp1_dupS,tp1_dupE);
					}					
				}else{
					if (tp1_dupS != 0 && tp2_dupS != 0 && (( tp2_dupS - tp1_dupE - 10000 ) * (tp1_dupS - tp2_dupE - 10000) >= 0)){
						sprintf(tempBP.extInfo,"ExtraDupCoord: %lld-%lld",min(tp1_dupS,tp2_dupS),max(tp1_dupE,tp2_dupE));
					}else{
						if (tp1_dupS!=0 && tp2_dupS!=0)
							sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld; %lld:%lld-%lld",tp1_chr,tp1_dupS,tp1_dupE,tp2_chr,tp2_dupS,tp2_dupE);
						else if (tp1_dupS!=0)
							sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tp1_chr,tp1_dupS,tp1_dupE);
						else
							sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tp1_chr,tp2_dupS,tp2_dupE);
					}
				}
				dup_flag = true;
			}
		}
                else if (occR_dis * 0.8 > occ_dis && occR_dis > (OM1.position[opt2St] - OM1.position[opt1End] - occInvR_dis)*0.3 && occR_dis > 5000){
			LL tp1_dupS = 0;
			LL tp1_dupE = 0;
			LL tp2_dupS = 0;
			LL tp2_dupE = 0;					
			for (LL ix = 0; ix <= occR_Ref1.size()-1; ix++){
				if (orient_sm1){
					if (occR_Ref1[ix]>0){								
						if (tp1_dupS == 0)
							tp1_dupS = OM1.position[max(opt1St+1,opt1End-gapNum2)+1+ix-1] - OM1.position[opt1End] + OM1.refEnd;
						tp1_dupE = OM1.position[max(opt1St+1,opt1End-gapNum2)+ix+1] - OM1.position[opt1End] + OM1.refEnd;
					}
				}else{
					if (occR_Ref1[ix]>0){
						if (tp1_dupE == 0)
							tp1_dupE = -OM1.position[max(opt1St+1,opt1End-gapNum2)+1+ix-1] + OM1.position[opt1End] + OM1.refStart;
						tp1_dupS = -OM1.position[max(opt1St+1,opt1End-gapNum2)+ix+1] + OM1.position[opt1End] + OM1.refStart;
					}
				}
			}
			for (LL ix = 0; ix <= occR_Ref2.size()-1; ix++){
				if (orient_sm2){
					if (occR_Ref2[ix]>0){								
						if (tp2_dupS == 0)
							tp2_dupS = OM1.position[opt2St+1+ix-1] - OM1.position[opt2St] + ref2St;
						tp2_dupE = OM1.position[opt2St+ix] - OM1.position[opt2St] + ref2St;
					}
				}else{
					if (occR_Ref2[ix]>0){
						if (tp2_dupE == 0)
							tp2_dupE = -OM1.position[opt2St+1+ix-1] + OM1.position[opt2St] + ref2End;
						tp2_dupS = -OM1.position[opt2St+1+ix] + OM1.position[opt2St] + ref2End;
					}
				}
			}
			if (tempBP.dupStart == 0){
				if (tp1_dupS == 0){
					tempBP.dupStart = tp2_dupS;
					tempBP.dupEnd = tp2_dupE;
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}else if (tp2_dupS == 0){
					tempBP.dupStart = tp1_dupS;
					tempBP.dupEnd = tp1_dupE;
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}else if ((max(tp1_dupS,tp2_dupS) - min(tp1_dupE,tp2_dupE)<10000)&&tp1_chr==tp2_chr){
					tempBP.dupStart = min(tp1_dupS,tp2_dupS);
					tempBP.dupEnd = max(tp1_dupE,tp2_dupE);
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}else if ((tp1_dupE - tp1_dupS)>(tp2_dupE-tp2_dupS)){
					tempBP.dupStart = tp1_dupS;
					tempBP.dupEnd = tp1_dupE;
					sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tp2_chr,tp2_dupS,tp2_dupE);
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}else{
					tempBP.dupStart = tp2_dupS;
					tempBP.dupEnd = tp2_dupE;
					sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tp2_chr,tp1_dupS,tp1_dupE);
					tempBP.dupStart_rag = tempBP.dupStart;//here can improve 
					tempBP.dupEnd_rag = tempBP.dupEnd;
				}
				dup_flag = true;
			}else{
				//merge the duplication segments if their distance < 10Kb
				if ((tp1_dupS!=0 && tp2_dupS!=0 && (( tp2_dupS - tp1_dupE - 10000 ) * (tp1_dupS - tp2_dupE - 10000) >= 0))&&(tp1_chr==tp2_chr)){
					sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld",tempBP.chr1,tempBP.dupStart,tempBP.dupEnd);
					tempBP.dupStart = min(tp1_dupS,tp2_dupS);
					tempBP.dupStart_rag = min(tp1_dupS,tp2_dupS);
					tempBP.dupEnd = max(tp1_dupE,tp2_dupE);
					tempBP.dupEnd_rag = max(tp1_dupE,tp2_dupE);
				}else{
					if (tp1_dupS!=0 && tp2_dupS!=0){
						if (tp1_dupE-tp1_dupS>tp2_dupE-tp2_dupS){
							tempBP.dupStart = tp1_dupS;
							tempBP.dupStart_rag = tp1_dupS;
							tempBP.dupEnd = tp1_dupE;
							tempBP.dupEnd_rag = tp1_dupE;
							sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld; %lld:%lld-%lld",tempBP.chr1,tempBP.dupStart,tempBP.dupEnd,tp2_chr,tp2_dupS,tp2_dupE);
						}else{
							tempBP.dupStart = tp2_dupS;
							tempBP.dupStart_rag = tp2_dupS;
							tempBP.dupEnd = tp2_dupE;
							tempBP.dupEnd_rag = tp2_dupE;
							sprintf(tempBP.extInfo,"ExtraDupCoord: %lld:%lld-%lld; %lld:%lld-%lld",tempBP.chr1,tempBP.dupStart,tempBP.dupEnd,tp1_chr,tp1_dupS,tp1_dupE);
						}
					}else if (tp1_dupS!=0){
						tempBP.dupStart = tp1_dupS;
						tempBP.dupStart_rag = tp1_dupS;
						tempBP.dupEnd = tp1_dupE;
						tempBP.dupEnd_rag = tp1_dupE;
					}else{
						tempBP.dupStart = tp2_dupS;
						tempBP.dupStart_rag = tp2_dupS;
						tempBP.dupEnd = tp2_dupE;
						tempBP.dupEnd_rag = tp2_dupE;
					}
				}
				dup_flag = true;
			}
		}
                if (occInv_dis > (OM1.position[opt2St] - OM1.position[opt1End] - occInvR_dis)*0.3 && occInv_dis > 5000){
			sprintf(tempBP.extInfo,"SplitInvCoord: %lld:%lld-%lld(%lld/%lld)",OM1.chrId,refStSeg,refEdSeg,occInv_dis,OM1.position[opt2St] -OM1.position[opt1End]);
			inv_flag = true;
		}
	}
	if (dup_flag&&!CNV_flag){
		tempBP.cnvSiz = tempBP.dupEnd - tempBP.dupStart;
		tempBP.cnvNum = round(tempBP.size * 1.0 / tempBP.cnvSiz);
	}
	if ((OM1.refStart<=OM2.refStart && OM1.refEnd >= OM2.refEnd)||(OM2.refStart<=OM1.refStart && OM2.refEnd >= OM1.refEnd))BP_flag = true;

}



vector<breakpoint> scanSplitedMap(vector<opticalMapType> opticalMap1, int miniM){
        vector<breakpoint> candidateSet;
        candidateSet.clear();
	sort(opticalMap1.begin(), opticalMap1.end());
	LL grpSt = 0, grpEnd;
	printf("There are %lld alignments being read\n",(LL)opticalMap1.size());
	for (LL i=0; i<opticalMap1.size()-1; i++){
		if (strcmp(opticalMap1[i].mapId, opticalMap1[i+1].mapId) != 0){
			grpEnd = i;
			if (grpEnd - grpSt > 0){
				for (LL k = grpSt; k <= grpEnd; k++){
					for (LL l = grpSt; l <= grpEnd; l++){
						if (k == l || strcmp(opticalMap1[l].mapId,"0")==0)continue;
						if (min(opticalMap1[k].optStartIndex,opticalMap1[k].optEndIndex)>=min(opticalMap1[l].optStartIndex,opticalMap1[l].optEndIndex) && max(opticalMap1[k].optStartIndex,opticalMap1[k].optEndIndex)<=max(opticalMap1[l].optStartIndex,opticalMap1[l].optEndIndex)){
							strcpy(opticalMap1[k].mapId,"0");
							break;
						}
					}
				}
				for (LL k = grpSt; k <= grpEnd; k++){
					if (strcmp(opticalMap1[k].mapId,"0")==0)continue;
					vector<int> indUniq(opticalMap1[k].position.size(),0);
					for (int xt = min(opticalMap1[k].optStartIndex,opticalMap1[k].optEndIndex); xt <= max(opticalMap1[k].optStartIndex,opticalMap1[k].optEndIndex); xt++)indUniq[xt] = 1;
					for (LL l = grpSt; l <= grpEnd; l++){
						if (l == k || strcmp(opticalMap1[l].mapId,"0")==0)continue;
						for (int xt = min(opticalMap1[l].optStartIndex,opticalMap1[l].optEndIndex); xt <= max(opticalMap1[l].optStartIndex,opticalMap1[l].optEndIndex); xt++)indUniq[xt] = 0;
					}
					if (accumulate(indUniq.begin(), indUniq.end(), 0) < 3 ){
						strcpy(opticalMap1[k].mapId,"0");
					}
				}
			}
			grpSt = i+1;
		}else{
			grpEnd = i+1;
			if (i == opticalMap1.size()-2){
				for (LL k = grpSt; k <= grpEnd; k++){
                                        vector<int> indUniq(opticalMap1[k].position.size(),0);
                                        for (int xt = min(opticalMap1[k].optStartIndex,opticalMap1[k].optEndIndex); xt <= max(opticalMap1[k].optStartIndex,opticalMap1[k].optEndIndex); xt++)indUniq[xt] = 1;
                                        for (LL l = grpSt; l <=grpEnd; l++){
                                                if (l == k)continue;
                                                for (int xt = min(opticalMap1[l].optStartIndex,opticalMap1[l].optEndIndex); xt <= max(opticalMap1[l].optStartIndex,opticalMap1[l].optEndIndex); xt++)indUniq[xt] = 0;
                                        }
                                        if (accumulate(indUniq.begin(), indUniq.end(), 0) < 3 ){
                                                strcpy(opticalMap1[k].mapId,"0");
                                        }
                                }
			}
		}
	}
	for (LL i=0; i<opticalMap1.size()-1; i++){
		// here if the alignment is overlapped with high-density region with > 50%, then remove it
		for (LL k = 0 ; k < highDens.size(); k++){
			LL hdSt = highDens[k].start;
			LL hdEnd = highDens[k].stop;
			if (highDens[k].chr == opticalMap1[i].chrId){
				if ((min(opticalMap1[i].refEnd,hdEnd) - max(opticalMap1[i].refStart,hdSt)) > (opticalMap1[i].refEnd - opticalMap1[i].refStart)/2 ){
					strcpy(opticalMap1[i].mapId,"0");
				}
			}
		}
	}
	sort(opticalMap1.begin(), opticalMap1.end());
	for (LL i = 0; i < opticalMap1.size(); i++){
		if (strcmp(opticalMap1[i].mapId,"0")!=0){
			if (i > 0){
				opticalMap1.erase(opticalMap1.begin(),opticalMap1.begin()+i);
			}
			break;
		}
	}


	vector<breakpoint> splitDup;
	splitDup.clear();
	for (LL i=0; i<opticalMap1.size()-1; i++){
		if (i==0 || strcmp(opticalMap1[i].mapId, opticalMap1[i-1].mapId) != 0){
			breakpoint tempBP;
			strcpy(tempBP.extInfo,"None");
			tempBP.supportMapId.clear();
			tempBP.supportScores.clear();
			tempBP.nonSupportMapId1.clear();
			tempBP.nonSupportMapId2.clear();
			tempBP.nonSupportMapId.clear();
			tempBP.nonSuppScores1.clear();
			tempBP.nonSuppScores2.clear();
			tempBP.nonSuppScores.clear();
			tempBP.chr1 = opticalMap1[i].chrId;
			tempBP.chr2 = opticalMap1[i].chrId;

			if (opticalMap1[i].orientation){
				tempBP.point1 = opticalMap1[i].refStart;
				tempBP.point1_rag = opticalMap1[i].refStart;
				tempBP.point1_bd = opticalMap1[i].refEnd;
				tempBP.pointIdx1 = opticalMap1[i].refStartIndex;
                	} else {
	                        tempBP.point1 = opticalMap1[i].refEnd;
        	                tempBP.point1_rag = opticalMap1[i].refEnd;
                	        tempBP.point1_bd = opticalMap1[i].refStart;
	                        tempBP.pointIdx1 = opticalMap1[i].refEndIndex;
        	        }
			tempBP.orient1 = opticalMap1[i].orientation?'+':'-';
			tempBP.orient2 = opticalMap1[i].orientation?'+':'-';
			tempBP.point2 = tempBP.point1;
			tempBP.point2_rag = tempBP.point1_rag;
			tempBP.point2_bd = tempBP.point1_bd;
			tempBP.pointIdx2 = tempBP.pointIdx1;
			if (calc_score(opticalMap1[i].hitEnum,miniM) > 0)
			{
				tempBP.supportMapId.push_back(atol(opticalMap1[i].mapId));
				tempBP.supportScores.push_back(calc_score(opticalMap1[i].hitEnum,miniM));
			}

			bool dup_flag = false;
			bool CNV_flag = false;
			bool inv_flag = false;
			bool tand_flag = true;
			bool BP_flag = false;

			clarify_insert(opticalMap1[i], opticalMap1[i], dup_flag, CNV_flag, inv_flag, BP_flag, tempBP, 1, 2);
			//***Here define the tandem as the duplications in adjacent regions: gap size between split map < 100Kb, the distance between duplicated regions and breakpoints < 100Kb***//
                        if (dup_flag)strcpy(tempBP.type,"Duplication");
                        else if (inv_flag) strcpy(tempBP.type,"Inversion");
                        else if (CNV_flag) strcpy(tempBP.type,"CNV");
			strcat(tempBP.type,"-Split");
			if (dup_flag || CNV_flag || inv_flag)candidateSet.push_back(tempBP);
		}
		if (strcmp(opticalMap1[i].mapId, opticalMap1[i+1].mapId) != 0){
			breakpoint tempBP;
			strcpy(tempBP.extInfo,"None");
			tempBP.supportMapId.clear();
			tempBP.supportScores.clear();
			tempBP.nonSupportMapId1.clear();
			tempBP.nonSupportMapId2.clear();
			tempBP.nonSupportMapId.clear();
			tempBP.nonSuppScores1.clear();
			tempBP.nonSuppScores2.clear();
			tempBP.nonSuppScores.clear();
			tempBP.chr1 = opticalMap1[i].chrId;
			tempBP.chr2 = opticalMap1[i].chrId;

			if (opticalMap1[i].orientation){
				tempBP.point1 = opticalMap1[i].refEnd;
				tempBP.point1_rag = opticalMap1[i].refEnd;
				tempBP.point1_bd = opticalMap1[i].refStart;
				tempBP.pointIdx1 = opticalMap1[i].refEndIndex;
        	        } else {
                	        tempBP.point1 = opticalMap1[i].refStart;
                        	tempBP.point1_rag = opticalMap1[i].refStart;
	                        tempBP.point1_bd = opticalMap1[i].refEnd;
        	                tempBP.pointIdx1 = opticalMap1[i].refStartIndex;
                	}
			tempBP.orient1 = opticalMap1[i].orientation?'+':'-';
			tempBP.orient2 = opticalMap1[i].orientation?'+':'-';
			tempBP.point2 = tempBP.point1;
			tempBP.point2_rag = tempBP.point1_rag;
			tempBP.point2_bd = tempBP.point1_bd;
			tempBP.pointIdx2 = tempBP.pointIdx1;
			if (calc_score(opticalMap1[i].hitEnum,miniM) > 0)
			{
				tempBP.supportMapId.push_back(atol(opticalMap1[i].mapId));
				tempBP.supportScores.push_back(calc_score(opticalMap1[i].hitEnum,miniM));
			}

			bool dup_flag = false;
			bool CNV_flag = false;
			bool inv_flag = false;
			bool tand_flag = true;
			bool BP_flag = false;
			clarify_insert(opticalMap1[i], opticalMap1[i], dup_flag, CNV_flag, inv_flag,BP_flag, tempBP, (LL)opticalMap1[i].position.size()-3 , 3);
                        if (dup_flag)strcpy(tempBP.type,"Duplication");
                        else if (inv_flag) strcpy(tempBP.type,"Inversion");
                        else if (CNV_flag) strcpy(tempBP.type,"CNV");
			strcat(tempBP.type,"-Split");
			if (dup_flag || CNV_flag || inv_flag)candidateSet.push_back(tempBP);

			continue;
		}


		breakpoint tempBP;
		strcpy(tempBP.extInfo,"None");
		tempBP.supportMapId.clear();
		tempBP.supportScores.clear();
		tempBP.nonSupportMapId1.clear();
		tempBP.nonSupportMapId2.clear();
		tempBP.nonSupportMapId.clear();
		tempBP.nonSuppScores1.clear();
		tempBP.nonSuppScores2.clear();
		tempBP.nonSuppScores.clear();
		tempBP.chr1 = opticalMap1[i].chrId;
		tempBP.chr2 = opticalMap1[i+1].chrId;

		if (opticalMap1[i].orientation){
			tempBP.point1 = opticalMap1[i].refEnd;
			tempBP.point1_rag = opticalMap1[i].refEnd;
			tempBP.point1_bd = opticalMap1[i].refStart;
			tempBP.pointIdx1 = opticalMap1[i].refEndIndex;
                } else {
                        tempBP.point1 = opticalMap1[i].refStart;
                        tempBP.point1_rag = opticalMap1[i].refStart;
                        tempBP.point1_bd = opticalMap1[i].refEnd;
                        tempBP.pointIdx1 = opticalMap1[i].refStartIndex;
                }
		tempBP.orient1 = opticalMap1[i].orientation?'+':'-';
                if (opticalMap1[i+1].orientation){
                        tempBP.point2 = opticalMap1[i+1].refStart;
                        tempBP.point2_rag = opticalMap1[i+1].refStart;
                        tempBP.point2_bd = opticalMap1[i+1].refEnd;
                        tempBP.pointIdx2 = opticalMap1[i+1].refStartIndex;
                } else {
                        tempBP.point2 = opticalMap1[i+1].refEnd;
                        tempBP.point2_rag = opticalMap1[i+1].refEnd;
                        tempBP.point2_bd = opticalMap1[i+1].refStart;
                        tempBP.pointIdx2 = opticalMap1[i+1].refEndIndex;
                }
		tempBP.orient2 = opticalMap1[i+1].orientation?'+':'-';

		double ovlpDist = max(opticalMap1[i+1].refStart,opticalMap1[i].refStart) - min(opticalMap1[i+1].refEnd,opticalMap1[i].refEnd);
//////////////////////////here consider to find out the overlapping matching sites between two split maps
                double ovlpMol = max( min(opticalMap1[i+1].optStart,opticalMap1[i+1].optEnd),min(opticalMap1[i].optStart,opticalMap1[i].optEnd)) - min(max(opticalMap1[i+1].optStart,opticalMap1[i+1].optEnd),max(opticalMap1[i].optStart,opticalMap1[i].optEnd));
		tempBP.size = ovlpMol - ovlpDist;
		tempBP.gapSiz = max(ovlpMol,0.0);
                


		LL forSt = opticalMap1[i].refStart;
		LL forEnd = opticalMap1[i].refEnd;
		LL latSt = opticalMap1[i+1].refStart;
		LL latEnd = opticalMap1[i+1].refEnd;

		if (min(forEnd,latEnd)>max(forSt,latSt)+10000){
			tempBP.dupStart_rag = max(forSt,latSt);
			tempBP.dupEnd_rag = min(forEnd,latEnd);
		}else{
			tempBP.dupStart_rag = 0;
			tempBP.dupEnd_rag = 0;
		}

		

		double ovlpRef = -min(ovlpDist,0.0) + min(ovlpMol, 0.0);
		//always trust the longer segment, and reduce the shorter one
		//longer segment as the reference part, and shorter segment as variant
		if (opticalMap1[i+1].refEnd-opticalMap1[i+1].refStart > opticalMap1[i].refEnd-opticalMap1[i].refStart){
			if (tempBP.orient1 == '+'){
				tempBP.point1 += min(ovlpMol,0.0);
				forEnd += min(ovlpMol,0.0);
			}
			else{
				tempBP.point1 -= min(ovlpMol,0.0);
				forSt -= min(ovlpMol,0.0);
			}
			if (tempBP.orient2 == '+'){
                                tempBP.point2_rag -= min(ovlpMol,0.0);
			}
                        else{
                                tempBP.point2_rag += min(ovlpMol,0.0);
			}
		}else{
			if (tempBP.orient1 == '+'){
				tempBP.point1_rag += min(ovlpMol,0.0);
			}
			else{
				tempBP.point1_rag -= min(ovlpMol,0.0);
			}
			if (tempBP.orient2 == '+'){
                                tempBP.point2 -= min(ovlpMol,0.0);
				latSt -= min(ovlpMol,0.0);
			}
                        else{
                                tempBP.point2 += min(ovlpMol,0.0);
				latEnd += min(ovlpMol,0.0);
			}
		}


		
		if (min(forEnd,latEnd)>max(forSt,latSt)+10000){
			tempBP.dupStart = max(forSt,latSt);
			tempBP.dupEnd = min(forEnd,latEnd);
		}else{
			tempBP.dupStart = 0;
			tempBP.dupEnd = 0;
		}

		bool cond1 = opticalMap1[i].chrId != opticalMap1[i+1].chrId;
		bool cond2 = (tempBP.size <= -1000000); // for the events span >1Mb, we call them as intra-translocation

			// for duplications, we need consider to remove the overlapping distance of molecules between alignments
		bool cond3 = (tempBP.size >= 10000);
		bool cond31 = (ovlpRef > 10000);
		bool cond4 = (opticalMap1[i].orientation != opticalMap1[i+1].orientation);
		bool cond5 = (opticalMap1[i].orientation != (opticalMap1[i].refStart < opticalMap1[i+1].refStart));

		LL dupRegS = min(tempBP.dupStart,tempBP.dupStart_rag);
		LL dupRegE = max(tempBP.dupEnd,tempBP.dupEnd_rag);
		LL bp1RegS = min(tempBP.point1, tempBP.point1_rag);
		LL bp1RegE = max(tempBP.point1, tempBP.point1_rag);
		LL bp2RegS = min(tempBP.point2, tempBP.point2_rag);
		LL bp2RegE = max(tempBP.point2, tempBP.point2_rag);
		bool tand_flag = ((tempBP.gapSiz < 100000) && ( max(dupRegS,bp1RegS) - min(dupRegE,bp1RegE) < 100000 ) && (max(dupRegS,bp2RegS) - min(dupRegE,bp2RegE) < 100000));

		if (tempBP.chr1 > tempBP.chr2){
			swapNums(tempBP.chr1,tempBP.chr2);
			swapNums(tempBP.point1,tempBP.point2);
			swapNums(tempBP.point1_rag,tempBP.point2_rag);
			swapNums(tempBP.point1_bd,tempBP.point2_bd);
			swapNums(tempBP.pointIdx1,tempBP.pointIdx2);
			tempBP.orient1 = '+'+'-'-tempBP.orient1;
			tempBP.orient2 = '+'+'-'-tempBP.orient2;			
			swapNums(tempBP.orient1,tempBP.orient2);
		}
		if (tempBP.chr1 == tempBP.chr2 && tempBP.point2 < tempBP.point1){
			swapNums(tempBP.point1,tempBP.point2);
			swapNums(tempBP.point1_rag,tempBP.point2_rag);
			swapNums(tempBP.point1_bd,tempBP.point2_bd);
			swapNums(tempBP.pointIdx1,tempBP.pointIdx2);
			tempBP.orient1 = '+'+'-'-tempBP.orient1;
			tempBP.orient2 = '+'+'-'-tempBP.orient2;			
			swapNums(tempBP.orient1,tempBP.orient2);
		}


		if (min(calc_score(opticalMap1[i].hitEnum,miniM),calc_score(opticalMap1[i+1].hitEnum,miniM)) > 0)
		{
			tempBP.supportMapId.push_back(atol(opticalMap1[i].mapId));
			tempBP.supportScores.push_back(min(calc_score(opticalMap1[i].hitEnum,miniM),calc_score(opticalMap1[i+1].hitEnum,miniM)));
		}
		
		if (opticalMap1[i].chrId == opticalMap1[i+1].chrId && opticalMap1[i].orientation == opticalMap1[i+1].orientation && (opticalMap1[i].orientation == (opticalMap1[i].refStart < opticalMap1[i+1].refStart)) && ovlpDist - ovlpMol < 100000){
			bool dup_flag = false;
			bool CNV_flag = false;
			bool inv_flag = false;
			bool BP_flag = false;
       			clarify_insert(opticalMap1[i], opticalMap1[i+1], dup_flag, CNV_flag, inv_flag, BP_flag, tempBP, 0, 1);
			if (tempBP.size > 10000 && abs(tempBP.size)>0.05*abs(tempBP.point2-tempBP.point1)){
				if (ovlpRef>10000)dup_flag = true;

				//***Here define the tandem as the duplications in adjacent regions: gap size between split map < 100Kb, the distance between duplicated regions and breakpoints < 100Kb***//
                                if (dup_flag)strcpy(tempBP.type,"Duplication-Split");
                                else if (inv_flag)strcpy(tempBP.type,"Inversion-Split");
                                else if (CNV_flag)strcpy(tempBP.type,"CNV-Split");
				else
					strcpy(tempBP.type,"Large-Insertion");

				candidateSet.push_back(tempBP);
			}
			else {
                                if (inv_flag)strcpy(tempBP.type,"Inversion-Split");
				else if (tempBP.size < -2000 && abs(tempBP.size)>0.05*abs(tempBP.point2-tempBP.point1)) {
					strcpy(tempBP.type,"Large-Deletion");// here call deletion or duplication
				}else {
					strcpy(tempBP.type,"Replacement");
				}
				candidateSet.push_back(tempBP);
                        }
			//it requires the opticalMaps are all + oriented
			//1. chr same, 2. orientation same, 3. gap is small, 4. gap orientation is same
		}
		else
		{
				auto tempBP1 = tempBP;
				auto tempBP2 = tempBP;
				if ((tempBP2.chr1 == opticalMap1[i].chrId) && (opticalMap1[i].refStart == min(tempBP2.point1,min(tempBP2.point1_rag,tempBP2.point1_bd))) && (opticalMap1[i].refEnd == max(tempBP2.point1,max(tempBP2.point1_rag,tempBP2.point1_bd)))){
					tempBP2.chr2 = tempBP2.chr1;
					tempBP2.point2 = tempBP2.point1;
					tempBP2.point2_rag = tempBP2.point1_rag;
					tempBP2.point2_bd = tempBP2.point1_bd;
					tempBP2.pointIdx2 = tempBP2.pointIdx1;
					tempBP1.chr1 = tempBP1.chr2;
					tempBP1.point1 = tempBP1.point2;
					tempBP1.point1_rag = tempBP1.point2_rag;
					tempBP1.point1_bd = tempBP1.point2_bd;
					tempBP1.pointIdx1 = tempBP1.pointIdx2;

				}else{
					tempBP2.chr1 = tempBP2.chr2;
					tempBP2.point1 = tempBP2.point2;
					tempBP2.point1_rag = tempBP2.point2_rag;
					tempBP2.point1_bd = tempBP2.point2_bd;
					tempBP2.pointIdx1 = tempBP2.pointIdx2;
					tempBP1.chr2 = tempBP1.chr1;
					tempBP1.point2 = tempBP1.point1;
					tempBP1.point2_rag = tempBP1.point1_rag;
					tempBP1.point2_bd = tempBP1.point1_bd;
					tempBP1.pointIdx2 = tempBP1.pointIdx1;
				}
				bool dup_flag = false;
				bool CNV_flag = false;
				bool inv_flag = false;
				bool tand_flag = true;
				bool BP_flag = false;
				tempBP2.orient1 = opticalMap1[i].orientation?'+':'-';
				tempBP2.orient2 = opticalMap1[i].orientation?'+':'-';
				tempBP1.orient1 = opticalMap1[i+1].orientation?'+':'-';
				tempBP1.orient2 = opticalMap1[i+1].orientation?'+':'-';
				clarify_insert(opticalMap1[i], opticalMap1[i], dup_flag, CNV_flag, inv_flag, BP_flag,tempBP2, min(min(opticalMap1[i+1].optStartIndex,opticalMap1[i+1].optEndIndex)-1, (LL)opticalMap1[i+1].position.size()-3) , 3);
                                if (dup_flag)strcpy(tempBP2.type,"Duplication");
                                else if (inv_flag)strcpy(tempBP2.type,"Inversion");
                                else if (CNV_flag)strcpy(tempBP2.type,"CNV");
				strcat(tempBP2.type,"-Split");

				if (dup_flag || CNV_flag || inv_flag)splitDup.push_back(tempBP2);
					
				dup_flag = false;
				CNV_flag = false;
				inv_flag = false;
				tand_flag = true;
				BP_flag = false;
				clarify_insert(opticalMap1[i+1], opticalMap1[i+1], dup_flag, CNV_flag, inv_flag,BP_flag, tempBP1, max((LL)1,max(opticalMap1[i].optEndIndex,opticalMap1[i].optStartIndex)+1) , 2);
                                if (dup_flag)strcpy(tempBP1.type,"Duplication");
                                else if (inv_flag)strcpy(tempBP1.type,"Inversion");
                                else if (CNV_flag)strcpy(tempBP1.type,"CNV");
				strcat(tempBP1.type,"-Split");
				if (dup_flag || CNV_flag || inv_flag)splitDup.push_back(tempBP1);

			int cct = 0;
			if (cond1){
				strcpy(tempBP.type,"Inter-TL-BP");
				tempBP.dupStart = 0;
				tempBP.dupEnd = 0;
				tempBP.dupStart_rag = 0;
				tempBP.dupEnd_rag = 0;
				candidateSet.push_back(tempBP);
				cct++;
				
			}else{
				if (cond2||((!cond4)&&cond5)){
					strcpy(tempBP.type,"Intra-TL-BP");
					if (((tempBP.point1 - tempBP.point1_bd)*(tempBP.point2_bd - tempBP.point2) >= 0) && ((tempBP.point1 - tempBP.point1_bd)*(tempBP.point2 - tempBP.point1) >= 0))strcat(tempBP.type,"/Deletion");
					candidateSet.push_back(tempBP);
					cct++;
				}
				if (cond4){
					strcpy(tempBP.type,"Inversion-BP");
					candidateSet.push_back(tempBP);
					cct++;
				}
			}
		}
	}

	printf("So far, there are %lld candidate complicated SV\n",(LL)candidateSet.size());
	for (auto tpBP:splitDup){
		candidateSet.push_back(tpBP);
	}
	printf("After adding the split dup scan, there are %lld candidate complicated SV\n",(LL)candidateSet.size());
        return candidateSet;
}

void printCandidate(vector<breakpoint> candidateSet){
	for(auto tempBP:candidateSet)
		printf("MapId:%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%g\t%lld\t%lld\n",tempBP.supportMapId[0],tempBP.chr1,tempBP.point1,tempBP.chr2,tempBP.point2,tempBP.type,tempBP.size,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size());
}


bool updateDupCoord(breakpoint &tempBP){
        LL p1min = min(tempBP.point1,tempBP.point1_bd);
        LL p1max = max(tempBP.point1,tempBP.point1_bd);
        LL p2min = min(tempBP.point2,tempBP.point2_bd);
        LL p2max = max(tempBP.point2,tempBP.point2_bd);
	if (((p2max-p1min)*(p1max-p2min)>=0)&&(p1min!=p2max)&&(p1max!=p2min)){
		tempBP.dupStart_rag=tempBP.dupStart=max(p1min,p2min);
		tempBP.dupEnd_rag=tempBP.dupEnd=min(p1max,p2max);
		tempBP.size = tempBP.dupEnd-tempBP.dupStart;
                return true;
	}else{
		tempBP.dupEnd_rag=tempBP.dupStart_rag=tempBP.dupStart=tempBP.dupEnd=0;
                return false;
	}
}

void printCandidateToFile(vector<breakpoint> completeSet, FILE* outputFile, LL minS, LL minSS){
	fprintf(outputFile,"#map id\tchr1\tpoint1\tchr2\tpoint2\ttype\taligned_region1\taligned_region2\tDupSeg_start_conf_int\tDupSeg_end_conf_int\tsize(span)\text1\text2\tCopy Size\tCopy Change\textra info\tSupport Map List\n");
        printf("Size of set: %lld\n",(LL)completeSet.size());
        int xi = 0;

        for(auto tempBP:completeSet){
                xi++;
                bool upded=false;
                if((tempBP.chr1==tempBP.chr2)&&(tempBP.point1!=tempBP.point1_bd)&&(tempBP.point2!=tempBP.point2_bd))upded = updateDupCoord(tempBP);
		if (strcmp(tempBP.type,"Replacement") == 0)continue;
		if (tempBP.supportScores.size() < minS)continue;
                char *output = NULL;
                if (strstr(tempBP.type,"Split") && (tempBP.supportScores.size() < minSS)) continue;
		if (strcmp(tempBP.type,"Inter-TL") == 0 || strcmp(tempBP.type,"Intra-TL-BP") == 0 || strcmp(tempBP.type,"Inter-TL-BP") == 0 || strcmp(tempBP.type,"Intra-TL") == 0)tempBP.size = 0;
                if (strcmp(tempBP.type,"Intra-TL") == 0){
                    if (upded){
                        if (strcmp(tempBP.type,"Inversion-Split") == 0 || strcmp(tempBP.type,"Duplication-Split") == 0)
                            strcpy(tempBP.type,"Duplication-Split");
                        else if (strcmp(tempBP.type,"Inversion-BP") == 0 || strcmp(tempBP.type,"Duplication-BP") == 0 || strcmp(tempBP.type,"Intra-TL-BP") == 0)
                            strcpy(tempBP.type,"Duplication-BP");
                        else
                            strcpy(tempBP.type,"Duplication");
                    }
                }
                if (strcmp(tempBP.type,"Intra-TL-BP") == 0){
                    if (abs(tempBP.point1-tempBP.point2)<500000){
                        if (tempBP.dupStart>1)strcpy(tempBP.type,"Duplication-BP");
                        else strcpy(tempBP.type,"Inversion-BP");
                    }
                }

		if (strcmp(tempBP.type,"Inter-TL")==0||strcmp(tempBP.type,"Intra-TL")==0){
			fprintf(outputFile,"%lld\t%lld\t%lld...%lld\t%lld\t%lld-%lld\t%s\t%lld%s",tempBP.supportMapId[0],tempBP.chr1,tempBP.point1,tempBP.point1_bd,tempBP.chr2,tempBP.point2,tempBP.point2_bd,tempBP.type,tempBP.invStart,tempBP.invStart<tempBP.dupStart?">>":"<<");
			if (tempBP.point1 == tempBP.dupStart)
				fprintf(outputFile,"%lld|...|",tempBP.point1);
			else
				if (tempBP.dupStart > tempBP.invStart)
					fprintf(outputFile,"[%lld-%lld]|...|",min(tempBP.dupStart,tempBP.point1),max(tempBP.dupStart,tempBP.point1));
				else
					fprintf(outputFile,"[%lld-%lld]|...|",max(tempBP.dupStart,tempBP.point1),min(tempBP.dupStart,tempBP.point1));
			if (tempBP.point1_bd == tempBP.dupStart_rag)
				fprintf(outputFile,"%lld%s%lld\t|",tempBP.point1_bd,tempBP.point1_bd<tempBP.invEnd?">>":"<<",tempBP.invEnd);
			else
				if (tempBP.point1_bd<tempBP.invEnd)
					fprintf(outputFile,"[%lld-%lld]%s%lld\t|",min(tempBP.point1_bd,tempBP.dupStart_rag),max(tempBP.point1_bd,tempBP.dupStart_rag),tempBP.point1_bd<tempBP.invEnd?">>":"<<",tempBP.invEnd);
				else
					fprintf(outputFile,"[%lld-%lld]%s%lld\t|",max(tempBP.point1_bd,tempBP.dupStart_rag),min(tempBP.point1_bd,tempBP.dupStart_rag),tempBP.point1_bd<tempBP.invEnd?">>":"<<",tempBP.invEnd);
			if (tempBP.point2 == tempBP.dupEnd)
				fprintf(outputFile,"%lld%s",tempBP.point2,tempBP.point2<tempBP.point2_bd?">>":"<<");
			else
				if (tempBP.point2<tempBP.point2_bd)
					fprintf(outputFile,"[%lld-%lld]%s",min(tempBP.point2,tempBP.dupEnd),max(tempBP.point2,tempBP.dupEnd),tempBP.point2<tempBP.point2_bd?">>":"<<");
				else
					fprintf(outputFile,"[%lld-%lld]%s",max(tempBP.point2,tempBP.dupEnd),min(tempBP.point2,tempBP.dupEnd),tempBP.point2<tempBP.point2_bd?">>":"<<");
			if (tempBP.point2_bd == tempBP.dupEnd_rag)
				fprintf(outputFile,"%lld|\t",tempBP.point2_bd);
			else
				if (tempBP.point2 < tempBP.point2_bd)
					fprintf(outputFile,"[%lld-%lld]|\t",min(tempBP.dupEnd_rag,tempBP.point2_bd),max(tempBP.dupEnd_rag,tempBP.point2_bd));
				else
					fprintf(outputFile,"[%lld-%lld]|\t",max(tempBP.dupEnd_rag,tempBP.point2_bd),min(tempBP.dupEnd_rag,tempBP.point2_bd));
			fprintf(outputFile,"0-0\t0-0\t%g\t%lld\t%lld\t%lld\t+%lld\t%s\t|",tempBP.size,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size(),tempBP.cnvSiz,tempBP.cnvNum,tempBP.extInfo);
			for (int t = 0; t < tempBP.supportMapId.size(); t++)fprintf(outputFile,"%lld|",tempBP.supportMapId[t]);
			fprintf(outputFile,"\n");
		}
		else{
			char ori1[5];
			char ori2[5];
			LL rag1S, rag1E, rag2S, rag2E;
			if (tempBP.point1_bd < tempBP.point1){
				strcpy(ori1,">>");
				rag1S = min(tempBP.point1,tempBP.point1_rag);
				rag1E = max(tempBP.point1,tempBP.point1_rag);
			}else{
				strcpy(ori1,"<<");
				rag1S = max(tempBP.point1,tempBP.point1_rag);
				rag1E = min(tempBP.point1,tempBP.point1_rag);
			}
			if (tempBP.point2 < tempBP.point2_bd){
				strcpy(ori2,">>");
				rag2S = min(tempBP.point2,tempBP.point2_rag);
				rag2E = max(tempBP.point2,tempBP.point2_rag);
			}else{
				strcpy(ori2,"<<");
				rag2S = max(tempBP.point2,tempBP.point2_rag);
				rag2E = min(tempBP.point2,tempBP.point2_rag);
			}
                        string typ(tempBP.type);
                        if (typ.find("-Split") != string::npos) {
                            tempBP.point1 = min(tempBP.point1,tempBP.point1_bd);
                            tempBP.point2 = max(tempBP.point2,tempBP.point2_bd);
                        }
			// Narrow duplication regions
                        if (tempBP.dupStart>1&&tempBP.dupEnd>1){tempBP.point1=tempBP.dupStart;tempBP.point1_rag=tempBP.dupStart_rag;tempBP.point2=tempBP.dupEnd;tempBP.point2_rag=tempBP.dupEnd_rag;}
                        tempBP.size = abs(tempBP.point2 - tempBP.point1);
			if (rag1S!=rag1E && rag2S!=rag2E)
				fprintf(outputFile,"%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%lld%s[%lld-%lld]|\t|[%lld-%lld]%s%lld\t%lld-%lld\t%lld-%lld\t%g\t%lld\t%lld\t%lld\t+%lld\t%s\t|",tempBP.supportMapId[0],tempBP.chr1,tempBP.point1,tempBP.chr2,tempBP.point2,tempBP.type,tempBP.point1_bd,ori1,rag1S,rag1E,rag2S,rag2E,ori2,tempBP.point2_bd,min(tempBP.dupStart,tempBP.dupStart_rag),max(tempBP.dupStart,tempBP.dupStart_rag),min(tempBP.dupEnd,tempBP.dupEnd_rag),max(tempBP.dupEnd,tempBP.dupEnd_rag),tempBP.size,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size(),tempBP.cnvSiz,tempBP.cnvNum,tempBP.extInfo);
			else if (rag2S!=rag2E)
				fprintf(outputFile,"%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%lld%s%lld|\t|[%lld-%lld]%s%lld\t%lld-%lld\t%lld-%lld\t%g\t%lld\t%lld\t%lld\t+%lld\t%s\t|",tempBP.supportMapId[0],tempBP.chr1,tempBP.point1,tempBP.chr2,tempBP.point2,tempBP.type,tempBP.point1_bd,ori1,rag1S,rag2S,rag2E,ori2,tempBP.point2_bd,min(tempBP.dupStart,tempBP.dupStart_rag),max(tempBP.dupStart,tempBP.dupStart_rag),min(tempBP.dupEnd,tempBP.dupEnd_rag),max(tempBP.dupEnd,tempBP.dupEnd_rag),tempBP.size,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size(),tempBP.cnvSiz,tempBP.cnvNum,tempBP.extInfo);
			else if (rag1S!=rag1E)
				fprintf(outputFile,"%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%lld%s[%lld-%lld]|\t|%lld%s%lld\t%lld-%lld\t%lld-%lld\t%g\t%lld\t%lld\t%lld\t+%lld\t%s\t|",tempBP.supportMapId[0],tempBP.chr1,tempBP.point1,tempBP.chr2,tempBP.point2,tempBP.type,tempBP.point1_bd,ori1,rag1S,rag1E,rag2S,ori2,tempBP.point2_bd,min(tempBP.dupStart,tempBP.dupStart_rag),max(tempBP.dupStart,tempBP.dupStart_rag),min(tempBP.dupEnd,tempBP.dupEnd_rag),max(tempBP.dupEnd,tempBP.dupEnd_rag),tempBP.size,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size(),tempBP.cnvSiz,tempBP.cnvNum,tempBP.extInfo);
			else
				fprintf(outputFile,"%lld\t%lld\t%lld\t%lld\t%lld\t%s\t%lld%s%lld|\t|%lld%s%lld\t%lld-%lld\t%lld-%lld\t%g\t%lld\t%lld\t%lld\t+%lld\t%s\t|",tempBP.supportMapId[0],tempBP.chr1,tempBP.point1,tempBP.chr2,tempBP.point2,tempBP.type,tempBP.point1_bd,ori1,rag1S,rag2S,ori2,tempBP.point2_bd,min(tempBP.dupStart,tempBP.dupStart_rag),max(tempBP.dupStart,tempBP.dupStart_rag),min(tempBP.dupEnd,tempBP.dupEnd_rag),max(tempBP.dupEnd,tempBP.dupEnd_rag),tempBP.size,(LL)tempBP.supportScores.size(),(LL)tempBP.nonSuppScores.size(),tempBP.cnvSiz,tempBP.cnvNum,tempBP.extInfo);
			for (int t = 0; t < tempBP.supportMapId.size(); t++)fprintf(outputFile,"%lld|",tempBP.supportMapId[t]);
			fprintf(outputFile,"\n");
		}
	}
	fclose(outputFile);
}


vector<breakpoint> mergeCandidateSVBreakpoint(vector<breakpoint> candidateSet, LL shift){
        vector<breakpoint> completeSet;
	completeSet.clear();
        printf("Let us check, candidateSet size: %lld\n",(LL)candidateSet.size());
	sort(candidateSet.begin(),candidateSet.end(),st);
	//printCandidate();
	LL pID1=-1,pID2=-1, pP1=-1, pP2=-1;
	bool orient1, orient2;
	char type[100] = "None";
	LL cnt = 0;
	LL shift_inv = 100000;
	LL shift_tl = 100000;
	shift = 100000;
	for (LL i = 0; i < candidateSet.size(); i++)
	{
                if (candidateSet[i].supportMapId.size()<=0)printf("How it could be!!!!!! %lld\n", i);
		if (candidateSet[i].supportScores.size()>0 && pID1 == candidateSet[i].chr1 && pID2 == candidateSet[i].chr2 && strcmp(candidateSet[i].type,type)==0 && (orient1 == (candidateSet[i].point1 > candidateSet[i].point1_bd)) && (orient2 == (candidateSet[i].point2_bd > candidateSet[i].point2) && abs(pP1-candidateSet[i].point1) <= shift && abs(pP2-candidateSet[i].point2) <= shift))
		{
			bool flag = true;
                        if (completeSet.size()<cnt)printf("What is it??????\n");
			for (LL k = 0; k < completeSet[cnt-1].supportMapId.size(); k++){
				if (completeSet[cnt-1].supportMapId[k]==candidateSet[i].supportMapId[0]){
					if (completeSet[cnt-1].supportScores[k] < candidateSet[i].supportScores[0]){
						completeSet[cnt-1].supportScores[k] = candidateSet[i].supportScores[0];
					}
					flag = false;
					break;
				}
			}
			if (flag){
				completeSet[cnt-1].supportScores.push_back(candidateSet[i].supportScores[0]);
				completeSet[cnt-1].supportMapId.push_back(candidateSet[i].supportMapId[0]);
			}
		}
		else
		{
			completeSet.push_back(candidateSet[i]);
			cnt++;
			pID1 = candidateSet[i].chr1;
			pID2 = candidateSet[i].chr2;
			pP1 = candidateSet[i].point1;
			pP2 = candidateSet[i].point2;
			orient1 = candidateSet[i].point1 > candidateSet[i].point1_bd;
			orient2 = candidateSet[i].point2_bd > candidateSet[i].point2;
			strcpy(type,candidateSet[i].type);
		}
	}
	printf("After merge the overlapping SVs/breakpoints with same type, there are %lld unique SV/breakpoints\n",(LL)completeSet.size());
	printf("After merge the overlapping SVs/breakpoints with same type, there are %lld unique SV/breakpoints\n",(LL)candidateSet.size());
        return completeSet;
}
vector<breakpoint> mergeCandidateBreakpoint(vector<breakpoint> candidateSet, LL shift){
        vector<breakpoint> completeSet;
        completeSet.clear();
	sort(candidateSet.begin(),candidateSet.end(),st);
	//printCandidate();
	LL pID1=-1,pID2=-1, pP1=-1, pP2=-1;
	bool orient1, orient2;
	char type[100] = "-";
	LL cnt = 0;
	LL shift_inv = 100000;
	LL shift_tl = 100000;
	shift = 100000;
	for (LL i = 0; i < candidateSet.size(); i++)
	{
		if (candidateSet[i].supportScores.size()>0 && pID1 == candidateSet[i].chr1 && pID2 == candidateSet[i].chr2 && strcmp(candidateSet[i].type,type)==0 && (orient1 == (candidateSet[i].point1 > candidateSet[i].point1_bd)) && (orient2 == (candidateSet[i].point2_bd > candidateSet[i].point2) && abs(pP1-candidateSet[i].point1) <= shift && abs(pP2-candidateSet[i].point2) <= shift))
		{
			bool flag = true;
			for (LL k = 0; k < completeSet[cnt-1].supportMapId.size(); k++){
				if (completeSet[cnt-1].supportMapId[k]==candidateSet[i].supportMapId[0]){
					if (completeSet[cnt-1].supportScores[k] < candidateSet[i].supportScores[0]){
						completeSet[cnt-1].supportScores[k] = candidateSet[i].supportScores[0];
					}
					flag = false;
					break;
				}
			}
			if (flag){
				completeSet[cnt-1].supportScores.push_back(candidateSet[i].supportScores[0]);
				completeSet[cnt-1].supportMapId.push_back(candidateSet[i].supportMapId[0]);
			}
		}
		else
		{
			completeSet.push_back(candidateSet[i]);
			cnt++;
			pID1 = candidateSet[i].chr1;
			pID2 = candidateSet[i].chr2;
			pP1 = candidateSet[i].point1;
			pP2 = candidateSet[i].point2;
			orient1 = candidateSet[i].point1 > candidateSet[i].point1_bd;
			orient2 = candidateSet[i].point2_bd > candidateSet[i].point2;
			strcpy(type,candidateSet[i].type);
		}
	}
	printf("After merge the overlapping breakpoints with same type, there are %lld unique breakpoints\n",(LL)completeSet.size());
        return completeSet;
}
vector<breakpoint> assembleBreakpoint(vector<breakpoint> candidateSet, vector<opticalMapType> opticalMap1, LL shift){
        vector<breakpoint> completeSet;
	completeSet.clear();
	sort(candidateSet.begin(),candidateSet.end(),st);
	//printCandidate();
	LL pID1=-1,pID2=-1, pP1=-1, pP2=-1;
	bool orient1, orient2;
	char type[100] = "None";
	LL cnt = 0;
	LL shift_inv = 100000;
	LL shift_tl = 100000;
	shift = shift_inv;
	for (LL i = 0; i < candidateSet.size(); )
	{
		if (i == candidateSet.size() - 1){
			completeSet.push_back(candidateSet[i]);
			break;
		}
		if ((candidateSet[i].chr1 == candidateSet[i+1].chr1) && (candidateSet[i].chr2 == candidateSet[i+1].chr2) && (strcmp(candidateSet[i].type,candidateSet[i+1].type)==0)){
			if ((abs(candidateSet[i].point1-candidateSet[i+1].point1)<shift) && (abs(candidateSet[i].point2-candidateSet[i+1].point2)<shift) && ((candidateSet[i].point1_bd - candidateSet[i].point1) != (candidateSet[i+1].point1_bd-candidateSet[i+1].point1)) && ((candidateSet[i].point2_bd-candidateSet[i].point2)!=(candidateSet[i+1].point2_bd-candidateSet[i+1].point2))){
				if ((strcmp(candidateSet[i].type,"Inversion-BP")==0) || ((strcmp(candidateSet[i].type,"Intra-TL-BP")==0)&&(candidateSet[i].orient1!=candidateSet[i].orient2))){
					auto tempBP = candidateSet[i];
					tempBP.point1 = max(candidateSet[i].point1,candidateSet[i+1].point1);
					tempBP.point1_rag = min(candidateSet[i].point1,candidateSet[i+1].point1);
					tempBP.point1_bd = min(min(candidateSet[i].point1_bd,candidateSet[i].point2_bd),min(candidateSet[i+1].point1_bd,candidateSet[i+1].point2_bd));
					tempBP.point2 = min(candidateSet[i].point2,candidateSet[i+1].point2);
					tempBP.point2_rag = max(candidateSet[i].point2,candidateSet[i+1].point2);
					tempBP.point2_bd = max(max(candidateSet[i].point1_bd,candidateSet[i].point2_bd),max(candidateSet[i+1].point1_bd,candidateSet[i+1].point2_bd));
					strcpy(tempBP.type,"Inversion");
					tempBP.orient1 = '+';
					tempBP.orient2 = '+';
					bool existMap;
					for (int s = 0; s < candidateSet[i+1].supportMapId.size(); s++){
						existMap = false;
						for (int t = 0; t < tempBP.supportMapId.size(); t++){
							if (candidateSet[i+1].supportMapId[s] == tempBP.supportMapId[t]){
								existMap = true;
								break;
							}
						}
						if (!existMap)tempBP.supportMapId.push_back(candidateSet[i+1].supportMapId[s]);
					}
					completeSet.push_back(tempBP);
					i++;
				}
				else{
					completeSet.push_back(candidateSet[i]);
				}
			}
			else{
				completeSet.push_back(candidateSet[i]);
			}
		}else{
			completeSet.push_back(candidateSet[i]);
		}
		i++;
	}
	printf("After breakpoint assembly(inversion and intra-translocation), there are %lld candidate SV/breakpoints\n",(LL)completeSet.size());
        return completeSet;
}
vector<breakpoint> mergeCandidateType(vector<breakpoint> candidateSet, vector<opticalMapType> opticalMap1,LL shift){
        vector<breakpoint> completeSet;
	completeSet.clear();
	LL pID1=-1,pID2=-1, pP1=-1, pP2=-1;
	char type[100] = "None";
	LL cnt = 0;
	shift = 0;
	for (LL i = 0; i < candidateSet.size(); )
	{
		if (i == candidateSet.size() - 1){
			completeSet.push_back(candidateSet[i]);
			break;
		}		
		if (candidateSet[i].chr1 == candidateSet[i+1].chr1 && candidateSet[i].chr2 == candidateSet[i+1].chr2 && candidateSet[i].supportMapId[0] == candidateSet[i+1].supportMapId[0] && strcmp(candidateSet[i].type,candidateSet[i+1].type)==0){
			if (strcmp(candidateSet[i].type,"Inversion-BP")==0){
				if ((candidateSet[i].point1-candidateSet[i+1].point2)*(candidateSet[i+1].point1-candidateSet[i].point2)>=0){
					
					auto tempBP = candidateSet[i];
					tempBP.point1 = max(candidateSet[i].point1,candidateSet[i+1].point1);
					tempBP.point1_rag = min(candidateSet[i].point1,candidateSet[i+1].point1);
					tempBP.point1_bd = min(min(candidateSet[i].point1_bd,candidateSet[i].point2_bd),min(candidateSet[i+1].point1_bd,candidateSet[i+1].point2_bd));
					tempBP.point2 = min(candidateSet[i].point2,candidateSet[i+1].point2);
					tempBP.point2_rag = max(candidateSet[i].point2,candidateSet[i+1].point2);
					tempBP.point2_bd = max(max(candidateSet[i].point1_bd,candidateSet[i].point2_bd),max(candidateSet[i+1].point1_bd,candidateSet[i+1].point2_bd));
					strcpy(tempBP.type,"Inversion");
					tempBP.orient1 = '+';
					tempBP.orient2 = '+';
                                        //the above segment is used to exclude the self-inverse segments
					completeSet.push_back(tempBP);
					i++;
				}
				else
					completeSet.push_back(candidateSet[i]);
			}
			else if (strcmp(candidateSet[i].type,"Inter-TL-BP")==0 || strcmp(candidateSet[i].type,"Intra-TL-BP")==0){
				bool swap1 = false, swap2 = false;
				if (candidateSet[i].chr1 == candidateSet[i+1].chr1 && (candidateSet[i].point1_bd == candidateSet[i+1].point1 || candidateSet[i].point1_bd == candidateSet[i+1].point1_rag) && (candidateSet[i+1].point1_bd == candidateSet[i].point1 || candidateSet[i+1].point1_bd == candidateSet[i].point1_rag)){
				// the segment in chromosome 1 has been translocated to chromosome 2
					swap1 = true;

				}
				else if (candidateSet[i].chr2 == candidateSet[i+1].chr2 && (candidateSet[i].point2_bd == candidateSet[i+1].point2 || candidateSet[i].point2_bd == candidateSet[i+1].point2_rag) && (candidateSet[i+1].point2_bd == candidateSet[i].point2 || candidateSet[i+1].point2_bd == candidateSet[i].point2_rag)){
				// the segment in chromosome 1 has been translocated to chromosome 2
					swap2 = true;
				}
				else if (candidateSet[i].chr1 == candidateSet[i+1].chr2 && (candidateSet[i].point1_bd == candidateSet[i+1].point2 || candidateSet[i].point1_bd == candidateSet[i+1].point2_rag) && (candidateSet[i+1].point2_bd == candidateSet[i].point1 || candidateSet[i+1].point2_bd == candidateSet[i].point1_rag)){
					swap1 = true;
					swap2 = true;
				}
				if (candidateSet[i].gapSiz + candidateSet[i+1].gapSiz<1000000){//a translocation requires the gaps (unaligned) between flanking regions < 1Mb
					if (swap1){
						swapNums(candidateSet[i].chr1,candidateSet[i].chr2);
						swapNums(candidateSet[i].point1,candidateSet[i].point2);
						swapNums(candidateSet[i].point1_rag,candidateSet[i].point2_rag);
						swapNums(candidateSet[i].point1_bd,candidateSet[i].point2_bd);
						swapNums(candidateSet[i].pointIdx1,candidateSet[i].pointIdx2);
						candidateSet[i].orient1 = '+'+'-'-candidateSet[i].orient1;
						candidateSet[i].orient2 = '+'+'-'-candidateSet[i].orient2;			
						swapNums(candidateSet[i].orient1,candidateSet[i].orient2);
					}
					if (swap2){
						swapNums(candidateSet[i+1].chr1,candidateSet[i+1].chr2);
						swapNums(candidateSet[i+1].point1,candidateSet[i+1].point2);
						swapNums(candidateSet[i+1].point1_rag,candidateSet[i+1].point2_rag);
						swapNums(candidateSet[i+1].point1_bd,candidateSet[i+1].point2_bd);
						swapNums(candidateSet[i+1].pointIdx1,candidateSet[i+1].pointIdx2);
						candidateSet[i+1].orient1 = '+'+'-'-candidateSet[i+1].orient1;
						candidateSet[i+1].orient2 = '+'+'-'-candidateSet[i+1].orient2;			
						swapNums(candidateSet[i+1].orient1,candidateSet[i+1].orient2);
					}
					auto tempBP = candidateSet[i];
					tempBP.invStart = candidateSet[i].point1_bd;
					tempBP.invEnd = candidateSet[i+1].point2_bd;
					tempBP.dupStart = candidateSet[i].point1_rag;
					tempBP.dupEnd = candidateSet[i].point2_rag;
					tempBP.dupStart_rag = candidateSet[i+1].point2_rag;
					tempBP.dupEnd_rag = candidateSet[i+1].point1_rag;
					tempBP.point1_bd = candidateSet[i+1].point2;
					tempBP.orient2 = candidateSet[i+1].orient2;
					tempBP.point2_bd = candidateSet[i+1].point1;
					tempBP.type[strlen(tempBP.type)-1] = '\0';
					tempBP.type[strlen(tempBP.type)-1] = '\0';
					tempBP.type[strlen(tempBP.type)-1] = '\0';
					i++;
					completeSet.push_back(tempBP);					
				}else{
					completeSet.push_back(candidateSet[i]);
				}
			}else{
				completeSet.push_back(candidateSet[i]);
			}	
		}else{
			completeSet.push_back(candidateSet[i]);
		}
		i++;
	}
	printf("After merge the breakpoint pairs supported by the same molecules, there are %lld candidate breakpoints\n",(LL)completeSet.size());
        return completeSet;
}

vector<breakpoint> extraDuplicationDetect(vector<breakpoint> candidateSet, vector<opticalMapType> opticalMap1,LL shift){
	//go first, detect extra duplications by combining multiple candidates in the same contig
	vector<breakpoint> completeSet = candidateSet;
	LL pID1=-1,pID2=-1, pP1=-1, pP2=-1;
	char type[100] = "None";
	LL cnt = 0;
	shift = 0;
	LL grpS = 0, grpE = -1;
	for (LL i = 0; i < candidateSet.size()-1; i++)
	{
		if (strstr(candidateSet[i].type,"_"))break;  // do not handle the extra dup found in the gaps between split map
		if (candidateSet[i].supportMapId[0] == candidateSet[i+1].supportMapId[0]){
			grpE = i+1;
			continue;
		}
		if (grpE - grpS >= 1){
			for (LL j = grpS; j <= grpE-1; j++){
				for (LL l = j+1; l <= grpE; l++){
				//don't worry about the multiply aligned contigs, there must be unique contigs in each split map
					if ((candidateSet[j].chr1 == candidateSet[l].chr1) && (candidateSet[j].chr2 == candidateSet[l].chr2) && (candidateSet[j].point1 == candidateSet[l].point1) && (candidateSet[j].point2 == candidateSet[l].point2))
						continue;
					LL v1s = -1, v2s = -1;
					int rt1, rt2, rt3, rt4;
					for (LL k = j+1; k <= l; k++){
						if ((candidateSet[k].chr1 == candidateSet[k-1].chr1) && (candidateSet[k].chr2 == candidateSet[k-1].chr2) && (candidateSet[k].point1 == candidateSet[k-1].point1) && (candidateSet[k].point2 == candidateSet[k-1].point2))
							continue;
						rt1 = 0;
						rt2 = 0;
						rt3 = 0;
						rt4 = 0;
						if (candidateSet[k-1].chr1 == candidateSet[k].chr1 && (candidateSet[k-1].point1_bd == candidateSet[k].point1 || candidateSet[k-1].point1_bd == candidateSet[k].point1_rag) && (candidateSet[k].point1_bd == candidateSet[k-1].point1 || candidateSet[k].point1_bd == candidateSet[k-1].point1_rag)){
							rt1 = 1;
						}
						if (candidateSet[k-1].chr2 == candidateSet[k].chr2 && (candidateSet[k-1].point2_bd == candidateSet[k].point2 || candidateSet[k-1].point2_bd == candidateSet[k].point2_rag) && (candidateSet[k].point2_bd ==candidateSet[k-1].point2 || candidateSet[k].point2_bd == candidateSet[k-1].point2_rag)){
							rt2 = 1;
						}
						if (candidateSet[k-1].chr1 == candidateSet[k].chr2 && (candidateSet[k-1].point1_bd == candidateSet[k].point2 || candidateSet[k-1].point1_bd == candidateSet[k].point2_rag) && (candidateSet[k].point2_bd == candidateSet[k-1].point1 || candidateSet[k].point2_bd == candidateSet[k-1].point1_rag)){
							rt3 = 1;
						}
						if (rt1+rt2+rt3==0){
							rt4 = 1;
						}
						if (v2s == 1){
							rt2 = 0;
							rt4 = 0;
						}else if (v2s == 2){
							rt1 = 0;
							rt3 = 0;
						}
						if (rt1==1){
							if (v1s == -1)v1s = 2;
							if (v2s == -1)v2s = 2;
							else if (v2s == 2)printf("Wrong linkage methods: route 1!\n");
							else if (v2s == 1)v2s = 2;
						}else if (rt2 == 1){
							if (v1s == -1)v1s = 1;
							if (v2s == -1)v2s = 1;
							else if (v2s == 1)printf("Wrong linkage methods: route 2!\n");
							else if (v2s == 2)v2s = 1;
						}else if (rt3 == 1){
							if (v1s == -1)v1s = 2;
							if (v2s == -1)v2s = 1;
							else if (v2s == 2)printf("Wrong linkage methods: route 3!\n");
							else if (v2s == 1)v2s = 1;
						}else{
							if (v1s == -1)v1s = 1;
							if (v2s == -1)v2s = 2;
							else if (v2s == 1)printf("Wrong linkage methods: route 4!\n");
							else if (v2s == 2)v2s = 2;
						}
						
					}
					auto tempBP = candidateSet[j];
					if (v1s == 2){
						tempBP.chr1 = tempBP.chr2;
						tempBP.point1 = tempBP.point2;
						tempBP.point1_bd = tempBP.point2_bd;
						tempBP.point1_rag = tempBP.point2_rag;
						tempBP.pointIdx1 = tempBP.pointIdx2;
					}
					if (v2s == 1){
						tempBP.chr2 = candidateSet[l].chr1;
						tempBP.point2 = candidateSet[l].point1;
						tempBP.point2_bd = candidateSet[l].point1_bd;
						tempBP.point2_rag = candidateSet[l].point1_rag;
						tempBP.pointIdx2 = candidateSet[l].pointIdx1;
					}else{
						tempBP.chr2 = candidateSet[l].chr2;
						tempBP.point2 = candidateSet[l].point2;
						tempBP.point2_bd = candidateSet[l].point2_bd;
						tempBP.point2_rag = candidateSet[l].point2_rag;
						tempBP.pointIdx2 = candidateSet[l].pointIdx2;
					}
					tempBP.gapSiz = 0;
					strcpy(tempBP.type,"Duplication");
					LL dupS1 = min(tempBP.point1,tempBP.point1_bd);
					LL dupE1 = max(tempBP.point1,tempBP.point1_bd);
					LL dupS2 = min(tempBP.point2,tempBP.point2_bd);
					LL dupE2 = max(tempBP.point2,tempBP.point2_bd);
					if (tempBP.chr1 == tempBP.chr2 && min(dupE1,dupE2) - max(dupS1,dupS2) > 10000){//10Kb duplications in the same chromosome
						tempBP.dupStart = max(dupS1,dupS2);
						tempBP.dupEnd = min(dupE1,dupE2);
						tempBP.dupStart_rag = tempBP.dupStart;
						tempBP.dupEnd_rag = tempBP.dupEnd;
						if (tempBP.point1 == tempBP.dupStart)tempBP.dupStart_rag = tempBP.point1_rag;
						else if (tempBP.point2 == tempBP.dupStart)tempBP.dupStart_rag = tempBP.point2_rag;
						if (tempBP.point1 == tempBP.dupEnd)tempBP.dupEnd_rag = tempBP.point1_rag;
						else if (tempBP.point2 == tempBP.dupEnd)tempBP.dupEnd_rag = tempBP.point2_rag;
						if ((tempBP.point1 - tempBP.point1_bd) * (tempBP.point2_bd - tempBP.point2) > 0) strcpy(tempBP.type,"Duplication");
						else strcpy(tempBP.type,"Invert-Duplication");
						completeSet.push_back(tempBP);
					}
				}
			}
		}
		grpS = i+1;
	}
	printf("After merge the breakpoint pairs supported by the same molecules, there are %lld candidate breakpoints\n",(LL)completeSet.size());
        return completeSet;
}

bool feq(double x, double y){
	//return fabs(x - y) <= eps;
	return fabs(x - y) <= eps*(min(fabs(x),fabs(y)));//should be normalized by the values of x and y
}

bool feq1(double x, double y){
        return fabs(x - y) <= eps;
       // return fabs(x - y) <= eps*fabs(min(x,y));//should be normalized by the values of x and y
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

char outputFolder[1000]; 

bool overlapped(LL x1, LL y1, LL x2, LL y2){
	if (x1 >= x2 && x1 <= y2) return true;
	if (y1 >= x2 && y1 <= y2) return true;
	if (x2 >= x1 && x2 <= y1) return true;
	if (y2 >= x1 && y2 <= y1) return true;
	return false;
}


int main(int argv, char *argc[]){
	printf("\n");
	//setDefault();
	LL miniM = 3;
	LL minSupp = 1;
	LL minSuppSplit = 1;
	LL shift = 100000;
        char chrMapFile[1000];
	char inputBN_CR[1000];
	char SVoutputFile2[1000];
	char inputHigh[1000]="nothing.txt";
        double minSegLen = 10000;
        int minSites = 3;
	if (argv == 1){
		printf("\nThe valid parameters are described as follows:\n");
		printf("\t-highDensityFile: \n\t\t  The file of inaccessible regions.\n");
		printf("\t-SVoutputFile_CR: \n\t\t The file name of SVs (.osv).\n");
		printf("\t-optCRAlignFile: \n\t\t  The file name of alignment file(.oma).\n");
		printf("\t-refMapFile: \n\t\t The file name of reference map file(.cmap).\n");
		printf("\t-minSegLen: \n\t\t Default value: %g. The minimum length of split map required to be kept.\n",minSegLen);
		printf("\t-minSites: \n\t\t Default value: %d. The minimum number of M required in HitEnum for an alignment to be kept.\n",minSites);
		printf("\t-minScore: \n\t\t  The minimum score of alignments to be kept.\n");
		printf("\t-minM: \n\t\t Default value: %lld. The minimum flanking Matches in HitEnum.\n",miniM);
		printf("\t-minS: \n\t\t Default value: %lld. The minimum number of molecules covering the break points.\n",minSupp);
		printf("\t-minSS: \n\t\t Default value: %lld. The minimum number of molecules covering the break points if they are detected by re-alignment (with flag -realigned).\n",minSuppSplit);
		return -1;
	}
	for (int i = 1; i < argv; i=i+2){
		string temp(argc[i]);
		if (temp.compare("-refMapFile")==0)
			strcpy(chrMapFile, argc[i+1]);
		else if (temp.compare("-optCRAlignFile")==0)
			strcpy(inputBN_CR, argc[i+1]);
		else if (temp.compare("-SVoutputFile_CR")==0)
			strcpy(SVoutputFile2,argc[i+1]);
		else if (temp.compare("-highDensityFile")==0)
			strcpy(inputHigh, argc[i+1]);
		else if (temp.compare("-minScore")==0)
			scoreLimit = atof(argc[i+1]);
		else if (temp.compare("-minSegLen")==0)
			minSegLen = atof(argc[i+1]);
		else if (temp.compare("-minSites")==0)
			minSites = atoi(argc[i+1]);
		else if (temp.compare("-minM")==0)
			miniM = atol(argc[i+1]);
		else if (temp.compare("-minS")==0)
			minSupp = atol(argc[i+1]);
		else if (temp.compare("-minSS")==0)
			minSuppSplit = atol(argc[i+1]);
		else
			printf("No such parameter or wrong : %s\n",argc[i]);
	}
	readHighDensity(inputHigh);


        genome.clear();
        int curChr = 0;
        getChromosomeList(chrMapFile);
        chrType eptChr = returnEmptyChr();
        for (LL tempChrId = 0; tempChrId < (LL)listOfChromosome.size(); tempChrId++){
            int chrId = listOfChromosome[tempChrId];
            for (int k = curChr; k < chrId-1; k++)
                genome.push_back(eptChr);
            curChr = chrId;
            genome.push_back(readChromosomeInfo(chrId,chrMapFile));
            cout << "The size of chromosome " << chrId << " is " << genome[genome.size()-1].numberOfSites << endl;
        }
	vector<opticalMapType> opticalMap3 = readSourceFile(inputBN_CR,minSegLen,minSites);
        printf("Done readSourceFile!\n");
        vector<breakpoint> candidateSet; 
	candidateSet.clear();

	candidateSet = scanSplitedMap(opticalMap3, miniM);
        printf("Done scanSplitedMap!\n");
	candidateSet = extraDuplicationDetect(candidateSet, opticalMap3,shift);
//        printCandidate(candidateSet);


	candidateSet = mergeCandidateType(candidateSet, opticalMap3,shift);
        sort(candidateSet.begin(),candidateSet.end(),st);
//        printCandidate(candidateSet);

	candidateSet = mergeCandidateBreakpoint(candidateSet, shift);
//        printCandidate(candidateSet);

	candidateSet = assembleBreakpoint(candidateSet, opticalMap3,shift);
	sort(candidateSet.begin(),candidateSet.end(),ss);
	FILE* outputCR = fopen(SVoutputFile2,"w");
	printCandidateToFile(candidateSet,outputCR,minSupp,minSuppSplit);
	return 0;
}
